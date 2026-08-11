nextflow.enable.dsl=2


include { VALIDATE_INPUTS } from './modules/local/validate_inputs'
include { CATALOG_ADAPTER } from './modules/local/catalog_adapter'
include { TRGT } from './modules/local/trgt'
include { LONGTR } from './modules/local/longtr'
include { ATARVA } from './modules/local/atarva'
include { STRDUST } from './modules/local/strdust'
include { NORMALIZE_TRGT } from './modules/local/normalize_trgt'
include { NORMALIZE_LONGTR } from './modules/local/normalize_longtr'
include { NORMALIZE_ATARVA } from './modules/local/normalize_atarva'
include { NORMALIZE_STRDUST } from './modules/local/normalize_strdust'
include { CONSENSUS } from './modules/local/consensus'


def normalizeCallers(value) {
    if (value == null) {
        error "--callers is required (for example: --callers trgt,longtr)"
    }
    def names = (value instanceof Collection ? value : value.toString().split(/,/)).collect { it.toString().trim().toLowerCase() }.findAll { it }
    if (!names) {
        error '--callers must not be empty'
    }
    def unknown = names.findAll { !(it ==~ /trgt|longtr|atarva|strdust/) }
    if (unknown) {
        error(/Unknown caller name; allowed: trgt,longtr,atarva,strdust/)
    }
    if (names.size() != names.toSet().size()) {
        error "Duplicate caller name(s) are not allowed: ${names.join(',')}"
    }
    return names
}

workflow {
    callers = normalizeCallers(params.callers)
    required = ['samplesheet', 'reference', 'reference_fai', 'locus_manifest']
    required.each { key -> if (!params[key]) error "--${key.replace('_','-')} is required" }

    reference = file(params.reference, checkIfExists: true)
    referenceFai = file(params.reference_fai, checkIfExists: true)
    manifest = file(params.locus_manifest, checkIfExists: true)
    sampleSheet = file(params.samplesheet, checkIfExists: true)
    trgtConfig = file(params.trgt_config, checkIfExists: true)
    longtrConfig = file(params.longtr_config, checkIfExists: true)
    atarvaConfig = file(params.atarva_config, checkIfExists: true)
    strdustConfig = file(params.strdust_config, checkIfExists: true)
    consensusConfig = file(params.consensus_config, checkIfExists: true)
    CATALOG_ADAPTER(Channel.of(tuple(manifest, referenceFai, 'trgt,longtr,atarva,strdust')))
    catalogAssets = CATALOG_ADAPTER.out.catalogs.map { catalogDir ->
        tuple(
            catalogDir.resolve('trgt.bed'),
            catalogDir.resolve('longtr.bed'),
            catalogDir.resolve('atarva.bed.gz'),
            catalogDir.resolve('atarva.bed.gz.tbi'),
            catalogDir.resolve('strdust.bed')
        )
    }

    validationInput = catalogAssets.map { trgtCatalog, longtrCatalog, atarvaCatalog, atarvaCatalogIndex, strdustCatalog ->
        tuple(
            sampleSheet, reference, referenceFai, manifest,
            trgtConfig, trgtCatalog, longtrConfig, longtrCatalog,
            atarvaConfig, atarvaCatalog, atarvaCatalogIndex,
            strdustConfig, strdustCatalog, consensusConfig, callers.join(',')
        )
    }
    alignmentAssets = Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> [file(row.alignment.toString(), checkIfExists: true), file(row.alignment_index.toString(), checkIfExists: true)] }
        .flatten()
        .collect()
    VALIDATE_INPUTS(validationInput, alignmentAssets)
    validationGate = VALIDATE_INPUTS.out.validated
    trgtCatalog = catalogAssets.map { it[0] }
    longtrCatalog = catalogAssets.map { it[1] }
    atarvaCatalog = catalogAssets.map { it[2] }
    atarvaCatalogIndex = catalogAssets.map { it[3] }
    strdustCatalog = catalogAssets.map { it[4] }

    samples = Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample_id.toString(), row.platform.toString().toLowerCase(),
                file(row.alignment.toString(), checkIfExists: true),
                file(row.alignment_index.toString(), checkIfExists: true),
                row.sex ? row.sex.toString() : '', row.ploidy ? row.ploidy.toString() : ''
            )
        }

    TRGT(samples.filter { it[1] == 'hifi' && callers.contains('trgt') }
        .combine(validationGate).combine(trgtCatalog)
        .map { sampleId, platform, alignment, alignmentIndex, sex, ploidy, gate, catalog ->
            tuple(sampleId, platform, alignment, alignmentIndex, sex, ploidy, reference, referenceFai, catalog, trgtConfig, gate)
        })
    LONGTR(samples.filter { it[1] in ['hifi', 'ont'] && callers.contains('longtr') }
        .combine(validationGate).combine(longtrCatalog)
        .map { sampleId, platform, alignment, alignmentIndex, sex, ploidy, gate, catalog ->
            tuple(sampleId, platform, alignment, alignmentIndex, sex, ploidy, reference, referenceFai, catalog, longtrConfig, gate)
        })
    ATARVA(samples.filter { it[1] in ['hifi', 'ont'] && callers.contains('atarva') }
        .combine(validationGate).combine(atarvaCatalog).combine(atarvaCatalogIndex)
        .map { sampleId, platform, alignment, alignmentIndex, sex, ploidy, gate, catalog, catalogIndex ->
            tuple(sampleId, platform, alignment, alignmentIndex, sex, ploidy, reference, referenceFai, catalog, catalogIndex, atarvaConfig, gate)
        })
    STRDUST(samples.filter { it[1] in ['hifi', 'ont'] && callers.contains('strdust') }
        .combine(validationGate).combine(strdustCatalog)
        .map { sampleId, platform, alignment, alignmentIndex, sex, ploidy, gate, catalog ->
            tuple(sampleId, platform, alignment, alignmentIndex, sex, ploidy, reference, referenceFai, catalog, strdustConfig, gate)
        })

    NORMALIZE_TRGT(TRGT.out.calls.map { s, p, native_dir -> tuple(s, p, native_dir, manifest) })
    NORMALIZE_LONGTR(LONGTR.out.calls.map { s, p, native_dir -> tuple(s, p, native_dir, manifest) })
    NORMALIZE_ATARVA(ATARVA.out.calls.map { s, p, native_dir -> tuple(s, p, native_dir, manifest) })
    NORMALIZE_STRDUST(STRDUST.out.calls.map { s, p, native_dir -> tuple(s, p, native_dir, manifest) })

    normalized = NORMALIZE_TRGT.out.calls
        .mix(NORMALIZE_LONGTR.out.calls)
        .mix(NORMALIZE_ATARVA.out.calls)
        .mix(NORMALIZE_STRDUST.out.calls)
    consensusAssets = Channel.of(tuple(manifest, sampleSheet, consensusConfig, callers.join(',')))
    CONSENSUS(normalized.map { sample_id, platform, normalized_file -> normalized_file }.collect(), consensusAssets)
}
