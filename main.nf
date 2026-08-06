nextflow.enable.dsl=2


include { VALIDATE_INPUTS } from './modules/local/validate_inputs'
include { TRGT } from './modules/local/trgt'
include { LONGTR } from './modules/local/longtr'
include { ATARVA } from './modules/local/atarva'
include { STRDUST } from './modules/local/strdust'
include { NORMALIZE_TRGT } from './modules/local/normalize_trgt'
include { NORMALIZE_LONGTR } from './modules/local/normalize_longtr'
include { NORMALIZE_ATARVA } from './modules/local/normalize_atarva'
include { NORMALIZE_STRDUST } from './modules/local/normalize_strdust'
include { CONSENSUS } from './modules/local/consensus'

def supportedCallers = ['trgt', 'longtr', 'atarva', 'strdust']

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

def loadCallerConfig(String configPath) {
    def cfgFile = new File(configPath)
    if (!cfgFile.exists()) error "Caller configuration does not exist: ${configPath}"
    return new groovy.yaml.YamlSlurper().parse(cfgFile)
}

def resolveConfigPath(String configPath, String childPath) {
    if (!childPath) error "native_catalog is required in ${configPath}"
    def candidate = new File(childPath)
    return (candidate.isAbsolute() ? candidate : new File(new File(configPath).parentFile, childPath)).canonicalPath
}

workflow {
    callers = normalizeCallers(params.callers)
    required = ['samplesheet', 'reference', 'reference_fai', 'locus_manifest']
    required.each { key -> if (!params[key]) error "--${key.replace('_','-')} is required" }

    trgtCfg = loadCallerConfig(params.trgt_config)
    longtrCfg = loadCallerConfig(params.longtr_config)
    atarvaCfg = loadCallerConfig(params.atarva_config)
    strdustCfg = loadCallerConfig(params.strdust_config)
    params.trgt_container = params.trgt_container ?: trgtCfg.container
    params.longtr_container = params.longtr_container ?: longtrCfg.container
    params.atarva_container = params.atarva_container ?: atarvaCfg.container
    params.strdust_container = params.strdust_container ?: strdustCfg.container

    reference = file(params.reference, checkIfExists: true)
    referenceFai = file(params.reference_fai, checkIfExists: true)
    manifest = file(params.locus_manifest, checkIfExists: true)
    sampleSheet = file(params.samplesheet, checkIfExists: true)
    trgtConfig = file(params.trgt_config, checkIfExists: true)
    longtrConfig = file(params.longtr_config, checkIfExists: true)
    atarvaConfig = file(params.atarva_config, checkIfExists: true)
    strdustConfig = file(params.strdust_config, checkIfExists: true)
    consensusConfig = file(params.consensus_config, checkIfExists: true)
    trgtCatalog = file(resolveConfigPath(params.trgt_config, trgtCfg.native_catalog as String), checkIfExists: true)
    longtrCatalog = file(resolveConfigPath(params.longtr_config, longtrCfg.native_catalog as String), checkIfExists: true)
    atarvaCatalog = file(resolveConfigPath(params.atarva_config, atarvaCfg.native_catalog as String), checkIfExists: true)
    atarvaCatalogIndex = file(resolveConfigPath(params.atarva_config,
        (atarvaCfg.native_catalog_index ?: "${atarvaCfg.native_catalog}.tbi") as String), checkIfExists: true)
    strdustCatalog = file(resolveConfigPath(params.strdust_config, strdustCfg.native_catalog as String), checkIfExists: true)

    validationInput = Channel.of(tuple(
        sampleSheet, reference, referenceFai, manifest,
        trgtConfig, trgtCatalog, longtrConfig, longtrCatalog,
        atarvaConfig, atarvaCatalog, atarvaCatalogIndex,
        strdustConfig, strdustCatalog, consensusConfig, callers.join(',')
    ))
    alignmentAssets = Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> [file(row.alignment.toString(), checkIfExists: true), file(row.alignment_index.toString(), checkIfExists: true)] }
        .flatten()
        .collect()
    VALIDATE_INPUTS(validationInput, alignmentAssets)
    validationGate = VALIDATE_INPUTS.out.validated

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

    TRGT(samples.filter { it[1] == 'hifi' && callers.contains('trgt') }.combine(validationGate)
        .map { s, p, a, i, sex, ploidy, gate -> tuple(s, p, a, i, sex, ploidy, reference, referenceFai, trgtCatalog, trgtConfig, gate) })
    LONGTR(samples.filter { it[1] in ['hifi', 'ont'] && callers.contains('longtr') }.combine(validationGate)
        .map { s, p, a, i, sex, ploidy, gate -> tuple(s, p, a, i, sex, ploidy, reference, referenceFai, longtrCatalog, longtrConfig, gate) })
    ATARVA(samples.filter { it[1] in ['hifi', 'ont'] && callers.contains('atarva') }.combine(validationGate)
        .map { s, p, a, i, sex, ploidy, gate -> tuple(s, p, a, i, sex, ploidy, reference, referenceFai, atarvaCatalog, atarvaCatalogIndex, atarvaConfig, gate) })
    STRDUST(samples.filter { it[1] in ['hifi', 'ont'] && callers.contains('strdust') }.combine(validationGate)
        .map { s, p, a, i, sex, ploidy, gate -> tuple(s, p, a, i, sex, ploidy, reference, referenceFai, strdustCatalog, strdustConfig, gate) })

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
