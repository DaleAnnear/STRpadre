nextflow.enable.dsl=2

include { CATALOG_ADAPTER } from './modules/local/catalog_adapter'

def normalizeCatalogCallers(value) {
    if (value == null) {
        error "--callers is required (for example: --callers trgt,longtr)"
    }
    def names = (value.toString().trim().toLowerCase() == 'all') ? ['trgt', 'longtr', 'atarva', 'strdust'] : value.toString().split(/,/).collect { it.trim().toLowerCase() }.findAll { it }
    if (!names) error '--callers must not be empty'
    def unknown = names.findAll { !(it ==~ /trgt|longtr|atarva|strdust/) }
    if (unknown) error(/Unknown caller name; allowed: trgt,longtr,atarva,strdust/)
    if (names.size() != names.toSet().size()) error "Duplicate caller name(s) are not allowed: ${names.join(',')}"
    return names
}

workflow {
    ['locus_manifest', 'reference_fai'].each { key ->
        if (!params[key]) error "--${key.replace('_', '-')} is required"
    }
    callers = normalizeCatalogCallers(params.callers)
    CATALOG_ADAPTER(Channel.of(tuple(
        file(params.locus_manifest, checkIfExists: true),
        file(params.reference_fai, checkIfExists: true),
        callers.join(',')
    )))
}
