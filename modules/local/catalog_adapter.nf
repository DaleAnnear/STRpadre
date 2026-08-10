process CATALOG_ADAPTER {
    tag "canonical-catalog:${callers}"
    label 'catalog_adapter'
    container(params.normalizer_container)
    publishDir "${params.outdir}/catalogs", mode: 'copy', overwrite: true

    input:
    tuple path(manifest), path(reference_fai), val(callers)

    output:
    path 'catalogs', emit: catalogs

    script:
    """
    python3 ${projectDir}/bin/catalog_adapter.py \\
      --locus-manifest '${manifest}' \\
      --reference-fai '${reference_fai}' \\
      --callers '${callers}' \\
      --output-dir catalogs
    """

    stub:
    """
    mkdir -p catalogs
    touch catalogs/trgt.bed catalogs/longtr.bed catalogs/atarva.bed.gz catalogs/atarva.bed.gz.tbi catalogs/strdust.bed
    printf 'canonical_locus_id\\tcaller\\tnative_locus_id\\tchr\\tstart\\tend\\tmotif\\tcatalog_file\\tcatalog_coordinate_system\\n' > catalogs/catalog_locus_mapping.tsv
    printf '{"schema_version":"1.0","stub":true}\\n' > catalogs/catalog_adaptation_report.json
    """
}
