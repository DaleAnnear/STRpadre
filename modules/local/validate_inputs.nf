process VALIDATE_INPUTS {
    tag 'preflight'
    label 'validation'
    container(params.normalizer_container)
    publishDir "${params.outdir}/validation", mode: 'copy', overwrite: true

    input:
    tuple path(samplesheet), path(reference), path(reference_fai), path(manifest),
          path(trgt_config, name: 'configs/trgt.yml'), path(trgt_catalog), path(longtr_config, name: 'configs/longtr.yml'), path(longtr_catalog),
          path(atarva_config, name: 'configs/atarva.yml'), path(atarva_catalog), path(atarva_catalog_index),
          path(strdust_config, name: 'configs/strdust.yml'), path(strdust_catalog), path(consensus_config), val(callers)
    path alignment_assets

    output:
    path 'validated', emit: validated
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${projectDir}/bin/validate_inputs.py \\
      --samplesheet '${samplesheet}' --reference '${reference}' --reference-fai '${reference_fai}' \\
      --locus-manifest '${manifest}' --callers '${callers}' \\
      --trgt-config '${trgt_config}' --trgt-catalog '${trgt_catalog}' \\
      --longtr-config '${longtr_config}' --longtr-catalog '${longtr_catalog}' \\
      --atarva-config '${atarva_config}' --atarva-catalog '${atarva_catalog}' --atarva-catalog-index '${atarva_catalog_index}' \\
      --strdust-config '${strdust_config}' --strdust-catalog '${strdust_catalog}' \\
      --consensus-config '${consensus_config}' --output-dir validated --staged-alignment-dir .
    cp validated/versions.yml versions.yml
    """

    stub:
    """
    mkdir -p validated
    printf 'schema_version: "1.0"\\nvalidation: stub\\n' > validated/versions.yml
    printf 'canonical_locus_id\\tcaller\\tnative_locus_id\\tchr\\tstart\\tend\\tmotif\\n' > validated/locus_mapping.tsv
    printf '{"stub": true}\\n' > validated/preflight_report.json
    cp validated/versions.yml versions.yml
    """
}
