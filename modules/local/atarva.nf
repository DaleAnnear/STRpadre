process ATARVA {
    tag "${sample_id}:atarva"
    label 'atarva'
    container(params.atarva_container)
    publishDir "${params.outdir}/raw/atarva", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), val(platform), path(alignment), path(alignment_index), val(sex), val(ploidy),
          path(reference), path(reference_fai), path(catalog), path(catalog_index), path(config), path(validation_gate)

    output:
    tuple val(sample_id), val(platform), path("${sample_id}.atarva.native"), emit: calls

    script:
    """
    python3 ${projectDir}/bin/run_caller.py --caller atarva --sample-id '${sample_id}' --platform '${platform}' \\
      --alignment '${alignment}' --alignment-index '${alignment_index}' --reference '${reference}' \\
      --catalog '${catalog}' --config '${config}' --sex '${sex}' --ploidy '${ploidy}' \\
      --output-dir '${sample_id}.atarva.native'
    """

    stub:
    """
    mkdir -p '${sample_id}.atarva.native'
    printf '##fileformat=VCFv4.3\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${sample_id}\\n' > '${sample_id}.atarva.native/calls.vcf'
    cat '${sample_id}.atarva.native/calls.vcf' > '${sample_id}.atarva.native/calls.vcf.gz'
    touch '${sample_id}.atarva.native/calls.vcf.gz.tbi'
    printf 'caller: atarva\\nversion: stub\\n' > '${sample_id}.atarva.native/versions.yml'
    """
}
