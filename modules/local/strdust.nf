process STRDUST {
    tag "${sample_id}:strdust"
    label 'strdust'
    container(params.strdust_container)
    publishDir "${params.outdir}/raw/strdust", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), val(platform), path(alignment), path(alignment_index), val(sex), val(ploidy),
          path(reference), path(reference_fai), path(catalog), path(config), path(validation_gate)

    output:
    tuple val(sample_id), val(platform), path("${sample_id}.strdust.native"), emit: calls

    script:
    """
    python3 ${projectDir}/bin/run_caller.py --caller strdust --sample-id '${sample_id}' --platform '${platform}' \\
      --alignment '${alignment}' --alignment-index '${alignment_index}' --reference '${reference}' \\
      --catalog '${catalog}' --config '${config}' --sex '${sex}' --ploidy '${ploidy}' \\
      --output-dir '${sample_id}.strdust.native'
    """

    stub:
    """
    mkdir -p '${sample_id}.strdust.native'
    printf '##fileformat=VCFv4.3\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${sample_id}\\n' > '${sample_id}.strdust.native/calls.vcf'
    cat '${sample_id}.strdust.native/calls.vcf' > '${sample_id}.strdust.native/calls.vcf.gz'
    touch '${sample_id}.strdust.native/calls.vcf.gz.tbi'
    printf 'caller: strdust\\nversion: stub\\n' > '${sample_id}.strdust.native/versions.yml'
    """
}
