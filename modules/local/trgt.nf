process TRGT {
    tag "${sample_id}:trgt"
    label 'trgt'
    container(params.trgt_container)
    publishDir "${params.outdir}/raw/trgt", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), val(platform), path(alignment), path(alignment_index), val(sex), val(ploidy),
          path(reference), path(reference_fai), path(catalog), path(config), path(validation_gate)

    output:
    tuple val(sample_id), val(platform), path("${sample_id}.trgt.native"), emit: calls

    script:
    """
    python3 ${projectDir}/bin/run_caller.py --caller trgt --sample-id '${sample_id}' --platform '${platform}' \\
      --alignment '${alignment}' --alignment-index '${alignment_index}' --reference '${reference}' \\
      --catalog '${catalog}' --config '${config}' --sex '${sex}' --ploidy '${ploidy}' \\
      --output-dir '${sample_id}.trgt.native'
    """

    stub:
    """
    mkdir -p '${sample_id}.trgt.native'
    printf '##fileformat=VCFv4.3\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${sample_id}\\n' > '${sample_id}.trgt.native/calls.vcf'
    cat '${sample_id}.trgt.native/calls.vcf' > '${sample_id}.trgt.native/calls.vcf.gz'
    touch '${sample_id}.trgt.native/calls.vcf.gz.tbi'
    printf 'caller: trgt\\nversion: stub\\n' > '${sample_id}.trgt.native/versions.yml'
    """
}
