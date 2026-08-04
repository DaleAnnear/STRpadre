process NORMALIZE_TRGT {
    tag "${sample_id}:normalize:trgt"
    label 'normalize'
    container(params.normalizer_container)
    publishDir "${params.outdir}/normalized/trgt", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), val(platform), path(native_dir), path(manifest)

    output:
    tuple val(sample_id), val(platform), path("${sample_id}.trgt.normalized.tsv.gz"), emit: calls
    path "${sample_id}.trgt.normalize.versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/normalize_calls.py --caller trgt --sample-id '${sample_id}' --platform '${platform}' \\
      --vcf '${native_dir}/calls.vcf.gz' --locus-manifest '${manifest}' \\
      --source-vcf 'raw/trgt/${sample_id}/calls.vcf.gz' --output '${sample_id}.trgt.normalized.tsv.gz' \\
      --versions '${sample_id}.trgt.normalize.versions.yml'
    """

    stub:
    """
    printf 'sample_id\\tplatform\\n' | cat > '${sample_id}.trgt.normalized.tsv.gz'
    printf 'caller: trgt\\nnormalizer: stub\\n' > '${sample_id}.trgt.normalize.versions.yml'
    """
}
