process NORMALIZE_LONGTR {
    tag "${sample_id}:normalize:longtr"
    label 'normalize'
    container(params.normalizer_container)
    publishDir "${params.outdir}/normalized/longtr", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), val(platform), path(native_dir), path(manifest)

    output:
    tuple val(sample_id), val(platform), path("${sample_id}.longtr.normalized.tsv.gz"), emit: calls
    path "${sample_id}.longtr.normalize.versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/normalize_calls.py --caller longtr --sample-id '${sample_id}' --platform '${platform}' \\
      --vcf '${native_dir}/calls.vcf.gz' --locus-manifest '${manifest}' \\
      --source-vcf 'raw/longtr/${sample_id}/calls.vcf.gz' --output '${sample_id}.longtr.normalized.tsv.gz' \\
      --versions '${sample_id}.longtr.normalize.versions.yml'
    """

    stub:
    """
    printf 'sample_id\\tplatform\\tcanonical_locus_id\\tchr\\tcanonical_start\\tcanonical_end\\treference_build\\tmotif\\tcaller\\tcall_status\\tallele1_length_bp\\tallele2_length_bp\\tallele1_copy_number\\tallele2_copy_number\\tallele1_sequence\\tallele2_sequence\\tgenotype_quality\\tread_support_allele1\\tread_support_allele2\\tphase_set\\thaplotype_labels\\tfilter_reasons\\tnative_record_id\\tsource_vcf\\tparser_schema_version\\n' | cat > '${sample_id}.longtr.normalized.tsv.gz'
    printf 'caller: longtr\\nnormalizer: stub\\n' > '${sample_id}.longtr.normalize.versions.yml'
    """
}
