process CONSENSUS {
    tag 'cohort-consensus'
    label 'consensus'
    container(params.normalizer_container)
    publishDir "${params.outdir}/consensus", mode: 'copy', overwrite: true

    input:
    path normalized_files
    tuple path(manifest), path(samplesheet), path(config), val(callers)

    output:
    path 'consensus', emit: results
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${projectDir}/bin/build_consensus.py --normalized-files ${normalized_files.collect { "'${it}'" }.join(' ')} \\
      --locus-manifest '${manifest}' --samplesheet '${samplesheet}' --config '${config}' \\
      --callers '${callers}' --output-dir consensus
    bgzip -f consensus/consensus.vcf
    tabix -f -p vcf consensus/consensus.vcf.gz
    cp consensus/versions.yml versions.yml
    """

    stub:
    """
    mkdir -p consensus
    printf '##fileformat=VCFv4.3\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' > consensus/consensus.vcf
    cat consensus/consensus.vcf > consensus/consensus.vcf.gz
    touch consensus/consensus.vcf.gz.tbi
    touch consensus/consensus.tsv.gz consensus/caller_matrix.tsv.gz consensus/discordant_calls.tsv.gz consensus/locus_summary.tsv.gz
    printf '{"stub": true}\\n' > consensus/provenance.json
    printf '# Stub QC report\\n' > consensus/qc_report.md
    printf 'consensus: stub\\n' > consensus/versions.yml
    cp consensus/versions.yml versions.yml
    """
}
