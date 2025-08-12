process COMBINE_QC_STATS {
    publishDir "${params.outdir}/table", mode: 'copy'
    
    input:
    val(all_stats)
    path(fail_reads_file)
    
    output:
    path("QC.xls")
    
    script:
    """
    combine_qc_stats.R
    """
    
    stub:
    """
    echo -e "sample_id\\tSize_GB\\tTotal_reads\\tpass_percent\\tMedian_pass_read_length\\tMean_pass_read_length\\tMapped_read\\tMapped_ratio" > QC.xls
    echo -e "test_sample\\t10.5\\t1000000\\t90.5\\t800\\t950\\t850000\\t85.0" >> QC.xls
    """
}