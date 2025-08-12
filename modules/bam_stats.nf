process BAM_FLAGSTAT {
    tag "$sample"
    publishDir "${params.qc_dir}/bam_map", mode: 'copy'
    
    input:
    tuple val(sample), path(bam)
    
    output:
    tuple val(sample), path("${sample}.map.stat"), emit: flagstat
    tuple val(sample), val(mapping_reads), emit: stats
    
    script:
    """
    samtools flagstat ${bam} > ${sample}.map.stat
    mapping_reads=\$(sed -n '5p' ${sample}.map.stat | cut -d ' ' -f 1)
    echo "${sample} \$mapping_reads" > ${sample}_mapping.txt
    """
    
    stub:
    """
    touch ${sample}.map.stat
    echo "${sample} 1000000" > ${sample}_mapping.txt
    """
}