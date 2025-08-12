process GET_FILE_SIZE {
    tag "$sample"
    
    input:
    tuple val(sample), path(fastq)
    
    output:
    tuple val(sample), val(size_gb), emit: sizes
    
    script:
    """
    size=\$(stat -c%s ${fastq})
    size_gb=\$(echo "scale=2; \$size / 1024 / 1024 / 1024" | bc)
    echo "${sample} \$size_gb" > ${sample}_size.txt
    """
    
    stub:
    """
    echo "${sample} 10.5" > ${sample}_size.txt
    """
}