#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// 定义项目目录
projectDir = file("$baseDir")

// Parameters
params.outdir = "/home/yvzeng/Project/mouseDRS/data/nextflow/testdata/result/"
params.qc_dir = "/home/yvzeng/Project/mouseDRS/data/nextflow/testdata/QC"
params.reference_fasta = "/home/public/reference/fa/mouse/GRCm39.primary_assembly.genomeV29.fa"
params.annotation_gtf = "/home/public/reference/gtf/mouse/gencode.vM29.primary_assembly.annotation.gtf"
params.ref_prot = "/home/yvzeng/Project/mouseDRS/data/refPep/uniprot_sprot_mouse.fasta"
params.manifest = "/home/yvzeng/Project/mouseDRS/data/nextflow/testdata/manifest.txt"
params.tpm_cutoff = 0.1
params.run_proteomics = false
params.proteomics_manifest = null
params.philosopher_path = "/path/to/MSPipe/philosopher/philosopher"
params.fragpipe_path = "/path/to/MSPipe/fragpipe/bin/fragpipe"
params.msfragger_jar = "/path/to/MSPipe/MSFragger-3.8/MSFragger-3.8.jar"
params.ionquant_jar = "/path/to/MSPipe/IonQuant-1.9.8/IonQuant-1.9.8.jar"

// Test module parameter
params.test_module = null

// Include processes
include { BAM_FLAGSTAT } from './modules/bam_stats.nf'
include { SEQKIT_LENGTH } from './modules/seqkit_stats.nf'
include { GET_FILE_SIZE } from './modules/file_size.nf'
include { COMBINE_QC_STATS } from './modules/combine_stats.nf'
include { READS_CHOP } from './modules/reads_chop.nf'
include { FLAIR_PIP } from './modules/flair_pip.nf'
include { GENERATE_REFERENCE } from './modules/generate_reference.nf'
include { REMOVE_HOMO } from './modules/remove_homo.nf'
include { MODIFICATION_DETECTION } from './modules/modification_detection.nf'
include { ADD_DECOY } from './modules/proteomics.nf'
include { FRAGPIPE } from './modules/proteomics.nf'
include { SUPPA2_ANALYSIS } from './modules/suppa2.nf'

workflow {
    // 从manifest文件中读取样本信息
    manifest_ch = Channel.fromPath(params.manifest)
    
    // 解析manifest文件，创建样本和FASTQ路径的通道
    samples_raw_fq_ch = Channel
      .fromPath(params.manifest)
      .splitCsv(sep: "\t", header: false)
      .map { row -> 
          def sample = row[0]
          def condition = row[1]
          def batch = row[2]
          def raw_fq = file(row[3].trim())
          def fast5_dir = file(row[4].trim())  // 新增FAST5目录列
          tuple(sample, condition, batch, raw_fq, fast5_dir)
      }
    
    // 预处理FASTQ文件
    READS_CHOP(samples_raw_fq_ch.map { sample, condition, batch, raw_fq -> tuple(sample, raw_fq) })
    clean_fq_ch = READS_CHOP.out.clean_fq
    
    // 创建包含清洁FASTQ路径的新manifest通道
    clean_manifest_ch = samples_raw_fq_ch
    .combine(clean_fq_ch)
    .map { sample, condition, batch, raw_fq, clean_fq -> 
        // 返回固定键 + 行内容
        tuple("manifest", "${sample}\t${condition}\t${batch}\t${clean_fq}")
    }
    .collectFile(
        name: "clean_manifest.txt",  // 输出文件名
        //seed: ["Sample\tCondition\tBatch\tClean_FASTQ"],
        storeDir: "${params.outdir}"
    ) { key, lines, result ->
        // 累积所有行
        result << lines
    }
    
    // 运行FLAIR管道（所有样本一起处理）
    FLAIR_PIP(clean_manifest_ch, file(params.annotation_gtf), file(params.reference_fasta), projectDir)
    filtered_gtf_ch = FLAIR_PIP.out.filtered_gtf
    quant_tsv_ch = FLAIR_PIP.out.quant_tsv
    bam_files_ch = FLAIR_PIP.out.bam_files  // 新增：获取BAM文件通道
    
    // 生成蛋白质参考序列
    GENERATE_REFERENCE(filtered_gtf_ch, file(params.reference_fasta), projectDir)
    orf_fa_ch = GENERATE_REFERENCE.out.orf_fa
    
    // 去除同源蛋白
    REMOVE_HOMO(orf_fa_ch, file(params.ref_prot), projectDir)
    filtered_orf_ch = REMOVE_HOMO.out.filtered_orf
    
    // 发布最终蛋白质参考序列
    filtered_orf_ch
        .map { orf -> 
            file(orf).copyTo("${params.outdir}/final_protein_reference.fa")
        }
    
    // 创建样本和BAM文件的通道
    bam_ch = bam_files_ch
        .flatMap { dir -> 
            fileTree(dir: dir, glob: "**/*.bam").collect { 
                def sample = it.getBaseName()
                tuple(sample, it)
            }
        }
    
    BAM_FLAGSTAT(bam_ch)
    SEQKIT_LENGTH(clean_fq_ch)
    GET_FILE_SIZE(clean_fq_ch)
    
    // 收集所有统计信息
    bam_stats = BAM_FLAGSTAT.out.stats
    length_stats = SEQKIT_LENGTH.out.lengths
    size_stats = GET_FILE_SIZE.out.sizes
    
    // 组合所有统计信息
    all_stats = bam_stats
        .join(length_stats)
        .join(size_stats)
    
    // 组合QC统计信息
    COMBINE_QC_STATS(
        all_stats.collect(),
        file("${params.qc_dir}/read_number_fail.txt")
    )
    
    // 获取FLAIR输出的参考转录本
    ref_tr_fa_ch = FLAIR_PIP.out.isoforms_fa.first()
    
    // 运行SUPPA2可变剪切分析
    SUPPA2_ANALYSIS(filtered_gtf_ch, quant_tsv_ch)
    
    // 准备修饰检测输入
    mod_input_ch = clean_fq_ch
      .join(samples_raw_fq_ch.map { s,c,b,fq,f5 -> tuple(s, f5) })
      .map { sample, clean_fq, sample2, fast5_dir -> 
          assert sample == sample2
          tuple(sample, fast5_dir, clean_fq) 
      }
      .combine(ref_tr_fa_ch)
      
    // 执行修饰检测
    MODIFICATION_DETECTION(mod_input_ch, projectDir)
    
    // 蛋白质谱分析模块
    if (params.run_proteomics) {
        // 检查manifest文件是否提供
        if (params.proteomics_manifest == null) {
            exit 1, "Proteomics manifest file not provided! Use --proteomics_manifest to specify."
        }
        
        // 添加诱饵处理
        ADD_DECOY(filtered_orf_ch.collect().first())
        decoy_db_ch = ADD_DECOY.out.decoy_db
        
        // 运行FragPipe
        FRAGPIPE(
            decoy_db_ch, 
            file(params.proteomics_manifest),
            file("$projectDir/bin/protein/workflow.workflow"),
            params.philosopher_path,
            params.fragpipe_path,
            params.msfragger_jar,
            params.ionquant_jar
        )
    }
}

// 测试入口点 - 移到主workflow外面
workflow test_module {
    // 1. 测试数据准备
    test_data_ch = Channel.fromPath("${projectDir}/test_data/*.fastq.gz")
        .map { file -> tuple(file.getBaseName(), file) }
    
    // 2. 模块选择器
    if (params.test_module == "reads_chop") {
        READS_CHOP(test_data_ch)
        READS_CHOP.out.clean_fq.view()
    }
    else if (params.test_module == "flair_pip") {
        // 需要额外测试数据
        test_gtf = file("${projectDir}/test_data/gencode.vM29.primary_assembly.annotation.gtf")
        test_fasta = file("${projectDir}/test_data/GRCm39.primary_assembly.genomeV29.fa")
        test_manifest = Channel.fromPath(file("${projectDir}/test_data/manifest.txt"))
        FLAIR_PIP(test_manifest, test_gtf, test_fasta, projectDir)
        FLAIR_PIP.out.isoforms_fa.view()
    }
    // 质量报告
    else if (params.test_module == "quality_report") {
        // 需要额外测试数据
        test_bam_ch = Channel.fromPath("${projectDir}/test_data/*.bam")
            .map { file -> tuple(file.getBaseName(), file) }
        BAM_FLAGSTAT(test_bam_ch.collect())
        SEQKIT_LENGTH(test_data_ch.collect())
        GET_FILE_SIZE(test_data_ch.collect())
        bam_stats = BAM_FLAGSTAT.out.stats
        length_stats = SEQKIT_LENGTH.out.lengths
        size_stats = GET_FILE_SIZE.out.sizes
        all_stats = bam_stats
            .join(length_stats)
            .join(size_stats)
        COMBINE_QC_STATS(
            all_stats.collect(),
            file("${projectDir}/test_data/read_number_fail.txt")
        )
    }
    // 蛋白质参考序列生成
    else if (params.test_module == "reference_generate") {
        test_gtf = file("${projectDir}/test_data/isoforms.gtf")
        test_fasta = file("${projectDir}/test_data/GRCm39.primary_assembly.genomeV29.fa")
        GENERATE_REFERENCE(test_gtf, test_fasta, file("${projectDir}"))
        orf_fa_ch = GENERATE_REFERENCE.out.orf_fa
        REMOVE_HOMO(orf_fa_ch, file("${projectDir}/test_data/uniprot_sprot_mouse.fasta"), file("${projectDir}"))
        filtered_orf_ch = REMOVE_HOMO.out.filtered_orf
        filtered_orf_ch
            .map { orf -> 
                file(orf).copyTo("${params.outdir}/final_protein_reference.fa")
            }
    }
    // 修饰检测
    else if (params.test_module == "modification") {  // 修正拼写错误
        test_mod_ch = Channel
          .fromPath("~/Project/mouseDRS/script/nextflow/v1/test_data/manifest.txt")
          .splitCsv(sep: "\t", header: false)
          .map { row -> 
              def sample = row[0]
              def condition = row[1]
              def batch = row[2]
              def raw_fq = file(row[3].trim())
              def fast5_dir = file(row[4].trim())  // 新增FAST5目录列
              tuple(sample, fast5_dir, raw_fq)
          }
        test_ref_tr_fa = file("${projectDir}/test_data/isoforms.fa")
        MODIFICATION_DETECTION(mod_input_ch, test_ref_tr_fa, projectDir)
    }
    else if (params.test_module == "proteomics") {
        test_ref = file("${projectDir}/test_data/test_protein.fa")
        test_manifest = file("${projectDir}/test_data/test_manifest.txt")
        ADD_DECOY(test_ref)
        FRAGPIPE(
            ADD_DECOY.out.decoy_db,
            test_manifest,
            file("${projectDir}/bin/protein/workflow.workflow"),
            params.philosopher_path,
            params.fragpipe_path,
            params.msfragger_jar,
            params.ionquant_jar
        )
    }
    else if (params.test_module == "suppa2") {
        test_gtf = file("${projectDir}/test_data/isoforms.gtf")
        test_tpm = file("${projectDir}/test_data/quant.tpm.tsv")
        SUPPA2_ANALYSIS(test_gtf, test_tpm)
        SUPPA2_ANALYSIS.out.ioe_files.view()
        SUPPA2_ANALYSIS.out.psi_files.view()
    }
    else {
        println "Unknown module: ${params.test_module}"
    }
}