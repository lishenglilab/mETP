import argparse
import os
import sys
import datetime
import pandas as pd
from flair_func import align,flair_correct,flair_collapse,flair_quant,transcript_filter
from preprocess import readsProcess
from dataCleaning import dataCleaning
from expAnalysis import expAnalysis
from pepFeature import pepFeature
from referenceGenerate import generateReference
from referenceGenerate import removeHomo
from fragpipe import fragpipe_pip,addDecoy

def main():
    parser = argparse.ArgumentParser(description="Long Reads RNA Sequencing and Protein Mass Spectrum Analysis Tool with Functions readsChop, flair_pip, generateReference, removeHomo, dataCleaning, pepFeature, addDecoy, fragpipe_pip, expAnalysis",
                                     usage='proteinAnalysis.py <function> [options]')
    subparsers = parser.add_subparsers(title='Functions', dest='function')

    # Subparser for Function readsChop
    parser_readsChop = subparsers.add_parser('readsChop', help='Remove adapter from long reads fastq files using porechop.')
    parser_readsChop.add_argument("-i","--in_fq", required=True,help="Path of fastq file.")
    parser_readsChop.add_argument("-o","--out_dir", required=True,help="Output path of processed fastq file.")
    parser_readsChop.add_argument("-n","--outname",default='clean.fastq.gz',type=str,help="Output name of processed fastq file (default: clean.fastq.gz).")
    parser_readsChop.add_argument("-t","--thread", type=int,default=4,help="number of thread (default: 4).")

    # Subparser for Function flair_pip
    parser_flair_pip = subparsers.add_parser('flair_pip', help='Long reads RNA sequencing assembly and quantification pipline using flair.')
    parser_flair_pip.add_argument("-m","--manifest", required=True,help="A tab delimiter file containing: sample id at first column, condition information at second column, batch information at third column, fastq file path at forth column.")
    parser_flair_pip.add_argument("-o","--out_dir", required=True,help="Output path of processed file.")
    parser_flair_pip.add_argument("-g","--annotation",required=True,help="Path of genome annotation GTF file.")
    parser_flair_pip.add_argument("-f","--fasta",required=True,help="Path of genome FASTA file.")
    parser_flair_pip.add_argument("-q","--tpm", type=int,default=0.1,help="TPM quantitative value cutoff (default: 0.1).")
    parser_flair_pip.add_argument("-t","--thread", type=int,default=4,help="Number of threads (default: 4)")

    # Subparser for Function generateReference
    parser_generateReference = subparsers.add_parser('generateReference', help='Calling ORF from transcript fasta file.')
    parser_generateReference.add_argument("-g","--gtf", required=True,help="path of gtf file")
    parser_generateReference.add_argument("-f","--fasta", required=True,help="path of genome fasta file")
    parser_generateReference.add_argument("-p","--protein",default=False,type=bool,help="whether remove protein coding transcript (default: False)")
    parser_generateReference.add_argument("-t","--thread", type=int,default=1,help="number of thread (default: 1)")

    # Subparser for Function removeHomo
    parser_removeHomo = subparsers.add_parser('removeHomo', help='Remove proteins that are homologous to the specified reference file')
    parser_removeHomo.add_argument("-r","--ORF",required=True,help="path of ORF fasta file")
    parser_removeHomo.add_argument("-b","--db",help="path of database to blast")
    parser_removeHomo.add_argument("-f","--ref",help="if there is not a database to blast, supply the path of reference protein fasta file")
    parser_removeHomo.add_argument("-t","--thread",default=10,type=int,help="number of threads (default:10)")

    # Subparser for Function dataCleaning
    parser_dataCleaning = subparsers.add_parser('dataCleaning', help='This function can only process the quants file of DIA-NN software, it will output the presurcor, peptide and protein expression matrix file.')
    parser_dataCleaning.add_argument("-q","--quantFile",required=True,help="path of DIA-NN report.tsv file")

    # Subparser for Function pepFeature
    parser_pepFeature = subparsers.add_parser('pepFeature', help='quantity transcript for each cell')
    parser_pepFeature.add_argument("-q","--quantFile",required=True,help="tab delimited file of protein expression matrix")
    parser_pepFeature.add_argument("-r","--pepRef",required=True,help="path of protein reference fasta file")
    parser_pepFeature.add_argument("-t","--transcriptRef",required=True,help="path of transcript reference fasta file which is used for ORF calling")
    parser_pepFeature.add_argument("-o","--outpath",help="directory of this step output")

    # Subparser for Function expAnalysis
    parser_expAnalysis = subparsers.add_parser('expAnalysis', help='allele specific transcript analysis')
    parser_expAnalysis.add_argument("-p","--pepQuant",required=True,help="tab delimited file of protein expression matrix")
    parser_expAnalysis.add_argument("-t","--transcriptQuant",required=True,help="tab delimited file of transcript expression matrix")
    parser_expAnalysis.add_argument("-c","--cor", default='pearson',type=str,help="method used for correlation analysis, default is 'pearson'")
    parser_expAnalysis.add_argument("-m","--manifest",help="tab delimited file containing the mapping between protein ids and transcript ids, optional only for protein reference generated by 'generateReference'")
    parser_expAnalysis.add_argument("-o","--outpath",help="directory of this step output")

    # Subparser for add decoy to protein reference
    parser_addDecoy = subparsers.add_parser('addDecoy', help='Add decoy to protein reference')
    parser_addDecoy.add_argument("-p","--philosopher",required=True,help="Path to binary file of philosopher.")
    parser_addDecoy.add_argument("-f","--ref_fa",required=True,help="Path to protein reference fasta file.")

    # Subparser for protein quantification of mass spectrum
    parser_fragpipe = subparsers.add_parser('fragpipe_pip', help='Protein quantification of mass spectrum using fragpipe')
    parser_fragpipe.add_argument("-f","--fragpipe",required=True,help="Path to binary file of fragpipe.")
    parser_fragpipe.add_argument("-m","--msfragger",required=True,help="Path to JAR file of msfragger.")
    parser_fragpipe.add_argument("-q","--ionquant",required=True,help="Path to JAR file of ionquant.")
    parser_fragpipe.add_argument("-p","--philosopher",required=True,help="Path to binary file of philosopher.")
    parser_fragpipe.add_argument("-w","--workflow",required=True,help="Work flow file which can be generated by fragpipe.")
    parser_fragpipe.add_argument("-i","--manifest",required=True,help="A tab delimiter file containing: mass spectrum file path at first column, sample id at second column, batch information at third column, 'DDA' or 'DIA' at forth column.")
    parser_fragpipe.add_argument("-o","--workdir",required=True,help="Output directory")
    parser_fragpipe.add_argument("-t","--thread",default=10,type=int,help="number of threads (default:10)")

    args = parser.parse_args()

    if args.function == 'readsChop':
        readsProcess(in_fq=args.in_fq,out_dir=args.out_dir,outname=args.outname,thread=args.thread)
    elif args.function == 'flair_pip':
        manifest=os.path.abspath(args.manifest)
        out_dir=os.path.abspath(args.out_dir)
        fasta=os.path.abspath(args.fasta)
        gtf=os.path.abspath(args.annotation)
        info=pd.read_table(manifest,header=None)
        fq=info.iloc[:,3].to_list()
        samples=info.iloc[:,0].to_list()
        align_path=os.path.join(out_dir,'align')
        if not os.path.exists(align_path):
            os.makedirs(align_path)
        for i in range(0,len(samples)):
            align_path_sub=os.path.join(align_path,samples[i])
            if not os.path.exists(align_path_sub):
                os.makedirs(align_path_sub)
            align(in_fq=fq[i],out_dir=align_path_sub,ref_fa=fasta,out_name = samples[i],thread = args.thread)
        merge_cmd='cat '+align_path+'/*/*.bed > '+align_path+'/merge.bed'
        subprocess.run(merge_cmd,shell=True,check=True)
        query_bed=os.path.join(align_path,'merge.bed')
        correct_path=os.path.join(out_dir,'correct')
        if not os.path.exists(correct_path):
            os.makedirs(correct_path)
        flair_correct(query_bed=query_bed,ref_gtf=gtf,ref_fa=fasta,out_dir=correct_path,out_name = 'flair',thread = args.thread)
        query_bed_corrected=os.path.join(correct_path,'flair_all_corrected.bed')
        collapse_path=os.path.join(out_dir,'collapse/')
        if not os.path.exists(collapse_path):
            os.makedirs(collapse_path)
        isoforms_path=os.path.join(collapse_path,'isoforms')
        flair_collapse(query_bed=query_bed_corrected,ref_gtf=gtf,ref_fa=fasta,out_put=isoforms_path,fq=','.join(fq),thread = args.thread)
        isoforms=os.path.join(collapse_path,'isoforms.gtf')
        quant_path=os.path.join(out_dir,'quant/')
        if not os.path.exists(quant_path):
            os.makedirs(quant_path)
        flair_quant(manifest=manifest,isoforms=isoforms,out_dir=quant_path,out_name='flair_quantify',thread=args.thread)
        quant_file=os.path.join(quant_path,'flair_quantify.tpm.tsv')
        transcript_filter(quant=quant_file,isoforms=isoforms,out_dir=out_dir,tpm = args.tpm)
    elif args.function == 'generateReference':
        generateReference(gtf=args.gtf,fasta=args.fasta,protein=args.protein,thread=args.thread)
    elif args.function == 'removeHomo':
        removeHomo(ORF=args.ORF,db=args.db,refFa=args.ref,thread=args.thread)
    elif args.function == 'dataCleaning':
        dataCleaning(quantFile=args.quantFile)
    elif args.function == 'pepFeature':
        pepFeature(quantFile=args.quantFile,pepRef=args.pepRef,transcriptRef=args.transcriptRef,outpath=args.outpath)
    elif args.function == 'expAnalysis':
        expAnalysis(pepQuant=args.pepQuant,transcriptQuant=args.transcriptQuant,cor=args.cor,manifest=args.manifest,outpath=args.outpath)
    elif args.function == 'addDecoy':
        addDecoy(ref_fa=args.ref_fa,philosopher=args.philosopher)
    elif args.function == 'fragpipe_pip':
        fragpipe_pip(fragpipe=args.fragpipe,workflow=args.workflow,msfragger=args.msfragger,ionquant=args.ionquant,philosopher=args.philosopher,manifest=args.manifest,workdir=args.workdir,thread=args.thread)
    else:
        print("Unknown function. Use -h for help.")


if __name__ == "__main__":
    main()
