import os
import subprocess
import sys

def align(in_fq,out_dir,ref_fa,out_name = 'flair.align',thread = 4):
    out_put = os.path.join(out_dir,out_name)
    align_cmd = f"flair align -g {ref_fa} -r {in_fq} -o {out_put} -t {thread}"
    subprocess.run(align_cmd,shell=True,check=True)


def flair_correct(query_bed,ref_gtf,ref_fa,out_dir,out_name = 'flair',thread = 4):
    out_put = os.path.join(out_dir,out_name)
    correct_cmd = f"flair correct -q {query_bed} -f {ref_gtf} -g {ref_fa} --output {out_put} --threads {thread}"
    subprocess.run(correct_cmd,shell=True,check=True)


def flair_collapse(query_bed,ref_gtf,ref_fa,out_put,fq,thread = 4):
    collapse_cmd=f"flair collapse -g {ref_fa} -q {query_bed} --gtf {ref_gtf} -r {fq} --output {out_put} --threads {thread} --temp_dir {out_dir}"
    subprocess.run(collapse_cmd,shell=True,check=True)


def flair_quant(manifest,isoforms,out_dir,out_name = 'flair_quantify',thread = 4):
    out_put = os.path.join(out_dir,out_name)
    quant_cmd=f"flair quantify -r {manifest} -i {isoforms} --output {out_put} --threads {thread} --temp_dir {out_dir} --sample_id_only --tpm"
    subprocess.run(quant_cmd,shell=True,check=True)


def transcript_filter(quant,isoforms,out_dir,tpm = 0.1):
    scriptPath=sys.path[0]+'/flair/transcriptFilter.R'
    filter_cmd='Rscript '+scriptPath+' -t '+quant+' -a '+isoforms+' -q '+str(tpm)+' -o '+out_dir
    subprocess.run(filter_cmd,shell=True,check=True)
