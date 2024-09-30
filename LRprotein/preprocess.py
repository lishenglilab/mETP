import os
import subprocess
import sys


def readsProcess(in_fq,out_dir,outname = 'clean.fastq.gz',thread = 4):
    out_fq = os.path.join(out_dir,outname)
    if not os.path.exists(out_fq):
        run_cmd = f"porechop -i {in_fq} -o {out_fq} -t {thread}"
        subprocess.run(run_cmd,shell=True,check=True)
    else:
        print('Clean data is existed, run alignment step.')

