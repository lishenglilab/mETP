import subprocess
import os
import argparse

def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Align Oxford Nanopore reads to genome fasta.")
    parser.add_argument("-i", "--input", help="Your input fastq file.")
    parser.add_argument("-o", "--outdir", help="Your output path.")
    parser.add_argument("-n", "--outname",default='align2genome', help="Your output file's name.")
    parser.add_argument("-g", "--genome", help="Fasta of reference genome")
    parser.add_argument("-t", "--thread",type=int,default=1,help="number of thread (default:1)")
    options = parser.parse_args(args)
    return options


options = getOptions(sys.argv[1:])


def align(in_fq = options.input,out_dir = options.outdir,out_name = options.outname,ref_fa = options.genome,thread = options.thread):
    out_put = os.path.join(out_dir,out_name)
    align_cmd = f"flair align -g {ref_fa} -r {in_fq} -o {out_put} -t {thread}"
    subprocess.run(align_cmd,shell=True,check=True)


if __name__ == '__main__':
  align()

