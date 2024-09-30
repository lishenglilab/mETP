import subprocess
import os
import argparse

def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Flair correct function.")
    parser.add_argument("-q", "--query", help="Your uncorrected bed12 file.")
    parser.add_argument("-a", "--annotation", help="GTF annotation file")
    parser.add_argument("-g", "--genome", help="Fasta of reference genome")
    parser.add_argument("-o", "--outdir", help="Your output path.")
    parser.add_argument("-n", "--outname",default='flair', help="Your output file's name.")
    parser.add_argument("-t", "--thread",type=int,default=4,help="number of thread (default:4)")
    options = parser.parse_args(args)
    return options


options = getOptions(sys.argv[1:])


def flair_correct(query_bed = options.query,ref_gtf = options.annotation,ref_fa = options.genome,out_dir = options.outdir,out_name = options.outname,thread = options.thread):
    out_put = os.path.join(out_dir,out_name)
    correct_cmd = f"flair correct -q {query_bed} -f {ref_gtf} -g {ref_fa} --output {out_put} --threads {thread}"
    subprocess.run(correct_cmd,shell=True,check=True)


if __name__ == '__main__':
  align()

