import os
import subprocess
import sys

def generateReference(gtf,fasta,protein=False,thread=1):
  scriptPath=sys.path[0]+'/referenceGenerate/'
  gtf=os.path.abspath(gtf)
  fasta=os.path.abspath(fasta)
  if protein:
    process_cmd='Rscript '+scriptPath+'0.gtf_process.R -g '+gtf
    subprocess.run(process_cmd,shell=True,check=True)
    gtf=gtf.replace('.gtf','_proteinCoding.gtf')
  outfa=gtf.replace('.gtf','.fa')
  process_cmd='gffread -w '+outfa+' -g '+fasta+' '+gtf
  subprocess.run(process_cmd,shell=True,check=True)
  outpath='/'
  outpath=outpath.join(outfa.split('/')[:-1])
  if thread == 1:
    process_cmd='python '+scriptPath+'1.ORF_generate.py -i '+outfa+' -o '+outpath
    subprocess.run(process_cmd,shell=True,check=True)
  if thread > 1:
    if not os.path.exists(outpath+'/split/trancript'):
      os.makedirs(outpath+'/split/trancript')
    process_cmd='Rscript '+scriptPath+'2.split_fasta.R -i '+outfa+' -s '+str(thread)+' -o '+outpath+'/split/trancript'
    subprocess.run(process_cmd,shell=True,check=True)
    if not os.path.exists(outpath+'/split/ORF'):
      os.makedirs(outpath+'/split/ORF')
    process_cmd='sed "s#1.ORF_generate.py#'+scriptPath+'1.ORF_generate.py#g" '+scriptPath+'ORF_generate.sh > '+scriptPath+'tmp_ORF_generate.sh'
    subprocess.run(process_cmd,shell=True,check=True)
    fa_split=outpath+'/split/trancript/'+outfa.split('/')[-1].replace('.fa','')
    process_cmd='sed -i "s#transcript#'+fa_split+'#g" '+scriptPath+'tmp_ORF_generate.sh'
    subprocess.run(process_cmd,shell=True,check=True)
    process_cmd='sed -i "s#path#'+outpath+'/split/ORF/'+'#g" '+scriptPath+'tmp_ORF_generate.sh'
    subprocess.run(process_cmd,shell=True,check=True)
    process_cmd='sed -i "s#number#'+str(thread)+'#g" '+scriptPath+'tmp_ORF_generate.sh'
    subprocess.run(process_cmd,shell=True,check=True)
    process_cmd='bash '+scriptPath+'tmp_ORF_generate.sh'
    subprocess.run(process_cmd,shell=True,check=True)
    process_cmd='cat '+outpath+'/split/ORF/* > '+outpath+'/ORF.fa'
    subprocess.run(process_cmd,shell=True,check=True)
  process_cmd='python '+scriptPath+'4.remove_ORF_redanduncy.py -i '+outpath+'/ORF.fa'+' -o '+outpath
  subprocess.run(process_cmd,shell=True,check=True)
  subprocess.run('rm -rf '+outpath+'/split',shell=True,check=True)


def removeHomo(ORF,db,refFa,thread=1):
  scriptPath=sys.path[0]+'/referenceGenerate/'
  ORF=os.path.abspath(ORF)
  orfPath='/'
  orfPath=orfPath.join(ORF.split('/')[:-1])
  if not db:
    if not os.path.exists(orfPath+'/db'):
      os.makedirs(orfPath+'/db')
    process_cmd='makeblastdb -in '+refFa+' -dbtype prot -parse_seqids -out '+orfPath+'/db/ref.db'
    subprocess.run(process_cmd,shell=True,check=True)
    db=orfPath+'/db/ref.db'
  if thread==1:
    process_cmd='blastp -query '+ORF+' -out '+ORF.replace('.fa','.blast')+' -db '+db+' -outfmt 6 -evalue 1e-5 -num_threads 1'
    subprocess.run(process_cmd,shell=True,check=True)
  if thread>1:
    if not os.path.exists(orfPath+'/split/ORF'):
      os.makedirs(orfPath+'/split/ORF')
    process_cmd='Rscript '+scriptPath+'2.split_fasta.R -i '+ORF+' -s '+str(thread)+' -o '+orfPath+'/split/ORF'
    subprocess.run(process_cmd,shell=True,check=True)
    if not os.path.exists(orfPath+'/split/blast'):
      os.makedirs(orfPath+'/split/blast')
    process_cmd='sed "s#number#'+str(thread)+'#g" '+scriptPath+'removeORFhomo.sh > '+scriptPath+'tmp_removeORFhomo.sh'
    subprocess.run(process_cmd,shell=True,check=True)
    fa_split=orfPath+'/split/ORF/'+ORF.split('/')[-1].replace('.fa','')
    process_cmd='sed -i "s#ORF#'+fa_split+'#g" '+scriptPath+'tmp_removeORFhomo.sh'
    subprocess.run(process_cmd,shell=True,check=True)
    process_cmd='sed -i "s#path#'+orfPath+'/split/blast/'+ORF.split('/')[-1].replace('.fa','')+'#g" '+scriptPath+'tmp_removeORFhomo.sh'
    subprocess.run(process_cmd,shell=True,check=True)
    process_cmd='sed -i "s#ref.db#'+db+'#g" '+scriptPath+'tmp_removeORFhomo.sh'
    subprocess.run(process_cmd,shell=True,check=True)
    process_cmd='bash '+scriptPath+'tmp_removeORFhomo.sh'
    subprocess.run(process_cmd,shell=True,check=True)
    process_cmd='cat '+orfPath+'/split/blast/* > '+ORF.replace('.fa','.blast')
    subprocess.run(process_cmd,shell=True,check=True)
  process_cmd='Rscript '+scriptPath+'3.remove_homo_ORF.R -f '+ORF+' -b '+ORF.replace('.fa','.blast')+' -n ORF_unmatch_ref.fa'
  subprocess.run(process_cmd,shell=True,check=True)
  subprocess.run('rm -rf '+orfPath+'/split',shell=True,check=True)
