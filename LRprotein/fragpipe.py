import os
import subprocess
import sys

def fragpipe_pip(fragpipe,msfragger,ionquant,philosopher,workflow,manifest,workdir,thread):
  fragpipe=os.path.abspath(fragpipe)
  workflow=os.path.abspath(fragpipe)
  msfragger=os.path.abspath(msfragger)
  ionquant=os.path.abspath(ionquant)
  philosopher=os.path.abspath(philosopher)
  manifest=os.path.abspath(manifest)
  workdir=os.path.abspath(workdir)
  if not os.path.exists(workdir):
    os.makedirs(workdir)
  fragpipe_cmd=fragpipe+' --headless --workflow '+workflow+' --manifest '+manifest+' --workdir '+workdir+' --config-msfragger '+msfragger+' --config-ionquant '+ionquant+' --config-philosopher '+philosopher+' --threads '+thread
  subprocess.run(fragpipe_cmd,shell=True,check=True)


def addDecoy(ref_fa,philosopher):
  philosopher=os.path.abspath(philosopher)
  ref_fa=os.path.abspath(ref_fa)
  work_dir=os.path.dirname()
  work_init_cmd=philosopher+' workspace --init --nocheck --temp '+work_dir
  addDecoy_cmd=philosopher+' database --custom '+ref_fa
  clean_cmd=philosopher+' workspace --clean --nocheck'
  subprocess.run(work_init_cmd,shell=True,check=True)
  subprocess.run(addDecoy_cmd,shell=True,check=True)
  subprocess.run(clean_cmd,shell=True,check=True)


