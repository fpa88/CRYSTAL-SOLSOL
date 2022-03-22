#!/usr/bin/env python3
# File name : generate_CRYSTAL_KFeF3_SOLSOL_8M.py
#
# Author : Dr Fabien PASCALE
# contact : fabien.pascale@univ-lorraine.fr

from itertools import product
from collections import Counter
import time
import os
import re

# This program generates automatically all 560 configurations for each composition ( 2/3/3 , 3/2/3 and 3/3/2 ) 
# of KFeF3 systems with 8 magnetic sites

# RUN variable allows to submit via SLURM each job. 
# A check is done to see if the output exists and contains a "OPT END" line. In this case, the job is not submit. 
RUN = False
# NEED 3 templates files in the same directory with each population of d orbitals. See the last part with EIGSHIFT.

def reverse_readline(filename, buf_size=8192):
    """A generator that returns the lines of a file in reverse order"""
    with open(filename,"rb") as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        #while remaining_size > 0:
        count=0
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split(br'\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first 
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
            count=count+1
            if count >8 :
              break
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment


#p={'V': [3, 4], 'B': [5, 6, 7], 'R': [0, 1, 2]}

def genere_config(p,n,CASE,RUN,listefile):
  template=LIST[CASE]["template"]

  with open(template, 'r',encoding='iso-8859-1') as filecry :
    filedata = filecry.readlines()
  if "233" in CASE :
    filedata[6] =L[p['V'][0]] +'#'+str(p['V'][0]+1) +' pop0 V\n'
    filedata[7] =L[p['V'][1]] +'#'+str(p['V'][1]+1) +' pop0 V\n'
    filedata[8] =L[p['B'][0]] +'#'+str(p['B'][0]+1) +' pop1 B\n'
    filedata[9] =L[p['B'][1]] +'#'+str(p['B'][1]+1) +' pop1 B\n'
    filedata[10]=L[p['B'][2]] +'#'+str(p['B'][2]+1) +' pop1 B\n'
    filedata[11]=L[p['R'][0]] +'#'+str(p['R'][0]+1) +' pop2 R\n'
    filedata[12]=L[p['R'][1]] +'#'+str(p['R'][1]+1) +' pop2 R\n'
    filedata[13]=L[p['R'][2]] +'#'+str(p['R'][2]+1) +' pop2 R\n'

  if "323" in CASE :
    filedata[6] =L[p['B'][0]] +'#'+str(p['B'][0]+1) +' pop0 B\n'
    filedata[7] =L[p['B'][1]] +'#'+str(p['B'][1]+1) +' pop0 B\n'
    filedata[8] =L[p['B'][2]] +'#'+str(p['B'][2]+1) +' pop0 B\n'
    filedata[9] =L[p['V'][0]] +'#'+str(p['V'][0]+1) +' pop1 V\n'
    filedata[10]=L[p['V'][1]] +'#'+str(p['V'][1]+1) +' pop1 V\n'
    filedata[11]=L[p['R'][0]] +'#'+str(p['R'][0]+1) +' pop2 R\n'
    filedata[12]=L[p['R'][1]] +'#'+str(p['R'][1]+1) +' pop2 R\n'
    filedata[13]=L[p['R'][2]] +'#'+str(p['R'][2]+1) +' pop2 R\n'

  if "332" in CASE :
    filedata[6] =L[p['B'][0]] +'#'+str(p['B'][0]+1) +' pop0 B\n'
    filedata[7] =L[p['B'][1]] +'#'+str(p['B'][1]+1) +' pop0 B\n'
    filedata[8] =L[p['B'][2]] +'#'+str(p['B'][2]+1) +' pop0 B\n'
    filedata[9] =L[p['R'][0]] +'#'+str(p['R'][0]+1) +' pop1 R\n'
    filedata[10]=L[p['R'][1]] +'#'+str(p['R'][1]+1) +' pop1 R\n'
    filedata[11]=L[p['R'][2]] +'#'+str(p['R'][2]+1) +' pop1 R\n'
    filedata[12]=L[p['V'][0]] +'#'+str(p['V'][0]+1) +' pop2 V\n'
    filedata[13]=L[p['V'][1]] +'#'+str(p['V'][1]+1) +' pop2 V\n'

  
  NEWFILE=LIST[CASE]["config"]+'_%03d.d12' % n
  listefile.write('# '+NEWFILE +str(p)+'\n')

  with open(NEWFILE, 'w',encoding='utf-8') as filecry:
            filecry.writelines(filedata)
  THEFILE=NEWFILE.replace('.d12','.out')
  if RUN :
    if not os.path.exists(THEFILE):
      listefile.write('# '+THEFILE + ' not found. New calculation.\n')
      os.system("sbatch submit.sh %s " % ( NEWFILE.replace('.d12','') ))
    else:
      with  open(THEFILE, "rb") as FILE:
        m=reverse_readline(THEFILE, buf_size=18192)
        text=list(m)[0]
        if not p_opt.search(text) :
          listefile.write('# '+THEFILE + ' no "OPT E" string found. New calculation.\n')
          os.system("sbatch submit.sh %s " % ( NEWFILE.replace('.d12','') ))
        else:
          listefile.write('# '+THEFILE + ' "OPT E" string found. No new calculation.\n')
  else :
    if not os.path.exists(THEFILE):
      listefile.write('# '+THEFILE + ' not found. New calculation.\n')
      listefile.write("sbatch submit.sh %s \n" % ( NEWFILE.replace('.d12','') ))
    else:
      with  open(THEFILE, "rb") as FILE:
        m=reverse_readline(THEFILE, buf_size=18192)
        text=list(m)[0]
        if not p_opt.search(text) :
          listefile.write('# '+THEFILE + ' no "OPT E" string found. New calculation.\n')
          listefile.write("sbatch submit.sh %s \n" % ( NEWFILE.replace('.d12','') ))
        else:
          listefile.write('# '+THEFILE + ' "OPT E" string found. No new calculation.\n')


def perm1(C,N):
  for r in product( C,repeat=N): # cartesian product, equivalent to a nested for-loop, 
    # product(A, repeat=4) means the same as product(A, A, A, A)
    # product(A, B)  returns the same as ((x,y) for x in A for y in B)
     v=list(r)
     if v.count('B')==3 and v.count('V')==2 :
       yield (v)

t0 = time.time()

TXT="""26   0.   0.   0.   
26   0.   0.   0.5  
26   0.   0.5  0.   
26   0.   0.5  0.5  
26   0.5  0.   0.   
26   0.5  0.   0.5  
26   0.5  0.5  0.   
26   0.5  0.5  0.5  """ 

L=TXT.split("\n")
p_opt  =re.compile(br'\s+\* OPT END \- CONVERGED \* E\(AU\):\s+(?P<ene>[-+\.E\d]+)\s+POINTS\s+(?P<nb_cyc>[\d]+)\s+\*')

LIST={
      "233":{"template":"template_233_optgeom.d12","config":"config_p233_optgeom"}, 
      "323":{"template":"template_323_optgeom.d12","config":"config_p323_optgeom"}, 
      "332":{"template":"template_332_optgeom.d12","config":"config_p332_optgeom"}, 
     }
N=8
C = ['R','V','B']
v=[]
# Generate all permutations 
v=perm1(C,N)

# Convert it to list
Lpermu=[i for i in v]

with open("combi_output.txt", 'w',encoding='utf-8') as listefile:
  for CASE in ["233","323","332"]:
    s=0
    for i in Lpermu:
      s=s+1
      print(CASE,' ',s)
      a=i
      p={'V':[],'B':[],'R':[]}
      for e in range(len(a)):
        p[a[e]].append(e)
      genere_config(p,s,CASE,RUN,listefile)
    t1 = time.time()
    total = t1-t0
  listefile.write("#\n# Temps total : %g secondes" % total)
  listefile.write("# Nb de combinaison : %d" % s)
