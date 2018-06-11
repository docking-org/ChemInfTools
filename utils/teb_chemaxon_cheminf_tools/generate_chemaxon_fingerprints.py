import sys, os
import tanimoto_tversky_cal_axon_lib as tccalc 
# writen by Trent Balius in 2016

def main():
  if (len(sys.argv) != 3): # if no input
     print "ERORR"
     print "syntexs: python generate_chemaxon_fingerprints.py smiles1 outputprefix"
     return

  smilesfile1   = sys.argv[1]
  outfileprefix = sys.argv[2]

  outfileF = outfileprefix +'.fp'

  pid = str(os.getpid()) # get the process idenifier so that we do not right over the same file. 

  # read in names from smiles file.  
  names = []
  file1 = open(smilesfile1,'r')
  for line in file1:
      print line
      name = line.split()[1]  
      #print line, name
      names.append(name)
  file1.close()

  
  # if the fp file already exists, just read in the fingerprints.  
  # otherwise caluclate the fingerprints.
  if (os.path.isfile(outfileF) ):
     print outfileF + "exists."
     exit()
  else:
     fpvec1 = tccalc.get_fp(smilesfile1,outfileF,pid)

main()

