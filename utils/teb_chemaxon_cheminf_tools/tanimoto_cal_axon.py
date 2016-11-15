
import sys,os
import tanimoto_tversky_cal_axon_lib as tccalc 

## Writen by Trent Balius in the Shoichet Group
## calculates the tanamoto matrics. 
## fingerprint are calculated with chemaxon
## this uses a simular chemaxon comand as sea.  
## bit comparisons are calculated in python

def main():
  if not (len(sys.argv) == 4 or len(sys.argv) == 5): # if no input
     print "ERORR"
     print "syntexs: python tanimoto_cal_axon.py -one smiles1 outputprefix"
     print "         this produces a squere symestric matrix of set1 with itself. "
     print "syntexs: python tanimoto_cal_axon.py -two smiles1 smiles2 outputprefix"
     print "         this produces a rectangular non-symestric matrix of set1 to set2"
     return

  pid = str(os.getpid()) # get the process idenifier so that we do not right over the same file. 
  print pid

  oneortwo    = sys.argv[1]
  smilesfile1 = sys.argv[2]
  if oneortwo == "-one":
    outfileprefix = sys.argv[3]
  elif  oneortwo == "-two":
    smilesfile2 = sys.argv[3]
    outfileprefix = sys.argv[4]
  else:
      print "the frist parameter must be -one or -two."
      exit()

  outfile1 = outfileprefix +'.1.fp'
  fpvec1 = tccalc.get_fp(smilesfile1,outfile1,pid)
  if oneortwo == "-one":
    fpvec2 = fpvec1
  if oneortwo == "-two":
    outfile2 = outfileprefix +'.2.fp'
    fpvec2 = tccalc.get_fp(smilesfile2,outfile2,pid)

  outfileM = outfileprefix +'.tanimoto.matrix'

  #print len(fpvec2)
  #print len(fpvec1)
  #exit(0)

  file1 = open(outfileM,'w')
  for fp1 in fpvec1:
     flag_frist = True 
     for fp2 in fpvec2:
        #print fp1
        #print fp2
        TC = tccalc.tanimoto(fp1,fp2)
        if (flag_frist):
           flag_frist = False
        else:
           file1.write(',')
        file1.write('%f' % TC )
     file1.write('\n' )
  file1.close()

  alpha = 0.2
  beta  = 0.2
  outfileM = outfileprefix +'.tversky.'+str(alpha)+'.'+str(beta)+'.matrix'
  file1 = open(outfileM,'w')
  for fp1 in fpvec1:
     flag_frist = True
     for fp2 in fpvec2:
        TV = tccalc.tversky_index(fp1,fp2,alpha,beta)
        if (flag_frist):
           flag_frist = False
        else:
           file1.write(',')
        file1.write('%f' % TV )
     file1.write('\n' )
  file1.close()

main()

