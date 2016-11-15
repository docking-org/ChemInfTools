
import sys,os

## Writen by Trent Balius in the Shoichet Group
## calculates the tanamoto matrics. 
## fingerprint are calculated with chemaxon
## this uses a simular chemaxon comand as sea.  
## bit comparisons are calculated in python

#this function converts a string of ones and zeros to a bit
def str_to_bit(s):
#    if len(s) != 8:
#       print "warning: sting not right length"
##       print "error: sting not right length"
##       exit()
    
    ## check that the string is all zerros and ones
    for c in s:
        if not (int(c)==1 or int(c)==0):
           print "error: string not zerros and ones"
           exit()

    # note the 2^7 = 10000000
    #          2^6 = 01000000
    #          2^5 = 00100000
    #          ... 
    #      1 = 2^0 = 00000001
    #            0 = 00000000

    #i = 7 # start at 7
    i = len(s) -1 # start at 7
    #print i
    b_int = int(0)
    for c in s:
        #print bin(b_int)
        x = int(c)
        if i >= 0: 
           b_int = b_int + ((2**i) * x)
        i = i -1

    b = bin(b_int)
    #print b_int, b
    return b,b_int

# this function counts the number of ones in the bitstring
def num_bit_ones(bitstring):
    count = 0
    for c in bitstring:
        if (c!='1'):
           continue
        count = count + 1
    return count

# this function computes the tanimoto between to fingerprints
def tanimoto(fp1,fp2):
    fp1_bits = fp1.split('|')
    fp2_bits = fp2.split('|')

    if len(fp1_bits) != len(fp2_bits):
       print "ERROR: bits do not agree in lenth"

    or_num_one  = 0
    and_num_one = 0
    for i in range(len(fp1_bits)):
        #print fp1_bits[i]
        bit1,int1 = str_to_bit(fp1_bits[i])
        bit2,int2 = str_to_bit(fp2_bits[i])
        and_bit = bin(int1 & int2)
        and_num_one = and_num_one + num_bit_ones(and_bit)
        #print str(bit1) + " AND " + str(bit2) + " = " + str(bin(int1 & int2)) + ', ones = ' + str(num_bit_ones(and_bit))
        or_bit  = bin(int1 | int2)
        or_num_one = or_num_one + num_bit_ones(or_bit)
        #print str(bit1) + " OR " + str(bit2) + " = " + str(bin(int1 | int2)) +', ones = ' + str(num_bit_ones(or_bit))
        #print str(bin(bit1)) + " AND " + str(bin(bit2))
    #print and_num_one, or_num_one   
    TC = float(and_num_one) / float(or_num_one)
    #print and_num_one, or_num_one, TC
    return TC

def tversky_index(fp1,fp2,alpha,beta):
  ## (A n B)/ [(A n B) + alpha * (A\B) + beta * (B\A)) where n is the intersection and \ is the inverse complement

    fp1_bits = fp1.split('|') # A
    fp2_bits = fp2.split('|') # B

    if len(fp1_bits) != len(fp2_bits):
       print "ERROR: bits do not agree in lenth"

    A_bs_B_num_one = 0 # inverse comlement 
    B_bs_A_num_one = 0 #
    and_num_one = 0 #interection
    for i in range(len(fp1_bits)):
        #print fp1_bits[i]
        bit1,int1 = str_to_bit(fp1_bits[i])
        bit2,int2 = str_to_bit(fp2_bits[i])
        and_bit = bin(int1 & int2)
        and_num_one = and_num_one + num_bit_ones(and_bit)

        B_bs_A  = bin(~int1 & int2)
        B_bs_A_num_one = B_bs_A_num_one + num_bit_ones(B_bs_A)

        A_bs_B  = bin(int1 & ~int2)
        A_bs_B_num_one = A_bs_B_num_one + num_bit_ones(A_bs_B)
    #print and_num_one, or_num_one   
    TV = float(and_num_one) / (float(and_num_one)+ alpha *float(A_bs_B_num_one) + beta * float(B_bs_A_num_one))
    #print and_num_one, or_num_one, TC
    return TV

# end tversky_index


## this function call chemaxon and computes a bunch of fingerprints
def fingerprint_vec(SmilesString_vec,pid):
    TMPDIR = "scratch"
    #TMPDIR = "tmp"
    # this will get the users name. need to write temp file
    #pid = str(os.getpid()) # get the process idenifier so that we do not right over the same file. 
    #print pid
    name = os.popen('whoami').readlines()[0].strip()
    #print name
    #Generatemd = "/nfs/software/jchem/5.10.3/bin/generatemd"
    #Generatemd = "/nfs/soft/jchem/jchem-5.10.3/bin/generatemd"
    Generatemd = "/nfs/soft/jchem/current/bin/generatemd"
    # write smiles to file
    fh = open("/"+TMPDIR+"/"+ name +"/temp"+pid+".smi",'w')
    for SmilesString in SmilesString_vec:
       fh.write(SmilesString+'\n')
    fh.close()

    #os.popen("/raid3/software/openbabel/openbabel-2.2.1-32/bin/babel -ismi /"+TMPDIR+"/tbalius/temp.smi -osdf /"+TMPDIR+"/tbalius/temp.sdf -d")
    #os.popen("/raid3/software/openbabel/openbabel-2.2.1-32/bin/babel -isdf /"+TMPDIR+"/tbalius/temp.sdf -osmi /"+TMPDIR+"/tbalius/temp2.smi -d")
 
    comand = Generatemd + " c /"+TMPDIR+"/" + name + "/temp"+pid+".smi -k ECFP -2"
    print "runing the comand:"+comand
    output = os.popen(comand).readlines()
    #print "output:"+str(output)
    fp_vec = []
    for line in output:
       fp = line.strip('\n')
       #print fp
       fp_vec.append(fp)
    # remove the temp file. 
    os.system("rm -fr "+"/"+TMPDIR+"/" + name + "/temp"+pid+".smi")
    #print fp
    return fp_vec

## this function call chemaxon and computes the molecular Mass of the molecule
def molecularMass(SmilesString,pid):
    TMPDIR = "scratch"
    #TMPDIR = "tmp"
    # this will get the users name. need to write temp file
    #pid = str(os.getpid())
    name = os.popen('whoami').readlines()[0].strip()
    #print name
    #Generatemd = "/nfs/software/jchem/5.10.3/bin/generatemd"
    #Generatemd = "/nfs/soft/jchem/jchem-5.10.3/bin/generatemd"
    Generatemd = "/nfs/soft/jchem/current/bin/generatemd"
    # write smiles to file
    fh = open("/"+TMPDIR+"/"+ name +"/temp"+pid+".smi",'w')
    fh.write(SmilesString+'\n')
    fh.close()

    #os.popen("/raid3/software/openbabel/openbabel-2.2.1-32/bin/babel -ismi /"+TMPDIR+"/tbalius/temp.smi -osdf /"+TMPDIR+"/tbalius/temp.sdf -d")
    #os.popen("/raid3/software/openbabel/openbabel-2.2.1-32/bin/babel -isdf /"+TMPDIR+"/tbalius/temp.sdf -osmi /"+TMPDIR+"/tbalius/temp2.smi -d")

    comand = Generatemd + " c /"+TMPDIR+"/" + name + "/temp"+pid+".smi -k Mass"
    print "runing the comand:"+comand
    output = os.popen(comand).readlines()
    #print output
    #print "output:"+str(output)
    #outlines = output.split('\n')
    #lastline = outlines[len(outlines)-1]
    lastline = output[len(output) - 1]
    #print lastline
    mass = lastline.split()[1]
    #print mass
    #print fp
    # remove the temp file. 
    os.system("rm -fr "+"/"+TMPDIR+"/" + name + "/temp"+pid+".smi")
    return mass

def heavyAtoms(SmilesString,pid):
    TMPDIR = "scratch"
    #TMPDIR = "tmp"
    # this will get the users name. need to write temp file
    #pid = str(os.getpid())
    name = os.popen('whoami').readlines()[0].strip()
    #print name
    #Generatemd = "/nfs/software/jchem/5.10.3/bin/generatemd"
    #Generatemd = "/nfs/soft/jchem/jchem-5.10.3/bin/generatemd"
    Generatemd = "/nfs/soft/jchem/current/bin/generatemd"
    # write smiles to file
    fh = open("/"+TMPDIR+"/"+ name +"/temp"+pid+".smi",'w')
    fh.write(SmilesString+'\n')
    fh.close()

    #os.popen("/raid3/software/openbabel/openbabel-2.2.1-32/bin/babel -ismi /"+TMPDIR+"/tbalius/temp.smi -osdf /"+TMPDIR+"/tbalius/temp.sdf -d")
    #os.popen("/raid3/software/openbabel/openbabel-2.2.1-32/bin/babel -isdf /"+TMPDIR+"/tbalius/temp.sdf -osmi /"+TMPDIR+"/tbalius/temp2.smi -d")

    comand = Generatemd + " c /"+TMPDIR+"/" + name + "/temp"+pid+".smi -k Heavy"
    print "runing the comand:"+comand
    output = os.popen(comand).readlines()
    #print "output:"+str(output)
    #outlines = output.split('\n')
    #lastline = outlines[len(outlines)-1]
    lastline = output[len(output) - 1]
    heavy = lastline.split()[1]

    print output,outlines,lastline,heavy
    # remove the temp file. 
    os.system("rm -fr "+"/"+TMPDIR+"/" + name + "/temp"+pid+".smi")
    #print fp
    return heavy


## this function reads in smiles and writes out the footprints.
## it returns a footprint vector
def get_fp(infile,outfile,pid):
  fpvec = []
  file = open(infile,'r')
  lines = file.readlines()
  file.close()
  file1 = open(outfile,'w')
  smiles_vec = []
  for line in lines:
     splitline = line.split()
     if len(splitline) > 2:
        print "ERROR:len(smiles) > 2"
        exit()
     print splitline

     smiles_vec.append(splitline[0])
     #print "simles = " + str(smiles);
  fp_vec = fingerprint_vec(smiles_vec,pid)
  #   fp = fingerprint(smiles,pid)
     #print "fingerprint = " + str(fp);
  for fp in fp_vec:
     file1.write("fingerprint = " + str(fp)+'\n')
     #fpvec.append(fp)
  file1.close()
  return fp_vec

#def main():
#  if not (len(sys.argv) == 4 or len(sys.argv) == 5): # if no input
#     print "ERORR"
#     print "syntexs: python tanimoto_cal_axon.py -one smiles1 outputprefix"
#     print "         this produces a squere symestric matrix of set1 with itself. "
#     print "syntexs: python tanimoto_cal_axon.py -two smiles1 smiles2 outputprefix"
#     print "         this produces a rectangular non-symestric matrix of set1 to set2"
#     return
#
#  oneortwo    = sys.argv[1]
#  smilesfile1 = sys.argv[2]
#  if oneortwo == "-one":
#    outfileprefix = sys.argv[3]
#  elif  oneortwo == "-two":
#    smilesfile2 = sys.argv[3]
#    outfileprefix = sys.argv[4]
#  else:
#      print "the frist parameter must be -one or -two."
#      exit()
#
#  outfile1 = outfileprefix +'.1.fp'
#  fpvec1 = get_fp(smilesfile1,outfile1)
#  if oneortwo == "-one":
#    fpvec2 = fpvec1
#  if oneortwo == "-two":
#    outfile2 = outfileprefix +'.2.fp'
#    fpvec2 = get_fp(smilesfile2,outfile2)
#
#  outfileM = outfileprefix +'.matrix'
#
#
#  file1 = open(outfileM,'w')
#  for fp1 in fpvec1:
#     flag_frist = True 
#     for fp2 in fpvec2:
#        TC = tanimoto(fp1,fp2)
#        if (flag_frist):
#           flag_frist = False
#        else:
#           file1.write(',')
#        file1.write('%f' % TC )
#     file1.write('\n' )
#  file1.close()
#main()

