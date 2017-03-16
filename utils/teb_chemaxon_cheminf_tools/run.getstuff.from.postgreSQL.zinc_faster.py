
# written by Trent E. Balius and Jiankun Lyu 

import sys,os,binascii
#import tanimoto_tversky_cal_axon_lib as tclib 
#  this function are takin from the web. 
# http://stackoverflow.com/questions/7396849/convert-binary-to-ascii-and-vice-versa
#def text_to_bits(text, encoding='utf-8', errors='surrogatepass'):
#    bits = bin(int(binascii.hexlify(text.encode(encoding, errors)), 16))[2:]
#    return bits.zfill(8 * ((len(bits) + 7) // 8))
#
#def text_from_bits(bits, encoding='utf-8', errors='surrogatepass'):
#    n = int(bits, 2)
#    return int2bytes(n).decode(encoding, errors)
#
#def int2bytes(i):
#    hex_string = '%x' % i
#    n = len(hex_string)
#    return binascii.unhexlify(hex_string.zfill(n + (n & 1)))

if len(sys.argv) != 3:
   print "wrong number of inputs:\n (1) input file with list of zinc id\n (2) output file name\n"
   exit()

inputfile  = sys.argv[1]
outputfile = sys.argv[2]

#machine = "samekh"
machine = "nun"

#file1 = open("zinc_id.5000.txt",'r')
file1 = open(inputfile,'r')
namelist = file1.readlines()
file1.close()
#file2 = open("zinc_id.5000.smi",'w')
file2 = open(outputfile,'w')

flush_size = 5000 # this is how meny zinc codes we want to quiery at one time. 
                 # and then write to file
count = 0
countnl = 0
#use a dict to covert the zid to 'ZINC0'
#make the generated smi file more readable
zinc_num_dict = {}
flag_start = True
command = 'echo "select smiles,sub_id from substance where ('
for line in namelist:
    splitline = line.split(',')
    zid = splitline[0].replace("ZINC","").lstrip("0")
    #print zid
    #print splitline[0]
    zinc_num_dict[str(zid)] = str(splitline[0])
    # get the simles
    #comand = 'echo "select smiles,sub_id from substance where sub_id='+zid+';" | psql -h samekh -U zincread zinc15'
    #po = os.popen('echo "select smiles,sub_id from substance where sub_id='+zid+';" | psql -h samekh -U zincread zinc15')
    if (flag_start):
       command = command + 'sub_id='+zid
       flag_start = False
    else: 
       command = command+'or sub_id='+zid

    count = count + 1

    if (count == flush_size or countnl >= len(namelist)-1):
        #print zinc_num_dict
        #po = os.popen(command + ');" | psql -h samekh -U zincread zinc15')
        po = os.popen(command + ');" | psql -h '+machine+' -U zincread zinc15')
        
        for linepo in po:
            if (len(linepo.split())!=3):
                continue
            print linepo
            smiles = linepo.split()[0]
            zinc_id = linepo.split()[2]+'\n'
            if smiles == 'smiles':
                continue
            else:
                file2.write('%s %s'%(smiles, zinc_num_dict[zinc_id]))
        count = 0
        command = 'echo "select smiles,sub_id from substance where ('
        flag_start = True
    countnl = countnl + 1 

print countnl,count 
    
#    #po[2]
##   for linepo in po:
##       print linepo
#    # get figureprints
#    po = os.popen('echo "select data,ecfp4_fk from fingerprints.ecfp4_new WHERE ecfp4_fk='+zid+';" | psql -h samekh -U zincread zinc15')
#    linepo = po.readlines()[2]
#    fingerprint = linepo.split()[0]
#    print fingerprint
#    if fingerprint == "(0":
#       #fplist.append("NA")
#       continue
#    #print fingerprint.lstrip("\\")
#    #bin_fp = binascii.a2b_uu(fingerprint.lstrip("\\"))
#    bin_fp = text_to_bits(fingerprint)
#    print bin_fp
#    fplist.append(bin_fp)
#    idlist.append(zid)
#    #print linepo
file2.close()
#fp1 = fplist[0]
#for i,fp2 in enumerate(fplist):
#    tc = tclib.tanimoto(fp1,fp2)
#    print idlist[0],",",idlist[i],":",tc
#    
