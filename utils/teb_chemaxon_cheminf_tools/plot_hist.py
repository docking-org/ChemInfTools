#! /usr/bin/python2.6
## this uses the version of python on sublime 

## This was adapted from the web:
## http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python

import sys
import copy
import math
import matplotlib
import scipy
import numpy
import pylab
#import scipy.cluster.hierarchy as sch
#from scipy.optimize import opt
import scipy.optimize as opt
import scipy.linalg   as la
import scipy.stats    as stats


def read_data_file(filename):

    print filename
    file = open(filename)
    lines = file.readlines()
    file.close()
    m = len(lines)
    data = scipy.zeros([m,1])
    i =0
    for line in lines:
        line = line.strip('\n')
        splitline = line.split()
        #print len(splitline)
        if (len(splitline) != 1): # data file has one column
            print line
            continue
        #print splitline[5] 
        data[i] = float(splitline[0])
        i = i+1
    return data



if len(sys.argv) != 3:
   print "error:  this program takes 1 input filename and a title  "    
   exit()

filename1     = sys.argv[1]
titlestr      = sys.argv[2]

print "filename1 = " + filename1
print "title     = " + titlestr

data = read_data_file(filename1)

n, bins, patches = pylab.hist(data,bins=50,range=[0,1])

midbin = scipy.zeros([len(n),1])
for i in range(0,len(bins)-1):
    midbin[i] = (bins[i] + bins[i+1])/2

print midbin
print bins
print len(n) 
print len(bins) 
#fig = pylab.figure(figsize=(16,8))
fig = pylab.figure(figsize=(8,8))
#axis = fig.add_axes([0.1,0.1,0.3,0.3])
axis = fig.add_axes([0.1,0.1,0.8,0.8])

#print bins, n
#print patches, n
#im = axis.plot(X,Y,'-',[0,100],[0,100],'--')
im = axis.plot(midbin, n,'b-') #,[0,100],[0,100],'--')
#im = axis.plot(X,Y,'.',[threshold1,threshold1],[Ymin,Ymax],'-',[Xmin,Xmax],[threshold2,threshold2],'-')
#axis.set_title('file='+xyfilename)

#axis.set_ylim(math.floor(Ymin),math.ceil(Ymax))
#axis.set_xlim(math.floor(Xmin),math.ceil(Xmax))
axis.set_ylim(min(n)-1,max(n)+1)
axis.set_xlim(-0.1, 1.1)
#axis.set_yticks([])
#axis.set_xticks([])
#axis = fig.add_axes([0.5,0.1,0.3,0.3])
#im = axis.plot(Xscore,Yscore,'o') #,[0,100],[0,100],'--')
pylab.title(titlestr)

fig.show()
fig.savefig('fig.png',dpi=600)



