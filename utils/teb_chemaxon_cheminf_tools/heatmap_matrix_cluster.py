#! /raid1/people/tbalius/zzz.virtualenvs/sgehead_python_env/bin/python

## Writen by Trent E. Balius in B. Shoichet group

## This was adapted from the web:
## http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python

import matplotlib  # must import first
matplotlib.use('Agg')  # allows you to not have an x-server running
#these lines must be first, if pylab is imported first it ruins this

import sys, os
import copy
import math
#import matplotlib
import scipy
import numpy
import pylab
import scipy.cluster.hierarchy as sch


def getlabel(labfilename):
  # import the label information
  
  file = open(labfilename)
  lines = file.readlines()
  file.close()
  
  #m_lab = len(lines)
  labels = []
  
  for line in lines:
      splitline = line.split(',')
      labels.append(splitline[0].strip())
      #print splitline[0], splitline[1]
  return labels


def mat_to_vector(Mat):
    m = len(Mat)
    n = len(Mat[0])

    if (m != n):
        print "inconsitancy in numbers of rows and columns in the matrix."
        sys.exit()

    print m,n

    X = scipy.zeros([m,n])
    Xvec = scipy.zeros(n*(n-1)/2)

    count2    = 0

    ## converts from a 2D array to Scipy Matrix 
    for i in range(0,n):
        for j in range(0,n):
               ## 1 - tc is more like a distance than tc.
               #X[i,j] = -Mat[i][j] + 1.0
               X[i,j] = Mat[i][j]

    for i in range(0,n):
        for j in range(i+1,n):
               ## 1 - tc is more like a distance than tc.
               #Xvec[count2] = -Mat[i][j] + 1.0
               Xvec[count2] = Mat[i][j]
               count2 = count2+1

    return X,Xvec

def get_cluster(X,labels,clusttypethreshold,dirname):
    print "in function get_cluster"
    if len(X) != len(labels):
       print "len(X) != len(labels)"
       print len(X), len(labels)
       exit()

    Xnew,Xvec = mat_to_vector(X)

    #Y = sch.linkage(Xvec, method='complete')
    #Y = sch.linkage(Xvec, method='single')
    Y = sch.linkage(Xvec, method=clusttype)
    clusters = sch.fcluster(Y, threshold, 'distance')
    print clusters
    cluster_list  = [] # list of pdb names in each cluster
    cluster_sizes = [] # list of the size of each cluster

    numOfClusters = clusters.max()

    ## intialize array that will store the labels for each cluster
    for i in range(numOfClusters):
        cluster_list.append('c'+ str(i+1)+' -- ')
        cluster_sizes.append(0)
    #
    ## fill array with labels by appending the string assosiated with each cluster
    for i in range(len(clusters)):
        cluster_list[clusters[i]-1] = cluster_list[clusters[i]-1] + labels[i] + ','
        cluster_sizes[clusters[i]-1] = cluster_sizes[clusters[i]-1] + 1

    ## write the cluster
    for i in range(numOfClusters):
        print cluster_list[i]

    ## write the cluster with more than 3 members
#    os.system('rm -rf   large_clusters'+dirname)
#    os.system('mkdir -p large_clusters'+dirname)
#    os.chdir('large_clusters'+dirname)


    filename = "cluster_rep.txt"
    fh_rep = open(filename,'w')
    print " larger clusters: "
    for i in range(numOfClusters):
       filename = "cluster" + str(i+1) + ".txt"
       fh = open(filename,'w')
       fh.write(cluster_list[i].replace(' ','').replace(',','\n').replace('-','\n'))
       fh.close()
       fh_rep.write(cluster_list[i].replace(' ','').split(',')[0].split('-')[2]+'\n')

       if cluster_sizes[i] > 3:
           print "  " + cluster_list[i]
#           name = cluster_list[i].split('--')[0].replace(' ','')
#           mols = cluster_list[i].split('--')[1].split(',')
           ## get images from zinc
#           os.system('mkdir -p '+ name)
#           os.chdir(name)
#           fout = open( name+"_info.txt",'w')
#           for mol in mols:
#               print mol
#               fout.write(mol+'\n')
#               os.system('wget http://zinc.docking.org/img/sub/' + mol.replace('C','').replace(' ','')+'.gif')
#               #os.system('wget http://zinc.docking.org/substance/' + mol.replace('C','').replace(' ','')+'')
#           fout.close()
#           os.chdir('../')
#    os.chdir('../')
    fh_rep.close()
              
    return Y

#def get_vec_2(mat,index):
#    vec = scipy.zeros([len(index),1])
#    if len(index) != len(mat):
#       print "warning"
#       print len(index), len(mat)
#       exit()
#    j = 0
#    for i in index:
#        vec[j] = mat[j,i]
#        j=j+1
#    return vec
#
#def get_vec_1(mat,index):
#    vec = scipy.zeros([len(index),1])
#    if len(index) != len(mat[0,:]):
#       print "warning"
#       print len(index), len(mat)
#       exit()
#    j = 0
#    for i in index:
#        vec[j] = mat[i,j]
#        j=j+1
#    return vec

#def get_min(Xorg,label1,label2,N):
#  X = copy.copy(Xorg)
#  print "In function get_min"
#  print "loop over matix", N, "times"
#  count = 0
#
#  os.system('rm -rf   min_pairs')
#  os.system('mkdir -p min_pairs')
#  os.chdir('min_pairs')
#  fout = open("links_info.txt",'w')
#  #os.system('cd       min_pairs')
#  while (count < N):
#    minname = 'minpair'+str(count)
#    os.system('mkdir -p minpair'+str(count))
#    os.chdir('minpair'+str(count))
#    indecies1 = numpy.argmin(X,axis=0)
#    indecies2 = numpy.argmin(X,axis=1)
#
#    vec1 = numpy.amin(X,axis=0)
#    vec2 = numpy.amin(X,axis=1)
#
#    j1 = numpy.argmin(vec1)
#    i1 = indecies1[j1]
#
#    i2 = numpy.argmin(vec2)
#    j2 = indecies2[i2]
#
#    print i1,j1, X[i1,j1]
#    print i2,j2,label1[j2],label2[i2], X[i2,j2]
#
#    fout.write('%s,%s,%s,%d,%d,%f\n' % (minname,label2[i2],label1[j2],i2,j2,X[i2,j2]))
#    # get images from zinc:
#    #os.system('wget http://zinc.docking.org/substance/' + label1[j2].replace('C',''))
#    #os.system('wget http://zinc.docking.org/substance/' + label2[i2].replace('C',''))
#    os.system('wget http://zinc.docking.org/img/sub/' + label2[i2].replace('C','') +'.gif')
#    os.system('wget http://zinc.docking.org/img/sub/' + label1[j2].replace('C','') +'.gif')
#
##    if ( i1 != i2 or j1 != j2): 
##        print "Error:: i1!=i2 or j1 == j2"
##        print j1, j2, i1, i2
##        continue
#    X[i2,j2] = 100
#    count = count + 1
#    os.chdir('../')
##    os.system('cd ../')
#  os.chdir('../')
#  fout.close()
#  #exit() 
#  return

def import_mat(matfilename):
     # Import data from matrix file:
     file = open(matfilename)
     lines = file.readlines()
     file.close()
     
     
     m = len(lines)
     n = len(lines[0].split(','))
     
     if (m != n):
         print "inconsitancy in numbers of rows and columns in the matrix."
     
     print m,n
     
     X = scipy.zeros([m,n])
     Xvec = scipy.zeros(n*(n-1)/2)
    
     countline = 0
     count2    = 0 
     for line in lines:
         line = line.strip('\n')
         splitline = line.split(',')
         if (n != (len(splitline))):
             print "ERROR: n != (len(splitline), inconsitancy in number of elements in rows"
             sys.exit()
     
         for i in range(0,n):
             val = float(splitline[i])
             X[countline,i] = 1-val 
             #X[countline,i] = 1-val ## 1-Tc is a metric
         countline = countline + 1
     return X ,n,m

pylab.matplotlib.use('Agg')

ZERRO = 0.0

matfilename  = sys.argv[1]
labfilename = sys.argv[2]
threshold   = float(sys.argv[3])
clustertype = sys.argv[4]

print "mat_filename = "+ matfilename  
print "lab_filename = "+ labfilename  
print "threshold = "+ str(threshold)
print "cluster type = "+ clustertype  

#Y = sch.linkage(Xvec, method='complete')
#Y = sch.linkage(Xvec, method='single')


if not ( clustertype == "complete" or clustertype == "single"):
    print "cluster type must be complete or single"

labels1 = getlabel(labfilename)
#labels2 = getlabel(lab2filename)

X,n,m           = import_mat(matfilename)
#dist_mat1,n1,m1 = import_mat(dist_mat1_fn)

#threshold =  0.7
#threshold =  0.47
#threshold = 0.51
#threshold =  0.4
## create a distance matrix --> dendogram by comparing all rows
Y1 = get_cluster(X,labels1,clustertype,threshold,'set2')

#get_min(X,labels1,labels2,50)

#
fig = pylab.figure(figsize=(8,8))
ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
##Z1 = sch.dendrogram(Y, orientation='right')
Z1 = sch.dendrogram(Y1, orientation='right',color_threshold=threshold)
matplotlib.pyplot.plot([threshold,threshold],[0,10*m],'k--') # draws a datshed line where dendogram is cut.
##help(sch.dendrogram)
#ax1.set_xticks([])
ax1.set_yticks([])
##exit()
##print ax1.get_ylim()
##ax1.set_ylim(-1, n)

## Plot distance matrix.
axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
idx1 = Z1['leaves']
#idx2 = Z2['leaves']
X = X[idx1,:]
X = X[:,idx1]
#X = X[:,idx2]

##labels_sort = labels[idx1]
## make sorted label list
#labels_sort = []
#for i in idx1:
#  labels_sort.append('c'+str(clusters[i]) + '-'+ labels[i])

cdict = {'red': ((0.0, 0.0, 0.0),
                  (0.0, 0.0, 0.0), 
                  (1.0, 1.0, 1.0)),
          'green': ((0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (1.0, 1.0, 1.0)),
          'blue': ((0.0, 0.0, 0.0),
                   (0.0, 0.0, 0.0),
                   (1.0, 1.0, 1.0))}


my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,20)

im = axmatrix.imshow(X, aspect='auto', origin='lower',interpolation='nearest', cmap=my_cmap)


#im.set_clim(0,threshold)
#im.set_clim(threshold,1)
im.set_clim(0.2,1)
axmatrix.set_ylim(-0.5, m-0.5)
axmatrix.set_xlim(-0.5, n-0.5)
axmatrix.set_yticks([])
axmatrix.set_xticks([])

# Plot colorbar.
axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
pylab.colorbar(im, cax=axcolor)
fig.show()
fig.savefig('dendrogram.png',dpi=600)

