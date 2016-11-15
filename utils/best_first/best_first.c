// modified by Trent Balius 2016


/* Copyright (C) 2008 Michael J Keiser, Brian K Shoichet, 
 * John J Irwin, Regents of the University of California
 *
 * Calculates cluster for docking hit lists using best first clustering
 * Start with best scoring molecule, find everything close in Tc, 
 * output to cluster, remove from list. on to best scoring molecule 
 * remaining and repet    
 */

#include <math.h>
#include "best_first.h"
#include "fast_tanimoto.h"
#include <stdio.h>

// fprintA is the bunch of list of binarys concatanated together.
// countA is the count of the number of ones in each fingerprint
// lenA is the is the number of fingerprints in fprintA
// db_len is the number of bit in one fingerprint (same for all fingerprint).  

void best_first(uint16_t *fprintA, int *countA, int lenA,  int db_len, char **names, double maxtc, int max_clusters) {

  int i, j;                       /* Generic counters         */
  //double tan,maxtc;                     /* Tanimoto coefficient     */
  double tan;                     /* Tanimoto coefficient     */
  uint16_t *fpA;                  /* Hold fpA[i] for speed    */ 

  //int indexlist[lenA];
  //int newlist[lenA];
  int * indexlist;
  int * newlist;
  int *swap;
  //int *temp;
  // this 3 pointers are used to swap arrays.
  //int *temp1;
  //int *temp2;

  //int temp[lenA];
  int cntA;
  int listsize = lenA;
  int count = 0;
  int clustercount = 0;
  char filename[16];
  FILE *fp;
  printf("I AM HERE\n");
  printf("tc threshold =%f\n",maxtc);

  indexlist = (int *)malloc(sizeof(int)*lenA);
  newlist = (int *)malloc(sizeof(int)*lenA);

  //printf("I AM HERE in best_first(1)\n");
  // intialize the arrays 
  for (i =0; i< lenA; i++){
       indexlist[i] = i;
       newlist[i] = -1;
  }  
  //printf("I AM HERE in best_first(2)\n");

  //MN = lenA*lenB;

  //maxtc = 0.6;
  //fp=fopen("clusters.txt", "w")
  /* Loop over fingerprints calculate tanimoto */

  while (listsize > 0 ){ // while there are still molecules
    if (clustercount > max_clusters) {
         printf("The maxium number of clusters allowed is reached. stoping clustering . . .");
         break;  // stop if we reach the max number of clusters.
    }
    //printf("I AM HERE in best_first(in while loop)\n");
    // open up file for cluser 
    sprintf(filename, "cluster_%05d.txt",clustercount);
    fp=fopen(filename, "w");
    fprintf(fp, "cluster%d\n",clustercount);
    
    count = 0;
    fpA = &(fprintA[indexlist[0]*db_len]); // index
    cntA = countA[indexlist[0]]; // number of 1 in the binary string. 
    for (i = 0; i < listsize; i++) {
       //fpA = &(fprintA[i*db_len]);
       //cntA = countA[i];
       //printf("I AM HERE in best_first(in for loop)\n");
       //printf("I AM HERE in best_first(in for loop):%d,%d,%d\n",i,listsize,db_len);
       //printf("tc compare: %d,%d\n",indexlist[0],indexlist[i]);
       tan = tanimoto_coeff(fpA, cntA, &(fprintA[indexlist[i]*db_len]),
                            countA[indexlist[i]], db_len);
 
       /* Score comparison */
       if (tan > maxtc){
           // output name to file
           printf("good tc:%d,%d,%s,%s\n",indexlist[0],indexlist[i],names[indexlist[0]],names[indexlist[i]]);
           fprintf(fp, "%d,%s,%f\n",indexlist[i],names[indexlist[i]],tan);
           // do not add to list
       }
       else if (tan <= maxtc) {
            // keep in list (add to newlist for next pass
            newlist[count] = indexlist[i];
            count++;
       }
    }
    // close file for cluster
    fclose(fp); 
    //printf("I AM HERE in best_first(3)\n");
    // swap the lists.
    //printf("pointers(1):%d,%d\n",*newlist,*indexlist);
    //printf("address(1):%d,%d\n",&newlist,&indexlist);
    // this 3 pointers are used to swap arrays.
    swap = newlist;
    newlist = indexlist;
    indexlist = swap;
    //printf("pointers(2):%d,%d\n",*newlist,*indexlist);
    //printf("address(2):%d,%d\n",&newlist,&indexlist);

    // intialize
    for (i = 0; i < listsize; i++) {
         //printf("nl=%d,il=%d\n",newlist[i],indexlist[i]);
         newlist[i] = -1;
    }
    listsize = count-1;
    clustercount++;
  }
  // deallocate memory for these arrays because malloc is used.
  free(indexlist);
  free(newlist);

  return;
}
