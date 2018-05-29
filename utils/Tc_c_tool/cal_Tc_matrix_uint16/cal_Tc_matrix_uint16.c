// modified by Trent Balius 2016
// modified by Jiankun Lyu 2016


/* Copyright (C) 2008 Michael J Keiser, Brian K Shoichet, 
 * John J Irwin, Regents of the University of California
 *
 * Calculates cluster for docking hit lists using best first clustering
 * Start with best scoring molecule, find everything close in Tc, 
 * output to cluster, remove from list. on to best scoring molecule 
 * remaining and repet    
 */

#include <math.h>
#include "cal_Tc_matrix_uint16.h"
#include "../common/fast_tanimoto.h"
#include <stdio.h>

// fprintA is the bunch of list of binarys concatanated together.
// countA is the count of the number of ones in each fingerprint
// lenA is the is the number of fingerprints in fprintA
// db_len is the number of bit in one fingerprint (same for all fingerprint).  

void cal_Tc_matrix(uint16_t *fprintA, int *countA, int lenA, char **namesA, uint16_t *fprintB, int *countB, int lenB, char **namesB, int db_len, char **prefix) {

  int i, j;                       /* Generic counters         */
  //double tan,maxtc;                     /* Tanimoto coefficient     */
  double tan;                     /* Tanimoto coefficient     */
  uint16_t *fpA;                  /* Hold fpA[i] for speed    */
  uint16_t *fpB;                  /* Hold fpB[i] for speed    */

  //int indexlist[lenA];
  //int newlist[lenA];
  int * indexlistA;
  int * indexlistB;

  int cntA;
  int cntB;
  int listsizeA = lenA;
  int listsizeB = lenB;
  char filename[16];
  FILE *fp,*fp_;
  char outputfile[50];
  printf("Start\n");
  //fflush( stdout );

  indexlistA = (int *)malloc(sizeof(int)*lenA);
  indexlistB = (int *)malloc(sizeof(int)*lenB);

  //printf("I AM HERE in best_first(1)\n");
  // intialize the arrays 
  for (i=0; i<lenA; i++){
       indexlistA[i] = i;
  }
  for (i=0; i<lenB; i++){
       indexlistB[i] = i;
  }  
  //printf("I AM HERE in best_first(2)\n");

  /* Loop over fingerprints calculate tanimoto */
  //fp=fopen("Tc_matrix", "w");
  sprintf(outputfile, "%s_max_TC.col",prefix);
  fp_=fopen(outputfile,"w");
  char *zinc_id_ori;
  for (i = 0; i < listsizeA; i++){ // while there are still molecules
    
    fpA = &(fprintA[indexlistA[i]*db_len]); // index
    cntA = countA[indexlistA[i]]; // number of 1 in the binary string.
    double tan_max = 0.0; //
    char *zinc_id;
    zinc_id_ori = namesA[indexlistA[i]];
 
    for (j = 0; j < listsizeB; j++) {

       fpB = &(fprintB[indexlistB[j]*db_len]);
       cntB = countB[indexlistB[j]];

       tan = tanimoto_coeff(fpA, cntA, fpB, cntB, db_len);
 
       if (tan > tan_max){
          tan_max = tan;
          zinc_id = namesB[indexlistB[j]]; 
          }

       /* Tc output */
       //printf("%d,%d,%s,%s\n",indexlistA[i],indexlistB[j],namesA[indexlistA[i]],namesB[indexlistB[j]]);
       //if (j == listsizeB-1){
          //fprintf(fp, "%f", tan);
          //}
       //else{
          //fprintf(fp, "%f,", tan);
          //}   
       }
    //fprintf(fp, "\n");
    fprintf(fp_, "%s,%f,%s\n",zinc_id_ori, tan_max, zinc_id);
    }
    // close file
    //fclose(fp); 
    fclose(fp_);

  // deallocate memory for these arrays because malloc is used.
  free(indexlistA);
  free(indexlistB);

  return;
}
