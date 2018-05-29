

#include "cal_Tc_matrix_uint16.h"
#include "../common/fast_tanimoto.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

// read in file, read in all footprints
// call Tc calculation routine. 
//

int line_count(FILE *n)
{
   char c;
   int lines = 0;
   
   while ((c = fgetc(n)) != EOF)
   {
   if (c == '\n') ++lines;
   }
   return lines;
}
// line is a string, len is the lenth of the input string
//char* process_str(char* line, int len){
void process_str(char* line, int len, char* newline, int newlen){
  //char outputline[len];
  bool flag =false;
  int  count = 0;
  int  countones = 0;

  for (int i=0; i<len; i++){
      //printf("%c::%d\n",line[i],count);
      if (count>newlen) printf("\n\ncount>newlen\n\n");
      if (line[i]=='|') continue;
      if (line[i]==' ') continue;
      if (line[i]=='\n') break;
      if (flag) {
         if (line[i] != '1' && line[i] != '0'){
             printf("\n\nERROR\n\n");
             exit(1);
         }
         newline[count] = line[i];
         count++;
      }
      if (line[i]=='1') countones++; // count how meny ones are there, for debuging
      if (line[i]=='=') flag = true;
  }
  //printf("count ones: %d\n",countones);
  //return outputline;
  return ;
}

///////////////////////

void process_smiles_str(char* line, int len, char* newline, int newlen){
  //char outputline[len];
  bool flag =false; // names in in the seconed column, after a space. 
  int  count = 0;

  for (int i=0; i<len; i++){
      //printf("%c::%d\n",line[i],count);
      if (count>newlen) printf("\n\ncount>newlen\n\n");
      if (line[i]=='\n') break;
      if (flag) {
         newline[count] = line[i];
         count++;
      }
      if (line[i]==' ' || line[i]=='\t') flag = true; // this will indecated when we are at the second column. 
  }
  return ;
}

////////////////////

void uint16tostr( uint16_t bits){

  int i;
  //uint16_t temp;
  int tempi;
  double tempbitsd,temp, tempd;
  uint16_t curentbit,tempbitsi;  
  
  tempbitsi = bits; 
  for(i=0;i<16;i++){
     //printf("%d\n",tempbitsi);
     tempbitsd  =  ((double) tempbitsi); // convert to double
     tempd =   pow(2.0, ( 15 - i ));
     temp  = floor(tempbitsd / tempd); 
     tempi = (int) temp; // convert to int
     curentbit = (uint16_t) temp;
     if (tempi ==  1) {
        //printf("1\n");
        printf("1");
        tempbitsd = tempbitsd - tempd;
        tempbitsi = (uint16_t) tempbitsd; // remove the left most part of the bits
     } 
     //else printf("0\n");
     else printf("0");
  } 
  printf("\n");
  return;
}

// this function converts a string of ones and zeros
// to a unit16_t integer. 
// str is a list of chartures, bits is a integer. 
void str2uint16(char* str, uint16_t* bits){
//
//     2^15 -> 10000000,00000000
// ...
//     2^8  -> 00000001,00000000
// 128 2^7  -> 00000000,10000000 
// ...
//   4 2^2  -> 00000000,00000100
//   2 2^1  -> 00000000,00000010
//   1 2^0  -> 00000000,00000001
//   0      -> 00000000,00000000
/*
*/
  int i;
  double tempd;
  uint16_t tempbits, temp;

  //printf("I AM HERE: %s\n",str);
  tempbits =  ((uint16_t) 0);
  for(i=0;i<16;i++){
     //if (str[i] == "1") {temp =  pow(2,(15-i));}
     if (str[i] == '1') {
         tempd =   pow(2.0, ( 15 - i ));
         temp  = (uint16_t) tempd;
         //printf("one at position: %d, value=%d\n",(15-i),temp);
     }
     else {
         temp = ((uint16_t) 0);
     }
     //printf("I AM HERE1.4:%d\n",i);
     tempbits = tempbits | temp;  
     //printf("I AM HERE1.6:\n");
     //tempbits = tempbits + temp;  
  }
  //printf("I AM HERE2:");
  //uint16tostr(tempbits);
  //printf("bits_int=%d\n",tempbits);
  //bits = (uint16_t) tempbits;
  *bits = (uint16_t*) tempbits;
  return; 
}

// wraper function that calls str2unt16.
void str2uint16_wraper(char* str,int len, uint16_t* bits_vec) {

   char temp_str[17]; 
   int temp_len = len/16;
   //int temp_len = len;
   int count = 0;
   uint16_t temp;
   //printf("%d\n",temp_len);
   for(int i=0; i<temp_len;i++){
      for(int j=0; j<16;j++){
         temp_str[j] = str[count];
         count++;
      }
      //printf("%s\n",temp_str);
      temp = bits_vec[i];
      //str2uint16(temp_str,&bits_vec[i]);
      str2uint16(temp_str,&temp);
      //printf("I AM HERE: str2uint16_wraper: %s -- %d\n",temp_str,temp);
      bits_vec[i] = temp;
   }
}


/*****************
 * main function.  opens up a file, reads in footprints and calculate Tc value
 *****************/
void main(int argc, char *argv[]) {

  char filename1[100];
  char filename_count1[100];
  char filename_smi1[100];
  char filename2[100];
  char filename_smi2[100];
  char filename_count2[100];
  char line[1000];
  FILE *fp1, *fps1, *fpc1;
  FILE *fp2, *fps2, *fpc2;
  int size1 = 0;
  int size2 = 0;
  int size_smi1 = 0;
  int size_smi2 = 0;  
  //char **fingerprints1;
  //char **fingerprints2;
  char **names1;
  char **names2;
  uint16_t *fingerprints_binary1;
  uint16_t *fingerprints_binary2;
  //char *fingerprints_binary_temp;
  //uint16_t *fingerprints_binary_temp1;
  //uint16_t *fingerprints_binary_temp2;
  int *fingerprints_one_count1;
  int *fingerprints_one_count2;
  int namesize = 100; // size of the name in the smiles file.
  int fpsize = 1024;
  int fpbsize = 1024/16; // lenth of bits in binary array when fingerprint is converted to binary. 
  int i,count,count2;
  char prefix[100];
  //double tc_threshold; 
  //int maxclusters; 

   if( argc == 8 ) {
      printf("The argument supplied are:\n (1) %s\n (2) %s \n (3) %s \n (4) %s\n (5) %s\n (6) %s\n", argv[1], argv[2], argv[3],argv[4],argv[5],argv[6]);
   }
   else if( argc > 8 ) {
      printf("Too many arguments supplied.\n");
      exit(1);
   }
   else {
      printf("7 argument expected.\n");
      printf("(1) first fingerprint file.\n");
      printf("(2) first smiles file.\n");
      printf("(3) first count file.\n");
      printf("(4) second fingerprint file.\n");
      printf("(5) second smiles file.\n");
      printf("(6) second count file.\n");
      printf("(7) prefix of the outputfile.\n");
      exit(1);
   }

  //filename = "tc0p6.fp";
  //strcpy(filename, "tc0p6.fp");
  strcpy(filename1, argv[1]);
  strcpy(filename_smi1, argv[2]);
  strcpy(filename_count1, argv[3]);
  strcpy(filename2, argv[4]);
  strcpy(filename_smi2, argv[5]);
  strcpy(filename_count2, argv[6]);
  strcpy(prefix, argv[7]);
  printf("figureprint file1 = %s\n",filename1);
  printf("smiles file1 = %s\n",filename_smi1);
  printf("count file1 = %s\n",filename_count1);
  printf("figureprint file2 = %s\n",filename2);
  printf("smiles file2 = %s\n",filename_smi2);
  printf("count file2 = %s\n",filename_count2);
  printf("prefix = %s\n",prefix);

  fps1=fopen(filename_smi1, "r");
  if (!fps1){
     printf("%s cannot be opened.\n",filename_smi1);
     exit(1);
  } 
  size_smi1 = line_count(fps1); // get the number of fingerprints
  rewind(fps1); // go back to the begining

  fps2=fopen(filename_smi2, "r");
  if (!fps2){
     printf("%s cannot be opened.\n",filename_smi2);
     exit(1);
  } 
  size_smi2 = line_count(fps2); // get the number of fingerprints
  rewind(fps2); // go back to the begining

  fp1=fopen(filename1, "r");
  if (!fp1){
     printf("%s cannot be opened.\n",filename1);
     exit(1);
  } 
  //fp=fopen("tc0p6.fp","r");
  //size1 = line_count(fp1); // get the number of fingerprints
  //rewind(fp1); // go back to the begining
  //if (size1 != size_smi1){
      //printf("files1 are not the same lenth.  exiting . . . \n");
      //exit(1);
  //} 

  fp2=fopen(filename2, "r");
  if (!fp2){
     printf("%s cannot be opened.\n",filename2);
     exit(1);
  } 
  //fp=fopen("tc0p6.fp","r");
  //size2 = line_count(fp2); // get the number of fingerprints
  //rewind(fp2); // go back to the begining
  //if (size2 != size_smi2){
      //printf("files2 are not the same lenth.  exiting . . . \n");
      //exit(1);
  //} 

  printf("number of lines in smi file1: %d\n",size_smi1);
  printf("number of lines in smi file2: %d\n",size_smi2);

  // alocate memory. 
  //fingerprints = malloc(sizeof(char) * fpsize * size);
  //char **fingerprints = (char **)malloc(size * sizeof(char *));
  //fingerprints1 = (char **)malloc(size1 * sizeof(char *));
  names1 = (char **)malloc(size_smi1 * sizeof(char *));
  for (i=0; i<size_smi1; i++){
       //fingerprints1[i] = (char *)malloc(fpsize * sizeof(char));
       names1[i] = (char *)malloc(namesize * sizeof(char));
  }
  // alocate memory. 
  //fingerprints = malloc(sizeof(char) * fpsize * size);
  //char **fingerprints = (char **)malloc(size * sizeof(char *));
  //fingerprints2 = (char **)malloc(size2 * sizeof(char *));
  names2 = (char **)malloc(size_smi2 * sizeof(char *));
  for (i=0; i<size_smi2; i++){
       //fingerprints2[i] = (char *)malloc(fpsize * sizeof(char));
       names2[i] = (char *)malloc(namesize * sizeof(char));
  }

  count = 0;
  // read in the names from file. 
  while(fgets(line, 2000, (FILE*) fps1)) { // 
      process_smiles_str(line, 2000, names1[count], namesize);
      //printf("%s\n",names1[count]);
      count++;
  }
  fclose(fps1);
 
  count = 0;
  // read in the names from file. 
  while(fgets(line, 2000, (FILE*) fps2)) { // 
      process_smiles_str(line, 2000, names2[count], namesize);
      //printf("%s\n",names2[count]);
      count++;
  }
  fclose(fps2);

  fingerprints_binary1 = (uint16_t *)malloc(fpbsize*size_smi1 * sizeof(uint16_t));
  fingerprints_binary2 = (uint16_t *)malloc(fpbsize*size_smi2 * sizeof(uint16_t));

  init_BitCountArray();
  //while ((line = fgetline(fp,2000)) != EOF){
  //      printf("%s\n",line);
  //} 
  count = 0;
  count2 = 0;// this is for fingerprints_binary
  fingerprints_one_count1 = (int *)malloc(size_smi1 * sizeof(int *));
  while(fgets(line, 2000, (FILE*) fp1)) { // read in fingerprints from file
     //printf("%s\n",line);
     fingerprints_binary1[count] = (uint16_t) atoi(line);
     //process_str2uint16(line,2000,fingerprints_binary1[count],fpsize);
     //printf("%d\n",fingerprints_binary1[count]);
     count++;
  }
  fclose(fp1);
  fpc1 = fopen(filename_count1,"r");
  if (!fpc1){
     printf("%s cannot be opened.\n",filename_count1);
     exit(1);
  } 
  while(fgets(line, 2000, (FILE*) fpc1)) { // read in fingerprints from file
     fingerprints_one_count1[count2] = atoi(line);
     //process_str2int(line,2000,fingerprints_one_count1[count2],fpsize);
     //printf("%d\n",fingerprints_one_count1[count2]);
     count2++;
  }
  fclose(fpc1);


  count = 0;
  count2 = 0;// this is for fingerprints_binary
  fingerprints_one_count2 = (int *)malloc(size_smi2 * sizeof(int *));
  while(fgets(line, 2000, (FILE*) fp2)) { // read in fingerprints from file
     fingerprints_binary2[count] = (uint16_t) atoi(line);
     //process_str2uint16(line,2000,fingerprints_binary2[count],fpsize);
     count++;
  }
  fclose(fp2);
  fpc2 = fopen(filename_count2,"r");
  if (!fpc2){
     printf("%s cannot be opened.\n",filename_count2);
     exit(1);
  } 
  while(fgets(line, 2000, (FILE*) fpc2)) { // read in fingerprints from file
     fingerprints_one_count2[count2] = atoi(line);
     //process_str2int(line,2000,fingerprints_one_count2[count2],fpsize);
     count2++;
  }
  fclose(fpc2);
  

  cal_Tc_matrix(fingerprints_binary1, fingerprints_one_count1, size_smi1, names1, fingerprints_binary2, fingerprints_one_count2, size_smi2, names2, fpbsize, prefix);
   
  // clean up arrays. 
  for (i=0; i<size1; i++){
       //free(fingerprints1[i]);
       free(names1[i]);
  }
  //free(fingerprints1);
  free(names1);
  free(fingerprints_binary1);
  //free(fingerprints_binary_temp1);
  free(fingerprints_one_count1);

  for (i=0; i<size2; i++){
       //free(fingerprints2[i]);
       free(names2[i]);
  }
  //free(fingerprints2);
  free(names2);
  free(fingerprints_binary2);
  //free(fingerprints_binary_temp2);
  free(fingerprints_one_count2);
}

//main();

