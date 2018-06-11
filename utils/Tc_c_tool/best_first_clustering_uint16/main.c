

#include "best_first.h"
#include "../common/fast_tanimoto.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

// read in file, read in all footprints
// call clustering routine. 
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
  printf("count ones: %d\n",countones);
  //return outputline;
  return ;
}

/////////////////////

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
 * main function.  opens up a file, reads in footprints and clusters them
 *****************/
void main(int argc, char *argv[]) {

  char filename[100];
  char filename_count[100];
  char filename_smi[100];
  char line[2000];
  FILE *fp, *fps, *fpc;
  int size = 0;  
  int size_smi = 0;  
  //char **fingerprints;
  char **names;
  uint16_t *fingerprints_binary;
  //char *fingerprints_binary_temp;
  //uint16_t *fingerprints_binary_temp2;
  int *fingerprints_one_count;
  int namesize = 100; // size of the name in the smiles file.
  int fpsize = 1024;
  int fpbsize = 1024/16; // lenth of bits in binary array when fingerprint is converted to binary. 
  int i,count,count2; 
  double tc_threshold; 
  int maxclusters; 

   if( argc == 6 ) {
      printf("The argument supplied are:\n (1) %s\n (2) %s \n (3) %s \n (4) %s\n (5) %s\n", argv[1], argv[2], argv[3], argv[4], argv[5]);
   }
   else if( argc > 6 ) {
      printf("Too many arguments supplied.\n");
      exit(1);
   }
   else {
      printf("five arguments expected.\n");
      printf("(1) fingerprint file.\n");
      printf("(2) cout file.\n");
      printf("(3) smiles file.\n");
      printf("(4) tanimoto coefficient threshold to define clustering. values between 0 and 1. smaller values will result in larger clusters. recommend value of 0.5\n");
      printf("(5) max number of clusters,  the program will terminate once this number is reached.\n  If you want to process the whole list than make this number grater than the number of molcules.\n");
      exit(1);
   }

  //filename = "tc0p6.fp";
  //strcpy(filename, "tc0p6.fp");
  strcpy(filename, argv[1]);
  strcpy(filename_count, argv[2]);
  strcpy(filename_smi, argv[3]);
  sscanf(argv[4], "%lf", &tc_threshold);
  sscanf(argv[5], "%d", &maxclusters);
  printf("figureprint file = %s\n",filename);
  printf("count file = %s\n",filename_count);
  printf("smiles file = %s\n",filename_smi);
  printf("tc threshold = %f\n",tc_threshold);
  printf("maxium clusters to be generated = %d\n",maxclusters);

  fps=fopen(filename_smi, "r");
  if (!fps){
     printf("%s cannot be opened.\n",filename_smi);
     exit(1);
  } 
  size_smi = line_count(fps); // get the number of fingerprints
  rewind(fps); // go back to the begining


  fp=fopen(filename, "r");
  //fp=fopen("tc0p6.fp","r");
  if (!fp){
     printf("%s cannot be opened.\n",filename);
     exit(1);
  } 
  size = line_count(fp); // get the number of fingerprints
  rewind(fp); // go back to the begining
  //if (size != size_smi){
      //printf("files are not the same lenth.  exiting . . . \n");
      //exit(1);
  //} 


  printf("number of lines in file: %d\n",size);

  // alocate memory. 
  //fingerprints = malloc(sizeof(char) * fpsize * size);
  //char **fingerprints = (char **)malloc(size * sizeof(char *));
  //fingerprints = (char **)malloc(size * sizeof(char *));
  names = (char **)malloc(size_smi * sizeof(char *));
  for (i=0; i<size_smi; i++){
       //fingerprints[i] = (char *)malloc(fpsize * sizeof(char));
       names[i] = (char *)malloc(namesize * sizeof(char));
  }
  count = 0;
  // read in the names from file. 
  while(fgets(line, 2000, (FILE*) fps)) { // 
      process_smiles_str(line, 2000, names[count], namesize);
      //printf("%s\n",names[count]);
      count++;
  }
  fclose(fps);
 
  fingerprints_binary = (uint16_t *)malloc(fpbsize*size_smi * sizeof(uint16_t));

  init_BitCountArray();
  //while ((line = fgetline(fp,2000)) != EOF){
  //      printf("%s\n",line);
  //} 
  count = 0;
  count2 = 0;// this is for fingerprints_binary
  //fingerprints_binary_temp = (char *)malloc(fpbsize * sizeof(char *));
  //fingerprints_binary_temp2 = (uint16_t *)malloc(fpbsize * sizeof(uint16_t *));
  fingerprints_one_count = (int *)malloc(size_smi * sizeof(int *));
  while(fgets(line, 2000, (FILE*) fp)) { // read in fingerprints from file
     //printf("%s\n", line);
     //process_str(line,2000,fingerprints[count],fpsize);
     //printf("%s\n", fingerprints[count]);
     //fingerprints_binary_temp = string2binary(fingerprints[count], fpbsize); 
     //fingerprints_binary_temp = binary2ascii(fingerprints[count], fpbsize); 
     //    fingerprints_binary_temp2[i] = (uint16_t) fingerprints_binary_temp[i];
     //str2uint16_wraper(fingerprints[count],fpsize, fingerprints_binary_temp2);
     //printf("\n");
     fingerprints_binary[count] = (uint16_t) atoi(line);
     //for(i=0;i<fpbsize;i++){
         //printf("%d-",fingerprints_binary_temp2[i]);
         //fingerprints_binary[count2] = fingerprints_binary_temp2[i];
         //count2++;
     //}
     //printf("\n");
     //init_BitCountArray();
     //fingerprints_one_count[count] = bitcount(fingerprints_binary_temp2, fpbsize);
     //printf("\n%d\n\n%s\n",fingerprints_one_count[count],fingerprints_binary_temp);
     //printf("\n%d\n",fingerprints_one_count[count]);
     //printf("%s\n", fingerprints[count]);
     count++;
  }
  
  fclose(fp);
  fpc = fopen(filename_count,"r");
  if (!fpc){
     printf("%s cannot be opened.\n",filename_count);
     exit(1);
  } 
  while(fgets(line, 2000, (FILE*) fpc)) { // read in fingerprints from file
     fingerprints_one_count[count2] = atoi(line);
     //process_str2int(line,2000,fingerprints_one_count1[count2],fpsize);
     //printf("%d\n",fingerprints_one_count1[count2]);
     count2++;
  }
  fclose(fpc);
  best_first(fingerprints_binary, fingerprints_one_count, size_smi, fpbsize, names, tc_threshold, maxclusters);
   
  // clean up arrays. 
  for (i=0; i<size_smi; i++){
       //free(fingerprints[i]);
       free(names[i]);
  }
  //free(fingerprints);
  free(names);
  free(fingerprints_binary);
  //free(fingerprints_binary_temp2);
  free(fingerprints_one_count);
}

//main();

