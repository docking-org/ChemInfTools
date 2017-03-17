/* fast_tanimoto.c: Performs fast tanimoto coefficient calculations in c
 *
 * Can be used in place of dt_fp_tanimoto on raw binary fingerprints
 * To use other fingerprint formats, use the provided routines to convert
 * to and from raw binary format.
 * 
 * Fast tanimoto coefficient code uses SUBSET bitvector.c code as a template 
 * SUBSET is found at http://cactus.nci.nih.gov/SUBSET/ 
 * Original code written by Dr. Bruno Bienfait
 * Daylight ASCII conversion code based on contrib c files, see below
 *
 * Revisions:
 * Mysinger 7/05 Created
 * Mysinger 7/05 Added conversion routines
 * Mysinger 5/07 Switched float to double
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fast_tanimoto.h"

#define MAXBUF 2048

/* Static array lookup table for on bitcounts of all possible 16 bit values */
/* Uses 64K memory */
static uint8_t BitCountArray[UINT16_MAX + 1];

/* Static flag to indicate when BitCountArray has been initialized */
static int INIT_FLAG = 0;

/* Initialize BitCountArray lookup table
 * There are UINT16_MAX + 1 different 16 bit unsigned integers.
 * For each integer, we keep track of the count of bits on.
 * Must be called before you use the bitcount or tanimoto_coeff functions
 */
void init_BitCountArray(void)
{
  int j, sum;
  uint32_t i = 0;
  /* Not a uint16_t value because that results in an infinite loop */
  /* due to last loop interation wrapping to 0 */

  for(i=0; i <= UINT16_MAX; i++) {
    sum = 0;
    for(j=0; j < 16; j++)
      sum += (((uint16_t) i & ( 1 << j)) != 0);
    BitCountArray[i] = sum;
  }

  INIT_FLAG = 1;
}

/* Slightly slower version of tanimoto_coeff. The advantages of this function 
 * are that you need not worry to initialize the BitCountArray nor precompute 
 * the bitcounts of the raw binary fingerprints.
 *
 * Calculates the tanimoto coefficent of fp1 with fp2
 * tanimoto coefficient = (number of AND bits)/(number of OR bits)
 * As arguments, takes a pointer to fp1, a pointer to fp2, and the length 
 * of the fingerprints in double bytes
 */
double tanimoto_coeff_with_sum(uint16_t *fp1, uint16_t *fp2, int db_len)
{
  int i;
  int common=0;
  int sum1, sum2;
  int db1, db2;

  /* Initialize the BitCountArray if it has not already been done */
  if(INIT_FLAG == 0)
    init_BitCountArray();

  sum1 = 0;
  sum2 = 0;

  for(i=0; i < db_len; i++) {
    db1 = fp1[i];
    db2 = fp2[i];

    /* Sum up the number of bits on in each incoming fingerprint */
    sum1 += BitCountArray[db1];
    sum2 += BitCountArray[db2];

    /* Count of common bits equals the count of bits in (db1 & db2) */
    common += BitCountArray[db1 & db2];
  }

  /* Prevent a division by zero */
  if ((sum1 == 0) && (sum2 == 0))
    return 1.0;

  return (((double) common)/(sum1 + sum2 - common));
}

/* WARNING: init_BitCountArray() must be called before using this function
 *
 * Calculates the number of bits on in a raw binary fingerprint.
 * As arguments, it takes a pointer to the first double byte of the 
 * fingerprint and the length of the fingerprint in double bytes
 */
int bitcount(uint16_t *fp, int db_len)
{
  int i, count = 0;
  uint16_t db;

  for(i=0; i < db_len; i++) {
    db = fp[i];
    count += BitCountArray[db];
  }

  return count;
}

/* WARNING: init_BitCountArray() must be called before using this function
 * 
 * Calculates the tanimoto coefficent of fp1 with fp2
 * tanimoto coefficient = (number of AND bits)/(number of OR bits)
 * As arguments, takes a pointer to fp1, the bitcount of fp1,
 * a pointer to fp2, the bitcount of fp2, and the length of the fingerprints 
 * in double bytes
 */
double tanimoto_coeff(uint16_t *fp1, int sum1, uint16_t *fp2, int sum2, int db_len)
{
  int i;
  int common=0;
  int db1, db2;

  /* Prevent a division by zero */
  if ((sum1 == 0) && (sum2 == 0))
    return 1.0;

  for(i=0; i < db_len; i++) {
    db1 = fp1[i];
    db2 = fp2[i];

    /* Count of common bits equals the count of bits in (db1 & db2) */
    common += BitCountArray[db1 & db2];
  }

  return (((double) common)/(sum1 + sum2 - common));
}

/* End of fast tanimoto code */
/***************************************************************************/
/* Start of format conversion code */

/* Convert from Daylight ASCII fingerprint format to raw binary fp format.
 *
 * ascii contains the input Daylight ASCII fingerprint
 * binary_length will be set to the length of the returned raw binary string
 * return value is the raw binary string
 *
 * $DY_ROOT/contrib/src/c/fingerprint/bits2ascii.c used as a guide
 * KEISER 06/05 Original conversion
 * MYSINGER 07/05 Extensive rewrite for library use
 */
char* ascii2binary(char *ascii, int *binary_length) 
{
  char binarray[MAXBUF];         /*** Binary array                 ***/
  char *q;                       /*** Ptr into ascii string        ***/
  char *p;                       /*** Ptr into binarray            ***/
  char c;                        /*** Terminal char of ascii       ***/
  int i;                         /*** Generic counter              ***/
  int len;                       /*** Length of ascii string       ***/
  int nbits = 0;                 /*** Number of bits in binarray   ***/
  char *binary;                  /*** Return string                ***/

  len = strlen(ascii);
  if ((len < 5) || ((len-1)%4)) {
    printf("Error: Length %d of ascii string must be 4n+1, n > 0!\n", len); 
    return NULL;
  }

  /* Convert ascii to binary in 4 byte swallows */
  p = binarray;
  q = ascii;
  for (i = 0; i < (len/4); i++) {
    a2b(q, p);
    p += 3;
    q += 4;
    nbits += 24;
    if (p >= binarray+MAXBUF) {
      printf("Error: Input ascii string too long for buffer!\n"); 
      return NULL;
    }
  }

  /* Use the terminal ascii char to modify length of binary string */
  c = ascii[len-1];
  if (c < '1' || c > '3') {
    printf("Error: Input ascii string has bad terminal char '%c'!\n", c); 
    return NULL;
  } 
  else {
    nbits -= 8 * ('3' - c);
  }
  
  /* Return the raw binary length */
  *binary_length = (nbits >> 3);

  /* Copy raw binary string to heap
   * Warning: cannot trust strlen/strncpy because of null chars
   */
  binary = (char *) malloc(*binary_length);
  memcpy(binary, binarray, *binary_length);

  return binary;
}

/* Convert from raw binary fingerprint format to Daylight ASCII format.
 *
 * binary contains the input raw binary string
 * binary_length contains the length of the raw binary string
 * return value is the Daylight ASCII string
 *
 * $DY_ROOT/contrib/src/c/fingerprint/bits2ascii.c used as a guide
 * MYSINGER 07/05 Original conversion for library use
 */
char* binary2ascii(char *binary, int binary_length)
{
  unsigned char ascarray[MAXBUF]; /*** ASCII array                  ***/
  unsigned char *p;               /*** Ptr into ascarray            ***/
  unsigned char *q;               /*** Ptr into ascarray            ***/
  unsigned char b3[3];            /*** 3 binary bytes (1 chunk)     ***/
  int i;                          /*** Generic counter              ***/
  int mod3;                       /*** binary_length % 3            ***/
  int len = 0;                    /*** Current ascii length         ***/
  char *ascii;                    /*** Return string                ***/

  if (binary_length <= 0) {
    printf("Error: Binary length %d is <= 0!\n", binary_length); 
    return NULL;
  }

  memset(ascarray, (unsigned char) 0x00, MAXBUF);

  p = ascarray;
  q = (unsigned char*) binary;
  for (i = 0; i < (binary_length/3); i++) {
    b2a(q, p);
    p += 4;
    q += 3;
    len += 4;
    if (p >= ascarray+MAXBUF) {
      printf("Error: Input binary string too long for buffer\n"); 
      return NULL;
    }
  }

  mod3 = (binary_length%3);
  if (mod3) {
    memset(b3, (unsigned char) 0x00, 3);
    for (i = 0; i < mod3; i++)
      b3[i] = *q++;
    b2a(b3, p);
    p += 4;
    len += 4;
  }
  if      (mod3 == 0) *p = '3';
  else if (mod3 == 1) *p = '1';
  else if (mod3 == 2) *p = '2';
  ++p;
  ++len;

  ascii = (char *) malloc((len+1) * sizeof(char));
  strcpy(ascii, (char *) ascarray);

  return ascii;
}

/* Convert from 0 or 1 binary string format to raw binary fp format.
 *
 * binstring contains the input 0 or 1 binary string
 * binary_length will be set to the length of the returned raw binary string
 * return value is the raw binary string
 *
 * MYSINGER 07/05 Created for library use
 */
char* string2binary(char *binstring, int *binary_length) {
  unsigned char binarray[MAXBUF];  /*** Binary array                 ***/
  char c;                          /*** Input char                   ***/
  int blen;                        /*** Length of binstring          ***/
  int i;                           /*** Generic counter              ***/
  int len = 0;                     /*** Length of raw binary         ***/
  char *binary;                    /*** Return string                ***/

  blen = strlen(binstring);
  if (blen == 0) {
    printf("Error: Input string is NULL!\n"); 
    return NULL;
  }
  memset(binarray, (unsigned char) 0x00, MAXBUF);
  
  /* Walk down binstring ORing bits into binarray */
  for (i = 0; i < blen; i++) {
    c = binstring[i];
    if (c == '1') 
      binarray[len] |= (1 << (7 - (i%8)));    
    if (!((i+1)%8))
      ++len;
  }

  /* Adjust length for extra bits */
  if (blen%8)
    ++len;

  /* Return the raw binary length */
  *binary_length = len;

  /* Copy raw binary string to heap
   * Warning: cannot trust strlen/strncpy because of null chars
   */
  binary = (char *) malloc(*binary_length);
  memcpy(binary, binarray, *binary_length);

  return binary;
}

/* Convert from raw binary fp format to 0 or 1 binary string format
 *
 * binary contains the input raw binary string
 * binary_length contains the length of the raw binary string
 * return value is the 0 or 1 binary string
 *
 * MYSINGER 07/05 Created for library use
 */ 
char* binary2string(char *binary, int binary_length) {
  unsigned char n;          /*** Input char                   ***/
  char c;                   /*** Binary char                  ***/
  char strarray[8*MAXBUF];  /*** String array                 ***/
  int i;                    /*** Generic counter              ***/
  int j;                    /*** Count bits inside byte       ***/
  int len = 0;              /*** Current string length        ***/
  char *binstring;          /*** Return string                ***/

  if (binary_length <= 0) {
    printf("Error: Binary length %d is <= 0!\n", binary_length); 
    return NULL;
  }

  for (i = 0; i < binary_length; i++) {
    n = binary[i];
    for (j = 7; j >= 0; j--) {
      c = '0' + ((n >> j) & 1);
      strarray[len++] = c;
      if (len >= 8*MAXBUF) {
	printf("Error: Binary input too long for buffer\n"); 
	return NULL;
      }
    }
  }
  strarray[len++] = '\0';

  binstring = (char *) malloc(len * sizeof(char));
  strcpy(binstring, strarray);

  return binstring;
}

/* Convert from 0 or 1 binary string format to Daylight ASCII format
 *
 * binstring contains the input 0 or 1 binary string
 * return value is the Daylight ASCII string
 *
 * MYSINGER 06/07 Added
 */
char* string2ascii(char *binstring) {
  int blen;
  char *binary;
  char *ascii;

  binary = string2binary(binstring, &blen);
  ascii = binary2ascii(binary, blen);
  free(binary);

  return ascii;
}

/* Convert from Daylight ASCII format to 0 or 1 binary string format
 *
 * binstring contains the input Daylight ASCII string
 * return value is the 0 or 1 binary string
 *
 * MYSINGER 06/07 Added
 */
char* ascii2string(char *ascii) {
  int blen;
  char *binary;
  char *binstring;

  binary = ascii2binary(ascii, &blen);
  binstring = binary2string(binary, blen);
  free(binary);

  return binstring;
}

/* Code verbatim from $DY_ROOT/contrib/src/c/fingerprint/ascii2bits.c
 * Helper function for ascii2binary()
 */
void a2b(char *a4, char *b3)
{
  int i;
  char b;
  char byte = 0x00;

  /*********************************************
  *** Use the Daylight mapping to convert each
  *** ascii char to its 6-bit code.
  ***
  *** a4: |xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx (printable)
  ***     |=======+=======+=======+=======|
  ***   becomes...
  *** a4: |00xxxxxx00xxxxxx00xxxxxx00xxxxxx
  ***     |=======+=======+=======+=======|
  *********************************************/
  for (i = 0; i < 4; ++i) {
    switch (a4[i]) {
      case '.': byte = 0x00; break;      /* 00 = __000000 */
      case '+': byte = 0x01; break;      /* 01 = __000001 */
      case '0': byte = 0x02; break;      /* 02 = __000010 */
      case '1': byte = 0x03; break;      /* 03 = __000011 */
      case '2': byte = 0x04; break;      /* 04 = __000100 */
      case '3': byte = 0x05; break;      /* 05 = __000101 */
      case '4': byte = 0x06; break;      /* 06 = __000110 */
      case '5': byte = 0x07; break;      /* 07 = __000111 */
      case '6': byte = 0x08; break;      /* 08 = __001000 */
      case '7': byte = 0x09; break;      /* 09 = __001001 */
      case '8': byte = 0x0a; break;      /* 10 = __001010 */
      case '9': byte = 0x0b; break;      /* 11 = __001011 */
      case 'A': byte = 0x0c; break;      /* 12 = __001100 */
      case 'B': byte = 0x0d; break;      /* 13 = __001101 */
      case 'C': byte = 0x0e; break;      /* 14 = __001110 */
      case 'D': byte = 0x0f; break;      /* 15 = __001111 */
      case 'E': byte = 0x10; break;      /* 16 = __010000 */
      case 'F': byte = 0x11; break;      /* 17 = __010001 */
      case 'G': byte = 0x12; break;      /* 18 = __010010 */
      case 'H': byte = 0x13; break;      /* 19 = __010011 */
      case 'I': byte = 0x14; break;      /* 20 = __010100 */
      case 'J': byte = 0x15; break;      /* 21 = __010101 */
      case 'K': byte = 0x16; break;      /* 22 = __010110 */
      case 'L': byte = 0x17; break;      /* 23 = __010111 */
      case 'M': byte = 0x18; break;      /* 24 = __011000 */
      case 'N': byte = 0x19; break;      /* 25 = __011001 */
      case 'O': byte = 0x1a; break;      /* 26 = __011010 */
      case 'P': byte = 0x1b; break;      /* 27 = __011011 */
      case 'Q': byte = 0x1c; break;      /* 28 = __011100 */
      case 'R': byte = 0x1d; break;      /* 29 = __011101 */
      case 'S': byte = 0x1e; break;      /* 30 = __011110 */
      case 'T': byte = 0x1f; break;      /* 31 = __011111 */
      case 'U': byte = 0x20; break;      /* 32 = __100000 */
      case 'V': byte = 0x21; break;      /* 33 = __100001 */
      case 'W': byte = 0x22; break;      /* 34 = __100010 */
      case 'X': byte = 0x23; break;      /* 35 = __100011 */
      case 'Y': byte = 0x24; break;      /* 36 = __100100 */
      case 'Z': byte = 0x25; break;      /* 37 = __100101 */
      case 'a': byte = 0x26; break;      /* 38 = __100110 */
      case 'b': byte = 0x27; break;      /* 39 = __100111 */
      case 'c': byte = 0x28; break;      /* 40 = __101000 */
      case 'd': byte = 0x29; break;      /* 41 = __101001 */
      case 'e': byte = 0x2a; break;      /* 42 = __101010 */
      case 'f': byte = 0x2b; break;      /* 43 = __101011 */
      case 'g': byte = 0x2c; break;      /* 44 = __101100 */
      case 'h': byte = 0x2d; break;      /* 45 = __101101 */
      case 'i': byte = 0x2e; break;      /* 46 = __101110 */
      case 'j': byte = 0x2f; break;      /* 47 = __101111 */
      case 'k': byte = 0x30; break;      /* 48 = __110000 */
      case 'l': byte = 0x31; break;      /* 49 = __110001 */
      case 'm': byte = 0x32; break;      /* 50 = __110010 */
      case 'n': byte = 0x33; break;      /* 51 = __110011 */
      case 'o': byte = 0x34; break;      /* 52 = __110100 */
      case 'p': byte = 0x35; break;      /* 53 = __110101 */
      case 'q': byte = 0x36; break;      /* 54 = __110110 */
      case 'r': byte = 0x37; break;      /* 55 = __110111 */
      case 's': byte = 0x38; break;      /* 56 = __111000 */
      case 't': byte = 0x39; break;      /* 57 = __111001 */
      case 'u': byte = 0x3a; break;      /* 58 = __111010 */
      case 'v': byte = 0x3b; break;      /* 59 = __111011 */
      case 'w': byte = 0x3c; break;      /* 60 = __111100 */
      case 'x': byte = 0x3d; break;      /* 61 = __111101 */
      case 'y': byte = 0x3e; break;      /* 62 = __111110 */
      case 'z': byte = 0x3f; break;      /* 63 = __111111 */
    }

    /*********************************************
    *** Now copy the 4x6=24 bits from a4 to b3. 
    ***
    *** a4: |--000000--111111--222222--333333
    ***     |=======+=======+=======+=======|
    ***
    *** b3: |000000111111222222333333
    ***     |=====+=====+=====+=====|
    *********************************************/
    if (i == 0)
      b3[0] = (byte << 2);             /*** 6 bits into 1st byte ***/
    else if (i == 1) {
      b3[0] |= ((b = byte) >> 4);      /*** 2 bits into 1st byte ***/
      b3[1] = ((b = byte) << 4);       /*** 4 bits into 2nd byte ***/
    } else if (i == 2) {
      b3[1] |= ((b = byte) >> 2);      /*** 4 bits into 2nd byte ***/
      b3[2] = ((b = byte) << 6);       /*** 2 bits into 3rd byte ***/
    } else if (i == 3)
      b3[2] |= byte;                   /*** 6 bits into 3rd byte ***/
  }
  return;
}

/* Code verbatim from $DY_ROOT/contrib/src/c/fingerprint/bits2ascii.c
 * Helper function for binary2ascii()
 */
void b2a(unsigned char *b3, unsigned char *a4)
{
  int i;

  /*********************************************
  *** Copy the 4 "6-bit bytes" from b3 to the
  *** 4 bytes of a4, right justified.
  ***
  *** b3: |000000111111222222333333            
  ***     |=====+=====+=====+=====| 
  ***
  *** a4: |--000000--111111--222222--333333    
  ***     |=======+=======+=======+=======| 
  *********************************************/
  a4[0] = b3[0] >> 2;
                   
  a4[1] = ((b3[0] & 0x03) << 4) | ((b3[1] & 0xf0) >> 4);

  a4[2] = ((b3[1] & 0x0f) << 2) | ((b3[2] & 0xc0) >> 6);

  a4[3] = b3[2] & 0x3f;

  /*********************************************
  *** Now use the Daylight mapping to convert
  *** each byte to its printable-ascii code.
  ***
  *** a4: |00xxxxxx00xxxxxx00xxxxxx00xxxxxx    
  ***     |=======+=======+=======+=======| 
  ***   becomes...
  *** a4: |xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx (printable)
  ***     |=======+=======+=======+=======| 
  *********************************************/
  for (i = 0; i < 4; ++i) {
    switch ((int) a4[i]) {
      case  0x00: a4[i] = '.'; break;      /* 00 = __000000 */
      case  0x01: a4[i] = '+'; break;      /* 01 = __000001 */
      case  0x02: a4[i] = '0'; break;      /* 02 = __000010 */
      case  0x03: a4[i] = '1'; break;      /* 03 = __000011 */
      case  0x04: a4[i] = '2'; break;      /* 04 = __000100 */
      case  0x05: a4[i] = '3'; break;      /* 05 = __000101 */
      case  0x06: a4[i] = '4'; break;      /* 06 = __000110 */
      case  0x07: a4[i] = '5'; break;      /* 07 = __000111 */
      case  0x08: a4[i] = '6'; break;      /* 08 = __001000 */
      case  0x09: a4[i] = '7'; break;      /* 09 = __001001 */
      case  0x0a: a4[i] = '8'; break;      /* 10 = __001010 */
      case  0x0b: a4[i] = '9'; break;      /* 11 = __001011 */
      case  0x0c: a4[i] = 'A'; break;      /* 12 = __001100 */
      case  0x0d: a4[i] = 'B'; break;      /* 13 = __001101 */
      case  0x0e: a4[i] = 'C'; break;      /* 14 = __001110 */
      case  0x0f: a4[i] = 'D'; break;      /* 15 = __001111 */
      case  0x10: a4[i] = 'E'; break;      /* 16 = __010000 */
      case  0x11: a4[i] = 'F'; break;      /* 17 = __010001 */
      case  0x12: a4[i] = 'G'; break;      /* 18 = __010010 */
      case  0x13: a4[i] = 'H'; break;      /* 19 = __010011 */
      case  0x14: a4[i] = 'I'; break;      /* 20 = __010100 */
      case  0x15: a4[i] = 'J'; break;      /* 21 = __010101 */
      case  0x16: a4[i] = 'K'; break;      /* 22 = __010110 */
      case  0x17: a4[i] = 'L'; break;      /* 23 = __010111 */
      case  0x18: a4[i] = 'M'; break;      /* 24 = __011000 */
      case  0x19: a4[i] = 'N'; break;      /* 25 = __011001 */
      case  0x1a: a4[i] = 'O'; break;      /* 26 = __011010 */
      case  0x1b: a4[i] = 'P'; break;      /* 27 = __011011 */
      case  0x1c: a4[i] = 'Q'; break;      /* 28 = __011100 */
      case  0x1d: a4[i] = 'R'; break;      /* 29 = __011101 */
      case  0x1e: a4[i] = 'S'; break;      /* 30 = __011110 */
      case  0x1f: a4[i] = 'T'; break;      /* 31 = __011111 */
      case  0x20: a4[i] = 'U'; break;      /* 32 = __100000 */
      case  0x21: a4[i] = 'V'; break;      /* 33 = __100001 */
      case  0x22: a4[i] = 'W'; break;      /* 34 = __100010 */
      case  0x23: a4[i] = 'X'; break;      /* 35 = __100011 */
      case  0x24: a4[i] = 'Y'; break;      /* 36 = __100100 */
      case  0x25: a4[i] = 'Z'; break;      /* 37 = __100101 */
      case  0x26: a4[i] = 'a'; break;      /* 38 = __100110 */
      case  0x27: a4[i] = 'b'; break;      /* 39 = __100111 */
      case  0x28: a4[i] = 'c'; break;      /* 40 = __101000 */
      case  0x29: a4[i] = 'd'; break;      /* 41 = __101001 */
      case  0x2a: a4[i] = 'e'; break;      /* 42 = __101010 */
      case  0x2b: a4[i] = 'f'; break;      /* 43 = __101011 */
      case  0x2c: a4[i] = 'g'; break;      /* 44 = __101100 */
      case  0x2d: a4[i] = 'h'; break;      /* 45 = __101101 */
      case  0x2e: a4[i] = 'i'; break;      /* 46 = __101110 */
      case  0x2f: a4[i] = 'j'; break;      /* 47 = __101111 */
      case  0x30: a4[i] = 'k'; break;      /* 48 = __110000 */
      case  0x31: a4[i] = 'l'; break;      /* 49 = __110001 */
      case  0x32: a4[i] = 'm'; break;      /* 50 = __110010 */
      case  0x33: a4[i] = 'n'; break;      /* 51 = __110011 */
      case  0x34: a4[i] = 'o'; break;      /* 52 = __110100 */
      case  0x35: a4[i] = 'p'; break;      /* 53 = __110101 */
      case  0x36: a4[i] = 'q'; break;      /* 54 = __110110 */
      case  0x37: a4[i] = 'r'; break;      /* 55 = __110111 */
      case  0x38: a4[i] = 's'; break;      /* 56 = __111000 */
      case  0x39: a4[i] = 't'; break;      /* 57 = __111001 */
      case  0x3a: a4[i] = 'u'; break;      /* 58 = __111010 */
      case  0x3b: a4[i] = 'v'; break;      /* 59 = __111011 */
      case  0x3c: a4[i] = 'w'; break;      /* 60 = __111100 */
      case  0x3d: a4[i] = 'x'; break;      /* 61 = __111101 */
      case  0x3e: a4[i] = 'y'; break;      /* 62 = __111110 */
      case  0x3f: a4[i] = 'z'; break;      /* 63 = __111111 */
      default: fprintf(stderr, "Aaack!\n");
    }
  }
  return;
}

/* End of format conversion code */
/***************************************************************************/
