/* fast_tanimoto.h: Performs fast tanimoto coefficient calculations in c
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
 * Mysinger 6/07 Use <stdint.h> and switch float to double
 */

#ifndef FAST_TANIMOTO
#define FAST_TANIMOTO

#include <stdint.h>

double tanimoto_coeff_with_sum(uint16_t *fp1, uint16_t *fp2, int db_len);

void init_BitCountArray(void);
int bitcount(uint16_t *fp, int db_len);
double tanimoto_coeff(uint16_t *fp1, int sum1, uint16_t *fp2, int sum2, int db_len);

char* ascii2binary(char *ascii, int *binary_length);
char* binary2ascii(char *binary, int binary_length);
void a2b(char *a4, char *b3);
void b2a(unsigned char *b3, unsigned char *a4);

char* string2binary(char *binstring, int *binary_length);
char* binary2string(char *binary, int binary_length);

#endif

