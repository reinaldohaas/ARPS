/*
*| Name:
*|      GouldToNative - converts Gould format floating point to native format.
*|
*| Interface:
*|      double
*|      GouldToNative(unsigned char *inVal)
*|
*| Input:
*|      inVal     - input value, pointer to 4 byte Gould format floating point
*|
*| Input and Output:
*|      none
*|
*| Output:
*|      none
*|
*| Return values:
*|      input value converted to native floating point format
*|
*| Remarks:
*|      Gould format description can be found in the NOAA Gvar Transmission
*| Format document, section 3.5.4.
*|
*| Categories:
*|      converter
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <netinet/in.h>

double GouldToNative (
  unsigned char *inVal_tmp
)

{
  double nativeVal;
  double temp;
  double dblMant;
  int compVal;
  int sign;
  int exponent;
  int need_swap;
  float mant;
  unsigned char tmpVal[4];
  unsigned char *inVal;
  union u_tag{
        unsigned char uchar[4];
        long l;
  } u;

 /*
  *  Make a local copy, as we might be doing byte swapping.
  */

  inVal = inVal_tmp;

/*
 *   Don't ask me why, but there is a byte order problem in the original
 *   data that must be corrected now!  This is independent of the endianess
 *   (is there a word?) of the processing machine.
 */

  swapword(inVal);

 /*
  *  Are we byte order dependent?  If so, this complicates things.
  */

  need_swap = 0;
  if (htonl(1) != 1) need_swap = 1;

 /*
  * an example conversion:
  *
  * input value (hex): BD DA 4D 07
  *  a) convert to binary: 
  *     1011 1110 1101 1010 0100 1101 0000 0111
  * BZZZ, wrong input!.  The above number is BE DA 4D 07!
  *  b) sign bit is set, so take twos complement:
  *     0100 0001 0010 0101 1011 0010 1111 1001
  *  c) convert this back to hex: 41 25 B2 F9
  *  d) mantissa = 2470649
  *  e) exponent = 65
  *  f) tempVal = mantissa / 16 exp (70 - 65)
  *  g) outputVal = tempVal * sign = -2.3561944
  *
  */

  /* load input into a temporary buffer for munging */
  memcpy(tmpVal, inVal, 4);

  /* determine sign of number; if negative, take 2's complement */
  sign = 1;
  if (tmpVal[0] & 0x80) {
    if (need_swap)
        swapword(inVal);
    sign = -1;
/*   compVal = TwosComp(* ((int *) inVal)); */
    memcpy(u.uchar, inVal, 4);
/*   printf(" long: %X\n",u.l); */
    compVal = ((unsigned long)u.l ^ (unsigned long)0xffffffff) +1;
/*  printf(" compVal: %X\n",compVal);  */
    memcpy(tmpVal, &compVal, 4);
    if (need_swap)
        swapword(tmpVal);
  }

  /* determine the exponent, biased at 64 */
  exponent = (tmpVal[0] & 0x7f);
  if (exponent == 0)
    exponent = 64;

  /* determine the value of the mantissa, load into a double */
  mant = (float) (tmpVal[3] + (tmpVal[2] * 256) + (tmpVal[1] * 65536));
  dblMant = (double) mant;

  /* now adjust the mantissa according to the exponent */
  temp = pow((double) 16, (double) (70 - exponent));
  nativeVal = dblMant / temp;
  nativeVal = nativeVal * sign;

/*printf(" sign:%d output:%lf\n",sign,nativeVal);*/
  /* return value as a double, caller can cast if needed */
  return nativeVal;

}
