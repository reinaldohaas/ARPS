#include <ctype.h>
#include <stdlib.h>
#include <string.h>

/**********************************************************************
 *
 * Trim a string with leading and trailing spaces
 *
 *********************************************************************/

char * trim(char *string) {
  size_t slen;
  char *tempstr;
  int i,j,k;

  slen = strlen(string);
  
  /* trim leading space */
  j = 0;
  while (string[j] == ' ' || string[j] == '\t') j++;

  /* trim trailing space */
  k = slen-1;
  while (string[k] == ' ' || string[k] == '\t' 
      || string[k] == '\n' || string[k] == '\r') k--;

  slen = k-j+2;
  tempstr = (char *)malloc(slen*sizeof(char));

  for(i=j;i<=k;i++) tempstr[i-j] = string[i];
  tempstr[i] = '\0';

  return tempstr;
}

/**********************************************************************
 *
 * Convert a string to lower case
 *
 *********************************************************************/

char * strtolower(char *str) {

  size_t slen;
  char *tempstr;
  int i;

  slen = strlen(str);
  tempstr = (char *)malloc((slen+1)*sizeof(char));

  for (i = 0; i < slen; i++) tempstr[i]=tolower(str[i]);
  tempstr[i] = '\0';

  return tempstr;
}

/**********************************************************************
 *
 * Extract the first n characters of a string
 *
 *    First trim the string;
 *    Then convert to lower case;
 *    At last extract the characters.
 *
 *********************************************************************/

char * strhead(char *str, int n) {
  char *tstr1,*tstr2;
  char *tem;
  int i;

  tstr1 = trim(str);
  tstr2 = strtolower(tstr1);
 
  tem = (char *) malloc((n+1)*sizeof(char));
  for (i=0;i<n;i++) 
  {
    tem[i] = tstr2[i];
  }
  tem[n] = '\0';

  free(tstr1);
  free(tstr2);

  return tem;
} 

/**********************************************************************
 *
 * Extract the last n characters of a string
 *
 *    extract the characters.
 *
 *********************************************************************/

char * strtail(char *str, int n) {
  char *tem;
  int slen,i;

  slen = strlen(str);
  if (slen < 4) perror("String is less than needed in strtail.\n");

  tem = (char *) malloc((n+1)*sizeof(char));
  for (i=0;i<n;i++) 
  {
    tem[i] = str[slen-n+i];
  }
  tem[n] = '\0';

  return tem;
} 
