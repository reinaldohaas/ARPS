/*  This c subroutine takes the place of the
    VAX system subroutine str$upcase 
    Keith Brewster March 1995 */
#ifdef UNDERSCORE
#define str_upcase str_upcase_
#endif
void str_upcase(char *str,char *strup,long *length)
{
  int i;
  for ( i = 0; i < *length; i++)
  {
    *(strup+i)=toupper(*(str+i));
  }
}
