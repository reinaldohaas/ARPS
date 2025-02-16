/*
   This program reads and decodes 1km topography data. It is similar to
   getgg.f program but implemented in C, to simply some of
   low level I/O to the topo files and for improve portability.

   xgetgg: Return an array of nlats by nlons dimension of terrain
   height values for a given lat/lon cell.

   programmer: Pedro Tsai,  NRL Monterey, Dec 1996
   rcs keywords: uvgrid.F,v
               2.0 1996/09/27 22:43:44
   SCCS IDENTIFICATION:  %W% %G%

*/

//#include <xgetgg.prol>

#define _POSIX_SOURCE 1

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>



#define ntypes 1
static char typetab[ntypes][32] = {
  "DTED1"
} ;

#define nparms 6
static char parmtab[nparms][32] = {
  "terr_ht","terr_ht_mean","terr_ht_median",
  "terr_ht_max","terr_ht_min","terr_ht_std_dev",
} ;

#define maxres 19

static int ndestab[ntypes] = { 4, } ;

static int nrestab[maxres] = { 19,} ;

static int restab[maxres]= {
  1,2,3,4,5,6,7,8,9,10,12,13,15,17,20,25,34,50,100
} ;

static int lapftab[maxres]= {
  1,1,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
} ;

static int lopftab[maxres]= {
  1,1,5,5,5,5,5,5,5,20,20,20,20,20,20,20,20,180,180
} ;

static int arestab[maxres]= {
  100.0 ,200.0 ,294.1, 400.0 ,500.0 ,588.2,666.7,759.2,
  833.3,1000.0,1111.1,1250.0,1428.6,1666.7,2000.0,
  2500.0,3333.3,5000.0,10000.0
} ;

static int nlattab[maxres]= {
  1200,600,408,300,240,204,180,156,144,120,108,
  96,84,72,60,48,36,24,12
} ;

/*
  Print out a bit pattern of an integer value
*/
static void bitprint(int v)
{
  int i , mask ;
  int wdsize ;

#ifdef DEBUG
  wdsize = sizeof(int) * 8 ;

  mask = 1 << (wdsize-1) ; /* mask = 100000000...0 */
  for (i=1; i <= wdsize ; i++) {
    putchar(((v & mask) == 0 ) ? '0' : '1' ) ;
    v <<= 1 ;
    if (i%8 == 0 && i != wdsize )
      putchar(' ') ;
  }
#endif
}

/*
  Find the length of a character string. String can consists of
  alphanumeric characters, dot . , underscore _ , slash / ,
  dash - , tilda ~ .
  */
static int string_len(char *f)
{
  int i ; i=0 ;
  /*  while ( ! isspace(f[i]) ) i++ ; */
  while ( isalnum(f[i])
	  || f[i] == '.'
	  || f[i] == '_'
	  || f[i] == '/'
	  || f[i] == '-'
          || f[i] == '~') {
    /* fprintf(stdout,">>>>> %c \n",f[i]) ; */
    i++ ;
  }


return i ;
}

/*
   Decode a 1-4 bytes value and return its integer value. Assuming
   byte[0] is the most significant byte (31-16 bits), and byte[1]
   is next significant byte (15-0 bits)
*/
static int bytes2int(char* byteary, int size_of_byteary)
{
  int val,i ,tmpval ;
  int x,negative ;

  if ( size_of_byteary > sizeof(int) ) {
    fprintf(stdout,
    "Error, the size of integer is to small to store the value of byte array.\n") ;
    fprintf(stdout,"Returned integer value is set to 0\n") ;
    return 0 ;
  }
  /*
    Check for negative sign on the most significant digit.
    If the most significant bit is 1, then this value is negative.
  */
  negative = 0 ;
  val = 0 ; tmpval = 0 ;
  if ( byteary[0] & 0x80 ) {
    negative = 1 ;
  }

  /* Store the first byte in the lower 8 bits */
  val =  byteary[0] & 0xff ;

  /*
     Process the remainder bytes.
     First copy the next byte and the tmpval by AND its bits to tmpval,
     Next, the lower 8 bits of val is shift left 8 position,
     and then the tmpval is OR into val variable.
  */
  for ( i = 1 ; i < size_of_byteary ; i++ ) {
    tmpval = byteary[i] & 0xff  ;
    val = val << 8 ;
    val = val | tmpval ;
    tmpval = 0 ;
  }

  /* Debug, print out bit pattern */
  /* bitprint(val) ; printf("\n") ; */

  /*
     If val is negative, we need to pad the remainding bits of val
     with 1. We do this by creating a bits mask of 1, for example,
     for two bytes values we want a bit mask of 0xffff0000, and OR
     this mask with val.
  */
  x = ~0 ; x = x << (size_of_byteary * 8) ;

  if ( negative ) {
    val = val | x  ;
  }
  return val ;
}

/*
  Set c function name so a fortarn program can call it.
  First, unset the macro HAS_FORTRAN_NAME to false, and
  chech for each system type. The Default name is c-function
  name with an under_score appended.
*/

#undef HAS_FORTRAN_NAME

#ifdef CRAY
void XGETGG(int *la, int *lo,
           int *nd, int *ides, int *lbuf, float *buf ,int *nla ,
	   int *nlo, int *istat)
#define HAS_FORTRAN_NAME
#endif

#if (defined HP || defined IBM_SP)
void xgetgg(int *la, int *lo,
           int *nd, int *ides, int *lbuf, float *buf ,int *nla ,
	   int *nlo, int *istat)
#define HAS_FORTRAN_NAME
#endif

#ifndef HAS_FORTRAN_NAME
void xgetgg_(int *la, int *lo,
           int *nd, int *ides, int *lbuf, float *buf ,int *nla ,
	   int *nlo, int *istat)
#undef HAS_FORTRAN_NAME
#endif

{
  int i, index, ndes, ires,nres, lat, lon, nlats, nlons , iskip ;
  int lat1,lon1 ;
  char type[32],parm[32],path[256], fname[256], cbuf[64] ;
  static char type_id[32];
  static char parm_id[32];
  static char path_id[256];
  int ihead[6] , jcell, icell, ncell , kchd ;
  int n, j, rc, ptsprec , celcnt, recl, n0dat ;
  int val ,nrec, currec,mcell, irec, len ;
  int x ;
  FILE *fdat;
  char byte[2] ;
  char rcd[248] ;
  long fpos ;
  long offset ;
  FILE *tmpf;
  static int doneflg ;
  int l_lat,l_lon;
  /* Get argument list from temp file */
  if ( doneflg != 100 ) {
    tmpf=fopen("./xgetgg_arg","r");
    if (tmpf==NULL) {
      fprintf(stdout,"XGETGG Error: Can't open temp file xgetgg_arg\n");
      *istat=-1;
      fclose(tmpf);
      return;
    } else {
      fscanf(tmpf,"%s",type_id);
      fscanf(tmpf,"%s",parm_id);
      fscanf(tmpf,"%s",path_id);
      doneflg=100;
      fclose(tmpf);
      /* print argument list */
      fprintf(stdout,"type_id = %s \n", type_id) ;
      fprintf(stdout,"parm_id = %s \n", parm_id) ;
      fprintf(stdout,"path_id = %s \n", path_id) ;
    }
  }

   /* print argument list
  fprintf(stdout,"type_id = %s \n", type_id) ;
  fprintf(stdout,"parm_id = %s \n", parm_id) ;
  fprintf(stdout,"path_id = %s \n", path_id) ;
  fprintf(stdout,"*la = %d \n",*la ) ;
  fprintf(stdout,"*lo = %d \n",*lo ) ;
  fprintf(stdout,"*nd = %d \n",*nd ) ;
  fprintf(stdout,"*ides = %d \n",*ides ) ;
  fprintf(stdout,"*lbuf = %d \n",*lbuf ) ;
   */

  /* copy the input argument parameters to local variables */
  lat = *la ;
  lon = *lo ;
  lat1 = lat ;
  lon1 = lon ;
  ndes = *nd ;

  /* check input values */
  if ( ( lon < -180 ) || ( lon > 180 ) ) {
    fprintf(stdout,"*** XGETGG ERROR - longitude out of range ***\n") ;
    fprintf(stdout,"value is %d, range is -180 to +180\n",lon) ;
    *istat = -4 ;
    return ;
  }
  if ( ( lat < -90 ) || ( lat > 90 ) ) {
    fprintf(stdout,"*** XGETGG ERROR - latitude out of range ***\n") ;
    fprintf(stdout,"value is %d, range is -90 to +90\n",lat) ;
    *istat = -4 ;
    return ;
  }

  /*
     Check the data base terrain type
     Note: currently only one typetab is supported: DTED1

     Find the location of '_' char symbol, the type name is the string
     before the '_' symbol

  */
  i=0 ;
  while ( type_id[i] != '_' ) {
    i++ ;
  }
  len=i ;
  strncpy(type,type_id,len) ;
  ires = atoi(&type_id[len+2]) ;

  /*
     Look up the type table, if we don't find the type_id in the table,
     issue error message and return .
  */
  for ( i = 0 ; i < ntypes ; i++ ) {
    if ( strncmp(type,typetab[i],len) == 0 ) {
      break ;
    } else {
      if ( i == (ntypes-1) ) {
	fprintf(stdout,
		"***XGETGG ERROR, type <%s> not recognized***\n",
		type);
	*istat=-3 ;
	return ;
      }
    }
  }
  /* index is the table location of the current data type */
  index = i ;
  if ( ndes < ndestab[index] ) {
    fprintf(stdout,"***XGETGG ERROR, ndes too small for type***\n");
    *istat=-5 ;
    return ;
  } else {
    ndes = ndestab[index] ;
  }

  /*
    Check the user request topo resolution with terrain database
    resolution store in the table (nrestab).

    ires is given in unit of 100m, so for resolution of 1km
    value of ires is 10
  */

#ifdef DEBUG_PRINT
   fprintf(stdout,">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n") ;
   fprintf(stdout,"Requested database resolution is %d\n",ires) ;
#endif
   /* process all 1 km resolution latitude and longitude circle
    */
 for ( l_lat=0; l_lat < 5; l_lat++) {
       lat=lat1+l_lat ;
   for ( l_lon=0; l_lon < 20; l_lon++) {
     if (lon < 0) {
       lon=lon1-l_lon ;
     } else {
       lon=lon1+l_lon ;
     }
#ifdef DEBUG_PRINT
   fprintf(stdout,"in xgetgg lat = %d\n",lat) ;
   fprintf(stdout,"in xgetgg lon = %d\n",lon) ;
   fprintf(stdout,"in xgetgg l_lon = %d\n",l_lon) ;
   fprintf(stdout,"in xgetgg l_lat = %d\n",l_lat) ;
   fprintf(stdout,"in xgetgg l_lon = %d\n",l_lon) ;
#endif

   nres = nrestab[index] ;
   for ( i = 0 ; i < nres ; i++ ) {
     if ( ires == restab[i]) {
       break ;
     } else {
       if ( i == (nres-1) ) {
	 fprintf(stdout,"***XGETGG ERROR - %d not a proper resolution\n***",
      	   ires) ;
	 *istat = -7 ;
	 return ;
       }
     }
   }

   /*
     compute parameters related to cell.
     */
   nlats = 1200/ires ;

#ifdef DEBUG_PRINT
   fprintf(stdout,"Number of latitude point is %d\n",nlats) ;
#endif

   /*
     set iskip to number of longitude columns in cell
   */
   if ( lat < -80 )
     iskip = 6 ;
   else if ( lat < -75)
     iskip = 4 ;
   else if ( lat < -70)
     iskip = 3 ;
   else if ( lat < -50)
     iskip = 2 ;
   else if ( lat < 50)
     iskip = 1 ;
   else if ( lat < 70)
     iskip = 2 ;
   else if ( lat < 75)
     iskip = 3 ;
   else if ( lat < 80)
     iskip = 4 ;
   else
     iskip = 6 ;

#ifdef DEBUG_PRINT
   fprintf(stdout,"longitudinal skip factor: %d \n",iskip) ;
   fprintf(stdout,"for lat/lon values of %d %d \n",lat, lon);
#endif

  /*
    compute parameters related to this file
    ihead(3) and ihead(4) are the size of the cell in degree latitude
    and longitude.
  */

  ihead[2] = lapftab[i] ;
  ihead[3] = lopftab[i] ;
  ihead[0] =  ((lat+90)/ihead[2])*ihead[2]-90 ;

#ifdef DEBUG_PRINT
   fprintf(stdout,"in xgetgg i = %d\n",lopftab[i]) ;
   fprintf(stdout,"in xgetgg lon = %d\n",lon) ;
#endif
  ihead[1] = ((lon+180)/ihead[3])*ihead[3]-180 ;
  jcell = lat - ihead[0] ;
  icell = lon - ihead[1] ;
  ncell = (jcell*ihead[3])+icell ;
  kchd  = (ncell*(ndes+1)) + 8 ;

#ifdef DEBUG_PRINT
  fprintf(stdout,"kchd = %d \n",kchd) ;
#endif

  /*
    Set the number of longitudinal and latitudinal points in the file
    ptsprec is the number of latitude value in this file
    recl is length (in bytes) of the a single latitude values record

    */

  nlons = nlats/iskip ;
  ptsprec = nlats + 1 ;
  celcnt = 0 ;
  ihead[4] = nlats ; *nla = nlats ;
  ihead[5] = nlons ; *nlo = nlons ;

  /*
     compute the size (bytes) of each record. Each record hold a nlats+1
     values.
     */
  recl = ptsprec*2 ;

  /*
    Round the recl size to be multiple of 8 bytes
    */

  if( (recl%8) != 0 ) recl = (recl/8 + 1) * 8 ; /* recl is 248 bytes */
  n0dat = ((((ndes+1)*(ihead[2]*ihead[3])) + 7)/ptsprec) + 1 ;

  /*
     Create path/filename
     */

 if( l_lat == 0 && l_lon == 0 ) {
  n=string_len(path_id) ;

  strncpy(path,path_id,n) ;
  if ( path[n-1] != '/' ) {
    path[n] = '/' ;
    path[n+1] = '\0' ;
    n=n+1 ;
  }
  /*
    file name is path/type.parm.latlon
    */
  fname[0]='\0';
  strncat(fname,path,n) ;
  n=string_len(type_id) ;
  strncat(fname,type_id,n) ;
  strcat(fname,".") ;
  n=string_len(parm_id) ;
  strncat(fname,parm_id,n) ;
  strcat(fname,".") ;
  if ( ihead[1] <  0 )
    strcat(fname,"W") ;
  else
    strcat(fname,"E") ;
  sprintf(cbuf,"%03d\0",abs(ihead[1])) ;
  strcat(fname,cbuf) ;
  if ( ihead[0] < 0 )
    strcat(fname,"S") ;
  else
    strcat(fname,"N") ;
  sprintf(cbuf,"%02d\0",abs(ihead[0])) ;
  strcat(fname,cbuf) ;


  if ( access(fname,F_OK|R_OK) != 0 ) {
    fprintf(stdout,"***XGETGG>> File does not exist: %s\n",fname) ;
    fprintf(stdout,"***      >> Assume all values in this cell are 0\n");
    *istat = 0 ;
    return ;
  }
  /*
     open datafile for reading
     */
  fdat=fopen(fname,"rb") ;
  if (fdat == NULL ) {
    fprintf(stdout,"***XGETGG>> Error: Can't open file: %s\n",fname) ;
    *istat = -1 ;
    fclose(fdat);
    return ;
  }
  fprintf(stdout,"Opened file %s for reading\n",fname) ;

  /*
     Decode the first 6 values in the record , each values is stored
     in two bytes
  */
  n=0 ;
  for (i=0 ; i < 6 ; i++ ) {
    rc=fread(byte,sizeof(byte),1,fdat) ;
    if (rc!=1) {
      fprintf(stdout,"***XGETGG ERROR>> reading file: %s\n",fname) ;
      *istat = -1 ;
      fclose(fdat);
      return ;
    }
    val = bytes2int(byte, sizeof(byte)) ;
    if ( val != ihead[i] ) {
      fprintf(stdout,"**XGETGG ERROR - header mismatch** %d %d %d\n",
	      i,val,ihead[i]) ;
      *istat = -6 ;
      fclose(fdat) ;
      return ;
    }

#ifdef DEBUG_PRINT
    fprintf(stdout,"ihead[%d]= %d\n",i,ihead[i]) ;
#endif

  }
  /* Get celcnt and ndes values */
  rc=fread(byte,sizeof(byte),1,fdat) ;
  val = bytes2int(byte, sizeof(byte)) ;
  if ( val != ndes ) {
    fprintf(stdout,"**XGETGG ERROR - header mismatch** %d %d %d\n",
	    i,val,ndes) ;
    *istat = -6 ;
    fclose(fdat) ;
    return ;
  }

  rc=fread(byte,sizeof(byte),1,fdat) ;
  celcnt = bytes2int(byte, sizeof(byte)) ;

#ifdef DEBUG_PRINT
  fprintf(stdout,"ndes= %d\n",ndes) ;
  fprintf(stdout,"celcnt= %d\n",celcnt) ;
#endif

  /* endif first loop check */
}
  /*
    general header processing finished
    check header for this specific cell

    move the file pointer to the new location and read the cell info
    and check to see if target cell present.

    */

  /*  each value is stored in 2 bytes */
  fpos = ( (kchd/ptsprec) * 248) + ( (kchd % ptsprec) * 2) ;
  fseek(fdat,fpos,SEEK_SET) ;
  rc=fread(byte,sizeof(byte),1,fdat) ;
  if (rc!=1) {
    fprintf(stdout,"***XGETGG ERROR>> Reading mcell info failed***: %s\n"
	    ,fname) ;
    *istat = -2 ;
    fclose(fdat) ;
    return ;
  }

  /*
     mcell is 0 , this mean the this cell is over ocean or not yet
     implemented.  Return to calling program.
     */
  mcell = bytes2int(byte, sizeof(byte)) ;
  offset=(nlons+(nlats+1)*(ptsprec-1)+1)*l_lon+(nlons+(nlats+1)*(ptsprec-1)+1)*20*l_lat ;
  if ( mcell == 0 ) {
    *istat = 0 ;
    for (i=0 ; i < nlons + 1 ; i++ ) {
      for (j=0 ; j < ptsprec ; j++ ) {
        x=i + (nlats +1) * j + offset ;
        buf[x]=0.;
      }
    }

#ifdef DEBUG_PRINT
  fprintf(stdout,"mcell is %d\n",mcell) ;
  fprintf(stdout,"in xgetgg offset = %d\n",offset) ;
  fprintf(stdout,"fpos %d\n",fpos) ;
    fprintf(stdout,"in xgetgg ptsprec = %d\n",ptsprec) ;
    fprintf(stdout,"in xgetgg nlats = %d\n",nlats) ;
    fprintf(stdout,"in xgetgg nlons = %d\n",nlons) ;
    fprintf(stdout,"mcell= %d\n",mcell) ;
    fprintf(stdout,"return to calling subroutine \n") ;
#endif
  } else {

#ifdef DEBUG_PRINT
  fprintf(stdout,"mcell= %d\n",mcell) ;
#endif

  for (i=0 ; i < ndes ; i++ ) {
    rc=fread(byte,sizeof(byte),1,fdat) ;
    if (rc!=1) {
      fprintf(stdout,"***XGETGG ERROR>> Reading ndes info failed***: %s\n"
	      ,fname) ;
      *istat = -2 ;
      fclose(fdat) ;
      return ;
    }
    ides[i]=bytes2int(byte, sizeof(byte)) ;

#ifdef DEBUG_PRINT
    fprintf(stdout,"ides[%d]= %d\n",i,ides[i]) ;
#endif

  }

  /*
    header processing complete, read the data
    */
  irec = n0dat + (mcell-1)*(nlons+1) ;

#ifdef DEBUG_PRINT
  fprintf(stdout,"in xgetgg irec = %d\n\n",irec) ;
  fprintf(stdout,">>>> Data starts at record %d\n",irec);
#endif

  n=0 ;
  for (i=0 ; i < nlons + 1 ; i++ ) {
    fpos = (irec+i) * 248 ;
    fseek(fdat,fpos,SEEK_SET) ;
    for (j=0 ; j < ptsprec ; j++ ) {
      rc=fread(byte,sizeof(byte),1,fdat) ;
      if (rc!=1) {
        fprintf(stdout,
		"***XGETGG ERROR>> Reading data file failed***: %s\n",
		fname) ;
	*istat = -2 ;
        fclose(fdat) ;
        return ;
      }
      x=i + (nlats +1) * j + offset ;
      buf[x]= bytes2int(byte, sizeof(byte)) ;
    }
  }
  /* endif mcell */
  }
  /*  close the for l_lat and l_lon loop
   */
  }
 }
  *istat = iskip ;
  fclose(fdat) ;

  return ;
}
