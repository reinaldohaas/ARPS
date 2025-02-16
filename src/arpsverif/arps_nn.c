#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define maxnin 100    /*max number of input nodes*/
#define maxnhd 100    /*max number of hidden nodes*/
#define maxnout 20    /*max number of output nodes */

static int nin,nh1,nout,act1,act2;
static double avg[maxnin],sigma[maxnin];
static double w[maxnout][maxnhd],wh1[maxnhd][maxnin];
static char wtfile[20],statfile[20];

void stats();
void weights();
void activation ( float *, float *);
double dotprod (int , double *, double *);
double act_func (double);
static void timetocode ( int, double *, double *, double *, double *);

/*****************************************************************/
arps_nn_( char station[5],float nn_input[maxnin], float nn_output[1] )
{
    strcpy(wtfile,"wts_");
    strcpy(statfile,"stat_");
    strncat(wtfile,station,4);
    strncat(statfile,station,4);

     stats();
     weights();
     activation(nn_input,&nn_output[0]);
}
/*****************************************************************/
void stats()
{
int i,n,nlines=19;          /* number of lines in stat file */
double min,max;
FILE *fs;

   if((fs=fopen(statfile,"r")) ==NULL) {printf("Missing stat file.\n"); exit(0);}

   for(i=0;i<nlines;i++)
   fscanf(fs,"%d %lf %lf %lf %lf",&n,&avg[i],&sigma[i],&min,&max);
   fclose(fs);
}
/******************************************************************/
void weights( )
{
int i,j;
FILE *fw;
char junk[150];

 if((fw=fopen(wtfile,"r")) ==NULL) {printf("Missing weight File.\n"); exit(0);}

     fscanf(fw,"%s",junk);    /* read in the specs */
     fscanf(fw,"%s",junk);
     fscanf(fw,"%s",junk);
     fscanf(fw,"%s %d",junk,&nin);
     fscanf(fw,"%s %d %d",junk,&nh1,&act1);
     fscanf(fw,"%s %d %d",junk,&nout,&act2);
     fscanf(fw,"%s",junk);
     fscanf(fw,"%s %s",junk,junk);
     fscanf(fw,"%s %s",junk,junk);
     fscanf(fw,"%s %s",junk,junk);
     fscanf(fw,"%s %s",junk,junk);
     fscanf(fw,"%s",junk);
     fscanf(fw,"%s %s",junk,junk);
     fscanf(fw,"%s %s",junk,junk);
     fscanf(fw,"%s",junk);
     fscanf(fw,"%s",junk);
     fscanf(fw,"%s",junk);
     fscanf(fw,"%s %s %s",junk,junk,junk);
     fscanf(fw,"%s %s",junk,junk);
     fscanf(fw,"%s %s",junk,junk);
     fscanf(fw,"%s %s %s",junk,junk,junk);
     fscanf(fw,"%s %s %s",junk,junk,junk);
     fscanf(fw,"%s %s %s",junk,junk,junk);

     if(!nh1){
        for(i=0;i<nout;i++){
        fscanf(fw,"%lf ",&w[i][nin]);
          for(j=0;j<nin;j++)
          fscanf(fw,"%lf ",&w[i][j]);
        }
     } else {
         for(i=0; i<nh1;i++){
         fscanf(fw,"%lf ",&wh1[i][nin]);
             for(j=0;j<nin;j++)
             fscanf(fw,"%lf ",&wh1[i][j]);
         }
         for(i=0; i<nout;i++){
         fscanf(fw,"%lf ",&w[i][nh1]);
             for(j=0;j<nh1;j++)
             fscanf(fw,"%lf ",&w[i][j]);
         }
     }

    fclose(fw);
}
/*****************************************************************/
double dotprod (
   int n ,
   double *vec1 ,
   double *vec2 )
{
   int k, m ;
   double summ ;

   summ = 0.0 ;
   k = n / 4 ;
   m = n % 4 ;

   while (k--) {
      summ += *vec1 * *vec2 ;
      summ += *(vec1+1) * *(vec2+1) ;
      summ += *(vec1+2) * *(vec2+2) ;
      summ += *(vec1+3) * *(vec2+3) ;
      vec1 += 4 ;
      vec2 += 4 ;
      }

   while (m--)
      summ += *vec1++ * *vec2++ ;

   return summ ;
}
/************************************************************************/
double act_func ( double x )
{
return  1./(1+exp(-x));
}
/**********************************************************************/
void activation ( float *input, float *output)
{
int i;
double in[maxnin],out[maxnout],outh1[maxnhd];

   timetocode(input[0],&in[0],&in[1],&in[2],&in[3]); /*convert time to input*/
   in[nin-1]=(input[1]-avg[2])/sigma[2]; /*scale arps_fore accordingly */
   in[nin]=1.;

          if(!nh1){              /*no hidden layers*/
            for(i=0;i<nout;i++){
              if(act1==0)
              out[i]=dotprod(nin+1,in,w[i]);
              if(act1==3)
              out[i]=act_func(dotprod(nin+1,in,w[i]));
            }
          }
          else{              /*one hidden layer*/
             for(i=0;i<nh1;i++){
               if(act1==0)
               outh1[i]=dotprod(nin+1,in,wh1[i]);
               if(act1==3)
               outh1[i]=act_func(dotprod(nin+1,in,wh1[i]));
             }
             outh1[nh1]=1.;
               for(i=0;i<nout;i++){
                  if(act2==0)
                  out[i]=dotprod(nh1+1,outh1,w[i]);
                  if(act2==3)
                  out[i]=act_func(dotprod(nh1+1,outh1,w[i]));
               }
          }
      *output= (float) out[0]*sigma[17]+avg[17];
}
/**************************************************************/
void timetocode ( int time, double *x, double *y, double *w, double *z)
{
int i,j,k,l,cnt=0;

     if(time>80){ printf("Error: Time problem.\n"); exit(0); }

          for(i=-1;i<=1;i++){
              for(j=-1;j<=1;j++){
                  for(k=-1;k<=1;k++){
                      for(l=-1;l<=1;l++){
                         if(cnt==time){
                         *x= (double) i;
                         *y= (double) j;
                         *w= (double) k;
                         *z= (double) l;
                         }
                      ++cnt;
                      }
                  }
              }
          }

}
