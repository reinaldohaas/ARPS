/*  calib.c 
    Calibration for IDD and gvar raw data from GOES-8,9,10

    Keith Brewster
    CAPS/Univ of Oklahoma
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mc_area.h"

int chtab[4] = {2, 3, 4, 5};
int fkidx[6] = {-1, -1, 0, 1, 2, 3};
int iridx[6] = {-1, -1, 3, 4, 1, 2};

/*
  These calibration data are obtained from the following URL
  http://cimss.ssec.wisc.edu/calibration/
  
*/

double g08fk1[4] = { 199980.00, 38792.00, 9737.80, 6944.50 };
double g08fk2[4] = {   3684.20,  2132.70, 1345.30, 1202.00 };
double g08tc1[4] = {   0.63568,  0.60603, 0.37351, 0.22171 };
double g08tc2[4] = {   0.99911,  0.99858, 0.99873, 0.99916 };

double g09fk1[4] = { 198810.00, 38732.00, 9717.20, 6899.50 };
double g09fk2[4] = {   3677.00,  2131.60, 1344.40, 1199.40 };
double g09tc1[4] = {   0.58637,  0.48409, 0.36225, 0.20144 };
double g09tc2[4] = {   0.99917,  0.99887, 0.99876, 0.99923 };

double g10fk1[4] = { 198400.00, 39086.00, 9774.30, 6828.50 };
double g10fk2[4] = {   3674.50,  2138.10, 1347.00, 1195.20 };
double g10tc1[4] = {   0.62223,  0.61437, 0.27790, 0.21145 };
double g10tc2[4] = {   0.99912,  0.99857, 0.99905, 0.99919 };

double g11fk1[4] = { 200180.00, 38789.00, 9653.40, 6877.80 };
double g11fk2[4] = {   3685.40,  2132.70, 1341.50, 1198.10 };
double g11tc1[4] = {   0.62864,  0.59339, 0.38284, 0.20258 };
double g11tc2[4] = {   0.99911,  0.99861, 0.99869, 0.99923 };

double g12fk1[4] = { 200960.00, 43702.00, 9685.90, 5047.10 };
double g12fk2[4] = {   3690.20,  2219.10, 1343.00, 1080.70 };
double g12tc1[4] = {   0.69703,  5.08315, 0.37554, 0.09537 };
double g12tc2[4] = {   0.99902,  0.98872, 0.99872, 0.99960 };

/* vis2radiance
   From raw GVAR visible data compute visible radiance */

int vis2radiance( gvcal, visim, vradiance, npixel )
     struct cal_gvar *gvcal;
     long *visim;
     float *vradiance;
     long npixel;

{
    static float vis1g;
    static float vis2g;
    static float vbias;
    static int avgd;
    long rawval;
    int i;

    if ( avgd != 1 ) {
        vis1g=0.;
        vis2g=0.;
        vbias=0.;
        for ( i=0 ; i < 8 ; i++) {
            vis1g+=gvcal->vis1gain[i];
            vis2g+=gvcal->vis2gain[i];
            vbias+=gvcal->visbias[i];
            }
        vis1g=0.125*vis1g;
        vis2g=0.125*vis2g;
        vbias=0.125*vbias;
        printf(" Average vis 1 gain: %10.4f\n",vis1g);
        printf(" Average vis 2 gain: %10.4f\n",vis2g);
        printf(" Average vis bias:   %10.4f\n",vbias);
        avgd=1;
        }

    for ( i=0 ; i < npixel ; i++ ){
        rawval=*(visim+i);
        *(vradiance+i)=
                     (vis2g*(rawval*rawval))+(vis1g*rawval)+vbias;
        }
    return(0);
}

/* vis2albedo
   From raw GVAR visible data compute visible albedo */

int vis2albedo( gvcal, visim, albedo, npixel)
     struct cal_gvar *gvcal;
     long *visim;
     float *albedo;
     long npixel;

{
    static float vis1g;
    static float vis2g;
    static float vbias;
    static int avgd;
    float albtem,albmin,albmax;
    long rawval;
    int i;

    if ( avgd != 1 ) {
        vis1g=0.;
        vis2g=0.;
        vbias=0.;
        for ( i=0 ; i < 8 ; i++) {
            vis1g+=gvcal->vis1gain[i];
            vis2g+=gvcal->vis2gain[i];
            vbias+=gvcal->visbias[i];
            }
        vis1g=0.125*vis1g;
        vis2g=0.125*vis2g;
        vbias=0.125*vbias;
        printf(" Average vis 1 gain: %10.4f\n",vis1g);
        printf(" Average vis 2 gain: %10.4f\n",vis2g);
        printf(" Average vis bias:   %10.4f\n",vbias);
        avgd=1;
        }


    albmin=1.0;
    albmax=0.0;
    for ( i=0 ; i < npixel ; i++ ){
        if( *(visim+i) > 0 ) {
            rawval=*(visim+i);
            albtem=gvcal->albedcon*
                   ( (vis2g*(rawval*rawval))+(vis1g*rawval)+vbias);
            if( albtem > 1.0)
                albtem = 1.0;
            if( albtem < 0.0 )
                albtem = 0.0;
            if( albtem > albmax)
                albmax = albtem;
            if( albtem < albmin)
                albmin = albtem;
            *(albedo+i)=albtem;
            }
        else
            *(albedo+i)=-999.0;
            
        if( (i % 500) == 0 ) {
            printf(" pixel %i raw %i albedo %10.4f\n",
                       i,(int)*(visim+i),*(albedo+i));
            }
        }
    printf(" vis2albedo: Minimum albedo: %8.4f   Maximum albedo: %8.4f\n",
             albmin,albmax);
    return(0);
}

/* ir2radiance
   From raw GVAR IR data compute IR radiance */

int ir2radiance( gvcal, irim, radiance, ichan, npixel )
     struct cal_gvar *gvcal;
     long *irim;
     float *radiance;
     int ichan;
     long npixel;

{
    double gaininv;
    double bias;
    double rawval;
    int i,index;

    index=iridx[ichan];

    printf( " gain is %10.4f   bias: %10.4f\n",
             gvcal->ir1gain[index],gvcal->ir1bias[index]);
    gaininv = 1./((double)(gvcal->ir1gain[index]));
    bias = (double)(gvcal->ir1bias[index]);
    for ( i=0 ; i < npixel ; i++ ){
        rawval=(double)(*(irim+i));
        *(radiance+i)=(float)(gaininv*(rawval-bias));
        }
    return(0);
}

/* ir2bright
   From raw GVAR IR data compute brightness temperature */

int ir2bright( gvcal, irim, bright, ibird, ichan, npixel )
     struct cal_gvar *gvcal;
     long *irim;
     float *bright;
     int ibird;
     int ichan;
     long npixel;

{
    double gaininv,bias;
    double rawval;
    double radiance;
    double expon,tt;
    float brttem,brtmin,brtmax;
    int i, index;

    index=iridx[ichan];
    gaininv = 1./((double) (gvcal->ir1gain[index]));
    bias = (double) (gvcal->ir1bias[index]);
    printf( " gain is %10.4f   bias: %10.4f\n",
             gvcal->ir1gain[index],bias);
    index=fkidx[ichan];
    if ( index > -1 ) {
       printf(" Using calibration for goes%2d Channel %d\n",
              ibird,chtab[index]);
       }
    else {
       printf(" Problem with table indexing for ichan %d\n",
              ichan);
       printf(" Stopping in ir2bright\n");
       exit(0);
       }

    brtmax=0.0;
    brtmin=999.0;
    if ( ibird == 8 ){
        printf( " fk1:%10.2f  fk2:%10.2f  tc1:%8.4f  tc2:%8.4f\n",
            g08fk1[index],g08fk2[index],g08tc1[index],g08tc2[index]);
        for ( i=0 ; i < npixel ; i++ ){
            if(*(irim+i) > 0 ) {
                rawval=(double)(*(irim+i));
                radiance=gaininv*(rawval-bias);
                expon=(g08fk1[index]/radiance) + 1.;
                tt=g08fk2[index]/log(expon);
                brttem=(tt-g08tc1[index])/
                           g08tc2[index];
                if( brttem < 0.0 ) 
                    brttem = 0.0;
                if( brttem > brtmax )
                    brtmax = brttem;
                if( brttem < brtmin )
                    brtmin = brttem;
                *(bright+i)=brttem;
                }
            else {
                *(bright+i)=-999.;
                }
            if( (i % 500) == 0 ) {
                printf(" pixel:%i raw:%i radiance:%10.2f bright:%10.2f\n",
                       i,(int)*(irim+i),radiance,*(bright+i));
                }
            }
        }
    else if ( ibird == 9 ){
        printf( " fk1:%10.2f  fk2:%10.2f  tc1:%8.4f  tc2:%8.4f\n",
            g09fk1[index],g09fk2[index],g09tc1[index],g09tc2[index]);
        for ( i=0 ; i < npixel ; i++ ){
            if(*(irim+i) > 0 ) {
                rawval=(double)(*(irim+i));
                radiance=gaininv*(rawval-bias);
                expon=(g09fk1[index]/radiance) + 1.;
                tt=g09fk2[index]/log(expon);
                brttem=(tt-g09tc1[index])/
                           g09tc2[index];
                if( brttem < 0.0 ) 
                    brttem = 0.0;
                if( brttem > brtmax )
                    brtmax = brttem;
                if( brttem < brtmin )
                    brtmin = brttem;
                *(bright+i)=brttem;
                }
            else {
                *(bright+i)=-999.;
                }
            if( (i % 500) == 0 )
                printf(" pixel:%i raw:%i radiance:%10.2f bright:%10.2f\n",
                       i,(int)*(irim+i),radiance,*(bright+i));
    
            }
        }
    else if ( ibird == 10 ){
        printf( " fk1:%10.2f  fk2:%10.2f  tc1:%8.4f  tc2:%8.4f\n",
            g10fk1[index],g10fk2[index],g10tc1[index],g10tc2[index]);
        for ( i=0 ; i < npixel ; i++ ){
            if(*(irim+i) > 0 ) {
                rawval=(double)(*(irim+i));
                radiance=gaininv*(rawval-bias);
                expon=(g10fk1[index]/radiance) + 1.;
                tt=g10fk2[index]/log(expon);
                brttem=(tt-g10tc1[index])/
                           g10tc2[index];
                if( brttem < 0.0 ) 
                    brttem = 0.0;
                if( brttem > brtmax )
                    brtmax = brttem;
                if( brttem < brtmin )
                    brtmin = brttem;
                *(bright+i)=brttem;
                }
            else {
                *(bright+i)=-999.;
                }
            if( (i % 500) == 0 )
                printf(" pixel:%i raw:%i radiance:%10.2f bright:%10.2f\n",
                       i,(int)*(irim+i),radiance,*(bright+i));
    
            }
        }
    else if ( ibird == 11 ){
        printf( " fk1:%10.2f  fk2:%10.2f  tc1:%8.4f  tc2:%8.4f\n",
            g11fk1[index],g11fk2[index],g11tc1[index],g11tc2[index]);
        for ( i=0 ; i < npixel ; i++ ){
            if(*(irim+i) > 0 ) {
                rawval=(double)(*(irim+i));
                radiance=gaininv*(rawval-bias);
                expon=(g11fk1[index]/radiance) + 1.;
                tt=g11fk2[index]/log(expon);
                brttem=(tt-g11tc1[index])/
                           g11tc2[index];
                if( brttem < 0.0 ) 
                    brttem = 0.0;
                if( brttem > brtmax )
                    brtmax = brttem;
                if( brttem < brtmin )
                    brtmin = brttem;
                *(bright+i)=brttem;
                }
            else {
                *(bright+i)=-999.;
                }
            if( (i % 500) == 0 )
                printf(" pixel:%i raw:%i radiance:%10.2f bright:%10.2f\n",
                       i,(int)*(irim+i),radiance,*(bright+i));
    
            }
        }
    else if ( ibird == 12 ){
        if ( index == 5 ) index = 4;
        printf( " fk1:%10.2f  fk2:%10.2f  tc1:%8.4f  tc2:%8.4f\n",
            g12fk1[index],g12fk2[index],g12tc1[index],g12tc2[index]);
        for ( i=0 ; i < npixel ; i++ ){
            if(*(irim+i) > 0 ) {
                rawval=(double)(*(irim+i));
                radiance=gaininv*(rawval-bias);
                expon=(g12fk1[index]/radiance) + 1.;
                tt=g12fk2[index]/log(expon);
                brttem=(tt-g12tc1[index])/
                           g12tc2[index];
                if( brttem < 0.0 ) 
                    brttem = 0.0;
                if( brttem > brtmax )
                    brtmax = brttem;
                if( brttem < brtmin )
                    brtmin = brttem;
                *(bright+i)=brttem;
                }
            else {
                *(bright+i)=-999.;
                }
            if( (i % 500) == 0 )
                printf(" pixel:%i raw:%i radiance:%10.2f bright:%10.2f\n",
                       i,(int)*(irim+i),radiance,*(bright+i));
    
            }
        }
    else {
        printf(" That satellite is unknown. ibird = %d\n",ibird);
        printf(" Time to update tables in calib.c \n");
        return(-1);
        }
    printf(" Minimum temp: %8.2f K   Maximum temp: %8.2f K\n",
             brtmin,brtmax);
    printf(" Minimum temp: %8.2f C   Maximum temp: %8.2f C\n",
             (brtmin-273.15),(brtmax-273.15));
    return(0);
}

/* cnt2albedo
   From 1-byte VIS counts estimate albedo */

int cnt2albedo( visim, albedo, npixel )
    long *visim;
    float *albedo;
    long npixel;
{
    int i;
    float scaleinv,albtem,albmin,albmax;
    long imgcnt,imgmin,imgmax;

    imgmin=9999;
    imgmax=-9999;
    albmin=999.;
    albmax=0.;
    scaleinv = 1./255.;
    for ( i = 0 ; i < npixel ; i++ ) {
        imgcnt=*(visim+i);
        albtem=scaleinv*(float)(*(visim+i));
        albtem=albtem*albtem;
        if ( albtem < 0. )
                albtem = 0.;
        if ( albtem > 1. )
                albtem = 1.0;
        if ( imgcnt < imgmin )
              imgmin = imgcnt;
        if ( imgcnt > imgmax )
              imgmax = imgcnt;
        if ( albtem < albmin )
              albmin = albtem;
        if ( albtem > albmax )
              albmax = albtem;
        *(albedo+i)=albtem;
        }

    printf(" cnt2albedo: Minimum imgcnt: %ld   Maximum imgcnt: %ld\n",
             imgmin,imgmax);
    printf(" cnt2albedo: Minimum albedo: %8.4f   Maximum albedo: %8.4f\n",
             albmin,albmax);
    return(0);
}

/* cnt2bright
   From 1-byte BRIT counts compute brightness temperature */

int cnt2bright( irim, bright, npixel )
    long *irim;
    float *bright;
    long npixel;
{
    static long c1    = 418;
    static long c2    = 660;
    static long irthr = 176;

    int i;
    float brttem,brtmin,brtmax;

    brtmin=999.;
    brtmax=0.;
    for ( i = 0 ; i < npixel ; i++ ) {
        if ( *(irim+i) > irthr ) {
           brttem=(float)(c1 - *(irim+i));
           if ( brttem < brtmin )
              brtmin = brttem;
           if ( brttem > brtmax )
              brtmax = brttem;
           *(bright+i)=brttem;
           }
        else if ( *(irim+i) > 0 ) {
           brttem=0.5 * (float)(c2 - *(irim+i));
           if ( brttem < 0. )
                brttem = 0.;
           if ( brttem < brtmin )
              brtmin = brttem;
           if ( brttem > brtmax )
              brtmax = brttem;
           *(bright+i)=brttem;
           }
        else
           *(bright+i)=-999.0;
        }

    printf(" Minimum temp: %8.2f K   Maximum temp: %8.2f K\n",
             brtmin,brtmax);
    printf(" Minimum temp: %8.2f C   Maximum temp: %8.2f C\n",
             (brtmin-273.15),(brtmax-273.15));
    return(0);
}
