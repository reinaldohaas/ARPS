/*######################################################################
  ######################################################################
  ########                                                      ########
  ########           C reader for ARPS HDF4 files               ########
  ########                                                      ########
  ######################################################################
  ######################################################################

!------------------------------------------------------------------------
!
!  PURPOSE:
!    To read ARPS HDF 4 history files.
!
!  Usage:
!
!    To compile:
!      $> h4cc -o readarps readarps.c
!
!    To run:
!      $> readarps [options] ARPS_HDF4_file_name
!
!      options are:
!        -n vname: print variable named as "vname"
!        -v      : verbose message
!        -h      : for help
!
!    To do meaningful job:
!      Please add your own code in the dataset loop below
!
!-----------------------------------------------------------------------
!
! Author:
!   Yunheng Wang (08/01/2005)
!
!----------------------------------------------------------------------
 */

#include "mfhdf.h"         // HDF header

/*#######################################################
 *
 * Read an integer attribute in the file
 *
 * #####################################################*/

int hdfrdi(int32 sd_id, char *attname, int32 *file_data) {
  int32  attr_index;
  intn   status;

  attr_index = SDfindattr(sd_id,attname);
  status     = SDreadattr(sd_id,attr_index,file_data);

  if (status < 0) {
    printf("\nERROR: attribute \"%s\" not found.\n",attname);
    *file_data = 1;
  }
  return status;
}

/*#######################################################
 *
 * Read an character string in the file
 *
 * #####################################################*/

int hdfrdc(int32 sd_id, char *attname, int8 *file_data) {
  int32  attr_index;
  intn   status;

  attr_index = SDfindattr(sd_id,attname);
  status     = SDreadattr(sd_id,attr_index,file_data);

  if (status < 0) {
    printf("\nERROR: attribute \"%s\" not found.\n",attname);
    sprintf(file_data,"Unknown");
  }
  return status;
}

/*#######################################################
 *
 * Read a float attribute in the file
 *
 * #####################################################*/

int hdfrdr(int32 sd_id, char *attname, float32 *file_data) {
  int32  attr_index;
  intn   status;

  attr_index = SDfindattr(sd_id,attname);
  status     = SDreadattr(sd_id,attr_index,file_data);

  if (status < 0) {
    printf("\nERROR: attribute \"%s\" not found.\n",attname);
    exit(0);
  }
  return status;
}

/*********************************************************
 *
 * Global variables for command line arguments
 *
 * ******************************************************/

extern char *optarg;
extern int optind;

/*==============================================================
 *
 * Main program begin here
 *
 * ============================================================= */

int main(int argc, char * argv[])
{

  /****************** Variable declaration *****************************/

  char    filename[128];
  char    runname[80],sizestr[80];
  int32   sd_id, sds_id;
  int32   n_datasets, n_file_attrs;
  int32   nx,ny,nz,nzsoil,nstyp,ns;

  int32   dindex;
  char    dname[MAX_NC_NAME];
  int32   rank, data_type, n_attrs;
  int32   dim_sizes[MAX_VAR_DIMS];
  intn    istatus;
  int     n,k,j,i,jd,kd,nd;
  int     ireturn;

  int32   start[MAX_VAR_DIMS],stride[MAX_VAR_DIMS];

  int32     bufsize;
  float32   amax, amin, fscale;
  float32 * buf_r;
  int32   * buf_i;
  int16   * bufc_i;
  float32 * buf_max, *buf_min;

  int   verbose;
  char  c, vname[40];

  /*@@@@@@@@@@@@@@@@@ Begin of executable code @@@@@@@@@@@@@@@@@@@@@@@*/

  ireturn = 0;
  verbose = 0;

  /****************** Get command line arguments **********************/

  while ((c = getopt(argc, argv, "n:hv")) != EOF) {
    switch (c) {
      case 'n': 
	sscanf(optarg,"%s",vname);
	break;
      case 'v':
	verbose = 1;
	break;
      case 'h':
      case '?':
	printf("\nUsage: %s [options] [ARPS HDF4 file name]\n",argv[0]);
	printf("\n    Options:\n");
	printf("       -n vname: Extract variable with name as \"vname\"\n");
	printf("       -v      : Print vebose message\n");
	printf("       -h      : Print usage and exit\n");
	printf("\n");
	return ireturn;
    }
  }
  
  if ( optind < argc ) {
    strcpy(filename,argv[optind]);
  } else {
    printf("Please enter the ARPS HDF4 history file name: ");
    scanf("%s",filename);
  }

  if(verbose) printf("\nREADARPS is about to read file: %s  ......\n",
                     filename);

  /****************** Open data file  ********************************/

  sd_id = SDstart(filename,DFACC_READ);
  if (sd_id == FAIL) {
    printf("ERROR: open HDF 4 file error. Program stopped.\n");
    return sd_id;
  }

  /****************** Get file attributes  ****************************/

  istatus = hdfrdc(sd_id,"runname", runname);
  if (verbose) printf("Runname in this file is: %s.\n",runname);

  istatus = hdfrdi(sd_id,"nx",    &nx);
  istatus = hdfrdi(sd_id,"ny",    &ny);
  istatus = hdfrdi(sd_id,"nz",    &nz);
  istatus = hdfrdi(sd_id,"nzsoil",&nzsoil);
  istatus = hdfrdi(sd_id,"nstyp", &nstyp);
  if (verbose) printf("Size of ARPS data are: nx = %d, ny = %d, nz = %d,"
         " nzsoil = %d, nstyp = %d.\n",nx,ny,nz,nzsoil,nstyp);

  ns = nstyp + 1;

  /****************** Allocate buffers  *******************************/

  bufsize = nx*ny*nz>nx*ny*nzsoil*ns?nx*ny*nz:nx*ny*nzsoil*ns;

  buf_r  = (float32 *) malloc(bufsize*sizeof(float32));
  buf_i  = (int32   *) malloc(bufsize*sizeof(int32));
  bufc_i = (int16   *) malloc(bufsize*sizeof(int16));

  buf_max = (float32 *) malloc(nz*sizeof(float32));
  buf_min = (float32 *) malloc(nz*sizeof(float32));

  /****************** Loop over all datasets  *************************/

  for (n=0;n<MAX_VAR_DIMS;n++) {
    start[n]  = 0;
    stride[n] = 1;
  }

  istatus = SDfileinfo(sd_id, &n_datasets, &n_file_attrs);
  if (verbose) printf("\nThere are totally %d datasets in the file\n",n_datasets);

  for(dindex=0;dindex<n_datasets;dindex++) {

    sds_id = SDselect(sd_id,dindex);

    istatus = SDgetinfo(sds_id, dname, &rank, dim_sizes, &data_type,
	                &n_attrs);
    if(verbose) {
      printf("\nREADARPS is reading variable: %s ",dname);
    
      if (rank == 1) {
        sprintf(sizestr,"%d",         dim_sizes[0]);
      } else if (rank == 2) {
        sprintf(sizestr,"%dx%d",      dim_sizes[0],dim_sizes[1]);
      } else if (rank == 3) {
        sprintf(sizestr,"%dx%dx%d",   dim_sizes[0],dim_sizes[1],
    	                              dim_sizes[2]);
      } else if (rank == 4) {
        sprintf(sizestr,"%dx%dx%dx%d",dim_sizes[0],dim_sizes[1],
  	                              dim_sizes[2],dim_sizes[3]);
      } else {
        printf(". Unsupport rank - %d --- Skipped.\n",rank);
        continue;
      }
    }

    /***************** 32-bit floating point without compression ******/

    if (data_type == DFNT_FLOAT32) {
      if(verbose) printf("(%dD, 32-bit floating point, %s) into buf_r",rank,sizestr);
      istatus = SDreaddata(sds_id,start,stride,dim_sizes,buf_r);

    /***************** 32-bit integer without compression ******/

    } else if (data_type == DFNT_INT32) {
      if(verbose) printf("(%dD, 32-bit signed integer, %s) into buf_i",rank,sizestr);
      istatus = SDreaddata(sds_id,start,stride,dim_sizes,buf_i);

    /***************** Compressed data *************************/

    } else if (data_type == DFNT_INT16) {
      if(verbose) printf("(%dD, 16-bit signed integer, %s) into bufc_i.\n",rank,sizestr);
      istatus = SDreaddata(sds_id,start,stride,dim_sizes,bufc_i);

      if (rank == 2) {
        bufsize =  dim_sizes[0]*dim_sizes[1];

	istatus = hdfrdr(sds_id,"max",&amax);
	istatus = hdfrdr(sds_id,"min",&amin);

	fscale = (amax - amin) / 65534.;
	for (i=0; i<bufsize;i++) {
	  buf_r[i] = fscale*(bufc_i[i]+32767) + amin;
	}

      } else if (rank == 3) {
        bufsize =  dim_sizes[1]*dim_sizes[2];

	istatus = hdfrdr(sds_id,"max",buf_max);
	istatus = hdfrdr(sds_id,"min",buf_min);
	for (k=0; k<dim_sizes[0];k++) {
	  fscale = (buf_max[k]-buf_min[k]) / 65534.;
	  kd = k*bufsize;
	  for (i=0; i<bufsize; i++) {
	    buf_r[i+kd] = fscale*(bufc_i[i+kd]+32767) + buf_min[k];
	  }
	}
      } else if (rank == 4) {
        bufsize =  dim_sizes[2]*dim_sizes[3];

	istatus = hdfrdr(sds_id,"max",buf_max);
	istatus = hdfrdr(sds_id,"min",buf_min);
	for (n=0; n < dim_sizes[0]; n++) {
	  nd = n*bufsize*dim_sizes[1];
	  for (k=0; k < dim_sizes[1]; k++) {
	    fscale = (buf_max[k]-buf_min[k]) / 65534.;
	    kd = k*bufsize;
	    for (i=0; i<bufsize; i++) {
	      buf_r[i+kd+nd] = fscale*(bufc_i[i+kd+nd]+32767) + buf_min[k];
	    }
	  }
	}

      }
      if(verbose) printf("         bufc_i was depacked into buf_r.");

    /***************** Unknown data type, Maybe error *********/

    } else {
      printf(". Unknown data type. Program stopped.\n");
      ireturn = -1;
      break;
    }

    if (istatus == FAIL) {
      if(verbose) printf(" --- ERROR!!!. Program stopped.\n");
      ireturn = istatus;
      break;
    } else {
      if(verbose) printf(" --- DONE.\n");
    }

    istatus = SDendaccess(sds_id);

    /*%%%%%%%%%%%%% Add your code below to do meaningful job %%%%%%%%%%*/

    /*=================================================================
     *
     * All data is now in eithe buf_r or buf_i
     *
     * soiltyp, vegtyp are in buf_i
     * all other data  is  in buf_r.
     *
     * You should add your own code to extract data from buf_r. The 
     * following example can be used as a template.
     *
     * ================================================================*/

    /*============= An example to print 3D 32-bit float  ==============
     *
     * Fortran ARPS array     C ARPS array    Data in buffer
     * ------------------     -------------   ---------------------
     *  var(i+1,j+1,k+1)   =  var[k][j][i]  = buf_r[i+j*nx+k*nx*ny]
     *
     *  where:
     *    i = 0 .. (nx-1)
     *    j = 0 .. (ny-1)
     *    k = 0 .. (nz-1)
     *  
     * =============================================================== */ 

    if (strcmp(dname,vname) == 0) {    // print one specific varaibe

      if (rank == 1) {
        nx = dim_sizes[0];
        ny = 1;
        nz = 1;
        ns = 1;
      } else if (rank == 2) {
        nx = dim_sizes[1];  // west-east
        ny = dim_sizes[0];  // south-north
        nz = 1;
        ns = 1;
      } else if (rank == 3) {
        nx = dim_sizes[2];  // west-east
        ny = dim_sizes[1];  // south-north
        nz = dim_sizes[0];  // bottom-top
        ns = 1;
      } else if (rank == 4) {
        nx = dim_sizes[3];
        ny = dim_sizes[2];
        nz = dim_sizes[1];  // soil depth
        ns = dim_sizes[0];  // number of soil type
      } else {
        continue;
      }
  
      printf("Dataset Name = %s\n",dname); 
      for (n=0;n<ns;n++) {
        nd = n*nx*ny*nz;
        for (k=0;k<nz;k++) {
          kd = k*nx*ny;
          for (j=0;j<ny;j++) {
            jd = j*nx;
            printf("%4d,%4d: ",k,j);
            for (i=0;i<nx;i++) {
              printf("%15.7f",buf_r[i+jd+kd+nd]);
            }
            printf("\n");
          }
          printf("\n");
        }
      }

    }

  } 

  /****************** Clear memory space before retrun **********************/

  free(buf_r);
  free(buf_i);
  free(bufc_i);
  free(buf_max);
  free(buf_min);

  istatus = SDend(sd_id);

  return ireturn;
}
