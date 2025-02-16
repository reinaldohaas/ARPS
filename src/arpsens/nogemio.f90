!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE initgemio                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE initgemio(nxlg,nylg,mapproj,trulat1,trulat2,trulon,          &
                     lat1,lon1,lat2,lon2,iret)

!  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! PURPOSE: Initial IO interface for GEMPAK package
!
!-----------------------------------------------------------------------
!
! Author: Yunheng Wang (04/02/2007)
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: nxlg, nylg
  INTEGER, INTENT(IN)  :: mapproj
  REAL,    INTENT(IN)  :: trulat1, trulat2, trulon
  REAL,    INTENT(IN)  :: lat1, lon1, lat2, lon2
  INTEGER, INTENT(OUT) :: iret
 
!-----------------------------------------------------------------------
!
!  GEMPAK variables
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=3), PARAMETER :: cproj(0:3) = (/'CED','STR','LCC','MER'/)

  CHARACTER (LEN=72) :: gdarea
  CHARACTER (LEN=72) :: chproj

!  REAL, PARAMETER :: deltan = 1.0
!  REAL, PARAMETER :: deltax = -9999., deltay = -9999.
  LOGICAL, PARAMETER :: angflg = .TRUE.

  REAL :: angle1,angle2,angle3
  REAL :: gbnds(4)

  LOGICAL :: notopn              ! should be deleted


  iret = 1

  WRITE(6,'(/1x,a,/)') 'GEMPAK library is not linked. Please use option "-io gempak" when compiling.'
  RETURN
END SUBROUTINE initgemio
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE enswrtgem                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE enswrtgem (nx,ny, filename, lens, nchannel,                  &
                    ctime,year,month,day,hour,minute,                   &
                    iprgem,nprgem,ihtgem,nhtgem,ilite,iltgci,icrtm,icitm,&
                    ibeg_offset,iend_offset,jbeg_offset,jend_offset,    &
                    hgtp,uwndp,vwndp,wwndp,tmpp,sphp,                   &
                    uwndh,vwndh,wwndh,                                  &
                    psf,mslp,vort500,temp2m,dewp2m,qv2m,u10m,v10m,      &
                    tempk2,dewpk2,qvk2,hgtsfc,refl1km,refl4km,cmpref,   &
                    accppt,convppt,cape,mcape,cin,mcin,lcl,             &
                    srh01,srh03,uh25,sh01,sh06,thck,li,brn,pw,vvelmax,  &
                    vvelmin,echotop,colqh,                              &
                    wspd10max,w_up_max,w_dn_max,refd_max,up_heli_max,grpl_max,&
                    wspd1kmmax,crefd_max,pblh,raingrpl,rainhail,up_heli16_max,&
                    refdm10c_max,wpbl,lfc,sat_bt1,sat_bt2, &
                    bku,bkv,bkspd,lllr,lr75,mlcape,mlcin,               &
                    snow_max,grpl05_max,up_heli_maxe,up_heli_maxp,      &
                    tmp1km,dwp1km,pre1km,                               &
                    ltg1_max,ltg2_max,ltg3_max)

  INTEGER :: nx,ny
  CHARACTER (LEN=*) :: filename
  CHARACTER (LEN=72) :: gdfile,gdarea
  INTEGER :: lens, nchout
  INTEGER :: imxgrd
  PARAMETER (imxgrd=9999)
  INTEGER :: ivsfc,ivprs,ivtheta,ivhgt
  PARAMETER (ivsfc=0,ivprs=1,ivtheta=2,ivhgt=3)
  INTEGER :: level(2)
  INTEGER :: ighdr(2)
  INTEGER :: ipktyp, nbits
  PARAMETER (ipktyp=1,nbits=16)
  CHARACTER (LEN=20) :: gemtime(2)
  CHARACTER (LEN=12) :: parm

  INTEGER :: nprgem,nhtgem
  INTEGER :: iprgem(nprgem),ihtgem(nhtgem)
  INTEGER :: ilite,iltgci,icrtm,icitm
  INTEGER :: ibeg_offset,iend_offset,jbeg_offset,jend_offset

  REAL :: hgtp(nx,ny,nprgem),uwndp(nx,ny,nprgem),vwndp(nx,ny,nprgem)
  REAL :: wwndp(nx,ny,nprgem),tmpp(nx,ny,nprgem),sphp(nx,ny,nprgem)
  REAL :: uwndh(nx,ny,nhtgem),vwndh(nx,ny,nhtgem),wwndh(nx,ny,nhtgem)
  REAL :: psf(nx,ny),mslp(nx,ny),vort500(nx,ny)
  REAL :: temp2m(nx,ny),dewp2m(nx,ny),qv2m(nx,ny)
  REAL :: tempk2(nx,ny),dewpk2(nx,ny),qvk2(nx,ny)
  REAL :: u10m(nx,ny),v10m(nx,ny),hgtsfc(nx,ny)
  REAL :: refl1km(nx,ny),refl4km(nx,ny),cmpref(nx,ny)
  REAL :: accppt(nx,ny),convppt(nx,ny)
  REAL :: cape(nx,ny),mcape(nx,ny),cin(nx,ny),mcin(nx,ny),lcl(nx,ny)
  REAL :: srh01(nx,ny),srh03(nx,ny),uh25(nx,ny),sh01(nx,ny),sh06(nx,ny)
  REAL :: thck(nx,ny),li(nx,ny),brn(nx,ny),pw(nx,ny),vvelmax(nx,ny)
  REAL :: vvelmin(nx,ny),echotop(nx,ny),colqh(nx,ny)
  REAL :: wspd10max(nx,ny),w_up_max(nx,ny),w_dn_max(nx,ny)
  REAL :: refd_max(nx,ny),up_heli_max(nx,ny),grpl_max(nx,ny)
  REAL :: ltg1_max(nx,ny),ltg2_max(nx,ny),ltg3_max(nx,ny)
  REAL :: wspd1kmmax(nx,ny),crefd_max(nx,ny),pblh(nx,ny)
  REAL :: raingrpl(nx,ny),rainhail(nx,ny),up_heli16_max(nx,ny)
  REAL :: wpbl(nx,ny),lfc(nx,ny)
  REAL :: refdm10c_max(nx,ny)
  REAL :: sat_bt1(nx,ny,nchannel),sat_bt2(nx,ny,nchannel)
  REAL :: bku(nx,ny),bkv(nx,ny),bkspd(nx,ny),lllr(nx,ny),lr75(nx,ny)
  REAL :: mlcape(nx,ny),mlcin(nx,ny),snow_max(nx,ny),grpl05_max(nx,ny)
  REAL :: up_heli_maxe(nx,ny),up_heli_maxp(nx,ny)
  REAL :: tmp1km(nx,ny),dwp1km(nx,ny),pre1km(nx,ny)

  INTEGER :: i,j,k,klev,iret,igdfil
  INTEGER :: year,month,day,hour,minute
  INTEGER :: iyr,ifhr,ifmin,ifile,ihd
  REAL :: ctime

  LOGICAL :: wrtflg,rewrite,notopn,gemExist


  WRITE(6,'(/1x,a,/)') 'GEMPAK library is not linked. Please use option "-io gempak" when compiling.'
  RETURN
END SUBROUTINE enswrtgem
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE enswrtgem2                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE enswrtgem2 (nx,ny, filename, lens,                              &
                    ctime,year,month,day,hour,minute,                      &
                    ilite1,icrtm,icitm,                                    &
                    ibeg_offset,iend_offset,jbeg_offset,jend_offset,       &
                    hgt250m,hgt500m,hgt850m,vort500m,                      &
                    u250m,u500m,u850m,v250m,v500m,v850m,                   &
                    mslpm,accpptm,acc03m,acc06m,temp2mm,dewp2mm,qv2mm,     &
                    pwatm,u10mm,v10mm,sbcapem,sbcinsm,sblclm,              &
                    mlcapem,mlcinsm,mucapem,mucinsm,lllrm,lr75m,           &
                    shr01m,shr06m,bkum,bkvm,bkspdm,                        &
                    accppt_max,acc03_max,acc06_max,uh_max_max,             &
                    wup_max_max,wdn_max_min,cqgmax_max,etop_max,           &
                    wspmax_max,refmax_max,                                 &
                    wsp1km_max,uh03mx_max,grpl01_max,grpl03_max,           &
                    wsp1km_3hmax,uh03mx_3hmax,                             &
                    uh_max_3hmax,wup_max_3hmax,cqgmax_3hmax,wspmax_3hmax,  &
                    llqg_max,mslpt,hgt500t,                                &
                    accppt_pm,acc03_pm,acc06_pm,                           &
                    refl1km_pm,refl4km_pm,cmpref_pm,                       &
                    p01m025_p,p01m050_p,p01m100_p,                         &
                    p06m050_p,p06m100_p,p06m200_p,                         &
                    sbcape0500_p,sbcape1000_p,sbcape2000_p,                &
                    mlcape0500_p,mlcape1000_p,mlcape2000_p,                &
                    mucape0500_p,mucape1000_p,mucape3000_p,                &
                    sbcins100_p,sbcins050_p,sbcins025_p,                   & 
                    mlcins100_p,mlcins050_p,mlcins025_p,                   &
                    mucins100_p,mucins050_p,mucins025_p,                   &
                    shr01_10_p,shr01_15_p,shr01_20_p,                      &
                    shr06_15_p,shr06_20_p,shr06_25_p,                      &
                    stp1_p,stp3_p,stp5_p,scp1_p,scp3_p,scp9_p,bkspd15_p,   &
                    sat_bt1_32_p,sat_bt1_52_p,sat_bt2_32_p,sat_bt2_52_p,   &
                    refl1km40_pn,refl4km40_pn,cmpref40_pn,                 &
                    uh_max025_pn,uh_max050_pn,uh_max100_pn,                &
                    wup_max10_pn,wup_max15_pn,wup_max20_pn,                &
                    wdn_max02_pn,wdn_max06_pn,wdn_max10_pn,                &
                    cqgmax20_pn,cqgmax30_pn,cqgmax40_pn,                   &
                    wspmax15_pn,wspmax20_pn,wspmax25_pn,                   &
                    llqg20_pn,refl1km3h40_pn,                              &
                    uh_3hmax025_pn,uh_3hmax050_pn,uh_3hmax100_pn,          &
                    wup_3hmax10_pn,wup_3hmax15_pn,wup_3hmax20_pn,          &
                    cqg3hmax20_pn,cqg3hmax30_pn,cqg3hmax40_pn,             &
                    wsp3hmax15_pn,wsp3hmax20_pn,wsp3hmax25_pn,             &
                    ltg3_05_pn,ltg3_3_pn,ltg3_6_pn,                        &
                    ltg3_002_pn,ltg3_1_pn,                                 &
                    etop25_pn,etop35_pn,etop45_pn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out ensemble product 2D fields for SPC/NSSL HWT in gempak format.
!  Called from ensemble post-processing program ens_ana.f90
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Fanyou Kong
!  04/04/2010
!-----------------------------------------------------------------------
!
!  IMPLICIT NONE

  INTEGER :: nx,ny
  CHARACTER (LEN=*)  :: filename
  CHARACTER (LEN=72) :: gdfile
  INTEGER :: lens, nchout
  INTEGER :: imxgrd
  PARAMETER (imxgrd=500)
  INTEGER :: ivsfc,ivprs,ivtheta,ivhgt
  PARAMETER (ivsfc=0,ivprs=1,ivtheta=2,ivhgt=3)
  INTEGER :: level(2)
  INTEGER :: ighdr(2)
  INTEGER :: ipktyp, nbits
  PARAMETER (ipktyp=1,nbits=16)
  CHARACTER (LEN=20) :: gemtime(2)
  CHARACTER (LEN=12) :: parm

  INTEGER :: ilite1,icrtm,icitm
  INTEGER :: ibeg_offset,iend_offset,jbeg_offset,jend_offset
  INTEGER :: isec

  REAL :: hgt250m(nx,ny),hgt500m(nx,ny),hgt850m(nx,ny),vort500m(nx,ny)
  REAL :: u250m(nx,ny),u850m(nx,ny),v250m(nx,ny),v850m(nx,ny)
  REAL :: u500m(nx,ny),v500m(nx,ny)
  REAL :: mslpm(nx,ny),accpptm(nx,ny),acc03m(nx,ny),acc06m(nx,ny)
  REAL :: temp2mm(nx,ny),dewp2mm(nx,ny),qv2mm(nx,ny),pwatm(nx,ny)
  REAL :: u10mm(nx,ny),v10mm(nx,ny)
  REAL :: shr01m(nx,ny),shr06m(nx,ny)
  REAL :: sbcapem(nx,ny),sbcinsm(nx,ny),sblclm(nx,ny)
  REAL :: mlcapem(nx,ny),mlcinsm(nx,ny)
  REAL :: mucapem(nx,ny),mucinsm(nx,ny)
  REAL :: lllrm(nx,ny),lr75m(nx,ny)
  REAL :: bkum(nx,ny),bkvm(nx,ny),bkspdm(nx,ny)
  REAL :: accppt_max(nx,ny),acc03_max(nx,ny),acc06_max(nx,ny)
  REAL :: uh_max_max(nx,ny),wup_max_max(nx,ny),wdn_max_min(nx,ny)
  REAL :: cqgmax_max(nx,ny),etop_max(nx,ny),wspmax_max(nx,ny)
  REAL :: refmax_max(nx,ny),wsp1km_max(nx,ny),uh03mx_max(nx,ny)
  REAL :: grpl01_max(nx,ny),grpl03_max(nx,ny)
  REAL :: wsp1km_3hmax(nx,ny),uh03mx_3hmax(nx,ny),uh_max_3hmax(nx,ny)
  REAL :: wup_max_3hmax(nx,ny),cqgmax_3hmax(nx,ny),wspmax_3hmax(nx,ny)
  REAL :: llqg_max(nx,ny),mslpt(nx,ny),hgt500t(nx,ny)
  REAL :: accppt_pm(nx,ny),acc03_pm(nx,ny),acc06_pm(nx,ny)
  REAL :: refl1km_pm(nx,ny),refl4km_pm(nx,ny),cmpref_pm(nx,ny)
  REAL :: p01m025_p(nx,ny),p01m050_p(nx,ny),p01m100_p(nx,ny)
  REAL :: p06m050_p(nx,ny),p06m100_p(nx,ny),p06m200_p(nx,ny)
  REAL :: sbcape0500_p(nx,ny),sbcape1000_p(nx,ny),sbcape2000_p(nx,ny)
  REAL :: mlcape0500_p(nx,ny),mlcape1000_p(nx,ny),mlcape2000_p(nx,ny)
  REAL :: mucape0500_p(nx,ny),mucape1000_p(nx,ny),mucape3000_p(nx,ny)
  REAL :: sbcins100_p(nx,ny),sbcins050_p(nx,ny),sbcins025_p(nx,ny)
  REAL :: mlcins100_p(nx,ny),mlcins050_p(nx,ny),mlcins025_p(nx,ny)
  REAL :: mucins100_p(nx,ny),mucins050_p(nx,ny),mucins025_p(nx,ny)
  REAL :: shr01_10_p(nx,ny),shr01_15_p(nx,ny),shr01_20_p(nx,ny)
  REAL :: shr06_15_p(nx,ny),shr06_20_p(nx,ny),shr06_25_p(nx,ny)
  REAL :: stp1_p(nx,ny),stp3_p(nx,ny),stp5_p(nx,ny)
  REAL :: scp1_p(nx,ny),scp3_p(nx,ny),scp9_p(nx,ny),bkspd15_p(nx,ny)
  REAL :: sat_bt1_32_p(nx,ny),sat_bt1_52_p(nx,ny)
  REAL :: sat_bt2_32_p(nx,ny),sat_bt2_52_p(nx,ny)
  REAL :: refl1km40_pn(nx,ny),refl4km40_pn(nx,ny),cmpref40_pn(nx,ny)
  REAL :: uh_max025_pn(nx,ny),uh_max050_pn(nx,ny),uh_max100_pn(nx,ny)
  REAL :: wup_max10_pn(nx,ny),wup_max15_pn(nx,ny),wup_max20_pn(nx,ny)
  REAL :: wdn_max02_pn(nx,ny),wdn_max06_pn(nx,ny),wdn_max10_pn(nx,ny)
  REAL :: cqgmax20_pn(nx,ny),cqgmax30_pn(nx,ny),cqgmax40_pn(nx,ny)
  REAL :: wspmax15_pn(nx,ny),wspmax20_pn(nx,ny),wspmax25_pn(nx,ny)
  REAL :: llqg20_pn(nx,ny),refl1km3h40_pn(nx,ny)
  REAL :: uh_3hmax025_pn(nx,ny),uh_3hmax050_pn(nx,ny),uh_3hmax100_pn(nx,ny)
  REAL :: wup_3hmax10_pn(nx,ny),wup_3hmax15_pn(nx,ny),wup_3hmax20_pn(nx,ny)
  REAL :: cqg3hmax20_pn(nx,ny),cqg3hmax30_pn(nx,ny),cqg3hmax40_pn(nx,ny)
  REAL :: wsp3hmax15_pn(nx,ny),wsp3hmax20_pn(nx,ny),wsp3hmax25_pn(nx,ny)
  REAL :: etop25_pn(nx,ny),etop35_pn(nx,ny),etop45_pn(nx,ny)
  REAL :: ltg3_05_pn(nx,ny),ltg3_3_pn(nx,ny),ltg3_6_pn(nx,ny)
  REAL :: ltg3_002_pn(nx,ny),ltg3_1_pn(nx,ny)

  INTEGER :: i,j,k,klev,iret,igdfil
  INTEGER :: year,month,day,hour,minute
  INTEGER :: iyr,ifhr,ifmin,ifile,ihd
  REAL :: ctime

  LOGICAL :: wrtflg,rewrite,notopn,gemExist

WRITE(6,'(/1x,a,/)') 'GEMPAK library is not linked. Please use option "-io gempak" when compiling.'

  RETURN
END SUBROUTINE enswrtgem2
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GEMDUMP2D                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE gemdump2d(igdfil,qfield, nx, ny,                          &
                     ibeg_offset,iend_offset,jbeg_offset,jend_offset,&
                     ighdr, gemtime, level, ivflag, parm,            &
                     rewrite, ipktyp, nbits, iret)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out one single 2D field in gempak format. This interface is
!  specially written for MPI dump to a joined GEMPAK file
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Fanyou Kong
!  02/20/2007
!-----------------------------------------------------------------------
!

  IMPLICIT NONE

  INCLUDE 'mp.inc'

  INTEGER ::  igdfil, iret, istat
  INTEGER ::  nx,ny
  INTEGER ::  ibeg_offset,iend_offset,jbeg_offset,jend_offset
  INTEGER :: ivflag
  INTEGER :: level(2)
  INTEGER :: ighdr(2)
  CHARACTER (LEN=20) :: gemtime(2)
  CHARACTER (LEN=12) :: parm
  REAL :: qfield(nx,ny)

  LOGICAL :: rewrite
  integer nbits, ipktyp

  INTEGER :: nxlg, nylg
  REAL,    ALLOCATABLE :: out2d(:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begining of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(6,'(/1x,a,/)') 'GEMPAK library is not linked. Please use option "-io gempak" when compiling.'

  RETURN
END SUBROUTINE gemdump2d
