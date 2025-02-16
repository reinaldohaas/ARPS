SUBROUTINE filgrd( mptr,msrc )
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'

  PARAMETER (lismax=200)
  INTEGER :: lisvar(lismax)
!
!-----------------------------------------------------------------------
!
!  This subroutine sets the data fields on the fine grid
!  mptr by interpolation from grid msrc.
!
!  The constants are set elsewhere.
!
!  Interpolation routines are called for each field that is to be
!  interpolated from coarse grid.
!
!  Whether a particular array needs to be initialized by interpolation,
!  is determined by parameter inixy and inixyz for 2-D and 3-D array
!  respectively.
!
!  Written by B. Skamarock and generalized by M. Xue.
!
!-----------------------------------------------------------------------
!
  IF(verbose6) WRITE(6,'(1x,a,i3,a,i3)')                                &
      'To call FILGRD to interpolate fields from grid ',msrc,           &
      ' to grid ',mptr
!
!-----------------------------------------------------------------------
!
!  Initialize 2d xy scalar arrays by interpolation
!
!-----------------------------------------------------------------------
!
  ndim = 2

  CALL inisvlst( lisvar, listno, ndim )

  IF(verbose6) THEN
    WRITE(6,'(1x,i4,i2,a)') listno,ndim,                                &
        '-D scalar variables selected to be initialized:'
    DO i = 1,listno
      WRITE(6,'(1x,a,i4)') 'ivar=',lisvar(i)
    END DO
  END IF

  IF(listno /= 0) CALL inivar( mptr,msrc,lisvar,listno,ndim )
!
!  set the t-dt fields to the values at time t
!  More accurate way is to do interpolation between t and t-dt fields.
!
  DO i=1,listno
    ivar = lisvar(i)
    ivar_dt = inixy(3,ivar)
    IF( ivar_dt /= 0 ) CALL cpyfld(mptr,ivar,ivar_dt,ndim)
  END DO
!
!-----------------------------------------------------------------------
!
!  Initialize 2d xy vector arrays by interpolation
!
!-----------------------------------------------------------------------
!
  CALL inivvlst( lisvar, listno, ndim )

  IF(verbose6) THEN
    WRITE(6,'(1x,i4,i2,a)') listno,ndim,                                &
        '-D vector variables selected to be initialized:'
    DO i = 1,listno
      WRITE(6,'(1x,a,i4)') 'ivar=',lisvar(i)
    END DO
  END IF


  IF(listno /= 0) CALL inivar( mptr,msrc,lisvar,listno,ndim )
!
!  set the t-dt fields to the values at time t
!  More accurate way is to do interpolation between t and t-dt fields.
!
  DO i=1,listno,2

    ivar = lisvar(i)
    ivar_dt = inixy(3,ivar)
    IF( ivar_dt /= 0 ) CALL cpyfld(mptr,ivar,ivar_dt,ndim)

    ivar = lisvar(i+1)
    ivar_dt = inixy(3,ivar)
    IF( ivar_dt /= 0 ) CALL cpyfld(mptr,ivar,ivar_dt,ndim)

  END DO

!
!-----------------------------------------------------------------------
!
!  Initialize 3d scalar arrays by interpolation
!
!-----------------------------------------------------------------------
!
  ndim = 3

  CALL inisvlst( lisvar, listno, ndim )

  IF(verbose6) THEN
    WRITE(6,'(1x,i4,i2,a)') listno,ndim,                                &
        '-D scalar variables selected to be initialized:'
    DO i = 1,listno
      WRITE(6,'(1x,a,i4)') 'ivar=',lisvar(i)
    END DO
  END IF

  IF(listno /= 0) CALL inivar( mptr,msrc,lisvar,listno,ndim )
!
!  set the t-dt fields to the values at time t
!  More accurate way is to do interpolation between t and t-dt fields.
!
  DO i=1,listno

    ivar = lisvar(i)
    ivar_dt = inixyz(3,ivar)
    IF( ivar_dt /= 0 ) CALL cpyfld(mptr,ivar,ivar_dt,ndim)

  END DO
!
!-----------------------------------------------------------------------
!
!  Initialize 3d vector arrays by interpolation
!
!-----------------------------------------------------------------------
!
  ndim = 3

  CALL inivvlst( lisvar, listno, ndim )

  IF(verbose6) THEN
    WRITE(6,'(1x,i4,i2,a)') listno,ndim,                                &
        '-D vector variables selected to be initialized:'
    DO i = 1,listno
      WRITE(6,'(1x,a,i4)') 'ivar=',lisvar(i)
    END DO
  END IF

  IF(listno /= 0) CALL inivar( mptr,msrc,lisvar,listno,ndim )
!
!  set the t-dt fields to the values at time t
!  More accurate way is to do interpolation between t and t-dt fields.
!
  DO i=1,listno,2

    ivar = lisvar(i)
    ivar_dt = inixyz(3,ivar)
    IF( ivar_dt /= 0 ) CALL cpyfld(mptr,ivar,ivar_dt,ndim)

    ivar = lisvar(i+1)
    ivar_dt = inixyz(3,ivar)
    IF( ivar_dt /= 0 ) CALL cpyfld(mptr,ivar,ivar_dt,ndim)

  END DO
!
!-----------------------------------------------------------------------
!
!  A number of vairables need to be initialized/defined by the user.
!  For ARPS, they include x,y,z,j1,j2,j3,wcont,rhostr and km.
!
!  Also many constants/parameters need to be set by the user for
!  the new grid.
!
!  These are done inside initngrd.
!
!-----------------------------------------------------------------------
!
  CALL initngrd(mptr,msrc)

  RETURN
END SUBROUTINE filgrd

SUBROUTINE inisvlst( lisvar, listno, ndim )
!
!-----------------------------------------------------------------------
!
!  Get a list of 2d xy or 3d xyz scalar arrays that
!  need initialization by interpolation in lisvar(listno)
!
!  Author: Ming Xue
!
!-----------------------------------------------------------------------
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'

  PARAMETER (lismax=200)
  INTEGER :: lisvar(lismax)

  IF( ndim == 2 ) THEN
!
    listno = 0
    DO ivar=1,nxy2d
      IF( inixy(1,ivar) == 1 .AND. inixy(2,ivar) == 0 ) THEN
        IF( listno+1 > lismax ) GO TO 999
        listno = listno + 1
        lisvar(listno) = ivar
      END IF
    END DO

  ELSE IF( ndim == 3 ) THEN

    listno = 0
    DO ivar=1,nxyz3d
      IF( inixyz(1,ivar) == 1 .AND. inixyz(2,ivar) == 0 ) THEN
        IF( listno+1 > lismax ) GO TO 999
        listno = listno + 1
        lisvar(listno) = ivar
      END IF
    END DO

  END IF

  RETURN

  999   WRITE(6,'(1x,a,a)') 'Array lisvar in INISVLST defined too small,', &
          'Job stopped.'
END SUBROUTINE inisvlst
!

SUBROUTINE inivvlst( lisvar, listno, ndim )
!
!-----------------------------------------------------------------------
!
!  Return a list of 2d xy or 3d xyz vector arrays that
!  need initialization by interpolation in lisvar(listno)
!
!-----------------------------------------------------------------------
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'

  PARAMETER (lismax=200)
  INTEGER :: lisvar(lismax)

  IF( ndim == 2 ) THEN
!
    listno = 0
    DO ivar=1,nxy2d
      IF( inixy(1,ivar) == 1 .AND. inixy(2,ivar) > 0 ) THEN
        IF( listno+2 > lismax ) GO TO 999
        listno = listno + 1
        lisvar(listno) = ivar
        listno = listno + 1
        lisvar(listno) = inixy(2,ivar)
      END IF
    END DO

  ELSE IF( ndim == 3 ) THEN

    listno = 0
    DO ivar=1,nxyz3d
      IF( inixyz(1,ivar) == 1 .AND. inixyz(2,ivar) > 0 ) THEN
        IF( listno+2 > lismax ) GO TO 999
        listno = listno + 1
        lisvar(listno) = ivar
        listno = listno + 1
        lisvar(listno) = inixyz(2,ivar)
      END IF
    END DO

  END IF


  RETURN

  999   WRITE(6,'(1x,a,a)') 'Array lisvar in INIVVLST defined too small,', &
          'Job stopped.'
  STOP

END SUBROUTINE inivvlst
