
c#######################################################################
c
c Update summary ( zzj, April, 1996 )
c   1. Including symbol font when call
c      PSfont(n)
c      n=5 corresponds to Symbol (Creek letter fond)
c   2. Modification was made on program setcolors to eliminate
C      some of unnecesary writing to PostCript io channel;
c      option 5 is available in setcolors
c#######################################################################

      SUBROUTINE xdevic
c-----------------------------------------------------------------------
c To be compatible with old programs. Newly created code should call
c xdevic_new directly.
c-----------------------------------------------------------------------
      CHARACTER(LEN=256) filename
      CALL xdevic_new(1,filename,0,13)
           ! filename will be ignored because len_filename = 0
      RETURN
      END

c#######################################################################
c
c Interface of ZXplot with Postscript.
c Created by Ming Xue at CIMMS/CAPS, Feburary, 1990-1996.
c
c#######################################################################
c
      SUBROUTINE xdevic_new(iwtype,filename,len_filename,iounit)
      INTEGER iwtype
      CHARACTER(LEN=256) filename  ! output graphical file name (base),
                                   ! the extension will be determined is always '.ps'.
      INTEGER len_filename         ! valid character length in filename
                                   ! = 0 to use the default "zxout.ps" or that set by xpsfn call
                                   !     It is for backward compatible.
      INTEGER iounit               ! A file unit for output, used only with PS interface.

      save nzxcal, mbordr
      data nzxcal /0/ ,mbordr/0/
      common /xpsize/ psize
      common /xpsd01/ xside, yside
      common /xoutch/ nch
      integer ncunique
      integer icoltable
      common /xcltbl/ ncunique,icoltable

      INTEGER outunit
c
c#######################################################################
c
c Set up device ( Paper for GHOST )
c
c#######################################################################
c
      nch = 6

      IF (len_filename > 0) THEN   ! else behavoir the same as before

        IF (iounit <= 0) THEN
          outunit = 13
        ELSE
          outunit = iounit
        END IF

        WRITE(nch,'(/,1x,3a,/)') '==== Creating Graphical file: ',
     +            filename(1:len_filename)//'.ps',' ===='

        CALL xpsfn(filename(1:len_filename)//'.ps',outunit)
      END IF

      CALL PSopn
      CALL xichar
      CALL xdspac( 1.0 )
      CALL xminit
      nzxcal=1

      CALL XSETCLRS(icoltable,0)
      CALL xcolor(1)
      CALL xthick(1)

      RETURN

c
c#######################################################################
c
c Advance plotting onto next frame
c
c#######################################################################
c
      ENTRY xframe
      CALL xqmap (  x1,x2,y1,y2 )
      CALL xlpnup(x1,y1)
      CALL PSfram
      RETURN
c
c#######################################################################
c
c Terminate graphic plotting
c
c#######################################################################
c
      ENTRY xgrend
      IF( nzxcal.eq.0)  then
        write(nch,'(1x,a)') 'XGREND called before device is set up.'
        RETURN
      ENDIF
      CALL xqmap (  x1,x2,y1,y2 )
      CALL xlpnup(x1,y1)

      CALL PScls
      RETURN

      ENTRY xfbord(mbord)
      mbordr=mbord
      RETURN
      END

      SUBROUTINE ppenup(x1,y1)
      CALL PSpnup(x1,y1)
      RETURN
      END

      SUBROUTINE ppendn(x,y)
      CALL PSpndn(x,y)
      RETURN
      END

c
c#######################################################################
c
c Define a normalized device (ND) space (0.0,XSIDE,0.0,YSIDE) on device
c provided. The size of this space is 'PSIZE' times of the max device
c space.  YSIDE should be 1.0, XSIDE can be bigger or smaller that 1.0
c depending on device.
c XSIDE, YSIDE should be passed to ZXPLOT through common block PSIDES.
c This routine sets up this space using GHOST.
c
c#######################################################################
c
      SUBROUTINE xdspac(psize)
      common /xpsize/ psize1
      common /xpsd01/ xside, yside
      common /xoutch/ nch
      data yside0 /1.0/
      save yside0

      psize1=psize

      xside=1.0
c
c#######################################################################
c
c To use only a squared plotting area, set yside=1.0:
c To make a full use of the US letter size plotting space, set yside=1.5:
c
c#######################################################################
c
      yside=yside0

      IF( psize.le.0.0) then
        write(nch,'(1x,a)')
     :  ' Zero or negative values for PSIZE not permitted!'
        write(nch,'(1x,a)')
     :  ' It should be in range of 0.1 to 1.0,please reset value'
        STOP
      ENDIF

      edgex=0.5*xside*(1.0/psize-1)
      edgey=edgex*yside/xside
      x1=-edgex
      x2= edgex+xside
      y1=-edgey
      y2= edgey+yside

      if(abs(yside-1.0).lt.0.001) then ! yside=1.0
c
c To use only a squared plotting area:
c
        call PSspac(50.0, 600.0, 50.0, 600.0, X1,X2,Y1,Y2)

      else ! yside=1.5
c
c To use full US letter size plotting area:
c
c       CALL PSspac(100.0,560.0,60.0,750.0,X1,X2,Y1,Y2)

        CALL PSspac(100.0,560.0,750.0-yside*460,750.0,X1,X2,Y1,Y2)

      endif

      return

      entry xpaprlnth( yside0a )

      yside0=yside0a

      RETURN
      END

      SUBROUTINE window(x1,x2,y1,y2)
      CALL PSwndw(x1,x2,y1,y2)
      RETURN
      END
c
c#######################################################################
c
c     Low level PostScript interface for ZXplot
c
c By Shian-Jiann Lin Feb.20, 1990 At CIMMS/CAPS, OU
c Modified and extended by Ming Xue to include font writing capablilities.
c
c This routine can only be called once.
c
c This routine opens the io & set up the min. header.
c
c#######################################################################
c
      SUBROUTINE PSopn
      common/psdef/io
      character char_io*132
      character*256 filename
      character psfn*(*)
      common /xoutch/ nch
      data filename /'zxout.ps'/
      data lfn /8/
      save filename,lfn
      integer iounit
      common /psbufferlines/ lines

      write(nch,'(/1x,a,a,a)')
     :  'PostScript output is in ',filename(1:lfn),'.'
      write(nch,'(1x,a,i2,a/)')
     :  'Data IO unit ',io,' to be used for writing the PS file.'

      open(unit=io,file=filename(1:lfn),
     :     status='unknown',form='formatted')

      lines = 0  ! set initial number of buffer lines to zero

      write(char_io,'(a)') '%!PS-Adobe-2.0'
      call write_ps(char_io)

      write(char_io,'(a)')'%%Title:ZX-PLOT'
      call write_ps(char_io)
      write(char_io,'(a)')'%%Pages:(atend)'
      call write_ps(char_io)
      write(char_io,'(a)')'%%DocumentFonts:(atend)'
      call write_ps(char_io)
      write(char_io,'(a)')'%%EndComments'
      call write_ps(char_io)

c
C Save the graphics state.
      write(char_io,'(a)')'gsave'
      call write_ps(char_io)
C
C Define Proc. here.......
C
      write(char_io,'(a)')'/w {setlinewidth} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/d {setdash} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/m {moveto} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/l {lineto} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/S {stroke} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/N {newpath} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/h {closepath} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/q {gsave} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/Q {grestore} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/g {setgray} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/ct 256 array def'
      call write_ps(char_io)
      write(char_io,'(a)') '/o {ct exch get aload pop setrgbcolor} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/f {fill} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/ff {findfont} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/scf {scalefont} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/stf {setfont} def'
      call write_ps(char_io)
      write(char_io,'(a)')'/rs { dup stringwidth pop'
      call write_ps(char_io)
      write(char_io,'(a)')'neg 0 rmoveto show } def'
      call write_ps(char_io)
      write(char_io,'(a)')'/cs { dup stringwidth pop'
      call write_ps(char_io)
      write(char_io,'(a)')'-0.5 mul 0 rmoveto show } def'
      call write_ps(char_io)
      write(char_io,'(a)')'/arrowdict 8 dict def'
      call write_ps(char_io)
      write(char_io,'(a)')'/arw %Procedure for an arrow'
      call write_ps(char_io)
      write(char_io,'(a)')'{arrowdict begin'
      call write_ps(char_io)
      write(char_io,'(a)')' /y4 exch def  /x4 exch def'
      call write_ps(char_io)
      write(char_io,'(a)')' /y3 exch def  /x3 exch def'
      call write_ps(char_io)
      write(char_io,'(a)')' /y2 exch def  /x2 exch def'
      call write_ps(char_io)
      write(char_io,'(a)')' /y1 exch def  /x1 exch def'
      call write_ps(char_io)
      write(char_io,'(a)')' S x1 y1 m x2 y2 l x3 y3 l S'
      call write_ps(char_io)
      write(char_io,'(a)')'   x4 y4 m x2 y2 l S'
      call write_ps(char_io)
      write(char_io,'(a)')' end } def '
      call write_ps(char_io)
      write(char_io,'(a)') '/landscape %transform to landscape layout '
      call write_ps(char_io)
      write(char_io,'(a)') '{792 0 translate 90 rotate'
      call write_ps(char_io)
      write(char_io,'(a)') '0 180 translate} def'
      call write_ps(char_io)
      write(char_io,'(a)')' %%EndProlog'
      call write_ps(char_io)
      write(char_io,'(a)') '%%Page: 1'
      call write_ps(char_io)
      write(char_io,'(a)') '%%set defaults:'
      call write_ps(char_io)
      write(char_io,'(a)') '%%landscape'
      call write_ps(char_io)
      write(char_io,'(a)') '/Helvetica ff 12 scf stf'
      call write_ps(char_io)

      return

      entry xpsfn(psfn,iounit )

c     write(6,'(a,1x,a)')'filename=', filename
c     write(6,'(a,1x,a)')'psfn    =', psfn

      filename = psfn
      lfn = 256
      call xstrlnth(filename,lfn)

      io = iounit

      return

      end


      SUBROUTINE PSfram
      common/psdef/io
      character char_io*132
      integer nofram
      common /pageno/ nofram
      nofram=nofram+1
      write(char_io,'(a,i5)') 'S %End of frame ',nofram
      call write_ps(char_io)
      write(char_io,'(a)') ' showpage'
      call write_ps(char_io)
      write(char_io,'(a,I4)') '%%Page: ',nofram+1
      call write_ps(char_io)
      write(char_io,'(a)') ' %%landscape %for landscape layout.'
      call write_ps(char_io)
      write(char_io,'(a)') ' 0.3 w %set default lw for new frame'
      call write_ps(char_io)
c
c write PS buffer constent to disk before moving onto next frame
c
      call flush_ps_buffer

      RETURN
      END

      BLOCK DATA
      integer nofram
      common /psdef/io
      common /pageno/ nofram
      data nofram/0/
      data io /13/
      END

C
      SUBROUTINE PScls
      common/psdef/io
      integer nofram
      common /pageno/nofram
      character char_io*132
      write(char_io,'(a)') 'S'
      call write_ps(char_io)
C Implement showpage and write the trailer here.
      write(char_io,'(a)')'showpage'
      call write_ps(char_io)
C Restore the original graphic state of the laser printer.
      write(char_io,'(a)') 'Q'
      write(char_io,'(a,I4)') '%%Pages: ', nofram+1
      call write_ps(char_io)

C write buffer content to disk
      call flush_ps_buffer

C Close the io.
      close(io)
      RETURN
      END
c
c#######################################################################
c
c (pl,pr,pt,pb) sets the area on the given device and map it with
c coordinate space (x1,x2,y1,y2).
c
c 11 X 8.5 " = 792 X 612 (1" = 72 points)
c
c Maximum allowable space in x -  792    (US paper size)
c
c#######################################################################
c
      SUBROUTINE PSspac(pl,pr,pb,pt, x1,x2,y1,y2)
      common /psscal/ p1,p2,p3,p4, xa,xb,ya,yb, xsca,ysca
      p1=pl
      p2=pr
      p3=pb
      p4=pt
      xa=x1
      xb=x2
      ya=y1
      yb=y2
      xsca = (p2-p1)/(xb-xa)
      ysca = (p4-p3)/(yb-ya)
      RETURN
      END
C
      SUBROUTINE PStran(x,y)
      common/psdef/io
      character char_io*132
      common /psscal/ pa,pb,pc,pd, xa,xb,ya,yb, xsca,ysca
      x = pa+ (x-xa)*xsca
      y = pc+ (y-ya)*ysca
      RETURN
      END

      SUBROUTINE PSpnup(x1,y1)
      common/psdef/io
      character char_io*132
      x=x1
      y=y1
      call PStran(x,y)
      write(char_io,100) x,y
      call write_ps(char_io)
100   format(1x,'S',1x,f10.2,1x,f10.2,1x,'m')
      RETURN

      ENTRY PSpndn(x2,y2)
      x=x2
      y=y2
      CALL PStran(x,y)
      write(char_io,110) x,y
      call write_ps(char_io)
110   format(3x,f10.2,1x,f10.2,1x,'l')
      RETURN
      END
C
      SUBROUTINE PSwndw(x1,x2,y1,y2)

      RETURN
      END
C

c#######################################################################
c
c Set gray degree to area filling routine xfilarea
c
c#######################################################################
c
      SUBROUTINE PSgray(gray)
      common/psdef/io
      character char_io*132
      write(char_io,'(1x,a)') 'S'
      call write_ps(char_io)
      write(char_io,100) gray
      call write_ps(char_io)
100   format(1x,f5.3,1x,'g')
      RETURN
      END

      SUBROUTINE PSlnwd(wd)
      common/psdef/io
      character char_io*132

      write(char_io,'(1x,a)') 'S'
      call write_ps(char_io)
      write(char_io,100) wd

      call write_ps(char_io)
100   format(1x,f5.3,1x,'w')
      RETURN
      END

c
c#######################################################################
c
c Re-set font by its number.
c
c#######################################################################
c
      SUBROUTINE PSfont(n)
      common/psdef/io
      character char_io*132
      character fontname*20,ftname*(*)
      save nftsiz,fontname
C Set default font name and font size:
      data nftsiz,fontname /12,'Helvetica'  /

      IF( n.eq.1) THEN
        fontname='Courier'
      ELSEIF(n.eq.2)THEN
        fontname='Helvetica'
      ELSEIF(n.eq.3)THEN
        fontname='Times-Roman'
      ELSEIF(n.eq.4)THEN
        fontname='Times-Oblique'
      ELSEIF(n.eq.5)THEN
        fontname='Symbol'
      ENDIF

      write(char_io,100)fontname
      call write_ps(char_io)
      write(char_io,101)nftsiz
      call write_ps(char_io)
 100  format(' /',a,' ff ')
 101  format(i6,' scf stf')
      RETURN

c
c#######################################################################
c
c Set postscript font size (height) in the number of points.
c
c#######################################################################
c
      ENTRY PSftsz(nsiz)
      nftsiz=nsiz
c     nftsiz=9+0.5*(nsiz-9)
      write(char_io,100)fontname
      call write_ps(char_io)
      write(char_io,101)nftsiz
      call write_ps(char_io)
      RETURN

c
c#######################################################################
c
c Re-set font option by its name, its size is set by PSftsz.
c
c#######################################################################
c
      ENTRY PSftnm(ftname)
      fontname=ftname
      write(char_io,100)fontname
      call write_ps(char_io)
      write(char_io,101)nftsiz
      call write_ps(char_io)
      RETURN
      END
c
c#######################################################################
c
C Write a character string CH with left, centered and right
C justification for just = -1,0,1.
c
c#######################################################################
c
      Subroutine PSstrg(x1,y1,ch,just)
      character ch*(*)
      common/psdef/io
      character char_io*256
      common /xags09/ da,ca,xangle,xsyman,srang,ksr,xx,yy
      x=x1
      y=y1
      CALL PStran(x,y)

c     print*,'mark 1'
      write(char_io,100) x,y
c     print*,'mark 2', char_io
      call write_ps(char_io)
100   format(1x,'S',1x,f7.2,1x,f7.2,1x,'m')
      CALL PSrot(xangle+xsyman)
      IF(just.lt.0)THEN
c       print*,'mark 3'
        write(char_io,'(a,a,a)') '(' ,ch, ') show'
c       print*,'mark 4'
        call write_ps(char_io)
      ELSEIF(just.eq.0) THEN
c       print*,'mark 5'
        write(char_io,'(a,a,a)') '(' ,ch, ') cs'
c       print*,'mark 6'
        call write_ps(char_io)
      ELSEIF(just.gt.0) THEN
c       print*,'mark 7'
        write(char_io,'(a,a,a)') '(' ,ch, ') rs'
c       print*,'mark 8'
        call write_ps(char_io)
      ENDIF
      CALL PSrot(-xangle-xsyman)
      RETURN
      END

c
c#######################################################################
c
c Rotate the coordinate frame anticlockwise through angle TangleU
c
c#######################################################################
c
      Subroutine PSrot(angle)
      common/psdef/io
      character char_io*132
      IF(angle.ne.0.0) then
        write(char_io,'(f8.2,'' rotate'')') angle
        call write_ps(char_io)
      endif
      RETURN
      END

      SUBROUTINE ZFILLN(X,Y,NP)
      integer np
      real x(np),y(np)
      CALL XFILAREA(x,y,np)
      RETURN
      END

      SUBROUTINE XFILAREA(X,Y,NP)

C This routine fills a polygon with gray pre-set by PSgray.
C The polygon is enclosed by curve (x(np), y(np)) defined in
C zxplot mathematical space. The polygon itself will be drawn
C in the width set by PSlnwd or by ZXplot routine Xthick.
      common/psdef/io
      character char_io*132
      common /psscal/ pa,pb,pc,pd, xa,xb,ya,yb, xsca,ysca
      real x(np),y(np)

      parameter (npmax=20000)
      REAL xra(npmax),yra(npmax)

      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      real x1,x2,y1,y2

      if(iwndon.eq.1)then
        nn=0
        do i=1,np-1
          x1=x(i)
          x2=x(i+1)
          y1=y(i)
          y2=y(i+1)
          call xtstwd(x1,y1,x2,y2,idispl)
          if(idispl.ne.0) then
            nn=nn+1
            IF( nn.gt.npmax) GOTO 999
            xra(nn)=x1
            yra(nn)=y1
            nn=nn+1
            IF( nn.gt.npmax) GOTO 999
            xra(nn)=x2
            yra(nn)=y2
          end if
        end do
      else
        IF( np.gt.npmax) then
          write(6,'(1x,a,/1x,a,i6)')
     :   'Work array xra and yra defined in XAREAFIL not large enough.',
     :   'Only ',npmax,' points were used.'
        endif
        do i=1,min(npmax,np)
          xra(i) = x(i)
          yra(i) = y(i)
        enddo
        nn = min(npmax,np)
      end if

      if(nn .ge. 3) then
        x1=xra(1)
        y1=yra(1)
        write(char_io,'(a)') ' N'
        call write_ps(char_io)
        call xtrans(x1,y1)
        call PStran(x1,y1)
        write(char_io,100) x1,y1
        call write_ps(char_io)
100     format(1x,f7.2,1x,f7.2,1x,'m')
        do 5 i=2,nn
          x1=xra(i)
          y1=yra(i)
          call xtrans(x1,y1)
          call PStran(x1,y1)
          write(char_io,110) x1,y1
          call write_ps(char_io)
110       format(1x,f7.2,1x,f7.2,1x,'l')
 5      continue
        write(char_io,'(a)') 'h f'
        call write_ps(char_io)
      endif

      return

999   write(6,'(1x,a,/1x,a,i6)')
     :'Work array xra and yra defined in XFILAREA were too small',
     :'Please increase the array size to more than ',npmax

      return
      end

      subroutine xichar
C This routine is not needed when characters are written using postscript
C facilitites.
      common /xchr30/ ncram(256)
      do 5 i=1,256
        ncram(i)=i
  5   continue

      return
      end

      subroutine xafstyl(nstyle)
      return
      end

      subroutine xafpatn(npat)
      return
      end

      subroutine xlncinx(ind)
      return
      end

      subroutine xafcinx(ind)
      return
      end


      subroutine set
      end
      subroutine xqset
      end
      subroutine xzx2ncar
      end

      subroutine strmln

        CALL PSstrg(0.5,0.8,
     :    'Sorry, no streamline for the Postscript version.',0)

        WRITE(6,'(a,/a)')
     :' Sorry, subroutine strmln is not available with the Postscript',
     :' version. No streamline field is plotted.'
      return
      end subroutine strmln

      subroutine strmln_new

        CALL PSstrg(0.5,0.8,
     :    'Sorry, no streamline for the Postscript version.',0)

        WRITE(6,'(a,/a)')
     :' Sorry, subroutine strmln is not available with the Postscript',
     :' version. No streamline field is plotted.'
      return
      end subroutine strmln_new

      SUBROUTINE XWRTCTBL(rgb_table,nc_max,exchange)

c#######################################################################
c
c     Activate pre-defined color table in rgb_table
c
c#######################################################################

      integer exchange
      integer nc_max,j,jj
      real rgb_table(3,nc_max)
      integer io
      common/psdef/io
      character char_io*132

      DO j=1,nc_max
       jj = j -1
       write(char_io,110)jj,rgb_table(1,j),rgb_table(2,j),rgb_table(3,j)
       call write_ps(char_io)
      END DO

110   format(' ct ',i3,' [',f6.3,1x,f6.3,1x,f6.3,'] put')

      RETURN
      END

      SUBROUTINE COLOR(icolor)
      integer icolor
        call xcolor(icolor)
      RETURN
      END

      SUBROUTINE XCOLOR(icolor)
c
c#######################################################################
c
c     Choose color by it index out of predefined color table.
c
c#######################################################################
c
c     INPUT:
c
c     icolor  index of the color
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none
      common /psdef/io
      character char_io*132
      integer io,icolor
      integer kolor, lcolor
      save kolor
      data kolor /1/

      kolor = icolor

      write(char_io,100) icolor
      call write_ps(char_io)
100   format(1x,'S ',i3,1x,'o')

      RETURN

      ENTRY XQCOLOR(lcolor)
      lcolor = kolor
      RETURN
      END

      subroutine xpscmnt(ch)
      character*(*) ch
      common /psdef/io
      character char_io*132
      write(char_io,'(a,a)') '%%PSCOMMENT:',ch
      call write_ps(char_io)
      return
      end

      subroutine write_ps( char_io )
      character char_io*(*)
      integer max_buffer
      parameter(max_buffer=500)
      character ps_buffer(max_buffer)*132
      common /psbuffer/ ps_buffer
      common /psbufferlines/ lines

      if( lines+1 .gt. max_buffer ) call flush_ps_buffer

      lines = lines + 1
      ps_buffer(lines)=char_io

      return
      end

      subroutine flush_ps_buffer
      common /psdef/io
      integer max_buffer
      parameter(max_buffer=500)
      character ps_buffer(max_buffer)*132
      common /psbuffer/ ps_buffer
      common /psbufferlines/ lines

      do i=1,lines
!        lps_buffer=132
!        call xstrlnth(ps_buffer(i), lps_buffer)
!        write(io,'(a)') ps_buffer(i)(1:lps_buffer)
        write(io,'(a)') trim(ps_buffer(i))
      enddo
      lines = 0

      return
      end


