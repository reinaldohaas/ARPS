*
*     ##################################################################
*     ##################################################################
*     ######                                                      ######
*     ######      Advanced Regional Prediction System (ARPS)      ######
*     ######                     Version 4.2                      ######
*     ######                                                      ######
*     ######             GrADS Display of ARPS Model              ######
*     ######                                                      ######
*     ######                     Developed by                     ######
*     ######     Center for Analysis and Prediction of Storms     ######
*     ######                University of Oklahoma                ######
*     ######                                                      ######
*     ##################################################################
*     ##################################################################
*
************************************************************************
*                                                                      *
* This GrADS script file includes standard scripts to display ARPS     *
* history outputs in GrADS format. The name of GrADS control file      *
* may be given as an argument, or otherwise be read in from            *
* arpsgrads.tmp.                                                       *
*                                                                      *
* AUTHOR: Yuhe Liu                                                     *
*         06/22/1995 beta 1.0                                          *
*         10/22/1996 1.0                                               *
*                                                                      *
* GrADS - Grid Analysis and Display System, Version 1.5.1.3            *
*         Copyright (c) 1988-1996 by Brian Doty                        *
*         Center for Ocean-Land-Atmosphere Studies                     *
*         Institute for Global Environment and Society                 *
*                                                                      *
************************************************************************

function main(args)

************************************************************************
*                                                                      *
* Get the data control file name from temporary file arpsgrads.tmp.    *
*                                                                      *
************************************************************************

    _ctlfile = subwrd(args,1)

    if ( _ctlfile = '' )
        say " "
        say "Usage: run sfcflx.gs sfcflx_control_filename"
        say " "
        return
    else
        'open ' _ctlfile
        dummy = sublin( result, 2 )
        dummy = subwrd( dummy, 2 )
        if ( dummy = 'Error:' )
            say " "
            say "Usage: run sfcflx.gs sfcflx_control_filename"
            say " "
            return
        endif
    endif

    _mbtns = 8
    _mbtn.1 = "X-Y"
    _mbtn.2 = "T-series"
    _mbtn.3 = "Zoom"
    _mbtn.4 = "Overlap"
    _mbtn.5 = "Shaded"
    _mbtn.6 = "Animation"
    _mbtn.7 = "Print"
    _mbtn.8 = "Quit"

    _tbtnlen = 10.6
    _tbtnhgt = 0.4

    _bbtnlen = 10.6
    _bbtnhgt = 0.8

    _rbtnhgt = 6.0
    _rbtnlen = 2.0

    _lbtnhgt = 6.0
    _lbtnlen = 2.0

    _tbtnx0  = 0.2
    _tbtny0  = 8.0

    _bbtnx0  = 0.2
    _bbtny0  = 0.1

    _lbtnx0  = 0.2
    _lbtny0  = _bbtny0 + _bbtnhgt + 0.3

    _rbtnx0  = 11 - _lbtnlen - _lbtnx0
    _rbtny0  = _lbtny0

    _lbtnmax = 30
    _rbtnmax = 30

    _pareax0 = 3.0
    _pareax1 = 8.0
    _pareay0 = 2.0
    _pareay1 = 7.0

    _plabx = 5.25
    _plaby = _bbtny0 + _bbtnhgt + 0.5

************************************************************************
*                                                                      *
* Get the run name from the data control file                          *
*                                                                      *
************************************************************************

    dummy = read(_ctlfile)
    runname = subwrd(dummy,8)
    dummy = close(_ctlfile)

************************************************************************
*                                                                      *
* Set the fore/back-ground colors and other environments               *
*                                                                      *
************************************************************************

    'set display color black'
    'clear'

    'set strsiz 0.50 0.5'
    'set string 7 bc 12'
    'draw string 5.5 7 GrADS'
    'set strsiz 0.50 0.50'
    'draw string 5.5 5.5 Display ARPS Surface Data'
    'set string 7 bc 6'
    'set strsiz 0.35 0.35'
    'draw string 5.5 3.5 Run name: 'runname
    'set strsiz 0.50 0.50'
    'set string 2 bc 12'
    'draw string 5.5 1.5 Click When Ready'
    'set string 1 l 1'
    'q bpos'

    dummy = initdis()
    dummy = sfcdisp()
    return

return

************************************************************************
*                                                                      *
* Initialization
*                                                                      *
************************************************************************

function initdis()

************************************************************************
*                                                                      *
* Get the grid and structure information for ARPS history data         *
*                                                                      *
* Define global variables for the dimension and grid.                  *
*                                                                      *
************************************************************************

    while (1)
        dummy = read(_ctlfile)
        if ( sublin( dummy, 1 ) != 0 ); break; endif
        xdim = sublin( dummy, 2 )
        defstr = subwrd(xdim,1)
        if ( defstr = 'XDEF' )
           _nx = subwrd(xdim,2)
           _xbgn = subwrd(xdim,4)
           _xint = subwrd(xdim,5)
           break
        endif
    endwhile
    dummy = close(_ctlfile)

    while (1)
        dummy = read(_ctlfile)
        if ( sublin( dummy, 1 ) != 0 ); break; endif
        ydim = sublin( dummy, 2 )
        defstr = subwrd(ydim,1)
        if ( defstr = 'YDEF' )
            _ny = subwrd(ydim,2)
            _ybgn = subwrd(ydim,4)
            _yint = subwrd(ydim,5)
            break
        endif
    endwhile
    dummy = close(_ctlfile)

    while (1)
        dummy = read(_ctlfile)
        if ( sublin( dummy, 1 ) != 0 ); break; endif
        zdim = sublin( dummy, 2 )
        defstr = subwrd(zdim,1)
        if ( defstr = 'ZDEF' )
            if ( subwrd(zdim,3) = 'LINEAR' )
                _nz = subwrd(zdim,2)
                _zbgn = subwrd(zdim,4)
                _zint = subwrd(zdim,5)

                i = 0
                while ( i < _nz )
                    _lev.i = _zbgn + (i-1)*_zint
                    i = i + 1
                endwhile
            else if ( subwrd(zdim,3) = 'LEVELS' )
                _nz = subwrd(zdim,2)
                i = 0
                dummy = read(_ctlfile)
                while ( i < _nz & sublin(dummy,1) = 0 )
                    j = 1
                    levels = sublin(dummy,2)
                    while ( subwrd(levels,j) != '' )
                        _lev.i = subwrd(levels,j)
                        i = i + 1
                        j = j + 1
                    endwhile
                    dummy = read(_ctlfile)
                endwhile
            endif
            break
        endif
    endwhile
    dummy = close(_ctlfile)

    while (1)
        dummy = read(_ctlfile)
        if ( sublin( dummy, 1 ) != 0 ); break; endif
        tdim = sublin( dummy, 2 )
        defstr = subwrd(tdim,1)
        if ( defstr = 'TDEF' )
            _nt = subwrd(tdim,2)
            _tbgn = subwrd(tdim,4)
            _tintstr = subwrd(tdim,5)

            _zhr0 = substr(_tbgn,1,2)
            if ( substr(_tbgn,3,3) = ':' )
                _min0 = substr(_tbgn,4,2)
                _day0 = substr(_tbgn,7,2)
                _mon0 = substr(_tbgn,9,3)
                _yr0  = substr(_tbgn,12,4)
            else
                _min0 = 0
                _day0 = substr(_tbgn,4,2)
                _mon0 = substr(_tbgn,6,3)
                _yr0  = substr(_tbgn,9,4)
            endif

            _tint = substr(_tintstr,1,2)
            _tunt = substr(_tintstr,3,2)

            break
        endif
    endwhile
    dummy = close(_ctlfile)

    while (1)
        dummy = read(_ctlfile)
        if ( sublin( dummy, 1 ) != 0 ); break; endif
        vars = sublin( dummy, 2 )
        defstr = subwrd(vars,1)
        if ( defstr = 'VARS' )
            _vnum = subwrd(vars,2)
            i = 1
            dummy = read(_ctlfile)
            while ( i <= _vnum & sublin(dummy,1) = 0 )
                var = sublin(dummy,2)
                if ( subwrd(var,1) != 'ENDVARS' )
                    _var.i = subwrd(var,1)
                endif
                i = i + 1
                dummy = read(_ctlfile)
            endwhile
            break
        endif
    endwhile
    dummy = close(_ctlfile)

    _vnpl = 14
    _varln = _vnum/_vnpl
    if ( _varln>4 )
      say 'Too many variables'
      return
    endif

    i = 1
    while ( i<=4 )
      if ( _varln>i-1 & _varln<=i )
        _varln = i
        break
      endif
      i = i + 1
    endwhile

    proj = 0
    while (1)
        dummy = read(_ctlfile)
        if ( sublin( dummy, 1 ) != 0 ); break; endif
        mprj = sublin( dummy, 2 )
        defstr = subwrd(mprj,1)
        if ( defstr = 'PDEF' )
            proj = 1
            break
        endif
    endwhile
    dummy = close(_ctlfile)
return

************************************************************************
*                                                                      *
* Function arpsdis is the most important function that does the main   *
* tasks of displaying ARPS GrADS history outputs.                      *
*                                                                      *
************************************************************************

function sfcdisp()

    'reinit' 

************************************************************************
*                                                                      *
* Open the data set.                                                   *
*                                                                      *
************************************************************************

    'open '_ctlfile

    dummy = sublin( result, 1 )
    dummy = sublin( result, 2 )

    if ( subwrd(dummy,2) = "Error" )
      return
    else
      _dfile = subwrd(dummy,3)
    endif

************************************************************************
*                                                                      *
* Set the fore/back-ground colors and other environments               * 
*                                                                      *
************************************************************************

    'set display color white'
    if ( proj != 1 )
        'set mproj off'
        _nx1 = _nx-1
        _ny1 = _ny-1
        _nz1 = _nz-1
    else
        'set mpdset hires'
        _nx1 = _nx
        _ny1 = _ny
        _nz1 = _nz
    endif
    'set grid off'
    'clear'

************************************************************************
*                                                                      *
* Set up the user's environment and draw an initial picture of pt      *
*                                                                      *
************************************************************************

    'set parea '_pareax0' '_pareax1' '_pareay0' '_pareay1
    'set string 7 bc 12'
    'set strsiz 0.4 0.4'

    'set x 1 '_nx1
    'set y 1 '_ny1
    'set z 2'
    'set t 1'

    'set csmooth on'
    'set clab forced'
    'set cterp off'

    'set strsiz 0.20 0.20'
    'set string 1 c 6'
    'set gxout contour'
*   'set gxout shaded'
    'set csmooth on'
    'set grads off'
    'd ts'
*   dummy = colorbars()
    'draw title Initial Ts'

************************************************************************
*                                                                      *
* Loop to display the environment.                                     *
*                                                                      *
************************************************************************

    i = 1
    while ( i<= _mbtns )
      _mchoi.i = 0
      i = i + 1
    endwhile

    _vchoi = 29
    _mchoi.1 = 1
    _time  = 1
    _level = 1

    while( _vchoi!=0 )

        dummy = colors() 
        dummy = qfield()

        if ( _mchoi._mbtns!=0 )
          return
        endif

        if ( _mchoi.2=0 & _mchoi.3=0 & _mchoi.4=0 & _mchoi.5=0 )
            'clear'
            dummy = colors() 
        endif

        drawflg = 1

        if ( _mchoi.2=2 | _mchoi.3=2 | _mchoi.4=2 | _mchoi.5=2 )
            drawflg = 0
        endif

        if ( _mchoi.1=2 )
            _mchoi.1 = 1
        endif
 
        if ( _mchoi.2=2 )
            _mchoi.2 = 1
        endif
 
        if ( _mchoi.3=2 )
            _mchoi.3 = 1
        endif
 
        if ( _mchoi.4=2 )
            _mchoi.4 = 1
        endif
 
        if ( _mchoi.5=2 )
            _mchoi.5 = 1
        endif
 
        if ( _mchoi.6=2 )
            _mchoi.6 = 1
        endif
 
        if ( drawflg = 1 )
            dummy = drawf()
        endif

    endwhile

return


************************************************************************
*                                                                      *
* Function drawf to display the picture that users selected.           *
*                                                                      *
************************************************************************

function drawf()

    'set display color white'
    'set rbcols auto'
    'set ccolor rainbow'
    'set clab forced'

    if ( _mchoi.5!=0 )
      gxtype = shaded
    else
      gxtype = contour
    endif

    if ( _mchoi.2=0 & _mchoi.3=0 & _mchoi.4=0 & _mchoi.5=0 )
        'set x 1 '_nx1
        'set y 1 '_ny1
        'set z '_level
        'set t '_time
        it = (_time-1)*_tint
        _plab = 't='it' '_tunt
    endif

    i = 1
    while ( i<=_vnum )
      if ( i=_vchoi )
        titl = _var.i
        vara = _var.i
        prnm = _var.i
      endif
      i = i + 1
    endwhile

    'set grads off'

    if ( _mchoi.7!=0 & _vchoi>0.5 )
        'gxprint 'prnm'.meta'
    endif

    if ( _vchoi>0.5 )
      if ( _mchoi.6!=0 )
        dummy = animate(_vchoi,vara,titl,_level)
        return
      else
        if ( _mchoi.4=1 )
          'set gxout contour'
          'set csmooth on'
          'set clab forced'
          'set ccolor 1'
          'set grads off'
          'd 'vara
          return
        endif

        'c'
        'set grads off'
        'set gxout 'gxtype
        'set csmooth on'

        'd 'vara

        'draw title 'titl

        if ( gxtype=shaded )
          dummy = colorbars()
        else
          'set string 1 c 5'
          'set strsiz .16 .16'
          'draw string '_plabx' '_plaby' '_plab
        endif
      endif
    endif

    if ( mchoi.7=1 )
        'print'
        'disable'
        '!gxps -i 'prnm'.meta -o 'prnm'.ps'
*       '!lpr -Pps 'prnm'.ps'
*       '!rm -f 'prnm'.meta 'prnm'.ps'
        _mchoi.7 = 0
    endif

return

************************************************************************
*                                                                      *
* This function implements the animation for selected variable. It     *
* calls function showtime                                              *
*                                                                      *
************************************************************************

function animate(ch,vara,titl,level)

    'c'
    'set dbuff on' 
    i=1
    'set t 'i''
    while( i <= _nt )

        'set grads off'
        'set csmooth on'
        'd 'vara
        'draw title 'titl

        dummy=showtime(i)
        'swap'
        i=i+1
        'set t 'i''
    endwhile
  
    while (1)
        'q bpos'
        x = subwrd(result,3)
        if (x != 0); break; endif;
    endwhile

    'set dbuff off'
    'set clab forced'

return

************************************************************************
*                                                                      *
* This function will get the current time and then prints it.          *
*                                                                      *
* One argument is needed for time index.                               *
*                                                                      *
************************************************************************

function showtime(i)

    'set t 'i
    'query time'
    time = subwrd(result,3)
    hour = substr(time,1,2)
    minu = substr(time,4,2)
    mnth = substr(time,6,3)
    year = substr(time,9,4)
    day  = substr(time,4,2)

    x = 5.25
    y = _bbtny0 + _bbtnhgt + 0.5
    inc=(i-1)*_tint
    'set string 1 c 5'
    'set strsiz .16 .16'
    'draw string '_plabx' '_plaby' 'mnth' 'day', 'year', 'hour':'minu'Z, t='inc' '_tunt

return
 
************************************************************************
*                                                                      *
* Function qfield draws the buttons of display options and determines  *
* what options the user chose.                                         *
*                                                                      *
************************************************************************

function qfield()

    'set string 1 c 3'
    push = 0
    'set strsiz 0.10 0.15'
    'set line 1 1 6'

    dx = _bbtnlen/_vnpl
    dy = _bbtnhgt/_varln

    y1 = _bbtny0 + _bbtnhgt
    y0 = y1

    vl = 1
    vn = 1
    while ( vl<=_varln )
      y1 = y0
      y0 = y1 - dy

      x0 = _bbtnx0
      x1 = x0

      i = 1
      while ( i<=_vnpl & vn<=_vnum )
        if (vn=_vchoi)
          'set string 2 c 6'
          push = 1
        else
          'set string 1 c 3'
          push = 0
        endif

        x0 = x1
        x1 = x0 + dx

        dummy = button(x0,y0,x1,y1,push,1,5,0.10,0.12,_var.vn)

        i = i + 1
        vn = vn + 1
      endwhile

      vl = vl + 1
    endwhile

    x0 = _lbtnx0 + _lbtnlen/2
    y0 = _lbtny0 + _lbtnhgt + 0.3

    'set strsiz 0.12 0.18'
    'set string 1 c 6'
    'draw string 'x0' 'y0' Forecast ('_tunt')'
    'set string 1 c 3'

    dx = _lbtnlen/2
    if ( _nt <= 30 )
        x = _lbtnx0
        dy = _lbtnhgt/_nt

        i = 0
        while ( i < _nt )
            y = i*dy + _lbtny0
            ai = i * _tint

            if ( i = _time-1 )
                'set string 2 c 6'
                push = 1
            else
                'set string 1 c 3'
                push = 0
            endif

            x0 = x + dx/2
            x1 = x0 + dx
            y0 = y
            y1 = y0 + dy
            dummy = button(x0,y0,x1,y1,push,1,5,0.12,0.08,ai)
            i = i + 1
        endwhile
    else
        nt1 = _nt/2
        x = _lbtnx0
        dy = _lbtnhgt/(nt1+1)

        i = 0
        while (i < nt1)
            it = 2*i
            itp1 = it + 1

            y = _lbtny0 + i*dy

            if ( itp1 < _nt )
                ai = itp1 * _tint

                if ( itp1 = _time-1 )
                    'set string 2 c 6'
                    push = 1
                else
                    'set string 1 c 3'
                    push = 0
                endif

                x0 = x + dx
                x1 = x0 + dx
                y0 = y + dy/2
                y1 = y0 + dy

                dummy = button(x0,y0,x1,y1,push,1,5,0.12,0.08,ai)
            endif

            ai = it * _tint

            if ( it = _time-1 )
                'set string 2 c 6'
                push = 1
            else
                'set string 1 c 3'
                push = 0
            endif

            x0 = x
            x1 = x0 + dx
            y0 = y
            y1 = y0 + dy

            dummy = button(x0,y0,x1,y1,push,1,5,0.12,0.08,ai)
            _time2 = _time2' 'ai
            i = i + 1
        endwhile

    endif

    x0 = _tbtnx0
    x1 = x0
    dx = _tbtnlen/_mbtns

    y0 = _tbtny0
    y1 = y0 + _tbtnhgt

    i = 1
    while ( i<=_mbtns )
      x0 = x1
      x1 = x0 + dx
      push = 0
      if ( _mchoi.i!=0 ); push=1; endif
      dummy = button(x0,y0,x1,y1,push,1,8,0.12,0.15,_mbtn.i)
      i = i + 1
    endwhile

    while (1)
        'q bpos'
        xms = subwrd(result,3)
        yms = subwrd(result,4)

* Check if the click point is on variable menu

        dx = _bbtnlen/_vnpl
        dy = _bbtnhgt/_varln
        y1 = _bbtny0 + _bbtnhgt
        y0 = y1

        vl = 1
        vn = 1
        brk = 0
        while ( vl<=_varln )
          y1 = y0
          y0 = y1 - dy

          x0 = _bbtnx0
          x1 = x0

          i = 1
          while ( i<=_vnpl & vn<=_vnum )
            x0 = x1
            x1 = x0 + dx

            if ( yms>y0 & yms<y1 & xms>x0 & xms<x1 )
                _vchoi=vn
                brk = 1
            endif

            i = i + 1
            vn = vn + 1
          endwhile

          vl = vl + 1
        endwhile
        if ( brk = 1 ); break; endif

* Check if the click point is on time menu

        dx = _lbtnlen/2
        if ( _nt <= 30 )
            brk = 0
            x = _lbtnx0
            dy = _lbtnhgt/_nt

            i = 0
            while ( i < _nt )
                y = i*dy + _lbtny0

                x0 = x + dx/2
                x1 = x0 + dx
                y0 = y
                y1 = y0 + dy
                if ( xms>x & xms<x1 & yms>y0 & yms<y1 )
                    _time=i+1
                    brk=1
                endif
                i = i + 1
            endwhile
            if ( brk = 1 ); break; endif
        else
            brk = 0
            nt1 = _nt/2
            x = _lbtnx0
            dy = _lbtnhgt/(nt1+1)

            i = 0
            while (i < nt1)
                it = 2*i
                itp1 = it + 1

                y = _lbtny0 + i*dy

                if ( itp1 < _nt )

                    x0 = x + dx
                    x1 = x0 + dx
                    y0 = y + dy/2
                    y1 = y0 + dy

                    if ( xms>x0 & xms<x1 & yms>y0 & yms<y1 )
                        _time=it+1
                        brk=1
                    endif

                endif

                x0 = x
                x1 = x0 + dx
                y0 = y
                y1 = y0 + dy

                if ( xms>x0 & xms<x1 & yms>y0 & yms<y1 )
                    _time=it+1
                    brk=1
                endif

                i = i + 1
            endwhile
            if ( brk = 1 ); break; endif

        endif

* Check if the click point is on main menu

        x0 = _tbtnx0
        x1 = x0
        dx = _tbtnlen/_mbtns

        y0 = _tbtny0
        y1 = y0 + _tbtnhgt

        brk = 0
        i = 1
        while ( i<=_mbtns )
          x0 = x1
          x1 = x0 + dx
          if ( yms>y0 & yms<y1 & xms>x0 & xms<x1 )
            brk = 1
            if ( _mchoi.i=0 )
               _mchoi.i = 2
            else
               _mchoi.i = 0
            endif

            if ( i=1 )
              _mchoi.2 = 0
            else
              if ( i=2 )
                _mchoi.1 = 0
              endif
            endif
            break
          endif
          i = i + 1
        endwhile

        if ( brk = 1 ); break; endif

* Check if the click point is within the plot area for zoom

        if ( xms>_pareax0 & xms<_pareax1 & yms>_pareay0 & yms<_pareay1 )
*         if ( _mchoi.1!=0 )
*           'set x 1 '_nx1
*           'set y 1 '_ny1
*           'set z '_level
*           'set t '_time
*         endif

          if ( _mchoi.2=1 )
            'q xy2w 'xms' 'yms
            cxlb = subwrd(result,1)
            cylb = subwrd(result,4)
            if ( cxlb = 'Lon' & cylb = 'Lat' )
              cx = subwrd(result,3)
              cy = subwrd(result,6)
              'set lon 'cx
              'set lat 'cy
              'set z '_level
              'set t 1 '_nt
              iz = _level-1
              iz = _lev.iz
              _plab = 'x='cx' km, y='cy' km'
            else
              _mchoi.2 = 0
            endif
            break
          endif
    
          if ( _mchoi.3!=0 )
            'q xy2w 'xms' 'yms
            cxlb = subwrd(result,1)
            cylb = subwrd(result,4)
            cx = subwrd(result,3)
            cy = subwrd(result,6)

            if ( cxlb = 'Lon' )
              xx = 1 + (cx-_xbgn)/_xint
              if ( xx <= _nx1/4 )
                xbgn = 1
                xend = _nx1/2
              else
                if ( xx > _nx1/4 & xx <= _nx1*3/4 )
                  xbgn = _nx1/4
                  xend = _nx1*3/4
                else
                  if ( xx > _nx1*3/4 )
                    xbgn = _nx1/2
                    xend = _nx1
                  endif
                endif
              endif

              'set x 'xbgn' 'xend

              if ( cylb = 'Lat' )
                yy = 1 + (cy-_ybgn)/_yint
                  if ( yy <= _ny1/4 )
                    ybgn = 1
                    yend = _ny1/2
                  else
                    if ( yy > _ny1/4 & yy <= _ny1*3/4 )
                      ybgn = _ny1/4
                      yend = _ny1*3/4
                    else
                      if ( yy > _ny1*3/4 )
                        ybgn = _ny1/2
                        yend = _ny1
                      endif
                    endif
                  endif
                  'set y 'ybgn' 'yend
                else
                  _mchoi.3 = 0
                endif
              endif
              _plab = 'Zoom'
              break
            endif
        endif
    endwhile

return


************************************************************************
*                                                                      *
* Function colors sets up the color table for shaded displays.         *
*                                                                      *
************************************************************************

function colors()

**  specifies some colors for snow **
    'set rgb 40 200 200 200'
    'set rgb 41 150 150 150'
    'set rgb 42 100 120 120'
    'set rgb 43 140 150 155'
    'set rgb 44 200 200 210'
    'set rgb 45 240 240 255'

**  specifies some colors for precip **
    'set rgb 50 0 50 0'
    'set rgb 51 0 100 0'
    'set rgb 52 0 150 0'
    'set rgb 53 0 200 0'
    'set rgb 54 0 250 0'
    'set rgb 55 0 255 0'

    'set rgb 56 100 100 100'
    'set rgb 57 100 100 74'
    'set rgb 58 73 77 50'
    'set rgb 59 50 50 50'

**  colors for temperature scale*
    'set rgb 64 233 233 255'
    'set rgb 65 200 200 255'
    'set rgb 66 166 166 255'
    'set rgb 67 133 133 255'
    'set rgb 68 100 100 255'
    'set rgb 69 50 50 255'
    'set rgb 70 0 180 255'
    'set rgb 71 0 255 175'
    'set rgb 72 0 255 0'
    'set rgb 73 175 255 0'
    'set rgb 74 255 255 0'
    'set rgb 75 255 233 0'
    'set rgb 76 255 200 0'
    'set rgb 77 255 166 0'
    'set rgb 78 255 133 0'
    'set rgb 79 255 105 0'
    'set rgb 80 250 0 0'
    'set rgb 81 230 0 0'
    'set rgb 82 215 0 0'
    'set rgb 83 175 0 0'
    'set rgb 84 125 0 0'
    'set rgb 85 255 225 225'

return 


************************************************************************
*                                                                      *
* Function tbar draws the bar scale for temperature                    *
*                                                                      *
************************************************************************

function tbar()

    i = 0.05
    c=64
    'set line 'c' 1'
    m = i + 0.425
    'draw recf 'i' 8.00 'm' 8.25'
    'set string 1 tc 4'
    'set strsiz 0.125'
    'draw string 0.265 7.95 -0'
    i= 0.85
    c=65
    t=5

    While (i < 9.05)
        'set line 'c' 1'
        m = i + 0.425
        'draw recf 'i' 8.00  'm' 8.25'
        'draw string 'm' 7.95 't
        i = i + 0.425
        c = c + 1
        t = t +5
    endwhile

    i = i + 0.4 
    m = i + 0.425
    t = t - 5
    'set line 'c' 1'
    'draw recf 'i' 8.00 'm' 8.25'
    n = (i + m)/2
    'draw string 'n' 7.95 't'+'

return

************************************************************************
*                                                                      *
* Function divscl sets up the scales for divergency or cinvergency.    *
*                                                                      *
************************************************************************

function divscl()

    'set line 9 1 1'
    'draw line 5.5 8.3 5.8 8.3'
    'set line 14 1'
    'draw line 5.8 8.3 6.1 8.3'
    'set line 4 1'
    'draw line 6.1 8.3 6.4 8.3'
    'draw line 6.35 8.25 6.4 8.3'
    'draw line 6.35 8.35 6.4 8.3'
    'set line 11 1'
    'draw line 6.4 8.3 6.7 8.3'
    'set line 5 1'
    'draw line 6.7 8.3 7.0 8.3'
    'set line 13 1'
    'draw line 7.0 8.3 7.3 8.3'
    'set line 3 1'
    'draw line 7.3 8.3 7.6 8.3'
    'draw line 7.55 8.25 7.6 8.3'
    'draw line 7.55 8.35 7.6 8.3'
    'set line 10 1'
    'draw line 7.6 8.3 7.9 8.3'
    'set line 7 1'
    'draw line 7.9 8.3 8.2 8.3'
    'set line 12 1'
    'draw line 8.2 8.3 8.5 8.3'
    'set line 8 1'
    'draw line 8.5 8.3 8.8 8.3'
    'draw line 8.75 8.25 8.8 8.3'
    'draw line 8.75 8.35 8.8 8.3'
    'set line 2 1'
    'draw line 8.8 8.3 9.1 8.3'
    'set line 6 1'
    'draw line 9.1 8.3 9.4 8.3'
    'draw string 5.8 7.95 Convergence'
    'draw string 9.1 7.95 Divergence'

return


************************************************************************
*                                                                      *
* Function button draw a filled box and write the given label with the *
* given string thickness and size.                                     *
*                                                                      *
* push = 0 out (grey), = 1 in (pink)
*                                                                      *
************************************************************************

function button(x1,y1,x2,y2,push,cstr,strthc,xsize,ysize,labl)

    y = (y1 + y2) * 0.5
    x = (x1 + x2) * 0.5

    'set rgb 90 100 100 100'
    'set rgb 91 150 150 150'
    'set rgb 92 200 200 200'
    'set rgb 93 230 170 170'

    'set line 90'
    if (push);'set line 92';endif
    'draw recf 'x1' 'y1' 'x2' 'y2

    y1 = y1 + 0.03
    x2 = x2 - 0.03

    'set line 92'
    if (push);'set line 90';endif
    'draw recf 'x1' 'y1' 'x2' 'y2

    y2 = y2 - 0.03
    x1 = x1 + 0.03

    'set line 91'
    if (push);'set line 93';endif
    'draw recf 'x1' 'y1' 'x2' 'y2

    'set string 'cstr' c 'strthc
    'set strsiz 'xsize' 'ysize
    'draw string 'x' 'y' 'labl

return

************************************************************************
*                                                                      *
* Color bars for shaded                                                *
*                                                                      *
* This function was modified from original distributed GrADS scripts,  *
* cbar.gs                                                              *
*                                                                      *
************************************************************************

*
*  Script to plot a colorbar
*
*  The script will assume a colorbar is wanted even if there is 
*  not room -- it will plot on the side or the bottom if there is
*  room in either place, otherwise it will plot along the bottom and
*  overlay labels there if any.  This can be dealt with via 
*  the 'set parea' command.  In version 2 the default parea will
*  be changed, but we want to guarantee upward compatibility in
*  sub-releases.
*
function colorbar (args)

*
*  Check shading information
*
  'query shades'
  shdinfo = result
  if (subwrd(shdinfo,1)='None') 
    say 'Cannot plot color bar: No shading information'
    return
  endif
* 
*  Get plot size info
*
  'query gxinfo'
  rec2 = sublin(result,2)
  rec3 = sublin(result,3)
  rec4 = sublin(result,4)
  xsiz = subwrd(rec2,4)
  ysiz = subwrd(rec2,6)
  ylo = subwrd(rec4,4)
  xhi = subwrd(rec3,6)
  xd = xsiz - xhi
*
*  Decide if horizontal or vertical color bar
*  and set up constants.
*
* if (ylo<0.6 & xd<1.0) 
*   say "Not enough room in plot for a colorbar"
*   return
* endif

  cnum = subwrd(shdinfo,5)
* if (ylo<0.6 | xd>1.5)
*   xl = xhi + xd/2 - 0.4
*   xr = xl + 0.2
*   xwid = 0.2
*   ywid = 0.5
*   if (ywid*cnum > ysiz*0.8) 
*     ywid = ysiz*0.8/cnum
*   endif
*   ymid = ysiz/2
*   yb = ymid - ywid*cnum/2
*   'set string 1 l 5'
*   vert = 1
* else
    ymid = _bbtny0+_bbtnhgt + 0.3
    yt = ymid + 0.2
    yb = ymid
    xmid = xsiz/2
    xwid = 0.5
    if (xwid*cnum > xsiz*0.8)
      xwid = xsiz*0.8/cnum
    endif
    xl = xmid - xwid*cnum/2
    'set string 1 tc 5'
    vert = 0
* endif
*
*  Plot colorbar
*
  'set strsiz 0.08 0.11'
  num = 0
  while (num<cnum) 
    rec = sublin(shdinfo,num+2)
    col = subwrd(rec,1)
    hi = subwrd(rec,3)
    'set line 'col
    if (vert) 
      yt = yb + ywid
    else 
      xr = xl + xwid
    endif
    'draw recf 'xl' 'yb' 'xr' 'yt
    if (num<cnum-1)
      if (vert) 
        'draw string '%(xr+0.05)%' 'yt' 'hi
      else
        'draw string 'xr' '%(yb-0.05)%' 'hi
      endif
    endif
    num = num + 1
    if (vert); yb = yt;
    else; xl = xr; endif;
  endwhile
