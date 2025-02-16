*
* Plot typhoon track
*
* Based on work by Yuanbing,2007,09,06
*
function main(args)
  file=subwrd(args,1)
  usage="usage: lujing.gs file legend ybsize tysize"
  if(file='')
    say usage
    return
  endif

* location of legend
* 1=>left-up; 2=>left-down; 3=>right-up; 4=>right-down
  legend=subwrd(args,2)
  if(legend=''); legend=1; endif

  ybsize=subwrd(args,3)
  if(ybsize=''); ybsize=0.14; endif
  tysize=subwrd(args,4)
  if(tysize=''); tysize=0.3; endif

*************************************************
  string=read(file)
  mask=sublin(string,1)
  data=sublin(string,2)

  if (mask=0)
    check=subwrd(data,1)
    check=substr(check,1,1)
    if(check!='#')
      say usage
      return
    endif
    runname=subwrd(data,2)
    lon1=subwrd(data,3)
    lon2=subwrd(data,4)
    lat1=subwrd(data,5)
    lat2=subwrd(data,6)
  else
    say usage
    return
  endif

  if(lat1=''); lat1=7.0;   endif
  if(lat2=''); lat2=43.0;  endif
  if(lon1=''); lon1=100.0; endif
  if(lon2=''); lon2=153.0; endif
*  say lon1' 'lon2' 'lat1' 'lat2

*****************************************
  'reinit'
  'open 'runname'_ght.ctl'

***********************************************
  'set lon 'lon1' 'lon2''
  'set lat 'lat1' 'lat2''
  'set grid off'
* 'set mpdset cnworld'
  'set mpdset mres'
* 'set map 2 1 2'
  'set map 1 1 2'
  'set vpage 0 11 0 8.5'
  'set parea 1 10 1 8'
  'set xlopts 1 5 0.15'
  'set ylopts 1 5 0.15'
  'set grads off'
  'set grid on'
  'set csmooth on'
  'set xlint 5'
  'set ylint 5'
* 'enable print lujing.gmf'
  'set gxout shaded'
  'd hgt'
*  must for a background

****************************************************************
  string=read(file)
  mask=sublin(string,1)
  data=sublin(string,2)

* Check file and find how many lines and store it in k
  k=0
  while(mask=0)
    check=subwrd(data,1)
    check=substr(check,1,1)
    if(check!='#' & check!=''); k=k+1; endif
    string=read(file)
    mask=sublin(string,1)
    data=sublin(string,2)
  endwhile
  myclose=close(file)
  k=k/2
* k is now the number of typhoon tracks

  i=1
  while(i<=k)
    string=read(file)
    lons=sublin(string,2)
    ttt=subwrd(lons,1)
    check=substr(ttt,1,1)
    while(check='#'|check='')
      if(check='#')
        sch.i=subwrd(lons,2)
        marker.i=subwrd(lons,3)
        linec.i=subwrd(lons,4)
        lines.i=subwrd(lons,5)
      endif
      string=read(file)
      lons=sublin(string,2)
      ttt=subwrd(lons,1)
      check=substr(ttt,1,1)
    endwhile
*   Track i with lons

    string=read(file)
    lats=sublin(string,2)
    ttt=subwrd(lats,1)
    while(substr(ttt,1,1)='#'|substr(ttt,1,1)='')
      string=read(file)
      lats=sublin(string,2)
      ttt=subwrd(lats,1)
    endwhile
*   Track i with lats

    'set line 'linec.i' 'lines.i

    x1=subwrd(lons,1)
    y1=subwrd(lats,1)
*   Convert lat/lon to x/y
    if(x1!=0&y1!=0)
      'query w2xy ' x1' 'y1''
      x1=subwrd(result,3)
      y1=subwrd(result,6)
*     Draw typhoon symbol 41，size is 0.3
      if (marker.i>0)
        'draw mark 'marker.i' 'x1' 'y1' 'tysize/10
      else
        'draw wxsym 41 'x1' 'y1' 'tysize
      endif
    endif

    num=2
    x2=subwrd(lons,num)
    y2=subwrd(lats,num)
    while(x2!=''&y2!=''&x2!=0&y2!=0)
*     say i'-'num': 'x2' 'y2''
      'query w2xy ' x2' 'y2''
      x2=subwrd(result,3)
      y2=subwrd(result,6)
*     say i'-'num': 'x2' 'y2''
      'draw line 'x1' 'y1' 'x2' 'y2
      x1=x2
      y1=y2
      if (marker.i>0)
        'draw mark 'marker.i' 'x1' 'y1' 'tysize/10
      else
        'draw wxsym 41 'x1' 'y1' 'tysize
      endif
      num=num+1
      x2=subwrd(lons,num)
      y2=subwrd(lats,num)
    endwhile

    i=i+1
  endwhile

  if(legend=1)
    'query w2xy ' lon1' 'lat2''
    xx=subwrd(result,3)
    yy=subwrd(result,6)
  endif
  if(legend=2)
    'query w2xy ' lon1' 'lat1''
    xx=subwrd(result,3)
    yy=subwrd(result,6)
  endif
  if(legend=3)
    'query w2xy ' lon2' 'lat2''
    xx=subwrd(result,3)
    yy=subwrd(result,6)
  endif
  if(legend=4)
    'query w2xy ' lon2' 'lat1''
    xx=subwrd(result,3)
    yy=subwrd(result,6)
  endif
  if(legend=3 | legend=4)
    i=1
    lmax=0
    while(i<=k)
      lstr.i=math_strlen(sch.i)
      if(lmax<lstr.i); lmax=lstr.i;endif
      i=i+1
    endwhile
  endif

* Starting legend

  strsize=0.13
  dstr=0.08
  dvet=0.3
  dvet1=0.18
  dhco=0.2
  dlam=0.05
  tysize=0.25
  ybsize=0.1
  'set strsiz 'strsize

  n=1
  while(n<=k)
    'set line 'linec.n' 'lines.n''

    if(legend=1)
      xx1=xx+dvet1
      yy1=yy-n*dvet
   endif
   if(legend=2)
     xx1=xx+dvet1
     yy1=yy+k*dvet-(n-1)*dvet
   endif
   if(legend=3)
     xx1=xx-3*dvet1-(lmax)*strsize
     yy1=yy-n*dvet
   endif
   if(legend=4)
     xx1=xx-3*dvet1-(lmax)*strsize
     yy1=yy+k*dvet-(n-1)*dvet
   endif

   xx2=xx1+dhco
   'draw line 'xx1' 'yy1' 'xx2' 'yy1
   xx3=xx2+dhco
   if (marker.n>0)
     'draw mark 'marker.n' 'xx2' 'yy1' 'tysize/10
   else
     'draw wxsym 41 'xx2' 'yy1' 'tysize
   endif
   xx4=xx3+dlam
   'draw line 'xx2' 'yy1' 'xx3' 'yy1
   yy2=yy1-dstr
   'draw string 'xx4' 'yy2' 'sch.n

   n=n+1
  endwhile

  say ' Plotted!'
* 'print'
* 'disable print'
  'printim 'runname'.png white x1200 y820'
  'quit'
