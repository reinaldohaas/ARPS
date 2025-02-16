reinit
open &0.gradscntl
enable print &0.meta
set mproj off
set grid off
*
set lon 0 16
set y 2
set lev 0 2400
*
* Plot for t=450 seconds
*
set t 4
set vpage 0.5 5.5 6.0 8.0
set grads off
set xaxis 0 16 2
set yaxis 0 24 4
set cmax -0.001
d ptprt
draw title ptprt at t=450s
*
set vpage 0.5 5.5 4.0 6.0
set grads off
set xaxis 0 16 2
set yaxis 0 24 4
d u
draw title u
*
set vpage 0.5 5.5 2.0 4.0
set grads off
set xaxis 0 16 2
set yaxis 0 24 4
d w
draw title w
*
set vpage 0.5 5.5 0.0 2.0
set grads off
set xaxis 0 16 2
set yaxis 0 24 4
d pprt
draw title pprt
*
* Plot for t=900 seconds
*
set t 7
*
set vpage 5.5 11.0 6.0 8.0
set grads off
set xaxis 0 16 2
set yaxis 0 24 4
set cmax -0.001
d ptprt
draw title ptprt at t=900s
*
set vpage 5.5 11.0 4.0 6.0
set grads off
set xaxis 0 16 2
set yaxis 0 24 4
d u
draw title u
*
set vpage 5.5 11.0 2.0 4.0
set grads off
set xaxis 0 16 2
set yaxis 0 24 4
d w
draw title w
*
set vpage 5.5 11.0 0.0 2.0
set grads off
set xaxis 0 16 2
set yaxis 0 24 4
d pprt
draw title pprt
*
set vpage off
set string 1 bc 6
set strsiz 0.2
draw string 5.8 8.2 Density Currecnt for &0
*
print
disable
*
!gxps -i &0.meta -o &0.ps
!rm -f &0.meta
quit
