reinit
open &0.sfcctl
set grid off
*
set x 2 2
set y 2 2
set z 2 2 
set t 1 289
*
d rn
d h
d le
d g
*
draw title &0 ARPS version &1
*
gxprint  &0.ps
*
quit
