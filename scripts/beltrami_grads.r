reinit
open &0.gradscntl
enable print &0_&1.meta
set mproj off
set grid off
*
* Define the exact solution
*
set x 1 91
set y 1 61
set z 1 47
define a = 2
define nu = 1
define lx = 267
define ly = 118
define lz = 44
define pi = 3.141596254
define kx = 4*pi/lx
define ky = 2*pi/ly
define kz = 2*pi/lz
define lambda = sqrt(kx*kx+ky*ky+kz*kz)
define ue0 = (-a/(kx*kx+ky*ky))*(lambda*ky*cos(kx*lon*1000)*sin(ky*lat*1000)*sin(kz*lev)+kz*kx*sin(kx*lon*1000)*cos(ky*lat*1000)*cos(kz*lev))
define ve0 = ( a/(kx*kx+ky*ky))*(lambda*kx*sin(kx*lon*1000)*cos(ky*lat*1000)*sin(kz*lev)-kz*ky*cos(kx*lon*1000)*sin(ky*lat*1000)*cos(kz*lev))
define we0 =  a*cos(kx*lon*1000)*cos(ky*lat*1000)*sin(kz*lev)
*
* Plot for t=&1
*
set t &1
*
define b = exp(-nu*lambda*lambda*(&1-1)*41)
define ue = ue0*b
define ve = ve0*b
define we = we0*b
*
set lev 0
set vpage 1.0 4.5 6.0 8.0
set grads off
set arrscl 0.4 2
d skip(ue,2,2);skip(ve,2,2)
draw title Exact u-v field at z = 0 m
*
set vpage 6.0 9.5 6.0 8.0
set grads off
set arrscl 0.4 1
d skip(u-ue,2,2);skip(v-ve,2,2)
draw title Difference u-v field (predicted - exact) at z = 0 m
*
set lev 5
set vpage 1.0 4.5 4.0 6.0
set grads off
set arrscl 0.4 2
d skip(ue,2,2);skip(ve,2,2)
draw title Exact u-v field at z = 5 m
*
set vpage 6.0 9.5 4.0 6.0
set grads off
set arrscl 0.4 1
d skip(u-ue,2,2);skip(v-ve,2,2)
draw title Difference u-v field (predicted - exact) at z = 5 m
*
set lev 10
set vpage 1.0 4.5 2.0 4.0
set grads off
set arrscl 0.4 2
d skip(ue,2,2);skip(ve,2,2)
draw title Exact u-v field at z = 10 m
*
set vpage 6.0 9.5 2.0 4.0
set grads off
set arrscl 0.4 1
d skip(u-ue,2,2);skip(v-ve,2,2)
draw title Difference u-v field (predicted - exact) at z = 10 m
*
set vpage 1.0 4.5 0.0 2.0
set grads off
d we
draw title Exact w field at z = 10 m
*
set vpage 6.0 9.5 0.0 2.0
set grads off
d w-we
draw title Difference w field (predicted - exact) at z = 10 m
*
set vpage off
set string 1 bc 6
set strsiz 0.2
draw string 5.5 8.2 Beltrami flow at t=(&1-1)*41s
*
print
disable
*
!gxps -c -i &0_&1.meta -o &0_&1.ps
!rm -f &0_&1.meta
*
quit
