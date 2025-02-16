TITLE   Soil initial data set
*
DSET    runname.soilvar
FILEHEADER 192
OPTIONS sequential big_endian
UNDEF   -9.e+33
*
XDEF      67  LINEAR       1    1
YDEF      67  LINEAR       1    1
ZDEF       1  LINEAR       0    1
TDEF       1  LINEAR  15:00z20May1977   1MN
*
VARS   05
tsfc    0      99   Ground surface temperature (K)
tsoil   0      99   Deep soil temperature (K)
wetsfc  0      99   Surface soil moisture
wetdp   0      99   Deep soil moisture
wetcnpy 0      99   Water amount on canopy
ENDVARS
