  INTEGER :: nprio   ! Maximum number of priority levels allowed
  PARAMETER (nprio=200)

  INTEGER :: max_verts, max_polys ! number of vertices, number of polygones
  PARAMETER (max_verts=20, max_polys=10)

  INTEGER :: maxarbvar   ! maximum of arbitrary variables allowed
  PARAMETER (maxarbvar=50)

  INTEGER :: maxunevm   ! maximum number of uneven contours
  INTEGER :: maxuneva   ! maximum number of uneven variables
  PARAMETER (maxunevm=32, maxuneva=50)


  INTEGER :: mxsfcobfl ! maximum number of surface obs files
  PARAMETER (mxsfcobfl=20)
  INTEGER :: mxsfcob   ! maximum number of surface observations
  PARAMETER (mxsfcob=9000)
  INTEGER :: mxstalo   ! maximum number of observation stations
  PARAMETER (mxstalo=9000)

  INTEGER :: nmax      ! maximum number of map data points to be
                       ! processed once by lltoxy
  INTEGER :: maxline   ! maximum number of line segments allowed in the
                       ! map data
  INTEGER :: maxmap    ! maximum number of map data files to be plotted
  PARAMETER (nmax = 100, maxline=50000,maxmap=10)
