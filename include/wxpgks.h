/*********************************************************************
**********************************************************************
*                                                                    *
*                        C O P Y R I G H T                           *
*                                                                    *
*                    Copyright 1988,1989,1990 (C) by                 *
*                    Purdue Research Foundation                      *
*                    Purdue University                               *
*                    West Lafayette, Indiana 47907                   *
*                                                                    *
* This software,  in whole or in part,  may be used  and copied only *
* with the written permission of the  Dept. of Earth and Atmospheric *
* Sciences  via  E.  M. Agee,  Purdue  University,  West  Lafayette, *
* Indiana,  and with the  inclusion  of the above  copyright notice. *
* This  software  or  any  copies  thereof  may  not be  provided or *
* otherwise made  available  to any  other person.  No title  to and *
* ownership of the software is hereby transferred.                   *
*                                                                    *
**********************************************************************

   HEADER FILE: wxpgks.h
   PROGRAMMER: Dan Vietor
   WXP VERSION: 4.8
   DATE: 931120

   COMPUTER:
      IBM RS/6000 running AIX 3.1 C Compiler
      IBM RT/PC running AIX 2.2.1 C Compiler
      IBM RT/PC running BSD 4.3 C Compiler
      IBM AT compatible running MSDOS 3.x using the
         Microsoft C 6.0 Compiler.
      SUN 3 and 4 running SUNOS 4.1
      DEC VaxStation running VMS 5.3 and GKS
      Interprocess communications based on either ATT UNIX System V
         Message queues or BSD Sockets.
      Asynchronous data stream function interface based on either
         System V ioctl, BSD 4.3 ioctl or Blaise Asynch Manager
         for MSDOS 3.x.
      Graphics calls based on ANSI C language Graphical Kernel
         System binding (revision dated July 1988)

   DESCRIPTION: This header file defines and declares structures
      needed by graphics calls within WXP.

******************************************************************/
#include "gks.h"

#define GDISPLAY 1
#define GMETAOUT 3
#define GPRINTER 4
#define GPLOTTER 5
#define GWINDOW  6
#ifdef GKS88
#define FONT1 1
#define FONT2 1
#define SPACING .1
#define LINEWIDTH 0.
#define M_SCALE .75
#endif
#ifdef XGKS
#define FONT1 1
#define FONT2 1
#define SPACING .25
#define M_SCALE 2.5
#define LINEWIDTH 1.
#define GKS86
#endif
#ifdef VMSGKS
#define FONT1 1
#define FONT2 1
#define SPACING .1
#define M_SCALE 2.5
#define LINEWIDTH 1.
#define GKS86
#endif

#ifdef X11GKS
#define GKS_EXTRA
#endif
#ifdef GSSGKS
#define GKS_EXTRA
#endif
#ifdef MSCGKS
#define GKS_EXTRA
#endif

#define MAX_CLR 75
#define MAX_FILL_CLR 60

#define COLOR_OFF -1
#define COLOR_END -2

#define DRAW_TEXT 0
#define DRAW_DATA 1
#define DRAW_MAP  2

#define DUMP_XWD  1
#define DUMP_GIF  2

/* Structure for graphics attributes */
typedef struct {
   Gint color;
   Gfloat width;
   Gfloat height;
   Gfloat expan;
   Gint style;
   Gint fill_pat;
   Gint font;
   } Gattribs;

/* Structure for the aspect ratio */
typedef struct {
   Gfloat x;
   Gfloat y;
   } Gaspect;

/* Structure for centering */
typedef struct {
   Gint x;
   Gint y;
   } Gcenter;

#ifdef X11GKS
#include <X11/Xlib.h>
Display *gget_xdisplay();
int gget_xscreen();
Window gget_xwindow();
Colormap gget_xcolormap();
GC gget_xgc();
int gget_xwindow_height();
int gget_xwindow_width();
FILE *gget_xoutfile();
#endif
#ifndef GKS_EXTRA
#define gset_text_width()
#define gset_expan()
#define gset_mouse_callback()
#endif

#ifdef GKS86
#define gclose_gks gclosegks
#define gemergency_close_gks gemergencyclosegks
#define gclose_ws gclosews
#define gactivate_ws gactivatews
#define gdeactivate_ws gdeactivatews
#define gclear_ws gclearws
#define gcreate_seg gcreateseg
#define gclose_seg gcloseseg
#define gcell_array gcellarray
#define gfill_area gfillarea
#define gset_linetype gsetlinetype
#define gset_linewidth gsetlinewidth
#define gset_line_colr_ind gsetlinecolourind
#define gset_marker_type gsetmarkertype
#define gset_marker_size gsetmarkersize
#define gset_marker_colr_ind gsetmarkercolourind
#define gset_fontprec gsettextfontprec
#define gset_space gsetcharspace
#define gset_text_colr_ind gsettextcolourind
#define gset_char_ht gsetcharheight
#define gset_char_up_vec gsetcharup
#define gset_text_path gsettextpath
#define gset_text_align gsettextalign
#define gset_fill_int_style gsetfillintstyle
#define gset_fill_style_ind gsetfillstyleind
#define gset_fill_colr_ind gsetfillcolourind
#define gset_colr_rep gsetcolourrep
#define gset_win gsetwindow
#define gset_vp gsetviewport
#define gsel_norm_tran gselntran
#define gset_ws_win gsetwswindow
#define gset_ws_vp gsetwsviewport
#define gcreate_seg gcreateseg
#define gclose_seg gcloseseg
#define gset_string_mode gsetstringmode
#define ginq_cur_norm_tran_num ginqcurntrannum
#endif
