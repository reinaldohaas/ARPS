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

   HEADER FILE: wxpresrc.h
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

   DESCRIPTION: This header file defines resrources needed by WXP.

******************************************************************/
/*
   WXP resource parameters
*/
#ifdef WXP_RESRC
#define RESRC_KEYED static WXP_deflist def_table = {
#define RESRC_POS ,
#define RESRC_END , NULL };
#define RESRC_INIT WXP_init( argv[0],"Wxp",argc,argv,def_table )

#define FILE_OPTIONS 
#define FILES_OPTIONS 
#define OPTIONAL 
#define OPTIONALS
#define NULL_OPTIONS NULL
#define NULL_COUNT 

#define WXPresrc_age          "age",        "wxpage",       POS1,  CV_int
#define WXPresrc_background   "background", "wxpback",      "bg",  CV_string
#define WXPresrc_batch        "batch",      "wxpbatch",     "ba",  CV_bool
#define WXPresrc_bin_path     "bin_path",   "wxpbin",       "bp",  CV_path
#define WXPresrc_bull_file    "bull_file",  "wxpbullfile",  "bf",  CV_string
#define WXPresrc_city_file    "city_file",  "wxpcity",      "cf",  CV_string
#define WXPresrc_color_cgrid  "color_cgrid","wxpclrcgrid",  "cocg",CV_string
#define WXPresrc_color_clabel "color_clabel","wxpclrclab",  "cocl",CV_string
#define WXPresrc_color_cloud  "color_cloud","wxpclrcld",    "cocd",CV_string
#define WXPresrc_color_cmap   "color_cmap", "wxpclrcmap",   "cocm",CV_string
#define WXPresrc_color_cont   "color_cont", "wxpclrcont",   "coco",CV_string
#define WXPresrc_color_conts  "color_cont", "wxpclrcont",   "coco",CV_string
#define WXPresrc_color_data   "color_data", "wxpclrdata",   "cod", CV_string
#define WXPresrc_color_echo   "color_echo", "wxpclrecho",   "coe", CV_string
#define WXPresrc_color_fill   "color_fill", "wxpclrfill",   "cof", CV_string
#define WXPresrc_color_front  "color_front","wxpclrfront",  "cofr",CV_string
#define WXPresrc_color_grid   "color_grid", "wxpclrgrid",   "cogr",CV_string
#define WXPresrc_color_hodo   "color_hodo","wxpclrhodo",    "coh", CV_string
#define WXPresrc_color_hodtxt "color_hodtxt","wxpclrhdtxt", "coht",CV_string
#define WXPresrc_color_label  "color_label","wxpclrlab",    "cola",CV_string
#define WXPresrc_color_labels "color_label","wxpclrlab",    "cola",CV_string
#define WXPresrc_color_line   "color_line", "wxpclrln",     "coln",CV_string
#define WXPresrc_color_map    "color_map",  "wxpclrmap",    "com", CV_string
#define WXPresrc_color_pop    "color_pop",  "wxpclrpop",    "copp",CV_string
#define WXPresrc_color_rep    "color_rep",  "wxpcolor",     "cor", CV_string
#define WXPresrc_color_sound  "color_sound","wxpclrsnd",    "cosd",CV_string
#define WXPresrc_color_table  "color_table","wxpclrtbl",    "ct",  CV_string
#define WXPresrc_color_text   "color_text", "wxpclrtext",   "cot", CV_string
#define WXPresrc_color_therm  "color_therm","wxpclrtherm",  "coth",CV_string
#define WXPresrc_color_thrmtxt "color_thrmtxt","wxpclrthtxt","cott",CV_string
#define WXPresrc_color_wind   "color_wind", "wxpclrwind",   "cowd",CV_string
#define WXPresrc_color_wx     "color_wx",   "wxpclrwx",     "cowx",CV_string
#define WXPresrc_command      "command",    "wxpcommand",   POS,   CV_strings
#define WXPresrc_commands     "command",    "wxpcommand",   "cd",  CV_string
#define WXPresrc_con_base     "con_base",   "wxpconbase",   "cb",  CV_float
#define WXPresrc_con_bases    "con_base",   "wxpconbase",   "cb",  CV_string
#define WXPresrc_con_interval "con_interval","wxpconint",   "in",  CV_float
#define WXPresrc_con_intervals "con_interval","wxpconint",  "in",  CV_string
#define WXPresrc_con_path     "con_path",   "wxpcon",       "cp",  CV_path
#define WXPresrc_current      "current",    "wxpcurrent",   "cu",  CV_current
#define WXPresrc_data_path    "data_path",  "wxpdata",      "dp",  CV_path
#define WXPresrc_data_src     "data_src",   "wxpdatasrc",   "ds",  CV_string
#define WXPresrc_default      "default",    "wxpdefault",   "df",  CV_string
#define WXPresrc_device       "device",     "wxpdevice",    "de",  CV_device
#define WXPresrc_dewpoint     "dewpoint",   "wxpdewpoint",  POS2,  CV_float
#define WXPresrc_draw         "draw",       "wxpdraw",      "dr",  CV_string
#define WXPresrc_file_param   "file_param", "wxpfileparam", "pf",  CV_string
#define WXPresrc_file_path    "file_path",  "wxpfile",      "fp",  CV_path
#define WXPresrc_filename     "filename",   "wxpfilename",  POS1,  CV_strcpy
#define WXPresrc_filenames    "filename",   "wxpfilename",  POS,   CV_strings
#define WXPresrc_font_list    "font_list",  "wxpfontlist",  "fl",  CV_string
#define WXPresrc_forecast     "forecast",   "wxpforecast",  "fc",  CV_string
#define WXPresrc_geometry     "geometry",   "wxpgeometry",  "ge",  CV_string
#define WXPresrc_grid_num     "grid_num",   "wxpgridnum",   "gn",  CV_string
#define WXPresrc_grid_path    "grid_path",  "wxpgrid",      "gp",  CV_path
#define WXPresrc_grid_type    "grid_type",  "wxpgridtype",  "gt",  CV_int
#define WXPresrc_header_file  "header_file","wxphdrfile",   "fh",  CV_string
#define WXPresrc_height       "height",     "wxpheight",    POS4,  CV_float
#define WXPresrc_help         "help",       "wxphelp",      "h",   CV_bool
#define WXPresrc_hour         "hour",       "wxphour",      "th",  CV_int
#define WXPresrc_icon         "icon",       "wxpicon",      "ic",  CV_bool
#define WXPresrc_identifier   "identifier", "wxpident",     "id",  CV_string
#define WXPresrc_identifiers  "identifier", "wxpident",     "id",  CV_string
#define WXPresrc_identifierp  "identifier", "wxpident",     POS,   CV_strings
#define WXPresrc_in_file      "in_file",    "wxpinfile",    "if",  CV_string
#define WXPresrc_in_files     "in_file",    "wxpinfile",    "if",  CV_string
#define WXPresrc_ingest_port  "ingest_port","wxpingport",   POS1,  CV_string
#define WXPresrc_ingestor     "ingestor",   "wxpingestor",  "di",  CV_ingestor
#define WXPresrc_input        "input",      "wxpinput",     "ip",  CV_string
#define WXPresrc_label        "label",      "wxplabel",     "la",  CV_string
#define WXPresrc_labels       "label",      "wxplabel",     POS,   CV_strings
#define WXPresrc_level        "level",      "wxplevel",     "le",  CV_string
#define WXPresrc_map_file     "map_file",   "wxpmapfile",   "mf",  CV_string
#define WXPresrc_message      "message",    "wxpmessage",   "me",  CV_mess
#define WXPresrc_model        "model",      "wxpmodel",     "mo",  CV_string
#define WXPresrc_name         "name",       "wxpname",      "na",  CV_string
#define WXPresrc_name_conv    "name_conv",  "wxpnameconv",  "nc",  CV_name
#define WXPresrc_num_grid     "num_grid",   "wxpnumgrid",   "ng",  CV_int
#define WXPresrc_num_hour     "num_hour",   "wxpnumhour",   "nh",  CV_string
#define WXPresrc_object_param "object_param","wxpobject",   "oa",  CV_string
#define WXPresrc_out_file     "out_file",   "wxpoutfile",   "of",  CV_string
#define WXPresrc_out_files    "out_file",   "wxpoutfile",   "of",  CV_string
#define WXPresrc_output       "output",     "wxpoutput",    "ou",  CV_string
#define WXPresrc_parameter    "parameter",  "wxpparam",     "pa",  CV_string
#define WXPresrc_plot_domain  "plot_domain","wxppltdomain", "pd",  CV_domain
#define WXPresrc_plot_format  "plot_format","wxppltformat", "fm",  CV_string
#define WXPresrc_plot_param   "plot_param", "wxppltparam",  "pp",  CV_string
#define WXPresrc_plot_range   "plot_range", "wxprange",     "ra",  CV_string
#define WXPresrc_plot_scale   "plot_scale", "wxpscale",     "sc",  CV_float
#define WXPresrc_plot_type    "plot_type",  "wxpplot",      "pl",  CV_string
#define WXPresrc_plot_types   "plot_type",  "wxpplot",      "pl",  CV_string
#define WXPresrc_position     "position",   "wxpposition",  "po",  CV_long
#define WXPresrc_pressure     "pressure",   "wxppressure",  POS3,  CV_float
#define WXPresrc_product      "product",    "wxpproduct",   "ph",  CV_string
#define WXPresrc_program      "program",    "wxpprogram",   "pg",  CV_string
#define WXPresrc_raw_path     "raw_path",   "wxpraw",       "rp",  CV_path
#define WXPresrc_read_conf    "read_conf",  "wxpreadconf",  "rc",  CV_bool
#define WXPresrc_region       "region",     "wxpregion",    "re",  CV_region
#define WXPresrc_save_param   "save_param", "wxpsaveparam", "sp",  CV_string
#define WXPresrc_speed        "speed",      "wxpspeed",     POS2,  CV_float
#define WXPresrc_stat_prior   "stat_prior", "wxpprior",     "pr",  CV_int
#define WXPresrc_state        "state",      "wxpstate",     "st",  CV_string
#define WXPresrc_table        "table",      "wxptable",     "ta",  CV_bool
#define WXPresrc_temp         "temp",       "wxptemp",      POS1,  CV_float
#define WXPresrc_time         "time",       "wxptime",      "ft",  CV_string
#define WXPresrc_title        "title",      "wxptitle",     "ti",  CV_string
#define WXPresrc_unit         "unit",       "wxpunit",      "un",  CV_string
#define WXPresrc_value        "value",      "wxpvalue",     POS1,  CV_float
#define WXPresrc_variable     "variable",   "wxpvariable",  "va",  CV_string
#define WXPresrc_window_num   "window_num", "wxpwindownum", "wn",  CV_int

#define CL_age          "age"
#define CL_background   "-bg clr"
#define CL_batch        "-ba"
#define CL_bin_path     "-bp path"
#define CL_bull_file    "-bf file"
#define CL_city_file    "-cf file"
#define CL_color_cgrid  "-cocg clr"
#define CL_color_clabel "-cocl clr"
#define CL_color_cloud  "-cocd clr"
#define CL_color_cmap   "-cocm clr"
#define CL_color_cont   "-coco clr"
#define CL_color_conts  "-coco clr[,...]"
#define CL_color_data   "-cod clr"
#define CL_color_echo   "-coe clr0,clr1,.."
#define CL_color_fill   "-cof clr[,...]"
#define CL_color_front  "-cofr hclr,lclr,cclr,wclr,oclr,tclr"
#define CL_color_grid   "-cogr clr"
#define CL_color_hodo   "-coh dclr,sclr,hclr"
#define CL_color_hodtxt "-coht dclr,sclr,hclr"
#define CL_color_label  "-cola clr"
#define CL_color_labels "-cola clr[,...]"
#define CL_color_line   "-coln clr"
#define CL_color_map    "-com clr"
#define CL_color_pop    "-copp clr"
#define CL_color_rep    "-cor rep"
#define CL_color_sound  "-cosd tclr,tdclr"
#define CL_color_table  "-ct file"
#define CL_color_text   "-cot clr"
#define CL_color_therm  "-coth pclr,tclr,daclr,maclr,mrclr"
#define CL_color_thrmtxt "-cott pclr,tclr,mrclr"
#define CL_color_wind   "-cowd clr"
#define CL_color_wx     "-cowx clr"
#define CL_command      "command"
#define CL_commands     "-cd cmd[;cmd]"
#define CL_con_base     "-cb base"
#define CL_con_bases    "-cb base[,...]"
#define CL_con_interval "-in intrvl"
#define CL_con_intervals "-in int[,...]"
#define CL_con_path     "-cp path"
#define CL_current      "-cu [hour]"
#define CL_data_path    "-dp path"
#define CL_data_src     "-ds src"
#define CL_default      "-df file"
#define CL_draw         "-dr param[,...]"
#define CL_device       "-de dev[,name]"
#define CL_file_param   "-pf param"
#define CL_file_path    "-fp path"
#define CL_filename     "file"
#define CL_filenames    "file [file...]"
#define CL_font_list    "-fl list"
#define CL_forecast     "-fc type"
#define CL_geometry     "-ge wxh+x+y"
#define CL_grid_num     "-gn number"
#define CL_grid_path    "-gp path"
#define CL_grid_type    "-gt type"
#define CL_help         "-h"
#define CL_hour         "-th hour"
#define CL_icon         "-ic"
#define CL_identifier   "-id ident"
#define CL_identifiers  "-id id[,...]"
#define CL_identifierp  "id [...]"
#define CL_in_file      "-if type"
#define CL_in_files     "-if type[,...]"
#define CL_ingest_port  "port"
#define CL_ingestor     "-di name"
#define CL_input        "-ip type"
#define CL_label        "-la label"
#define CL_labels       "label [label]"
#define CL_level        "-le level"
#define CL_map_file     "-mf file[,file]"
#define CL_message      "-me level"
#define CL_model        "-mo model"
#define CL_name         "-na name"
#define CL_name_conv    "-nc type"
#define CL_num_grid     "-ng number"
#define CL_num_hour     "-nh hour[,step]"
#define CL_object_param "-oa smooth[,radinf,pass,conv,min]"
#define CL_out_file     "-of type"
#define CL_out_files    "-of type[,...]"
#define CL_output       "-ou type"
#define CL_parameter    "-pa param[,...]"
#define CL_plot_domain  "-pd proj,lat1,lon1,nx,ny,dx,dy"
#define CL_plot_format  "-fm format"
#define CL_plot_param   "-pp param"
#define CL_plot_range   "-ra min,max"
#define CL_plot_scale   "-sc scale"
#define CL_plot_type    "-pl type"
#define CL_plot_types   "-pl type[,...]"
#define CL_position     "-po pos"
#define CL_product      "-ph header"
#define CL_raw_path     "-rp path"
#define CL_program      "-pg program"
#define CL_read_conf    "-rc"
#define CL_region       "-re region[,clat,clon][,nx,ny][,size]"
#define CL_save_param   "-sp param[,...]"
#define CL_stat_prior   "-pr prior"
#define CL_state        "-st state"
#define CL_table        "-ta"
#define CL_time         "-ft time"
#define CL_title        "-ti title"
#define CL_unit         "-un unit"  
#define CL_variable     "-va var"
#define CL_window_num   "-wn num"
#endif
#ifdef GUISE_RESRC
/*
   Declare GUISE conversion functions
*/
int Cvt_mess();
int Cvt_current();
int Cvt_name();
int Cvt_ingestor();
int Cvt_device();
int Cvt_region();
int Cvt_domain();
int Cvt_strcpy();
int Cvt_path();

#define RESRC_KEYED static int resrc_itemp; static char *resrc_ptr; static UdKeyTab KeyTab = {
#define RESRC_POS , NULL }; static UdPosTab PosTab = {
#define RESRC_END , NULL };
#define RESRC_INIT if (!udinit(KeyTab, PosTab, &argc, argv, argv[0], "Wxp")) WXP_exit( WNOINIT );

#define FILE_OPTIONS &resrc_itemp, UDOPT
#define FILES_OPTIONS UDOPT
#define OPTIONAL &resrc_itemp, UDOPT
#define OPTIONALS UDOPT
#define NULL_OPTIONS &resrc_ptr, &resrc_itemp
#define NULL_COUNT &resrc_itemp

#define WXPresrc_age          "age",          0,      UdInt
#define WXPresrc_background   "background",   1,      UdString
#define WXPresrc_batch        "batch",        0,      UdBool
#define WXPresrc_bin_path     "bin_path",     1,      Cvt_path
#define WXPresrc_bull_file    "bull_file",    1,      UdString
#define WXPresrc_city_file    "city_file",    1,      UdString
#define WXPresrc_color_cgrid  "color_cgrid",  1,      UdString
#define WXPresrc_color_clabel "color_clabel", 1,      UdString
#define WXPresrc_color_cloud  "color_cloud",  1,      UdString
#define WXPresrc_color_cmap   "color_cmap",   1,      UdString
#define WXPresrc_color_cont   "color_cont",   1,      UdString
#define WXPresrc_color_conts  "color_cont",   UDARB,  UdString
#define WXPresrc_color_data   "color_data",   1,      UdString
#define WXPresrc_color_echo   "color_echo",   UDARB,  UdString
#define WXPresrc_color_fill   "color_fill",   UDARB,  UdString
#define WXPresrc_color_front  "color_front",  UDARB,  UdString
#define WXPresrc_color_grid   "color_grid",   1,      UdString
#define WXPresrc_color_hodo   "color_hodo",   UDARB,  UdString
#define WXPresrc_color_hodtxt "color_hodtxt", UDARB,  UdString
#define WXPresrc_color_label  "color_label",  1,      UdString
#define WXPresrc_color_labels "color_label",  UDARB,  UdString
#define WXPresrc_color_line   "color_line",   1,      UdString
#define WXPresrc_color_map    "color_map",    1,      UdString
#define WXPresrc_color_pop    "color_pop",    1,      UdString
#define WXPresrc_color_rep    "color_rep",    1,      UdString
#define WXPresrc_color_sound  "color_sound",  UDARB,  UdString
#define WXPresrc_color_table  "color_table",  1,      UdString
#define WXPresrc_color_text   "color_text",   1,      UdString
#define WXPresrc_color_therm  "color_therm",  UDARB,  UdString
#define WXPresrc_color_thrmtxt "color_thrmtxt",UDARB, UdString
#define WXPresrc_color_wind   "color_wind",   1,      UdString
#define WXPresrc_color_wx     "color_wx",     1,      UdString
#define WXPresrc_command      "command",      UDARB,  UdString
#define WXPresrc_commands     "command",      1,      UdString
#define WXPresrc_con_base     "con_base",     1,      UdFloat
#define WXPresrc_con_bases    "con_base",     UDARB,  UdFloat
#define WXPresrc_con_interval "con_interval", 1,      UdFloat
#define WXPresrc_con_intervals "con_interval",UDARB,  UdFloat
#define WXPresrc_con_path     "con_path",     1,      Cvt_path
#define WXPresrc_current      "current",      UDARB,  Cvt_current
#define WXPresrc_data_path    "data_path",    1,      Cvt_path
#define WXPresrc_data_src     "data_src",     1,      UdString
#define WXPresrc_device       "device",       UDARB,  UdString
#define WXPresrc_dewpoint     "dewpoint",     1,      UdFloat
#define WXPresrc_draw         "draw",         UDARB,  UdString
#define WXPresrc_file_param   "file_param",   1,      UdString
#define WXPresrc_file_path    "file_path",    1,      Cvt_path
#define WXPresrc_filename     "filename",     1,      Cvt_strcpy
#define WXPresrc_filenames    "filename",     UDARB,  UdString
#define WXPresrc_font_list    "font_list",    UDARB,  UdString
#define WXPresrc_forecast     "forecast",     1,      UdString
#define WXPresrc_geometry     "geometry",     1,      UdString
#define WXPresrc_grid_num     "grid_num",     1,      UdString
#define WXPresrc_grid_path    "grid_path",    1,      Cvt_path
#define WXPresrc_grid_type    "grid_type",    1,      UdInt
#define WXPresrc_height       "height",       1,      UdFloat
#define WXPresrc_help         "help",         0,      UdBool
#define WXPresrc_hour         "hour",         1,      UdInt
#define WXPresrc_icon         "icon",         1,      UdBool
#define WXPresrc_identifier   "identifier",   1,      UdString
#define WXPresrc_identifiers  "identifier",   UDARB,  UdString
#define WXPresrc_identifierp  "identifier",   UDARB,  UdString
#define WXPresrc_in_file      "in_file",      1,      UdString
#define WXPresrc_in_files     "in_file",      UDARB,  UdString
#define WXPresrc_ingest_port  "ingest_port",  1,      UdString
#define WXPresrc_ingestor     "ingestor",     1,      Cvt_ingestor
#define WXPresrc_input        "input",        1,      UdString
#define WXPresrc_label        "label",        UDARB,  UdString
#define WXPresrc_labels       "label",        UDARB,  UdString
#define WXPresrc_level        "level",        1,      UdString
#define WXPresrc_map_file     "map_file",     UDARB,  UdString
#define WXPresrc_message      "message",      1,      Cvt_mess
#define WXPresrc_model        "model",        1,      UdString
#define WXPresrc_name_conv    "name_conv",    1,      Cvt_name
#define WXPresrc_num_grid     "num_grid",     1,      UdInt
#define WXPresrc_num_hour     "num_hour",     UDARB,  UdString
#define WXPresrc_object_param "object_param", UDARB,  UdString
#define WXPresrc_out_file     "out_file",     1,      UdString
#define WXPresrc_out_files    "out_file",     UDARB,  UdString
#define WXPresrc_output       "output",       1,      UdString
#define WXPresrc_parameter    "parameter",    UDARB,  UdString
#define WXPresrc_plot_domain  "plot_domain",  UDARB,  UdString
#define WXPresrc_plot_format  "plot_format",  1,      UdString
#define WXPresrc_plot_param   "plot_param",   UDARB,  UdString
#define WXPresrc_plot_range   "plot_range",   UDARB,  UdFloat
#define WXPresrc_plot_scale   "plot_scale",   1,      UdFloat
#define WXPresrc_plot_type    "plot_type",    1,      UdString
#define WXPresrc_plot_types   "plot_type",    UDARB,  UdString
#define WXPresrc_position     "position",     1,      UdInt
#define WXPresrc_pressure     "pressure",     1,      UdFloat
#define WXPresrc_product      "product",      1,      UdString
#define WXPresrc_program      "program",      1,      UdString
#define WXPresrc_raw_path     "raw_path",     1,      Cvt_path
#define WXPresrc_read_conf    "read_conf",    1,      UdBool
#define WXPresrc_region       "region",       UDARB,  UdString
#define WXPresrc_save_param   "save_param",   UDARB,  UdString
#define WXPresrc_speed        "speed",        1,      UdFloat
#define WXPresrc_stat_prior   "stat_prior",   1,      UdInt
#define WXPresrc_state        "state",        1,      UdString
#define WXPresrc_table        "table",        0,      UdBool
#define WXPresrc_temp         "temp",         1,      UdFloat
#define WXPresrc_time         "time",         1,      UdString
#define WXPresrc_title        "title",        1,      UdString
#define WXPresrc_unit         "unit",         1,      UdString
#define WXPresrc_value        "value",        1,      UdFloat
#define WXPresrc_variable     "variable",     1,      UdString
#define WXPresrc_window_num   "window_num",   1,      UdInt

#define CL_age          "age"
#define CL_background   "-background clr"
#define CL_batch        "-batch"
#define CL_bin_path     "-bin_path path"
#define CL_bull_file    "-bull_file file"
#define CL_city_file    "-city_file file"
#define CL_color_cgrid  "-color_cgrid clr"
#define CL_color_clabel "-color_clabel clr"
#define CL_color_cloud  "-color_cloud clr"
#define CL_color_cmap   "-color_cmap clr"
#define CL_color_cont   "-color_cont clr"
#define CL_color_conts  "-color_cont clr[,...]"
#define CL_color_data   "-color_data clr"
#define CL_color_echo   "-color_echo clr0,clr1,.."
#define CL_color_fill   "-color_fill clr[,...]"
#define CL_color_front  "-color_front hclr,lclr,cclr,wclr,oclr,tclr"
#define CL_color_grid   "-color_grid clr"
#define CL_color_hodo   "-color_hodo dclr,sclr,hclr"
#define CL_color_hodtxt "-color_hodtxt dclr,sclr,hclr"
#define CL_color_label  "-color_label clr"
#define CL_color_labels "-color_label clr[,...]"
#define CL_color_line   "-color_line clr"
#define CL_color_map    "-color_map clr"
#define CL_color_pop    "-color_pop clr"
#define CL_color_rep    "-color_rep rep"
#define CL_color_sound  "-color_sound tclr,tdclr"
#define CL_color_table  "-color_table clr"
#define CL_color_text   "-color_text clr"
#define CL_color_therm  "-color_therm pclr,tclr,daclr,maclr,mrclr"
#define CL_color_thrmtxt "-color_thrmtxt pclr,tclr,mrclr"
#define CL_color_wind   "-color_wind clr"
#define CL_color_wx     "-color_wx clr"
#define CL_command      "command"
#define CL_commands     "-command cmd[;cmd]"
#define CL_con_base     "-con_base base"
#define CL_con_bases    "-con_base base[,...]"
#define CL_con_interval "-con_interval intrvl"
#define CL_con_intervals "-con_interval int[,...]"
#define CL_con_path     "-con_path path"
#define CL_current      "-current [hour]"
#define CL_data_path    "-data_path path"
#define CL_data_src     "-data_src src"
#define CL_default      "-default file"
#define CL_draw         "-draw param[,...]"
#define CL_device       "-device dev[,name]"
#define CL_file_param   "-file_param param"
#define CL_file_path    "-file_path path"
#define CL_filename     "file"
#define CL_filenames    "file [file...]"
#define CL_font_list    "-font_list list"
#define CL_forecast     "-forecast type"
#define CL_geometry     "-geometry wxh+x+y"
#define CL_grid_num     "-grid_num num"
#define CL_grid_path    "-grid_path path"
#define CL_grid_type    "-grid_type type"
#define CL_help         "-help"
#define CL_hour         "-hour hour"
#define CL_icon         "-icon"
#define CL_identifier   "-identifier ident"
#define CL_identifiers  "-identifier ident[,...]"
#define CL_identifierp  "ident [...]"
#define CL_in_file      "-in_file type"
#define CL_in_files     "-in_file type[,...]"
#define CL_ingest_port  "port"
#define CL_ingestor     "-ingestor name"
#define CL_input        "-input type"
#define CL_label        "-label label"
#define CL_labels       "label [label]"
#define CL_level        "-level level"
#define CL_map_file     "-map_file file[,file]"
#define CL_message      "-message level"
#define CL_model        "-model model"
#define CL_name         "-name name" 
#define CL_name_conv    "-name_conv type"
#define CL_num_grid     "-num_grid num"
#define CL_num_hour     "-num_hour hour[,step]"
#define CL_object_param "-object_param smooth[,radinf,pass,conv,min]"
#define CL_out_file     "-out_file type"
#define CL_out_files    "-out_file type[,...]"
#define CL_output       "-output type"
#define CL_parameter    "-parameter param[,...]"
#define CL_plot_domain  "-plot_domain proj,lat1,lon1,nx,ny,dx,dy"
#define CL_plot_format  "-plot_format format"
#define CL_plot_param   "-plot_param param"
#define CL_plot_range   "-plot_range min,max"
#define CL_plot_scale   "-plot_scale scale"
#define CL_plot_type    "-plot_type type"
#define CL_plot_types   "-plot_type type[,...]"
#define CL_position     "-position pos"
#define CL_product      "-product header"
#define CL_program      "-program program"
#define CL_raw_path     "-raw_path path"
#define CL_read_conf    "-read_conf"
#define CL_region       "-region region[,clat,clon][,nx,ny][,size]"
#define CL_save_param   "-save_param param[,...]"
#define CL_stat_prior   "-stat_prior prior"
#define CL_state        "-state state"
#define CL_table        "-table"  
#define CL_time         "-time time"  
#define CL_title        "-title title"  
#define CL_unit         "-unit unit"  
#define CL_variable     "-variable var"
#define CL_window_num   "-window_num num"
#endif
