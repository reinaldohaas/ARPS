/*
 * dorade.c --
 *
 *        This file defines functions that access and manage data in DORADE sweep
 *        files.
 *
 *        The struct's are based on those found in the soloii source distribution
 *        from National Center for Atmospheric Research (NCAR), Earth Observing
 *        Laboratory, Research Data Program (RDP), downloaded from
 *        ftp.atd.ucar.edu.
 *
 * Copyright (c) 2007 Gordon D. Carrie
 *
 * Licensed under the Open Software License version 2.1
 *
 * Please send feedback to user0@tkgeomap.org
 *
 * @(#) $Id: dorade_lib.c,v 4.0 2013/02/08 18:31:44 ywang Exp $
 *
 **********************************************************************
 *
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include "errMsg.h"
#include "getVals.h"
#include "dorade.h"

/*
 * This macro adjusts a float value so that when it truncates, value rounds
 * to the nearest integer.
 */

#define N_INT(f) ((f) > 0.0 ? (f) + 0.5 : (f) - 0.5)

/*
 * Macros needed for decompression
 */

#define MASK15 0x7fff
#define SIGN16 0x8000

/*
 * Declarations for local functions.
 */

static int angle_cmp(const void *a1, const void *a2);

/*
 * This array temporarily stores information for new command sets.
 */

static struct Hash_Info *info;

/*
 * Input buffer
 */

static char *inBufP;                /* Input buffer */
static size_t in_sz;                /* Size of allocation at inBufP */

char *strhead(char *,int);
/*
 *------------------------------------------------------------------------
 *
 * Dorade_InitSweep --
 *
 *         This function initializes a new DORADE sweep structure.
 *
 * Arguments:
 *        struct Dorade_Sweep *swpP        - Sweep structure.  The structure
 *                                          should already be allocated.  Its
 *                                          contents should be garbage.
 *
 * Results:
 *         None
 *
 * Side effects:
 *         The structure is filled in with bogus values.
 *
 *------------------------------------------------------------------------
 */

void Dorade_InitSweep(struct Dorade_Sweep *swpP)
{
    char blank[22] = "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";

    strncpy(swpP->radar_name, blank, 9);
    strncpy(swpP->proj_name, blank, 21);
    swpP->n_parms = 0;
    swpP->n_use_parms = 0;
    swpP->vol_num = 0;
    swpP->time.day = 0;
    swpP->time.second = 0;
    swpP->n_sensors = 0;
    swpP->radar_const = 0.0;
    swpP->peak_power = 0.0;
    swpP->noise_power = 0.0;
    swpP->receiver_gain = 0.0;
    swpP->antenna_gain = 0.0;
    swpP->system_gain = 0.0;
    swpP->horz_beam_width = 0;
    swpP->vert_beam_width = 0;
    swpP->radar_type = 0;
    swpP->scan_mode = DORADE_CALIBRATION;
    swpP->req_rotat_vel = 0.0;
    swpP->data_reduction = 1;
    swpP->data_red_parm0 = 0.0;
    swpP->data_red_parm1 = 0.0;
    swpP->radar_location = GeoPtNowhere();
    swpP->radar_altitude = 0.0;
    swpP->eff_unamb_vel = 0.0;
    swpP->eff_unamb_range = 0.0;
    swpP->num_freq_trans = 1;
    swpP->num_ipps_trans = 1;
    swpP->freq1 = 0.0;
    swpP->freq2 = 0.0;
    swpP->freq3 = 0.0;
    swpP->freq4 = 0.0;
    swpP->freq5 = 0.0;
    swpP->interpulse_per1 = 0.0;
    swpP->interpulse_per2 = 0.0;
    swpP->interpulse_per3 = 0.0;
    swpP->interpulse_per4 = 0.0;
    swpP->interpulse_per5 = 0.0;
    swpP->n_cells = 0;
    swpP->distP = NULL;
    swpP->parmP = NULL;
    Hash_InitTable(&swpP->parm_tbl);
    swpP->azimuth_corr = AngleFmDeg(0.0);
    swpP->elevation_corr = AngleFmDeg(0.0);
    swpP->range_delay_corr = 0.0;
    swpP->longitude_corr = AngleFmDeg(0.0);
    swpP->latitude_corr = AngleFmDeg(0.0);
    swpP->pressure_alt_corr = 0.0;
    swpP->radar_alt_corr = 0.0;
    swpP->ew_gndspd_corr = 0.0;
    swpP->ns_gndspd_corr = 0.0;
    swpP->vert_vel_corr = 0.0;
    swpP->heading_corr = AngleFmDeg(0.0);
    swpP->roll_corr = AngleFmDeg(0.0);
    swpP->pitch_corr = AngleFmDeg(0.0);
    swpP->drift_corr = AngleFmDeg(0.0);
    swpP->rot_angle_corr = AngleFmDeg(0.0);
    swpP->tilt_corr = AngleFmDeg(0.0);
    swpP->sweep_num = 0;
    swpP->n_rays = 0;
    swpP->n_good_rays = 0;
    swpP->start_angle = 0;
    swpP->stop_angle = 0;
    swpP->fixed_angle = 0;
    swpP->filter_flag = 0;
    swpP->rayHdrP = NULL;
    swpP->dat = NULL;
    swpP->interpolated = 0;
}

/*
 *------------------------------------------------------------------------
 *
 * Dorade_AllocDat --
 *
 *         This function allocates memory for the dat member of
 *         a Dorade_Sweep structure.
 *
 * Arguments:
 *        long n_parms                - number of parameters
 *        long n_rays                - number of rays
 *        unsigned n_cells        - number of cells per ray
 *
 * Results:
 *         Return value is a three dimensional array of float values.
 *
 * Side effects:
 *         Memory is allocated.  It will be freed with a call to Dorade_FreeSweep.
 *         If something goes wrong, an error message is stored with a call to
 *         ErrMsg_Append.
 *
 *------------------------------------------------------------------------
 */

float ***Dorade_AllocDat(long n_parms, long n_rays, unsigned n_cells)
{
    float ***dat = NULL;
    unsigned p, r;
    size_t ray_psz, cell_psz;

    dat = (float ***)MALLOC(n_parms * sizeof(float **));
    for (p = 0; p < n_parms; p++) {
        dat[p] = NULL;
    }
    ray_psz = n_rays * sizeof(float *);
    cell_psz = n_rays * n_cells * sizeof(float);
    for (p = 0; p < n_parms; p++) {
        dat[p] = (float **)MALLOC(ray_psz);
        dat[p][0] = (float *)MALLOC(cell_psz);
        for (r = 1; r < n_rays; r++) {
            dat[p][r] = dat[p][r - 1] + n_cells;
        }
    }
    return dat;
}

/*
 *------------------------------------------------------------------------
 *
 * Dorade_FreeSweep --
 *
 *         This function frees internal allocations in  a DORADE sweep structure.
 *
 * Arguments:
 *        struct Dorade_Sweep *swpP        - Sweep structure initialized with
 *                                          Dorade_InitSweep
 * Results:
 *         None
 *
 * Side effects:
 *         Dynamic memory associated with the structure, if any, is freed.
 *         The space for structure itself remains allocated.
 *         The sweep structure is filled in with bogus values.
 *
 *------------------------------------------------------------------------
 */

void Dorade_FreeSweep(struct Dorade_Sweep *swpP)
{
    int p;

    if ( !swpP ) {
        return;
    }
    FREE(swpP->distP);
    FREE(swpP->parmP);
    FREE(swpP->rayHdrP);
    if (swpP->dat) {
        for (p = 0; p < swpP->n_parms; p++) {
            if (swpP->dat[p]) {
                FREE(swpP->dat[p][0]);
            }
            FREE(swpP->dat[p]);
        }
        FREE(swpP->dat);
    }
    Hash_ClearTable(&swpP->parm_tbl);
    Dorade_InitSweep(swpP);
}

/*
 *------------------------------------------------------------------------
 *
 * Dorade_InitParm --
 *
 *         This function initializes a new Dorade_ParmDesc structure.
 *
 * Arguments:
 *        struct Dorade_ParmDesc *parmP        - Parameter descriptor structure.
 *                                          The structure should already be
 *                                          allocated.  Its contents should
 *                                          be garbage.
 * Results:
 *         None
 *
 * Side effects:
 *         The structure is filled in with bogus values.
 *
 *------------------------------------------------------------------------
 */

void Dorade_InitParm(struct Dorade_ParmDesc *parmP)
{
    char blank[42] = "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
        "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";

    parmP->use = 0;
    strncpy(parmP->name, blank, 9);
    strncpy(parmP->description, blank, 41);
    strncpy(parmP->units, blank, 9);
    parmP->interpulse_time = 0;
    parmP->xmitted_freq = 0;
    parmP->recvr_bandwidth = 0.0;
    parmP->pulse_width = 0;
    parmP->polarization = 0;
    parmP->num_samples = 0;
    parmP->binary_format = 0;
    strncpy(parmP->threshold_field, "NONE", 9);
    parmP->threshold_value = 0.0;
    parmP->scale = 0.0;
    parmP->bias = 0.0;
    parmP->bad_data = 0;
    parmP->min = 0.0;
}

/*
 *------------------------------------------------------------------------
 *
 * Dorade_InitRayHdr --
 *
 *         This function initializes a new Dorade ray header structure.
 *
 * Arguments:
 *        struct Dorade_RayHdr *rayHdrP        - ray header to initialize.
 *                                          The structure should already be
 *                                          allocated.  Its contents should
 *                                          be garbage.
 * Results:
 *         None
 *
 * Side effects:
 *         The structure is filled in with bogus values.
 *
 *------------------------------------------------------------------------
 */

void Dorade_InitRayHdr(struct Dorade_RayHdr *rayHdrP)
{
    rayHdrP->good = 0;
    rayHdrP->time.day = 0;
    rayHdrP->time.second = 0.0;
    rayHdrP->azimuth = AngleFmDeg(0.0);;
    rayHdrP->elevation = AngleFmDeg(0.0);;
    rayHdrP->peak_power = 0.0;;
    rayHdrP->true_scan_rate = AngleFmDeg(0.0);;
    rayHdrP->ray_status = 0;
    rayHdrP->longitude = AngleFmDeg(0.0);;
    rayHdrP->latitude = AngleFmDeg(0.0);;
    rayHdrP->altitude_msl = 0.0;;
    rayHdrP->altitude_agl = 0.0;;
    rayHdrP->ew_velocity = 0.0;;
    rayHdrP->ns_velocity = 0.0;;
    rayHdrP->vert_velocity = 0.0;;
    rayHdrP->heading = AngleFmDeg(0.0);;
    rayHdrP->roll = AngleFmDeg(0.0);;
    rayHdrP->pitch = AngleFmDeg(0.0);;
    rayHdrP->drift_angle = AngleFmDeg(0.0);;
    rayHdrP->rotation_angle = AngleFmDeg(0.0);;
    rayHdrP->tilt = AngleFmDeg(0.0);;
    rayHdrP->ew_horiz_wind = 0.0;;
    rayHdrP->ns_horiz_wind = 0.0;;
    rayHdrP->vert_wind = 0.0;;
    rayHdrP->heading_change = AngleFmDeg(0.0);;
    rayHdrP->pitch_change = AngleFmDeg(0.0);;
}

/*
 *------------------------------------------------------------------------
 *
 * Dorade_HashParms --
 *
 *         This function sets up a hash table that maps parameter names
 *         to their structures in the parameter array.
 *
 * Arguments:
 *         swpP        - a DORADE sweep.
 *
 * Results:
 *         Return value is a success or failure code.  The parm_tbl member of
 *         the sweep structure is populated.  It will be cleared when the
 *         sweep is freed.
 *
 *------------------------------------------------------------------------
 */

int Dorade_HashParms(struct Dorade_Sweep *swpP)
{
    struct Hash_Info *ip;
    struct Dorade_ParmDesc *pp;
    struct Dorade_ParmDesc *parmP;
    long n_use_parms, n_parms;
    size_t sz;

    if ( !swpP ) {
        ErrMsg_Append("Tried to create parameter table for nonexistent sweep.");
        return 0;
    }
    parmP = swpP->parmP;
    n_parms = swpP->n_parms;
    n_use_parms = swpP->n_use_parms;
    sz = swpP->n_use_parms * sizeof(struct Hash_Info);
    info = (struct Hash_Info *)REALLOC(info, sz);
    for (pp = parmP, ip = info; pp < parmP + n_parms; pp++) {
        if (pp->use) {
            ip->key = pp->name;
            ip->val = pp - parmP;
            ip++;
        }
    }
    if ( !Hash_FillTable(&swpP->parm_tbl, n_use_parms, info) ) {
        ErrMsg_Append("Failed to create sweep parameters table.");
        return 0;
    }
    return 1;
}

/*
 * Dorade_GetParm --
 *
 *        This access function returns the parameter structure for a sweep
 *        parameter given the parameter name.
 *
 * Arguments:
 *        struct Dorade_Sweep *swpP        - DORADE sweep
 *        char *parm                        - name of parameter
 *
 * Results:
 *        Return value is the address of the parameter in the sweep parameter
 *        array, or NULL if there is no parameter named parm in the sweep.
 */

struct Dorade_ParmDesc * Dorade_GetParm(struct Dorade_Sweep *swpP, char *parm)
{
    struct Dorade_ParmDesc *parmP;
    long l;

    if ( !swpP ) {
        return NULL;
    }
    if ( Hash_GetVal(&swpP->parm_tbl, parm, &l) ) {
        parmP = swpP->parmP + l;
        if (parmP->use) {
            return parmP;
        } else {
            ErrMsg_Append(parm);
            ErrMsg_Append(" not in use.\n");
            return NULL;
        }
    } else {
        ErrMsg_Append("No parameter named ");
        ErrMsg_Append(parm);
        ErrMsg_Append(" in sweep.\n");
    }
    return NULL;
}

/*
 *------------------------------------------------------------------------
 *
 * Dorade_ReadSweepHeaders --
 *
 *         This function reads and stores the header information from
 *        a DORADE sweep file.  It skips ray headers and data.
 *
 * Arguments:
 *        FILE *in                        - Input file or stream. It should be
 *                                          positioned at the start of a DORADE
 *                                          sweep file.
 *        struct Dorade_Sweep *swpP        - Place to store values, should have
 *                                          been initialized with a call to
 *                                          Dorade_InitSweep.
 *
 * Results:
 *         Return value is non-zero if successful, zero if something goes wrong.
 *
 * Side effects:
 *         The sweep structure is updated.  Memory might be allocated, which
 *         should be eventually freed with a call to Dorade_FreeSweep.
 *         The input stream remains open, although it will be repositioned to
 *         the end of the sweep.
 *
 *------------------------------------------------------------------------
 */

int Dorade_ReadSweepHeaders(FILE *in, struct Dorade_Sweep *swpP)
{
    char desc_id[4];                        /* Descriptor identifier */
    size_t desc_len;                        /* Descriptor length */
    size_t front = 4 + sizeof(long);        /* Size of id and size specifier */
    int start;                                /* If true, stream is at start
                                         * of sweep */

    /*
     * Convenience variables
     */

    int r;
    struct Dorade_ParmDesc *parmP;
    size_t sz;

    /*
     * Loop parameters
     */

    float *fInP;
    float *fP, *f1P;
    char *c, *e, *end;

    swpP->n_parms = 0;
    r = 0;

    if (front > in_sz) {
        inBufP = (char *)REALLOC(inBufP, front);
        in_sz = front;
    }
    start = 1;
    while ( !feof(in) ) {
        /*
         * Determine type and size of block.
         */

        if (fread(inBufP, 1, front, in) != front) {
            if ( !feof(in) ) {
                ErrMsg_Append("End of file before end of block\n");
                goto error;
            }
            return 1;
        }
        strncpy(desc_id, inBufP, 4);
        desc_len = GetVals_GetSINT32(inBufP + 4);
        if (desc_len < 0) {
            ErrMsg_Append("Negative descriptor length.\n");
            goto error;
        }
        if (start) {
            if (strncmp(desc_id, "SSWB", 4) != 0
                    && strncmp(desc_id, "COMM", 4) != 0
                    && strncmp(desc_id, "VOLD", 4) != 0) {
                ErrMsg_Append("Not a Dorade file.\n");
                goto error;
            }
            start = !start;
        }
        if (desc_len > 0) {
            if (desc_len > in_sz) {
                inBufP = (char *)REALLOC(inBufP, desc_len);
                in_sz = desc_len;
            }
            desc_len -= front;
            if (fread(inBufP + front, 1, desc_len, in) != desc_len) {
                ErrMsg_Append("Failed to read descriptor block\n");
                goto error;
            }
        }

        /*
         * Read and store block contents.
         */

        if (strncmp(desc_id, "SSWB", 4) == 0) {
            time_t start_time = GetVals_GetSINT32(inBufP + 12);
            struct tm *tm = gmtime(&start_time);
            swpP->time = GeoTime_CalToJul(
                    GeoTime_CalSet(tm->tm_year + 1900, tm->tm_mon + 1,
                        tm->tm_mday, tm->tm_hour, tm->tm_min, tm->tm_sec));
        } else if (strncmp(desc_id, "VOLD", 4) == 0) {
            swpP->vol_num = GetVals_GetSINT16(inBufP + 10);
            strncpy(swpP->proj_name, inBufP + 16, 20);
            swpP->n_sensors = GetVals_GetSINT16(inBufP + 70);
        } else if (strncmp(desc_id, "RADD", 4) == 0) {
            strncpy(swpP->radar_name, inBufP + 8, 8);

            /*
             * Remove trailing spaces from radar_name.
             */

            for (c = e = swpP->radar_name, end = swpP->radar_name + 8;
                    c < end; c++) {
                if ( !isspace((int)*c) ) {
                    e = c;
                }
            }
            if (e + 1 != end) {
                *(e + 1) = '\0';
            }

            swpP->radar_const = GetVals_GetFLOAT32(inBufP + 16);
            swpP->peak_power = GetVals_GetFLOAT32(inBufP + 20);
            swpP->noise_power = GetVals_GetFLOAT32(inBufP + 24);
            swpP->receiver_gain = GetVals_GetFLOAT32(inBufP + 28);
            swpP->antenna_gain = GetVals_GetFLOAT32(inBufP + 32);
            swpP->system_gain = GetVals_GetFLOAT32(inBufP + 36);
            swpP->horz_beam_width
                = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 40));
            swpP->vert_beam_width
                = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 44));
            swpP->radar_type = GetVals_GetSINT16(inBufP + 48);
            swpP->scan_mode
                = (enum Dorade_ScanMode)GetVals_GetUINT16(inBufP + 50);
            swpP->req_rotat_vel = GetVals_GetFLOAT32(inBufP + 52);
            swpP->compression = GetVals_GetSINT16(inBufP + 68);
            swpP->data_reduction = GetVals_GetSINT16(inBufP + 70);
            swpP->data_red_parm0 = GetVals_GetFLOAT32(inBufP + 72);
            swpP->data_red_parm1 = GetVals_GetFLOAT32(inBufP + 76);
            swpP->radar_location = GeoPtFmDeg(
                    GetVals_GetFLOAT32(inBufP + 84),
                    GetVals_GetFLOAT32(inBufP + 80));
            swpP->radar_altitude = GetVals_GetFLOAT32(inBufP + 88);
            swpP->eff_unamb_vel = GetVals_GetFLOAT32(inBufP + 92);
            swpP->eff_unamb_range = GetVals_GetFLOAT32(inBufP + 96);
            swpP->num_freq_trans = GetVals_GetSINT16(inBufP + 100);
            swpP->num_ipps_trans = GetVals_GetSINT16(inBufP + 102);
            swpP->freq1 = GetVals_GetFLOAT32(inBufP + 104);
            swpP->freq2 = GetVals_GetFLOAT32(inBufP + 108);
            swpP->freq3 = GetVals_GetFLOAT32(inBufP + 112);
            swpP->freq4 = GetVals_GetFLOAT32(inBufP + 116);
            swpP->freq5 = GetVals_GetFLOAT32(inBufP + 120);
            swpP->interpulse_per1 = GetVals_GetFLOAT32(inBufP + 124);
            swpP->interpulse_per2 = GetVals_GetFLOAT32(inBufP + 128);
            swpP->interpulse_per3 = GetVals_GetFLOAT32(inBufP + 132);
            swpP->interpulse_per4 = GetVals_GetFLOAT32(inBufP + 136);
            swpP->interpulse_per5 = GetVals_GetFLOAT32(inBufP + 140);
        } else if (strncmp(desc_id, "PARM", 4) == 0) {
            sz = (swpP->n_parms + 1) * sizeof(struct Dorade_ParmDesc);
            swpP->parmP
                = (struct Dorade_ParmDesc *)REALLOC(swpP->parmP, sz);
            parmP = swpP->parmP + swpP->n_parms;
            Dorade_InitParm(parmP);
            swpP->n_parms++;
            parmP->use = 1;
            swpP->n_use_parms++;
            strncpy(parmP->name, inBufP + 8, 8);
            strncpy(parmP->description, inBufP + 16, 40);
            strncpy(parmP->units, inBufP + 56, 8);
            parmP->interpulse_time = GetVals_GetSINT16(inBufP + 64);
            parmP->xmitted_freq = GetVals_GetSINT16(inBufP + 66);
            parmP->recvr_bandwidth = GetVals_GetFLOAT32(inBufP + 68);
            parmP->pulse_width = GetVals_GetSINT16(inBufP + 72);
            parmP->polarization = GetVals_GetSINT16(inBufP + 74);
            parmP->num_samples = GetVals_GetSINT16(inBufP + 76);
            parmP->binary_format = GetVals_GetSINT16(inBufP + 78);
            strncpy(parmP->threshold_field, inBufP + 80, 8);
            parmP->threshold_value = GetVals_GetFLOAT32(inBufP + 88);
            parmP->scale = GetVals_GetFLOAT32(inBufP + 92);
            parmP->bias = GetVals_GetFLOAT32(inBufP + 96);
            parmP->bad_data = GetVals_GetSINT32(inBufP + 100);
        } else if (strncmp(desc_id, "CELV", 4) == 0) {
            swpP->n_cells = GetVals_GetUINT32(inBufP + 8);
            swpP->distP = (float *)MALLOC(swpP->n_cells * sizeof(float));
            for (fInP = (float *)(inBufP + 12),
                    fP = swpP->distP,
                    f1P = fP + swpP->n_cells;
                    fP < f1P; fInP++, fP++) {
                *fP = *fInP;
            }
        } else if (strncmp(desc_id, "CFAC", 4) == 0) {
            swpP->azimuth_corr
                = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 8));
            swpP->elevation_corr
                = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 12));
            swpP->range_delay_corr = GetVals_GetFLOAT32(inBufP + 16);
            swpP->longitude_corr
                = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 20));
            swpP->latitude_corr
                = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 24));
            swpP->pressure_alt_corr = GetVals_GetFLOAT32(inBufP + 28);
            swpP->radar_alt_corr = GetVals_GetFLOAT32(inBufP + 32);
            swpP->ew_gndspd_corr = GetVals_GetFLOAT32(inBufP + 36);
            swpP->ns_gndspd_corr = GetVals_GetFLOAT32(inBufP + 40);
            swpP->vert_vel_corr = GetVals_GetFLOAT32(inBufP + 44);
            swpP->heading_corr = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 48));
            swpP->roll_corr = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 52));
            swpP->pitch_corr = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 56));
            swpP->drift_corr = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 60));
            swpP->rot_angle_corr = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 64));
            swpP->tilt_corr = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 68));
        } else if (strncmp(desc_id, "SWIB", 4) == 0) {
            struct Dorade_RayHdr *rhp, *rhp1;

            strncpy(swpP->swib_comment, inBufP + 8, 8);
            swpP->swib_comment[8] = '\0';
            swpP->sweep_num = GetVals_GetSINT32(inBufP + 16);
            swpP->n_rays = GetVals_GetSINT32(inBufP + 20);
            swpP->n_good_rays = swpP->n_rays;
            sz = swpP->n_rays * sizeof(struct Dorade_RayHdr);
            swpP->rayHdrP = (struct Dorade_RayHdr *)MALLOC(sz);
            for (rhp = swpP->rayHdrP, rhp1 = rhp + swpP->n_rays;
                    rhp < rhp1;
                    rhp++) {
                Dorade_InitRayHdr(rhp);
            }
            swpP->start_angle
                = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 24));
            swpP->stop_angle
                = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 28));
            swpP->fixed_angle
                = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 32));
            swpP->filter_flag = GetVals_GetSINT32(inBufP + 36);
            break;
        }
    }
    Dorade_HashParms(swpP);
    return 1;

error:
    ErrMsg_Append("Dorade_ReadSweep failed\n");
    Dorade_FreeSweep(swpP);
    return 0;
}

/*
 *------------------------------------------------------------------------
 *
 * Dorade_ReadSweep --
 *
 *         This function reads and stores the contents of a DORADE sweep file.
 *
 * Arguments:
 *        FILE *in                        - Input file or stream. It should be
 *                                          positioned at the start of a DORADE
 *                                          sweep file.
 *        struct Dorade_Sweep *swpP        - Place to store values, should have
 *                                          been initialized with a call to
 *                                          Dorade_InitSweep.
 *
 * Results:
 *         Return value is non-zero if successful, zero if something goes wrong.
 *
 * Side effects:
 *         The sweep structure is updated.  Memory might be allocated, which
 *         should be eventually freed with a call to Dorade_FreeSweep.
 *         The input stream remains open, although it will be repositioned to
 *         the end of the sweep.
 *
 *------------------------------------------------------------------------
 */

int Dorade_ReadSweep(FILE *in, struct Dorade_Sweep *swpP)
{
    char desc_id[4];                        /* Descriptor identifier */
    size_t desc_len;                        /* Descriptor length */
    size_t front = 8;                        /* Size of id + size specifier */
    struct GeoTime_Jul jul0;                /* Julian day of start of year in
                                         * which this sweep was taken.
                                         * This is needed for ray times,
                                         * which use "Julian day" relative
                                         * to start of year. */

    /*
     * Convenience variables
     */

    struct GeoTime_Cal cal;
    struct Dorade_RayHdr *rayHdrP;
    long n_parms;
    unsigned n_cells;
    long n_rays;
    int r, p;
    struct Dorade_ParmDesc *parmP;
    double scale_inv, bias;
    long julian_day;
    short hour, minute;
    float second, millisecond;

    /*
     * Loop parameters
     */

    char *cp;
    short *sp, *sp1;
    float *fp, *fp1, *fp2;
    int cnt;

    swpP->n_parms = 0;
    r = 0;

    if (front > in_sz) {
        inBufP = (char *)REALLOC(inBufP, front);
        in_sz = front;
    }
    in_sz = front;
    if ( !Dorade_ReadSweepHeaders(in, swpP) ) {
        ErrMsg_Append("Allocation failed.\n");
        goto error;
    }
    cal = GeoTime_JulToCal(swpP->time);
    jul0 = GeoTime_CalToJul(GeoTime_CalSet(cal.year, 1, 1, 0, 0, 0.0));

    /*
     * Allocate data array.
     */

    n_parms = swpP->n_parms;
    n_cells = swpP->n_cells;
    n_rays = swpP->n_rays;
    if (n_parms == LONG_MIN || n_cells == UINT_MAX || n_rays == LONG_MIN) {
        ErrMsg_Append("Ray data found before sweep "
                "dimensions known.\n");
        goto error;
    }
    if (n_parms == 0 || n_cells == 0 || n_rays == 0) {
        ErrMsg_Append("Ray data found for sweep "
                "of zero size.\n");
        goto error;
    }

    if ( !(swpP->dat = Dorade_AllocDat(n_parms, n_rays, n_cells)) ) {
        goto error;
    }

    /*
     * Read in rays.
     */

    for (r = 0; r < swpP->n_rays; r++) {
        double incr;

        if (fread(inBufP, 1, front, in) != front) {
            ErrMsg_Append("Failed to read id and length for ray "
                    "info descriptor\n");
            goto error;
        }
        strncpy(desc_id, inBufP, 4);
        desc_len = GetVals_GetSINT32(inBufP + 4);
        if (desc_len < 0) {
            ErrMsg_Append("Negative descriptor length.\n");
            goto error;
        }
        if (desc_len > 0) {
            if (desc_len > in_sz) {
                inBufP = (char *)REALLOC(inBufP, desc_len);
                in_sz = desc_len;
            }
            desc_len -= front;
            if (fread(inBufP + front, 1, desc_len, in) != desc_len) {
                ErrMsg_Append("Failed to read ray info descriptor\n");
                goto error;
            }
        } else {
            ErrMsg_Append("Unexpected empty block\n");
            goto error;
        }
        if (strncmp(desc_id, "RYIB", 4) != 0) {
            ErrMsg_Append("Missing ray info block\n");
            goto error;
        }
        rayHdrP = swpP->rayHdrP + r;
        rayHdrP->time = jul0;
        julian_day = GetVals_GetSINT32(inBufP + 12);
        hour = GetVals_GetSINT16(inBufP + 16);
        minute = GetVals_GetSINT16(inBufP + 18);
        second = GetVals_GetSINT16(inBufP + 20);
        millisecond = GetVals_GetSINT16(inBufP + 22);
        incr = 86400.0 * (julian_day - 1) + 3600.0 * hour
            + 60.0 * minute + second + 0.001 * millisecond;
        GeoTime_Incr(&rayHdrP->time, incr);
        rayHdrP->azimuth = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 24));
        rayHdrP->elevation = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 28));
        rayHdrP->peak_power = GetVals_GetFLOAT32(inBufP + 32);
        rayHdrP->true_scan_rate
            = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 36));

        /*
         * Get the platform info, which should come immediately
         * after the ray header.
         */

        if (fread(inBufP, 1, front, in) != front) {
            ErrMsg_Append("Failed to read id and length for "
                    "platform info descriptor\n");
            goto error;
        }
        strncpy(desc_id, inBufP, 4);
        desc_len = GetVals_GetSINT32(inBufP + 4);
        if (desc_len < 0) {
            ErrMsg_Append("Negative descriptor length.\n");
            goto error;
        }
        if (desc_len > 0) {
            if (desc_len > in_sz) {
                inBufP = (char *)REALLOC(inBufP, desc_len);
                in_sz = desc_len;
            }
            desc_len -= front;
            if (fread(inBufP + front, 1, desc_len, in) != desc_len) {
                ErrMsg_Append("Failed to read platform info descriptor\n");
                goto error;
            }
        } else {
            ErrMsg_Append("Unexpected empty block\n");
            goto error;
        }
        rayHdrP->ray_status = GetVals_GetSINT32(inBufP + 40);
        if (strncmp(desc_id, "ASIB", 4) != 0) {
            ErrMsg_Append("Missing platform info block\n");
            goto error;
        }
        rayHdrP = swpP->rayHdrP + r;
        rayHdrP->longitude = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 8));
        rayHdrP->latitude = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 12));
        rayHdrP->altitude_msl = GetVals_GetFLOAT32(inBufP + 16);
        rayHdrP->altitude_agl = GetVals_GetFLOAT32(inBufP + 20);
        rayHdrP->ew_velocity = GetVals_GetFLOAT32(inBufP + 24);
        rayHdrP->ns_velocity = GetVals_GetFLOAT32(inBufP + 28);
        rayHdrP->vert_velocity = GetVals_GetFLOAT32(inBufP + 32);
        rayHdrP->heading = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 36));
        rayHdrP->roll = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 40));
        rayHdrP->pitch = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 44));
        rayHdrP->drift_angle = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 48));
        rayHdrP->rotation_angle
            = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 52));
        rayHdrP->tilt = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 56));
        rayHdrP->ew_horiz_wind = GetVals_GetFLOAT32(inBufP + 60);
        rayHdrP->ns_horiz_wind = GetVals_GetFLOAT32(inBufP + 64);
        rayHdrP->vert_wind = GetVals_GetFLOAT32(inBufP + 68);
        rayHdrP->heading_change
            = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 72));
        rayHdrP->pitch_change = AngleFmDeg(GetVals_GetFLOAT32(inBufP + 76));
        rayHdrP->good = 1;

        for (p = 0; p < n_parms; p++) {
            if (fread(inBufP, 1, front, in) != front) {
                ErrMsg_Append("Failed to read id and length "
                        "for data descriptor\n");
                goto error;
            }
            strncpy(desc_id, inBufP, 4);
            desc_len = GetVals_GetSINT32(inBufP + 4);
            if (desc_len < 0) {
                ErrMsg_Append("Negative descriptor length.\n");
                goto error;
            }
            if (desc_len > 0) {
                if (desc_len > in_sz) {
                    inBufP = (char *)REALLOC(inBufP, desc_len);
                    in_sz = desc_len;
                }
                desc_len -= front;
                if (fread(inBufP + front, 1, desc_len, in)
                        != desc_len) {
                    ErrMsg_Append("Failed to read data descriptor\n");
                    goto error;
                }
            } else {
                ErrMsg_Append("Unexpected empty block\n");
                goto error;
            }
            if (strncmp(desc_id, "XSTF", 4) == 0) {
                p--;
                continue;
            }
            if (strncmp(desc_id, "RDAT", 4) != 0) {
                ErrMsg_Append("Missing data block.\n");
                goto error;
            }
            parmP = swpP->parmP + p;
            if (strncmp(inBufP + front, parmP->name, 8) != 0) {
                ErrMsg_Append("Parameters out of order in ray.\n");
                goto error;
            }
            scale_inv = 1.0 / parmP->scale;
            bias = parmP->bias;
            switch (parmP->binary_format) {
                case DD_8_BITS:
                    for (cp = (char *)(inBufP + 16),
                            fp  = swpP->dat[p][r],
                            fp1 = swpP->dat[p][r] + swpP->n_cells;
                            fp < fp1; cp++, fp++) {
                        if (*cp == parmP->bad_data) {
                            *fp = Radar_NoData();
                        } else {
                            *fp = *cp * scale_inv - bias;
                            if (*fp < parmP->min) {
                                parmP->min = *fp;
                            }
                        }
                    }
                    break;
                case DD_16_BITS:
                    if (swpP->compression) {
                        sp = (short *)(inBufP + 16);
                        fp  = swpP->dat[p][r];
                        fp1 = fp + swpP->n_cells;
                        while (*sp != 1) {
                            cnt = *sp & 0x7fff;
                            if (fp + cnt > fp1) {
                                ErrMsg_Append("Pointer went "
                                        "out of data array while "
                                        "decompressing ray.  ");
                                goto error;
                            }
                            if (*sp & 0x8000) {
                                /*
                                 * Run of good data.  Transfer cnt
                                 * values from input buffer to
                                 * data array in sweep.
                                 */

                                sp++;
                                for (sp1 = sp + cnt; sp < sp1; sp++, fp++) {
                                    if (*sp == parmP->bad_data) {
                                        *fp = Radar_NoData();
                                    } else {
                                        *fp = *sp * scale_inv - bias;
                                        if (*fp < parmP->min) {
                                            parmP->min = *fp;
                                        }
                                    }
                                }
                            } else {
                                /*
                                 * Run of no data or bad data.  Put
                                 * cnt no-data values into data array
                                 * in sweep.  Ignore input buffer.
                                 */

                                sp++;
                                for (fp2 = fp + cnt; fp < fp2; fp++) {
                                    *fp = Radar_NoData();
                                }
                            }
                        }
                        if (fp != fp1) {
                            ErrMsg_Append("Decompression "
                                    "finished before end of ray.  ");
                            goto error;
                        }
                    } else {
                        for (sp = (short *)(inBufP + 16),
                                fp = swpP->dat[p][r],
                                fp1 = fp + swpP->n_cells;
                                fp < fp1;
                                sp++, fp++) {
                            if (*sp == parmP->bad_data) {
                                *fp = Radar_NoData();
                            } else {
                                *fp = *sp * scale_inv - bias;
                                if (*fp < parmP->min) {
                                    parmP->min = *fp;
                                }
                            }
                        }
                    }
                    break;
                case DD_24_BITS:
                    ErrMsg_Append("No support for 24 bit "
                            "integers.\n");
                    goto error;
                case DD_32_BIT_FP:
                    ErrMsg_Append("No support for 32 bit "
                            "floats.\n");
                    goto error;
                case DD_16_BIT_FP:
                    ErrMsg_Append("No support for 16 bit "
                            "floats.\n");
                    goto error;
            }
        }
    }
    return 1;

error:
    ErrMsg_Append("Dorade_ReadSweep failed\n");
    Dorade_FreeSweep(swpP);
    return 0;
}

/*
 *------------------------------------------------------------------------
 *
 * Dorade_RmParm --
 *
 *         This function removes a parameter from a sweep.
 *
 * Arguments:
 *         swpP        - DORADE sweep
 *         parm        - name of parameter to remove.
 *
 * Results:
 *         A parameter's use flag is set to false.
 *         Parameter allocations remain, but their contents are ignored.
 *         This function quietly does nothing if the parameter does not exist
 *         or has already been removed.
 *
 *------------------------------------------------------------------------
 */

void Dorade_RmParm(struct Dorade_Sweep *swpP, char *parm)
{
    struct Dorade_ParmDesc *parmP;

    if ( !swpP ) {
        ErrMsg_Append("Parameter copy function given "
                "non-existent DORADE structure.\n");
        return;
    }
    if ( (parmP = Dorade_GetParm(swpP, parm) ) ) {
        parmP->use = 0;
        --swpP->n_use_parms;
    }
}

/*
 *------------------------------------------------------------------------
 *
 * Dorade_SweepFName --
 *
 *         This function generates the name of a sweep file for a DORADE structure.
 *
 * Arguments:
 *         char *nm                        - if not NULL, a memory address at
 *                                           which to write the name.
 *        size_t *l;                        - allocation at nm, if nm is not NULL.
 *        struct Dorade_Sweep *swpP        - a DORADE sweep
 *
 * Results:
 *         Return value is a character string with the name of a sweep file.
 *
 * Side effects:
 *         If nm is not NULL, it is assumed to point to allocated memory.
 *         Previous contents will be clobbered by a string.  Upon return
 *         the allocation might have been modified by REALLOC, and l will
 *         contain the size of the final allocation.  Return value is nm, or
 *         the address of the new allocation.
 *
 *         If nm is NULL, memory is allocated with MALLOC.  The address of the
 *         new memory is returned, and the size of the allocation is stored into l.
 *         The returned memory should eventually be freed with a call to FREE.
 *
 *------------------------------------------------------------------------
 */

char * Dorade_SweepFName(char *nm, size_t *l, struct Dorade_Sweep *swp)
{
    struct GeoTime_Cal cal;
    char *scan_mode;
    size_t dflt_sz;
    char *solodir, *p;
    struct stat st;

    //dflt_sz = strlen("swp.YYYMMDDhhmmss.RADARXXX.10.0.4_PPI_v1_ref") + 1;
    //dflt_sz = strlen("swp.YYYMMDDhhmmss.RADARXXX.10.0.4_PPI_v10_ref") + 1;
    dflt_sz = strlen("RRRRSOLOQC/DBZ/swp.YYYMMDDhhmmss.RADARXXX.10.0.4_PPI_v10_ref") + 1;
    if ( !nm ) {
        nm = MALLOC(dflt_sz);
    } else if (nm && *l < dflt_sz) {
        nm = REALLOC(nm, dflt_sz);
    }
    *l = dflt_sz;
    switch (swp->scan_mode) {
        case DORADE_CALIBRATION:
            scan_mode = "CAL";
            break;
        case DORADE_PPI:
            scan_mode = "PPI";
            break;
        case DORADE_COPLANE:
            scan_mode = "CPL";
            break;
        case DORADE_RHI:
            scan_mode = "RHI";
            break;
        case DORADE_VERTICAL_POINTING:
            scan_mode = "VPT";
            break;
        case DORADE_TARGET_MANUAL:
            scan_mode = "TGM";
            break;
    }
    cal = GeoTime_JulToCal(swp->time);

    solodir = MALLOC(16);

    strcpy(solodir,swp->radar_name);
    p = strchr(swp->proj_name,',');
    strncat(solodir,swp->proj_name,(p-swp->proj_name));
    if(stat(solodir,&st) != 0) mkdir(solodir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    strcat(solodir,"/");

    strncat(solodir,swp->parmP->name,3);
    if(stat(solodir,&st) != 0) mkdir(solodir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    char *tstr =  strhead(swp->parmP->name,3);

    sprintf(nm,"%s/swp.%0d%02d%02d%02d%02d%02d.%.8s.0.%01.1f_%s_v%d_%s",
        solodir,
        cal.year - 1900, cal.month, cal.day, cal.hour, cal.minute,
        (int)(cal.second + 0.5), swp->radar_name, AngleToDeg(swp->fixed_angle),
        scan_mode,swp->sweep_num,tstr);

    free(tstr);
    free(solodir);

    return nm;
}

/*
 *------------------------------------------------------------------------
 *
 * Dorade_WriteSweep --
 *
 *         This function writes the contents of a sweep structure to a DORADE
 *         sweep file.
 *
 * Arguments:
 *        FILE *out                        - Output file or stream
 *        struct Dorade_Sweep *swpP        - Sweep to send to output
 *
 * Results:
 *         Return value is non-zero if successful, zero if something goes wrong.
 *
 * Side effects:
 *         Output is generated.
 *
 *------------------------------------------------------------------------
 */

int Dorade_WriteSweep(FILE *out, struct Dorade_Sweep *swpP)
{
    char *outBufP;                        /* Output buffer */
    struct GeoTime_Jul jul0;                /* Julian day of start of year in
                                         * which this sweep was taken.
                                         * This is needed for ray times,
                                         * which use "Julian day" relative
                                         * to start of year. */
    struct GeoTime_Jul sweep_time;
    struct GeoTime_Jul time0;                /* Start of Unix epoch */
    time_t time;
    size_t file_sz;                        /* File size for SSWB block */

    /*
     * Output block sizes
     */

    size_t sswb_sz = 196;
    size_t rot_ang_table_sz = 32;
    size_t vol_descr_sz = 72;
    size_t rad_descr_sz = 144;
    size_t parm_descr_sz = 104;
    size_t cfac_sz = 72;
    size_t swib_sz = 40;
    size_t ryib_sz = 44;
    size_t platfm_info_blk_sz = 80;
    size_t null_sz = 8;
    size_t celv_sz;
    size_t rdat_sz;
    size_t rktb_sz;
    size_t max_sz;

    /*
     * Values for RKTB block.
     */

    struct rotation_table_entry {
        float angle;
        int offset;
        int size;
    };

    /*
     * Convenience variables.
     */

    long n_parms, n_use_parms;
    long n_rays;
    long n_good_rays;
    unsigned n_cells;
    struct Dorade_ParmDesc *parmP, *pe;
    double scale, bias;
    struct GeoTime_Cal cal;

    /*
     * Loop parameters
     */

    int r, p;                                /* Ray and parameter indeces*/
    char *cp;
    short *sp;
    float *fOutP;
    float *fp, *fp1;

    /*
     * Variables from dd_rotang_table function
     * in soloii:translate/ddb_common.c
     */

    int angle_ndx_size = 480;
    float angle2ndx;
    int ndx_que_size;
    int first_key_offset;
    int angle_table_offset;
    int *angle_table;
    struct rotation_table_entry *rotation_table, entry;
    int i;

    if ( !out ) {
        ErrMsg_Append("No output file.\n");
        return 0;
    }
    if ( !swpP ) {
        ErrMsg_Append("Bogus sweep.\n");
        return 0;
    }
    n_parms = swpP->n_parms;
    n_use_parms = swpP->n_use_parms;
    n_rays = swpP->n_rays;
    n_good_rays = swpP->n_good_rays;
    n_cells = swpP->n_cells;
    cal = GeoTime_JulToCal(swpP->time);
    jul0 = GeoTime_CalToJul(GeoTime_CalSet(cal.year, 1, 1, 0, 0, 0.0));
    file_sz = 0;

    /*
     * Use scale and bias seen in sample sweep files.  Ignore value in
     * parameter block.
     */

    scale = 100.0;
    bias = 0.0;

    /*
     * Figure out the sizes of the output blocks and the size of the
     * largest block.  Make sure block sizes are multiples of 4.
     */

    max_sz = vol_descr_sz;
    max_sz = (rad_descr_sz > max_sz) ? rad_descr_sz : max_sz;
    max_sz = (n_use_parms * parm_descr_sz > max_sz)
        ? n_use_parms * parm_descr_sz : max_sz;
    max_sz = (cfac_sz > max_sz) ? cfac_sz : max_sz;
    max_sz = (swib_sz > max_sz) ? swib_sz : max_sz;
    max_sz = (ryib_sz > max_sz) ? ryib_sz : max_sz;
    max_sz = (platfm_info_blk_sz > max_sz) ? platfm_info_blk_sz : max_sz;
    celv_sz = 16 + 4 * n_cells;
    max_sz = (celv_sz > max_sz) ? celv_sz : max_sz;
    rdat_sz = 16 + ((2 * n_cells + 3) / 4) * 4;
    max_sz = (rdat_sz > max_sz) ? rdat_sz : max_sz;
    rktb_sz = rot_ang_table_sz
        + ((angle_ndx_size > n_good_rays ? angle_ndx_size : n_good_rays))
        * (4 + sizeof(struct rotation_table_entry));
    max_sz = (rktb_sz > max_sz) ? rktb_sz : max_sz;
    outBufP = (char *)MALLOC(max_sz);

    /*
     * Send the "SSWB" block
     */

    memset(outBufP, 0, sswb_sz);
    strncpy(outBufP, "SSWB", 4);
    GetVals_PutSINT32(sswb_sz, outBufP + 4);
    time0 = GeoTime_CalToJul(GeoTime_CalSet(1970, 1, 1, 0, 0, 0.0));
    sweep_time = swpP->time;
    time = (int)((double)(sweep_time.day - time0.day) * 86400.0
            + sweep_time.second - time0.second + 0.5);
    GetVals_PutSINT32(time, outBufP + 12);
    GetVals_PutSINT32(1, outBufP + 64);
    GetVals_PutSINT32(sswb_sz + vol_descr_sz + rad_descr_sz
            + n_use_parms * parm_descr_sz + celv_sz + cfac_sz + swib_sz
            + n_good_rays * (ryib_sz + platfm_info_blk_sz
                + n_use_parms * rdat_sz) + null_sz, outBufP + 100);
    GetVals_PutSINT32(rktb_sz, outBufP + 104);
    GetVals_PutSINT32(2, outBufP + 108);
    if (fwrite(outBufP, 1, sswb_sz, out) != sswb_sz) {
        ErrMsg_Append("Could not write outBufP descriptor.\n");
        goto error;
    }
    file_sz += sswb_sz;

    /*
     * Send the volume descriptor.
     */

    memset(outBufP, 0, vol_descr_sz);
    strncpy(outBufP, "VOLD", 4);
    GetVals_PutSINT32(vol_descr_sz, outBufP +4);
    GetVals_PutSINT16(1, outBufP + 10);
    GetVals_PutSINT32(65500, outBufP + 12);
    strncpy(outBufP + 16, swpP->proj_name, 20);
    cal = GeoTime_JulToCal(swpP->time);
    GetVals_PutSINT16(cal.year, outBufP + 36);
    GetVals_PutSINT16(cal.month, outBufP + 38);
    GetVals_PutSINT16(cal.day, outBufP + 40);
    GetVals_PutSINT16(cal.hour, outBufP + 42);
    GetVals_PutSINT16(cal.minute, outBufP + 44);
    GetVals_PutSINT16(cal.second, outBufP + 46);
    strncpy(outBufP + 56, swpP->radar_name, 8);
    GetVals_PutSINT16(swpP->n_sensors, outBufP + 70);
    if (fwrite(outBufP, 1, vol_descr_sz, out) != vol_descr_sz) {
        ErrMsg_Append("Could not write volume descriptor.\n");
        goto error;
    }
    file_sz += vol_descr_sz;

    /*
     * Fill in the radar descriptor.
     */

    memset(outBufP, 0, rad_descr_sz);
    strncpy(outBufP, "RADD", 4);
    GetVals_PutSINT32(rad_descr_sz, outBufP + 4);
    strncpy(outBufP + 8, swpP->radar_name, 8);
    GetVals_PutFLOAT32(swpP->radar_const, outBufP + 16);
    GetVals_PutFLOAT32(swpP->peak_power, outBufP + 16);
    GetVals_PutFLOAT32(swpP->noise_power, outBufP + 16);
    GetVals_PutFLOAT32(swpP->receiver_gain, outBufP + 16);
    GetVals_PutFLOAT32(swpP->antenna_gain, outBufP + 16);
    GetVals_PutFLOAT32(swpP->system_gain, outBufP + 16);
    GetVals_PutFLOAT32(AngleToDeg(swpP->horz_beam_width), outBufP + 40);
    GetVals_PutFLOAT32(AngleToDeg(swpP->vert_beam_width), outBufP + 44);
    GetVals_PutSINT16(swpP->radar_type, outBufP + 48);
    GetVals_PutSINT16(swpP->scan_mode, outBufP + 50);
    GetVals_PutSINT16(swpP->n_use_parms, outBufP + 64);
    GetVals_PutSINT16(0, outBufP + 68);
    GetVals_PutSINT16(swpP->data_reduction, outBufP + 70);
    GetVals_PutFLOAT32(swpP->data_red_parm0, outBufP + 72);
    GetVals_PutFLOAT32(swpP->data_red_parm1, outBufP + 76);
    GetVals_PutFLOAT32(AngleToDeg(swpP->radar_location.lon),
            outBufP + 80);
    GetVals_PutFLOAT32(AngleToDeg(swpP->radar_location.lat),
            outBufP + 84);
    GetVals_PutFLOAT32(swpP->radar_altitude, outBufP + 88);
    GetVals_PutFLOAT32(swpP->eff_unamb_vel, outBufP + 92);
    GetVals_PutFLOAT32(swpP->eff_unamb_range, outBufP + 96);
    GetVals_PutSINT16(swpP->num_freq_trans, outBufP + 100);
    GetVals_PutSINT16(swpP->num_ipps_trans, outBufP + 102);
    GetVals_PutFLOAT32(swpP->freq1, outBufP + 104);
    GetVals_PutFLOAT32(swpP->freq2, outBufP + 108);
    GetVals_PutFLOAT32(swpP->freq3, outBufP + 112);
    GetVals_PutFLOAT32(swpP->freq4, outBufP + 116);
    GetVals_PutFLOAT32(swpP->freq5, outBufP + 120);
    GetVals_PutFLOAT32(swpP->interpulse_per1, outBufP + 124);
    GetVals_PutFLOAT32(swpP->interpulse_per2, outBufP + 128);
    GetVals_PutFLOAT32(swpP->interpulse_per3, outBufP + 132);
    GetVals_PutFLOAT32(swpP->interpulse_per4, outBufP + 136);
    GetVals_PutFLOAT32(swpP->interpulse_per5, outBufP + 140);
    if (fwrite(outBufP, 1, rad_descr_sz, out) != rad_descr_sz) {
        ErrMsg_Append("Could not write radar descriptor.\n");
        goto error;
    }
    file_sz += rad_descr_sz;

    /*
     * Send the parameter descriptors.
     */

    for (parmP = swpP->parmP, pe = parmP + swpP->n_parms;
            parmP < pe;
            parmP++) {
        if ( !parmP->use ) {
            continue;
        }
        memset(outBufP, 0, parm_descr_sz);
        strncpy(outBufP, "PARM", 4);
        GetVals_PutSINT32(parm_descr_sz, outBufP + 4);
        strncpy(outBufP + 8, parmP->name, 8);
        strncpy(outBufP + 16, parmP->description, 40);
        strncpy(outBufP + 56, parmP->units, 56);
        GetVals_PutSINT16(parmP->interpulse_time, outBufP + 64);
        GetVals_PutSINT16(parmP->xmitted_freq, outBufP + 66);
        GetVals_PutFLOAT32(parmP->recvr_bandwidth, outBufP + 68);
        GetVals_PutSINT16(parmP->pulse_width, outBufP + 72);
        GetVals_PutSINT16(parmP->polarization, outBufP + 74);
        GetVals_PutSINT16(parmP->num_samples, outBufP + 76);
        GetVals_PutSINT16(parmP->binary_format, outBufP + 78);
        strncpy(outBufP + 80, parmP->threshold_field, 8);
        GetVals_PutFLOAT32(parmP->threshold_value, outBufP + 88);
        GetVals_PutFLOAT32(scale, outBufP + 92);
        GetVals_PutFLOAT32(bias, outBufP + 96);
        GetVals_PutSINT32(parmP->bad_data, outBufP + 100);
        if (fwrite(outBufP, 1, parm_descr_sz, out) != parm_descr_sz) {
            ErrMsg_Append("Could not write parameter descriptor.\n");
            goto error;
        }
        file_sz += parm_descr_sz;
    }

    /*
     * Send the cell range vector.
     */

    memset(outBufP, 0, celv_sz);
    strncpy(outBufP, "CELV", 4);
    GetVals_PutSINT32(celv_sz, outBufP + 4);
    GetVals_PutSINT32(swpP->n_cells, outBufP + 8);
    for (fp = swpP->distP, fp1 = fp + swpP->n_cells,
            fOutP = (float *)(outBufP + 12);
            fp < fp1; fp++, fOutP++) {
        *fOutP = *fp;
    }
    if (fwrite(outBufP, 1, celv_sz, out) != celv_sz) {
        ErrMsg_Append("Could not write cell range vector.\n");
        goto error;
    }
    file_sz += celv_sz;

    /*
     * Send the correction factor descriptor.
     */

    memset(outBufP, 0, cfac_sz);
    strncpy(outBufP, "CFAC", 4);
    GetVals_PutSINT32(cfac_sz, outBufP + 4);
    GetVals_PutFLOAT32(swpP->azimuth_corr, outBufP + 8);
    GetVals_PutFLOAT32(swpP->elevation_corr, outBufP + 12);
    GetVals_PutFLOAT32(swpP->range_delay_corr, outBufP + 16);
    GetVals_PutFLOAT32(swpP->longitude_corr, outBufP + 20);
    GetVals_PutFLOAT32(swpP->latitude_corr, outBufP + 24);
    GetVals_PutFLOAT32(swpP->pressure_alt_corr, outBufP + 28);
    GetVals_PutFLOAT32(swpP->radar_alt_corr, outBufP + 32);
    GetVals_PutFLOAT32(swpP->ew_gndspd_corr, outBufP + 36);
    GetVals_PutFLOAT32(swpP->ns_gndspd_corr, outBufP + 40);
    GetVals_PutFLOAT32(swpP->vert_vel_corr, outBufP + 44);
    GetVals_PutFLOAT32(swpP->heading_corr, outBufP + 48);
    GetVals_PutFLOAT32(swpP->roll_corr, outBufP + 52);
    GetVals_PutFLOAT32(swpP->pitch_corr, outBufP + 56);
    GetVals_PutFLOAT32(swpP->drift_corr, outBufP + 60);
    GetVals_PutFLOAT32(swpP->rot_angle_corr, outBufP + 64);
    GetVals_PutFLOAT32(swpP->tilt_corr, outBufP + 68);
    if (fwrite(outBufP, 1, cfac_sz, out) != cfac_sz) {
        ErrMsg_Append("Could not write correction factor descriptor.\n");
        goto error;
    }
    file_sz += cfac_sz;

    /*
     * Send the sweep info block
     */

    memset(outBufP, 0, swib_sz);
    strncpy(outBufP, "SWIB", 4);
    GetVals_PutSINT32(swib_sz, outBufP + 4);
    GetVals_PutSINT32(swpP->sweep_num, outBufP + 16);
    GetVals_PutSINT32(swpP->n_good_rays, outBufP + 20);
    GetVals_PutFLOAT32(AngleToDeg(swpP->start_angle), outBufP + 24);
    GetVals_PutFLOAT32(AngleToDeg(swpP->stop_angle), outBufP + 28);
    GetVals_PutFLOAT32(AngleToDeg(swpP->fixed_angle), outBufP + 32);
    GetVals_PutSINT32(swpP->filter_flag, outBufP + 36);
    if (fwrite(outBufP, 1, swib_sz, out) != swib_sz) {
        ErrMsg_Append("Could not write sweep info block.\n");
        goto error;
    }
    file_sz += swib_sz;

    /*
     * Loop through the rays.
     */

    for (r = 0; r < n_rays; r++) {
        struct Dorade_RayHdr *rayHdrP = swpP->rayHdrP + r;

        if ( !rayHdrP->good ) {
            continue;
        }

        /*
         * Send the ray info block.
         */

        memset(outBufP, 0, ryib_sz);
        strncpy(outBufP, "RYIB", 4);
        GetVals_PutSINT32(ryib_sz, outBufP +4);
        GetVals_PutSINT32(swpP->sweep_num, outBufP + 8);
        GetVals_PutSINT32(rayHdrP->time.day - jul0.day + 1, outBufP + 12);
        cal = GeoTime_JulToCal(rayHdrP->time);
        GetVals_PutSINT16(cal.hour, outBufP + 16);
        GetVals_PutSINT16(cal.minute, outBufP + 18);
        GetVals_PutSINT16(cal.second, outBufP + 20);
        GetVals_PutSINT16((cal.second - (short)cal.second) * 1000,
                outBufP + 22);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->azimuth), outBufP + 24);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->elevation), outBufP + 28);
        GetVals_PutFLOAT32(rayHdrP->peak_power, outBufP + 32);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->true_scan_rate),
                outBufP + 36);
        GetVals_PutSINT32(rayHdrP->ray_status, outBufP + 40);
        if (fwrite(outBufP, 1, ryib_sz, out) != ryib_sz) {
            ErrMsg_Append("Could not write ray info block.\n");
            goto error;
        }
        file_sz += ryib_sz;

        /*
         * Send the platform info block.
         */

        memset(outBufP, 0, platfm_info_blk_sz);
        strncpy(outBufP, "ASIB", 4);
        GetVals_PutSINT32(platfm_info_blk_sz, outBufP +4);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->longitude), outBufP + 8);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->latitude), outBufP + 12);
        GetVals_PutFLOAT32(rayHdrP->altitude_msl, outBufP + 16);
        GetVals_PutFLOAT32(rayHdrP->altitude_agl, outBufP + 20);
        GetVals_PutFLOAT32(rayHdrP->ew_velocity, outBufP + 24);
        GetVals_PutFLOAT32(rayHdrP->ns_velocity, outBufP + 28);
        GetVals_PutFLOAT32(rayHdrP->vert_velocity, outBufP + 32);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->heading), outBufP + 36);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->roll), outBufP + 40);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->pitch), outBufP + 44);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->drift_angle), outBufP + 48);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->rotation_angle),
                outBufP + 52);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->tilt), outBufP + 56);
        GetVals_PutFLOAT32(rayHdrP->ew_horiz_wind, outBufP + 60);
        GetVals_PutFLOAT32(rayHdrP->ns_horiz_wind, outBufP + 64);
        GetVals_PutFLOAT32(rayHdrP->vert_wind, outBufP + 68);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->heading_change),
                outBufP + 72);
        GetVals_PutFLOAT32(AngleToDeg(rayHdrP->pitch_change), outBufP + 76);
        if (fwrite(outBufP, 1, platfm_info_blk_sz, out)
                != platfm_info_blk_sz) {
            ErrMsg_Append("Could not write platform info block.\n");
            goto error;
        }
        file_sz += platfm_info_blk_sz;

        /*
         * Send data blocks.
         */

        for (p = 0; p < n_parms; p++) {
            parmP = swpP->parmP + p;

            if ( !parmP->use ) {
                continue;
            }
            memset(outBufP, 0, rdat_sz);
            strncpy(outBufP, "RDAT", 4);
            GetVals_PutSINT32(rdat_sz, outBufP + 4);
            strncpy(outBufP + 8, parmP->name, 8);
            switch (parmP->binary_format) {
                case DD_8_BITS:
                    for (cp = (char *)(outBufP + 16),
                            fp = swpP->dat[p][r],
                            fp1 = fp + swpP->n_cells;
                            fp < fp1; cp++, fp++) {
                        *cp = Radar_ValIsData(*fp)
                            ? N_INT(scale * (*fp + bias))
                            : parmP->bad_data;
                    }
                    break;
                case DD_16_BITS:
                    for (sp = (short *)(outBufP + 16),
                            fp = swpP->dat[p][r],
                            fp1 = fp + swpP->n_cells;
                            fp < fp1; sp++, fp++) {
                        *sp = Radar_ValIsData(*fp)
                            ? N_INT(scale * (*fp + bias))
                            : parmP->bad_data;
                    }
                    break;
                case DD_24_BITS:
                    ErrMsg_Append("No support for 24 bit integers.\n");
                    goto error;
                case DD_32_BIT_FP:
                    ErrMsg_Append("No support for 32 bit floats.\n");
                    goto error;
                case DD_16_BIT_FP:
                    ErrMsg_Append("No support for 16 bit floats.\n");
                    goto error;
            }
            if (fwrite(outBufP, 1, rdat_sz, out) != rdat_sz) {
                ErrMsg_Append("Could not write ray data.\n");
                goto error;
            }
            file_sz += rdat_sz;
        }
    }

    /*
     * Write NULL and RKTB blocks
     */

    memset(outBufP, 0, null_sz);
    strncpy(outBufP, "NULL", 4);
    GetVals_PutSINT32(null_sz, outBufP + 4);
    if (fwrite(outBufP, 1, null_sz, out) != null_sz) {
        ErrMsg_Append("Could not write outBufP block.\n");
        goto error;
    }
    file_sz += null_sz;

    /*
     * Fill in the RKTB block.  This code is derived from the
     * dd_rotang_table function in soloii:translate/ddb_common.c.
     * I'm pretty clueless about why it's needed or what it's supposed to
     * do.  GDC
     */

    memset(outBufP, 0, rktb_sz);
    angle_table_offset = rot_ang_table_sz;
    angle2ndx = angle_ndx_size / 360.0;
    ndx_que_size = angle_ndx_size;
    first_key_offset = angle_table_offset + 4 * ndx_que_size;
    angle_table = (int *)(outBufP + angle_table_offset);
    rotation_table
        = (struct rotation_table_entry *)(outBufP + first_key_offset);
    strncpy(outBufP, "RKTB", 4);
    GetVals_PutSINT32(rktb_sz, outBufP + 4);
    GetVals_PutFLOAT32(angle2ndx, outBufP + 8);
    GetVals_PutSINT32(ndx_que_size, outBufP + 12);
    GetVals_PutSINT32(first_key_offset, outBufP + 16);
    GetVals_PutSINT32(angle_table_offset, outBufP + 20);
    GetVals_PutSINT32(n_good_rays, outBufP + 24);

    /*
     * Initialize angle table for RKTB block.
     */

    for(i = 0; i < angle_ndx_size; i++) {
        angle_table[i] = -1;
    }

    /*
     * Put entries in the RKTB tables
     */

    for (r = 0; r < n_rays; r++) {
        float az;

        if ( !swpP->rayHdrP[r].good ) {
            continue;
        }

        /*
         * RKTB angle table entry
         */

        az = AngleToDeg(swpP->rayHdrP[r].azimuth);
        if(az < 0.0) {
            az += 360.0;
            i = (az < 0.0) ? 0.0 : az * angle2ndx;
        } else if (az >= 360.0) {
            i = angle_ndx_size - 1;
        } else {
            i = az * angle2ndx;
        }
        if(i >= angle_ndx_size || i < 0) {
            i = 0;
        }
        angle_table[i] = r;

        /*
         * RKTB rotation table entry.
         */

        entry.angle = az;
        entry.size = ryib_sz + platfm_info_blk_sz + n_use_parms * rdat_sz;
        entry.offset = sswb_sz + vol_descr_sz + rad_descr_sz
            + n_use_parms * parm_descr_sz + celv_sz + cfac_sz + swib_sz
            + r * entry.size;
        rotation_table[r] = entry;
    }

    if (fwrite(outBufP, 1, rktb_sz, out) != rktb_sz) {
        ErrMsg_Append("Could not write outBufP block.\n");
        goto error;
    }
    file_sz += rktb_sz;

    /*
     * Put file size into SSWB block
     */

    fseek(out, 20, SEEK_SET);
    GetVals_PutUINT32(file_sz, outBufP);
    if (fwrite(outBufP, 1, 4, out) != 4) {
        ErrMsg_Append("Could not write file size.\n");
        goto error;
    }

    /*
     * Give back memory.
     */

    FREE(outBufP);
    return 1;

error:
    ErrMsg_Append("Dorade_WriteSweep failed\n");
    FREE(outBufP);
    return 0;
}

/*
 *------------------------------------------------------------------------
 *
 * Dorade_ToASCII --
 *
 *         This function dumps a sweep structure to an ASCII stream.
 *
 * Arguments:
 *        FILE *out                        - Output file or stream
 *        struct Dorade_Sweep *swpP        - Sweep to send to output
 *
 * Results:
 *         Return value is non-zero if successful, zero if something goes wrong.
 *
 * Side effects:
 *         Sweep headers and data are written to stdout in ascii format.
 *
 *------------------------------------------------------------------------
 */

int Dorade_ToASCII(FILE *out, struct Dorade_Sweep *swpP)
{
    struct GeoTime_Cal cal;
    double lat, lon;
    unsigned b, p, r;

    fprintf(out, "radar_name %s\n", swpP->radar_name);
    fprintf(out, "proj_name %s\n", swpP->proj_name);
    fprintf(out, "n_use_parms %ld\n", swpP->n_use_parms);
    fprintf(out, "vol_num %d\n", swpP->vol_num);
    cal = GeoTime_JulToCal(swpP->time);
    fprintf(out, "time %d/%d/%d %d:%d:%.1f\n",
            cal.year, cal.month, cal.day, cal.hour, cal.minute, cal.second);
    fprintf(out, "n_sensors %d\n", swpP->n_sensors);
    fprintf(out, "radar_const %.2f\n", swpP->radar_const);
    fprintf(out, "peak_power %.2f\n", swpP->peak_power);
    fprintf(out, "noise_power %.2f\n", swpP->noise_power);
    fprintf(out, "receiver_gain %.2f\n", swpP->receiver_gain);
    fprintf(out, "antenna_gain %.2f\n", swpP->antenna_gain);
    fprintf(out, "system_gain %.2f\n", swpP->system_gain);
    AnglePrintVal(out, swpP->horz_beam_width, "horz_beam_width %.2f\n",
            "horz_beam_width bad_angle\n");
    AnglePrintVal(out, swpP->vert_beam_width, "vert_beam_width %.2f\n",
            "vert_beam_width bad_angle\n");
    fprintf(out, "radar_type %d\n", swpP->radar_type);
    fprintf(out, "swpP->scan_mode = %s\n",
            (swpP->scan_mode == 0) ? "Calibration" :
            (swpP->scan_mode == 1) ? "PPI (constant Elevation)" :
            (swpP->scan_mode == 2) ? "Coplane" :
            (swpP->scan_mode == 3) ? "RHI (Constant Azimuth)" :
            (swpP->scan_mode == 4) ? "Vertical Pointing" :
            (swpP->scan_mode == 5) ? "Target (Stationary)" :
            (swpP->scan_mode == 6) ? "Manual" : "Idle (Out of Control)");
    fprintf(out, "req_rotat_vel %.2f\n", swpP->req_rotat_vel);
    fprintf(out, "compression %d\n", swpP->compression);
    fprintf(out, "data_reduction %d\n", swpP->data_reduction);
    fprintf(out, "data_red_parm0 %.2f\n", swpP->data_red_parm0);
    fprintf(out, "data_red_parm1 %.2f\n", swpP->data_red_parm1);
    GeoPtGetDeg(swpP->radar_location, &lat, &lon);
    fprintf(out, "radar_altitude %.2f\n", swpP->radar_altitude);
    fprintf(out, "eff_unamb_vel %.2f\n", swpP->eff_unamb_vel);
    fprintf(out, "eff_unamb_range %.2f\n", swpP->eff_unamb_range);
    fprintf(out, "num_freq_trans %d\n", swpP->num_freq_trans);
    fprintf(out, "num_ipps_trans %d\n", swpP->num_ipps_trans);
    fprintf(out, "freq1 %.2f\n", swpP->freq1);
    fprintf(out, "freq2 %.2f\n", swpP->freq2);
    fprintf(out, "freq3 %.2f\n", swpP->freq3);
    fprintf(out, "freq4 %.2f\n", swpP->freq4);
    fprintf(out, "freq5 %.2f\n", swpP->freq5);
    fprintf(out, "interpulse_per1 %.2f\n", swpP->interpulse_per1);
    fprintf(out, "interpulse_per2 %.2f\n", swpP->interpulse_per2);
    fprintf(out, "interpulse_per3 %.2f\n", swpP->interpulse_per3);
    fprintf(out, "interpulse_per4 %.2f\n", swpP->interpulse_per4);
    fprintf(out, "interpulse_per5 %.2f\n", swpP->interpulse_per5);
    fprintf(out, "n_cells %d\n", swpP->n_cells);
    fprintf(out, "Cell vector: ");
    for (b = 0; b < swpP->n_cells; b++) {
        fprintf(out, "%.2f ", swpP->distP[b]);
    }
    fprintf(out, "\n");
    for (p = 0; p < swpP->n_parms; p++) {
        struct Dorade_ParmDesc parm = swpP->parmP[p];

        if ( !parm.use ) {
            continue;
        }
        fprintf(out, "Parameter: name %s.  description %s.  units %s.  "
                "interpulse_time %d.  xmitted_freq %d.  ",
                parm.name, parm.description, parm.units,
                parm.interpulse_time, parm.xmitted_freq);
        Radar_PrintVal(out, parm.recvr_bandwidth, "recvr_bandwidth %.4g.  ",
                "recvr_bandwidth noval.  ");
        fprintf(out, "pulse_width %d.  polarization %d.  num_samples %d.  "
                "binary_format %d.  threshold_field %s.  ",
                parm.pulse_width, parm.polarization, parm.num_samples,
                parm.binary_format, parm.threshold_field);
        Radar_PrintVal(out, parm.threshold_value, "threshold_value %.4g.  ",
                "threshold_value noval.  ");
        fprintf(out, "scale %.2f.  bias %.2f.  bad_data %ld.  ",
                parm.scale, parm.bias, parm.bad_data);
        Radar_PrintVal(out, parm.min, "min %.4g.  ", "min noval.  ");
        fprintf(out, "\n");
    }
    AnglePrintVal(out, swpP->azimuth_corr, "azimuth_corr %.2f\n",
            "azimuth_corr bad_angle\n");
    AnglePrintVal(out, swpP->elevation_corr, "elevation_corr %.2f\n",
            "elevation_corr bad_angle\n");
    fprintf(out, "range_delay_corr %.2f\n", swpP->range_delay_corr);
    AnglePrintVal(out, swpP->longitude_corr, "longitude_corr %.2f\n",
            "longitude_corr bad_angle\n");
    AnglePrintVal(out, swpP->latitude_corr, "latitude_corr %.2f\n",
            "latitude_corr bad_angle\n");
    fprintf(out, "pressure_alt_corr %.2f\n", swpP->pressure_alt_corr);
    fprintf(out, "radar_alt_corr %.2f\n", swpP->radar_alt_corr);
    fprintf(out, "ew_gndspd_corr %.2f\n", swpP->ew_gndspd_corr);
    fprintf(out, "ns_gndspd_corr %.2f\n", swpP->ns_gndspd_corr);
    fprintf(out, "vert_vel_corr %.2f\n", swpP->vert_vel_corr);
    AnglePrintVal(out, swpP->heading_corr, "heading_corr %.2f\n",
            "heading_corr bad_angle\n");
    AnglePrintVal(out, swpP->roll_corr, "roll_corr %.2f\n",
            "roll_corr bad_angle\n");
    AnglePrintVal(out, swpP->pitch_corr, "pitch_corr %.2f\n",
            "pitch_corr bad_angle\n");
    AnglePrintVal(out, swpP->drift_corr, "drift_corr %.2f\n",
            "drift_corr bad_angle\n");
    AnglePrintVal(out, swpP->rot_angle_corr, "rot_angle_corr %.2f\n",
            "rot_angle_corr bad_angle\n");
    AnglePrintVal(out, swpP->tilt_corr, "tilt_corr %.2f\n",
            "tilt_corr bad_angle\n");
    fprintf(out, "swib_comment %s\n", swpP->swib_comment);
    fprintf(out, "sweep_num %ld\n", swpP->sweep_num);
    fprintf(out, "n_rays %ld\n", swpP->n_rays);
    fprintf(out, "n_good_rays %ld\n", swpP->n_good_rays);
    AnglePrintVal(out, swpP->start_angle, "start_angle %.2f\n",
            "start_angle bad_angle\n");
    AnglePrintVal(out, swpP->stop_angle, "stop_angle %.2f\n",
            "stop_angle bad_angle\n");
    AnglePrintVal(out, swpP->fixed_angle, "fixed_angle %.2f\n",
            "fixed_angle bad_angle\n");
    fprintf(out, "filter_flag %ld\n", swpP->filter_flag);
    for (r = 0; r < swpP->n_rays; r++) {
        struct Dorade_RayHdr ray_hdr;

        ray_hdr = swpP->rayHdrP[r];
        cal = GeoTime_JulToCal(ray_hdr.time);
        fprintf(out, "ray %d:  good %d.  time %d/%d/%d %d:%d:%.1f.  ",
                r, ray_hdr.good, cal.year, cal.month, cal.day, cal.hour,
                cal.minute, cal.second);
        AnglePrintVal(out, ray_hdr.azimuth, "azimuth %f.  ",
                "azimuth bad_angle.  ");
        AnglePrintVal(out, ray_hdr.elevation, "elevation %f.  ",
                "elevation bad_angle.  ");
        Radar_PrintVal(out, ray_hdr.peak_power, "peak_power %.4g.  ",
                "peak_power noval.  ");
        fprintf(out, "ray_status %ld.  altitude_msl %f.  altitude_agl %f.  ",
                ray_hdr.ray_status, ray_hdr.altitude_msl, ray_hdr.altitude_agl);
        AnglePrintVal(out, ray_hdr.true_scan_rate, "true_scan_rate %f.  ",
                "true_scan_rate bad_angle.  ");
        AnglePrintVal(out, ray_hdr.longitude, "longitude %f.  ",
                "longitude bad_angle.  ");
        AnglePrintVal(out, ray_hdr.latitude, "latitude %f.  ",
                "latitude bad_angle.  ");
        Radar_PrintVal(out, ray_hdr.ew_velocity, "ew_velocity %.4g.  ",
                "ew_velocity noval.  ");
        Radar_PrintVal(out, ray_hdr.ns_velocity, "ns_velocity %.4g.  ",
                "ns_velocity noval.  ");
        Radar_PrintVal(out, ray_hdr.vert_velocity, "vert_velocity %.4g.  ",
                "vert_velocity noval.  ");
        AnglePrintVal(out, ray_hdr.heading, "heading %f.  ",
                "heading bad_angle.  ");
        AnglePrintVal(out, ray_hdr.roll, "roll %f.  ", "roll bad_angle.  ");
        AnglePrintVal(out, ray_hdr.pitch, "pitch %f.  ", "pitch bad_angle.  ");
        AnglePrintVal(out, ray_hdr.drift_angle, "drift_angle %f.  ",
                "drift_angle bad_angle.  ");
        AnglePrintVal(out, ray_hdr.rotation_angle, "rotation_angle %f.  ",
                "rotation_angle bad_angle.  ");
        AnglePrintVal(out, ray_hdr.tilt, "tilt %f.  ", "tilt bad_angle.  ");
        Radar_PrintVal(out, ray_hdr.ew_horiz_wind, "ew_horiz_wind %.4g.  ",
                "ew_horiz_wind noval.  ");
        Radar_PrintVal(out, ray_hdr.ns_horiz_wind, "ns_horiz_wind %.4g.  ",
                "ns_horiz_wind noval.  ");
        Radar_PrintVal(out, ray_hdr.vert_wind, "vert_wind %.4g.  ",
                "vert_wind noval.  ");
        AnglePrintVal(out, ray_hdr.heading_change, "heading_change %f.  ",
                "heading_change bad_angle.  ");
        AnglePrintVal(out, ray_hdr.pitch_change, "pitch_change %f.  ",
                "pitch_change bad_angle.  ");
        fprintf(out, "\n");
    }
    for (p = 0; p < swpP->n_parms; p++) {
        if ( !swpP->parmP[p].use ) {
            continue;
        }
        for (r = 0; r < swpP->n_rays; r++) {
            fprintf(out, "%s %4d: ", swpP->parmP[p].name, r);
            for (b = 0; b < swpP->n_cells; b++) {
                Radar_PrintVal(out, swpP->dat[p][r][b], "%.4g ", "noval ");
            }
            fprintf(out, "\n");
        }
    }
    fprintf(out, "interpolated %d\n", swpP->interpolated);
    return 1;
}

/*
 * This comparison function is needed for the qsort call in Dorade_SweepElev.
 */

static int angle_cmp(const void *a1, const void *a2)
{
    return AngleCmp(*(Angle *)a1, *(Angle *)a2);
}

/*
 *------------------------------------------------------------------------
 *
 * Dorade_SweepElev --
 *
 *         This function computes the tilt of a sweep.
 *
 * Arguments:
 *        struct Dorade_Sweep *swpP        - a DORADE sweep.
 *
 * Results:
 *         Return value is the median tilt of rays in the sweep, or BadAngle()
 *         if something goes wrong.
 *
 * Side effects:
 *         None.
 *
 *------------------------------------------------------------------------
 */

Angle Dorade_SweepElev(struct Dorade_Sweep *swpP)
{
    Angle *elevP, *eP;
    struct Dorade_RayHdr *rayHdrP, *endP;
    float elev;

    elevP = (Angle *)MALLOC(swpP->n_rays * sizeof(float));
    for (rayHdrP = swpP->rayHdrP,
            endP = rayHdrP + swpP->n_rays,
            eP = elevP;
            rayHdrP < endP;
            rayHdrP++) {
        if ( rayHdrP->good ) {
            *eP++ = rayHdrP->elevation;
        }
    }
    qsort(elevP, swpP->n_good_rays, sizeof(Angle), angle_cmp);
    elev = elevP[swpP->n_good_rays / 2];
    FREE(elevP);
    return elev;
}

/*
 *-------------------------------------------------------------------------
 *
 * Dorade_ShiftAz --
 *
 *         This function increments the ray azimuths in a sweep structure.
 *
 * Arguments:
 *        struct Dorade_Sweep *swpP        - Sweep to modify
 *        Angle d_az                        - Amount to add to each ray azimuth
 *
 * Results:
 *         Return value is true if function succeeds, otherwise false.
 *         Values in sweep and ray headers are modified.
 *
 * Side effects:
 *         If the function fails, information is added to the radar error message.
 *
 *-------------------------------------------------------------------------
 */

int Dorade_ShiftAz(struct Dorade_Sweep *swpP, Angle d_az)
{
    struct Dorade_RayHdr *rP, *eP;

    if ( !swpP ) {
        ErrMsg_Append("Dorade_ShiftAz called with bogus sweep.\n");
        return 0;
    }

    switch (swpP->scan_mode) {
        case DORADE_PPI:
            swpP->start_angle += d_az;
            swpP->stop_angle += d_az;
            for (rP = swpP->rayHdrP, eP = rP + swpP->n_rays;
                    rP < eP;
                    rP++) {
                if (rP->good) {
                    rP->azimuth += d_az;
                }
            }
            break;
        case DORADE_RHI:
            swpP->fixed_angle += d_az;
            break;
        case DORADE_CALIBRATION:
            ErrMsg_Append("Unable to shift azimuths on "
                    "calibration scan.\n");
            return 0;
        case DORADE_COPLANE:
            ErrMsg_Append("Unable to shift azimuths on "
                    "coplane scan.\n");
            return 0;
        case DORADE_VERTICAL_POINTING:
            ErrMsg_Append("Unable to shift azimuths on "
                    "vertical pointing scan.\n");
            return 0;
        case DORADE_TARGET_MANUAL:
            ErrMsg_Append("Unable to shift azimuths on "
                    "target manual scan.\n");
            return 0;
    }
    return 1;
}

/*
 *-------------------------------------------------------------------------
 *
 * Dorade_IncrTime --
 *
 *         This function increments all of the times in a sweep structure.
 *
 * Arguments:
 *        struct Dorade_Sweep *swpP        - Sweep to modify
 *        double dt                        - Amount to add to each ray azimuth
 *
 * Results:
 *         Return value is true if function succeeds, otherwise false.
 *
 * Side effects:
 *         Times in the sweep structure are modified.
 *         If the function fails, information is added to the radar error message.
 *
 *-------------------------------------------------------------------------
 */

int Dorade_IncrTime(struct Dorade_Sweep *swpP, double dt)
{
    struct Dorade_RayHdr *rayHdrP, *endP;

    if ( !swpP ) {
        ErrMsg_Append("Dorade_IncrTime called with bogus sweep.\n");
        return 0;
    }
    GeoTime_Incr(&swpP->time, dt);
    for (rayHdrP = swpP->rayHdrP,
            endP = rayHdrP + swpP->n_rays;
            rayHdrP < endP;
            rayHdrP++) {
        if (rayHdrP->good) {
            GeoTime_Incr(&rayHdrP->time, dt);
        }
    }
    return 1;
}

/*
 *-------------------------------------------------------------------------
 *
 * Dorade_CopyParm --
 *
 *        This function creates a new parameter that is a copy of an existing one.
 *
 * Arguments:
 *        swpP                - sweep
 *        parm1                - name of the parameter to copy.
 *                          Up to 9 characters, 8 readable
 *                          characters and a terminating nul
 *        parm2                - Name for the new parameter.  Up
 *                          to 9 characters, 8 readable
 *                          characters and a terminating nul
 *        description        - Description for the new parameter.
 *                          Up to 41 characters, 40 readable
 *                          characters and a terminating nul.
 *                          If NULL, use description from
 *                          original parameter.
 *
 * Results:
 *         Return value is true if function succeeds or false if something goes
 *         wrong.  The parameter and data arrays in the sweep structure will be
 *         enlarged to accommodate the new parameter.  If a parameter with
 *         the given name already exists, the previously existing parameter
 *         will be renamed.
 *         The memory will be freed when Dorade_FreeSweep is called.
 *
 *-------------------------------------------------------------------------
 */

int Dorade_CopyParm(struct Dorade_Sweep *swpP, char parm1[], char parm2[],
        char description[])
{
    struct Dorade_ParmDesc *parmP1, *parmP2;
    size_t sz;
    long n_parms_old, n_parms_new;
    float *fp, *gp, *gp1;                /* Loop parameters */
    int p1, p2;
    float ***dat_new;                        /* New data array */

    if ( !swpP ) {
        ErrMsg_Append("Parameter copy function given "
                "non-existent DORADE structure.\n");
        return 0;
    }
    if (strncmp(parm1, parm2, 8) == 0) {
        ErrMsg_Append("Copying a parameter onto itself is forbidden.\n");
        return 0;
    }
    if ( !(parmP1 = Dorade_GetParm(swpP, parm1) ) ) {
        ErrMsg_Append("No parameter named ");
        ErrMsg_Append(parm1);
        ErrMsg_Append(".\n");
        return 0;
    }
    if ( !parmP1->use ) {
        ErrMsg_Append(parm1);
        ErrMsg_Append(" not in use.\n");
        return 0;
    }
    p1 = parmP1 - swpP->parmP;
    n_parms_old = swpP->n_parms;
    n_parms_new = n_parms_old + 1;
    p2 = n_parms_old;
    sz = n_parms_new * sizeof(struct Dorade_ParmDesc);
    swpP->parmP = (struct Dorade_ParmDesc *)REALLOC(swpP->parmP, sz);
    parmP1 = swpP->parmP + p1;
    parmP2 = swpP->parmP + p2;
    Dorade_InitParm(parmP2);
    swpP->n_parms = n_parms_new;
    swpP->n_use_parms++;
    strncpy(parmP2->name, parm2, 8);
    if (description) {
        strncpy(parmP2->description, description, 40);
    } else {
        strncpy(parmP2->description, parmP1->description, 40);
    }
    strncpy(parmP2->units, parmP1->units, 8);
    parmP2->use = 1;
    parmP2->interpulse_time = parmP1->interpulse_time;
    parmP2->xmitted_freq = parmP1->xmitted_freq;
    parmP2->recvr_bandwidth = parmP1->recvr_bandwidth;
    parmP2->pulse_width = parmP1->pulse_width;
    parmP2->polarization = parmP1->polarization;
    parmP2->num_samples = parmP1->num_samples;
    parmP2->binary_format = parmP1->binary_format;
    strncpy(parmP2->threshold_field, parmP1->threshold_field, 8);
    parmP2->threshold_value = parmP1->threshold_value;
    parmP2->scale = parmP1->scale;
    parmP2->bias = parmP1->bias;
    parmP2->bad_data = parmP1->bad_data;
    parmP2->min = parmP1->min;
    Dorade_HashParms(swpP);
    if (swpP->dat) {
        int p;

        /*
         * Create a new data array.
         */

        dat_new = Dorade_AllocDat(n_parms_new, swpP->n_rays, swpP->n_cells);
        if ( !dat_new ) {
            ErrMsg_Append("Could not reallocate data array "
                    "when adding new parameter to sweep.\n");
            return 0;
        }

        /*
         * Copy values from current sweep data array to the new array.
         * Copy data from source parameter to the new parameter.
         */

        for (p = 0; p < n_parms_old; p++) {
            for (fp = swpP->dat[p][0],
                    gp = dat_new[p][0],
                    gp1 = gp + swpP->n_rays * swpP->n_cells;
                    gp < gp1;
                    fp++, gp++) {
                *gp = *fp;
            }
        }
        for (fp = swpP->dat[p1][0],
                gp = dat_new[p2][0],
                gp1 = gp + swpP->n_rays * swpP->n_cells;
                gp < gp1;
                fp++, gp++) {
            *gp = *fp;
        }

        /*
         * Delete the old data array.
         */

        for (p = 0; p < n_parms_old; p++) {
            if (swpP->dat[p]) {
                FREE(swpP->dat[p][0]);
            }
            FREE(swpP->dat[p]);
        }
        FREE(swpP->dat);

        swpP->dat = dat_new;
    }
    return 1;
}

/*
 *-------------------------------------------------------------------------
 *
 * Dorade_ParmPlus --
 *
 *         This function adds one field to a constant or another field.
 *
 * Arguments:
 *         swpP        - a DORADE sweep
 *         parm1        - field to which an amount will be added. Sum will replace
 *                   values in this field.
 *         a_s        - string representation of a floating point value.  If not a
 *                   floating point value, then name of a field to add to parm1.
 *
 * Results:
 *         Data for one of the fields in the sweep is replaced with the sum.
 *         Return value is true is function is successful.  Otherwise, return
 *         value is false.
 *
 *-------------------------------------------------------------------------
 */

int Dorade_ParmPlus(struct Dorade_Sweep *swpP, char *parm1, char *a_s)
{
    double a;
    long p1, p2;                        /* Offsets to parameters */
    float *f1, *f1e, *f2;                /* Loop parameters */
    float ***dat;                        /* Convenience variable */
    unsigned n_cells;

    if ( !swpP ) {
        ErrMsg_Append("Parameter addition function given "
                "non-existent DORADE structure.\n");
        return 0;
    }
    n_cells = swpP->n_rays * swpP->n_cells;
    dat = swpP->dat;
    if ( !Hash_GetVal(&swpP->parm_tbl, parm1, &p1) ) {
        ErrMsg_Append("No parameter named ");
        ErrMsg_Append(parm1);
        ErrMsg_Append(" in sweep.\n");
        return 0;
    }
    if ( !swpP->parmP[p1].use ) {
        ErrMsg_Append(parm1);
        ErrMsg_Append(" not in use.\n");
        return 0;
    }
    if (sscanf(a_s, "%lf", &a) == 1) {
        for (f1 = dat[p1][0], f1e = dat[p1][0] + n_cells; f1 < f1e; f1++) {
            if (Radar_ValIsData(*f1)) {
                *f1 += a;
            }
        }
    } else if ( Hash_GetVal(&swpP->parm_tbl, a_s, &p2) ) {
        if ( !swpP->parmP[p2].use ) {
            ErrMsg_Append(a_s);
            ErrMsg_Append(" not in use.\n");
            return 0;
        }
        for (f1 = dat[p1][0], f2 = dat[p2][0], f1e = dat[p1][0] + n_cells;
                f1 < f1e;
                f1++, f2++) {
            if (Radar_ValIsData(*f1)) {
                if (Radar_ValIsData(*f2)) {
                    *f1 += *f2;
                } else {
                    *f1 = Radar_NoData();
                }
            }
        }
    } else {
        ErrMsg_Append(a_s);
        ErrMsg_Append(" is not a float value or parameter name.\n");
        return 0;
    }
    return 1;
}

/*
 *-------------------------------------------------------------------------
 *
 * Dorade_ParmMinus --
 *
 *         This function subtracts one field to a constant or another field.
 *
 * Arguments:
 *         swpP        - a DORADE sweep
 *         parm1        - field to which an amount will be subtracted. Difference will
 *                   replace values in this field.
 *         a_s        - string representation of a floating point value.  If not a
 *                   floating point value, then name of a field to subtract from
 *                   parm1.
 *
 * Results:
 *         Data for one of the fields in the sweep is replaced with the difference.
 *         Return value is true is function is successful.  Otherwise, return
 *         value is false.
 *
 *-------------------------------------------------------------------------
 */

int Dorade_ParmMinus(struct Dorade_Sweep *swpP, char *parm1, char *a_s)
{
    double a;
    long p1, p2;                        /* Offsets to parameters */
    float *f1, *f1e, *f2;                /* Loop parameters */
    float ***dat;                        /* Convenience variable */
    unsigned n_cells;

    if ( !swpP ) {
        ErrMsg_Append("Parameter subtraction function given "
                "non-existent DORADE structure.\n");
        return 0;
    }
    n_cells = swpP->n_rays * swpP->n_cells;
    dat = swpP->dat;
    if ( !(Hash_GetVal(&swpP->parm_tbl, parm1, &p1)) ) {
        ErrMsg_Append("No parameter named ");
        ErrMsg_Append(parm1);
        ErrMsg_Append(".\n");
        return 0;
    }
    if ( !swpP->parmP[p1].use ) {
        ErrMsg_Append(parm1);
        ErrMsg_Append(" not in use.\n");
        return 0;
    }
    if (sscanf(a_s, "%lf", &a) == 1) {
        for (f1 = dat[p1][0], f1e = dat[p1][0] + n_cells; f1 < f1e; f1++) {
            if (Radar_ValIsData(*f1)) {
                *f1 -= a;
            }
        }
    } else if ( Hash_GetVal(&swpP->parm_tbl, a_s, &p2) ) {
        for (f1 = dat[p1][0], f2 = dat[p2][0], f1e = dat[p1][0] + n_cells;
                f1 < f1e;
                f1++, f2++) {
            if (Radar_ValIsData(*f1)) {
                if (Radar_ValIsData(*f2)) {
                    *f1 -= *f2;
                } else {
                    *f1 = Radar_NoData();
                }
            }
        }
    } else {
        ErrMsg_Append(a_s);
        ErrMsg_Append(" is not a float value or parameter name.\n");
        return 0;
    }
    return 1;
}

/*
 *-------------------------------------------------------------------------
 *
 * Dorade_ParmMultiply --
 *
 *         This function multiplies a field by a scalar.
 *
 * Arguments:
 *         swpP        - a DORADE sweep
 *         parm        - field to scale.  Product will replace values in this field.
 *         a        - scale factor
 *
 * Results:
 *         Data for one of the fields in the sweep is replaced with the sum.
 *         Return value is true is function is successful.  Otherwise, return
 *         value is false.
 *
 *-------------------------------------------------------------------------
 */

int Dorade_ParmMultiply(struct Dorade_Sweep *swpP, char *parm, double a)
{
    long p;                                /* Offsets to parameter */
    float *f, *fe;                        /* Loop parameters */
    float ***dat;                        /* Convenience variable */
    unsigned n_cells;

    if ( !swpP ) {
        ErrMsg_Append("Parameter multiplication function given "
                "non-existent DORADE structure.\n");
        return 0;
    }
    n_cells = swpP->n_rays * swpP->n_cells;
    dat = swpP->dat;
    if ( !Hash_GetVal(&swpP->parm_tbl, parm, &p) ) {
        ErrMsg_Append("No parameter named ");
        ErrMsg_Append(parm);
        ErrMsg_Append(".\n");
        return 0;
    }
    if ( !swpP->parmP[p].use ) {
        ErrMsg_Append(parm);
        ErrMsg_Append(" not in use.\n");
        return 0;
    }
    for (f = dat[p][0], fe = dat[p][0] + n_cells; f < fe; f++) {
        if (Radar_ValIsData(*f)) {
            *f *= a;
        }
    }
    return 1;
}

/*
 *-------------------------------------------------------------------------
 *
 * Dorade_ParmLog10 --
 *
 *         This function replaces field data values with their logarithm.
 *
 * Arguments:
 *         swpP        - a DORADE sweep
 *         parm        - field to modify.
 *
 * Results:
 *         Data for one of the fields in the sweep is replaced with its base 10
 *         logarithm.  Return value is true is function is successful.  Otherwise,
 *         return value is false and information is stored in the global
 *         error message.
 *
 *-------------------------------------------------------------------------
 */

int Dorade_ParmLog10(struct Dorade_Sweep *swpP, char *parm)
{
    long p;                                /* Offsets to parameter */
    float *f, *fe;                        /* Loop parameters */
    float ***dat;                        /* Convenience variable */
    unsigned n_cells;

    if ( !swpP ) {
        ErrMsg_Append("Parameter multiplication function given "
                "non-existent DORADE structure.\n");
        return 0;
    }
    n_cells = swpP->n_rays * swpP->n_cells;
    dat = swpP->dat;
    if ( !Hash_GetVal(&swpP->parm_tbl, parm, &p) ) {
        ErrMsg_Append("No parameter named ");
        ErrMsg_Append(parm);
        ErrMsg_Append(".\n");
        return 0;
    }
    if ( !swpP->parmP[p].use ) {
        ErrMsg_Append(parm);
        ErrMsg_Append(" not in use.\n");
        return 0;
    }
    for (f = dat[p][0], fe = dat[p][0] + n_cells; f < fe; f++) {
        if (Radar_ValIsData(*f)) {
            *f = log10(*f);
        }
    }
    return 1;
}

/*
 *-------------------------------------------------------------------------
 *
 * Dorade_SweepLatLonToBin --
 *
 *        This function computes ray and bin index for the gate above a
 *        given location.
 *
 * Arguments:
 *        struct Dorade_Sweep *swpP        - (in) Sweep
 *        GeoPt geoPt                        - (in) Lat-lon for which we want gate
 *        int *rP                        - (out) Ray
 *        int *bP                        - (out) Bin
 *
 * Results:
 *        If the sweep has a gate over the given point, the ray and bin
 *        indeces are filled in with values for the gate and the function
 *        returns true.  Otherwise, return value is false.
 *
 * Side effects:
 *        None.
 *
 *-------------------------------------------------------------------------
 */

int Dorade_SweepLatLonToBin(struct Dorade_Sweep *swpP, GeoPt geoPt,
        int *rP, int *bP)
{
    GeoPt ctr;                        /* Radar location */
    double R;                        /* Distance from center of Earth to antenna */
    Angle tilt;                        /* Sweep tilt */
    Angle delta;                /* Distance along ground from radar to point */
    double d;                        /* Distance along beam from radar to point */
    Angle az;                        /* Azimuth from radar to point */
    Angle az0, az1;                /* Sweep azimuths */
    unsigned b, r;                /* Candidates for *bP and *rP */
    int n_rays, b_min;                /* Loop parameters */
    double d_min;                /* Smallest distance seen so far */
    int increasing;                /* 1 if azimuth increases from ray to ray */

    if ( !swpP || GeoPtIsNowhere(geoPt) ) {
        return 0;
    }
    ctr = swpP->radar_location;
    n_rays = swpP->n_rays;

    /*
     * Compute the bin index for the gate over the point.
     */

    R = REarth() + 1000.0 * swpP->radar_altitude;
    tilt = swpP->fixed_angle;
    delta = GeoDistance(ctr, geoPt);
    if (AngleCmp(delta + tilt, AngleFmDeg(90.0)) == 1) {
        return 0;
    }
    d = R * sin(AngleToRad(delta)) / cos(AngleToRad(delta + tilt));
    for (b = 0, b_min = -1, d_min = DBL_MAX; b < swpP->n_cells; b++) {
        double d_curr = fabs(swpP->distP[b] - d);
        if (d_curr < d_min) {
            d_min = d_curr;
            b_min = b;
        }
    }
    if (b_min < 0 || b_min >= swpP->n_cells) {
        return 0;
    }

    /*
     * Compute the ray index for the gate over the point.
     */

    az0 = swpP->rayHdrP[0].azimuth;
    az1 = swpP->rayHdrP[1].azimuth;
    increasing = AngleCmp(az1, az0) == 1;
    az = Azimuth(ctr, geoPt);
    for (r = 0; r < swpP->n_rays; r++) {
        if (swpP->rayHdrP[r].good) {
            az0 = swpP->rayHdrP[r].azimuth;
            az1 = increasing
                ? az0 + swpP->horz_beam_width
                : az0 - swpP->horz_beam_width;
            if (LonBtwn(az, az0, az1)) {
                break;
            }
        }
    }
    if (r >= swpP->n_rays) {
        return 0;
    }

    *bP = b_min;
    *rP = r;
    return 1;
}
#ifdef GEOLN

/*
 *-------------------------------------------------------------------------
 *
 * Dorade_SweepBinOutline --
 *
 *        This function computes the corners of a bin.
 *        It assumes the ray azimuth points to the center of the bin.
 *
 * Arguments:
 *        struct Dorade_Sweep *swpP        - Dorade sweep
 *        unsigned p                        - Parameter
 *        unsigned r                        - Ray
 *        unsigned b                        - Bin
 *
 * Results:
 *        A geoline (see geoLn (3)) with four points corresponding to the
 *        corners of a bin.
 *
 * Side effects:
 *        A geoline object is created.  The user should eventually destroy it
 *        with a call to GeoLnDestroy.
 *
 *-------------------------------------------------------------------------
 */

struct GeoLn * Dorade_SweepBinOutline(struct Dorade_Sweep *swpP,
        unsigned p, unsigned r, unsigned b)
{
    struct GeoLn *ln;                /* Result */
    GeoPt ctr, cnr[4];                /* Radar location.  Corners of bin. */
    float R;                        /* Distance from center of Earth to antenna */
    Angle az0, az1;                /* Ray azimuth at start and end of ray */
    Angle len;                        /* Length of bin over ground */
    Angle tilt;                        /* Ray tilt */
    double cos_tilt;                /* Cosine of sweep tilt (used several times) */
    Angle delta;                /* Distance to start of bin */
    double step;                /* Bin size */
    double d;                        /* Distance to bin */
    int n;

    if ( !swpP || p >= swpP->n_parms || r >= swpP->n_rays
            || b >= swpP->n_cells) {
        return NULL;
    }
    if ( !swpP->rayHdrP[r].good ) {
        return NULL;
    }
    if ( !(ln = GeoLnCreate(5)) ) {
        return NULL;
    }
    ctr = swpP->radar_location;
    R = REarth() + swpP->radar_altitude;
    tilt = swpP->fixed_angle + 0.5 * swpP->vert_beam_width;
    cos_tilt = ICos(tilt);
    az0 = swpP->rayHdrP[r].azimuth - swpP->horz_beam_width;
    az1 = swpP->rayHdrP[r].azimuth + swpP->horz_beam_width;
    az1 = DomainLon(az1, az0);
    d = swpP->distP[b];
    delta = AngleFmRad(asin(d * cos_tilt
                / sqrt(R * R + d * d + 2 * R * d * ISin(tilt))));
    if (b < swpP->n_cells - 1) {
        step = swpP->distP[b + 1] - swpP->distP[b];
    } else {
        step = swpP->distP[b] - swpP->distP[b - 1];
    }
    len = AngleFmRad(step * cos_tilt / REarth());
    cnr[0] = GeoStep(ctr, az1, delta);
    cnr[1] = GeoStep(ctr, az1, delta + len);
    cnr[2] = GeoStep(ctr, az0, delta + len);
    cnr[3] = GeoStep(ctr, az0, delta);
    for (n = 0; n <= 4; n++) {
        GeoLnAddPt(cnr[n % 4], ln);
    }
    return ln;
}
#endif

/*
 *------------------------------------------------------------------------
 *
 * Dorade_SweepToGrnd --
 *
 *        This function interpolates data from a sweep to a polar grid on the
 *        ground.
 *
 * Arguments:
 *        struct Dorade_Sweep *swpP                - (in) sweep to interpolate
 *        unsigned p                                - (in) param to interpolate
 *        const struct Radar_PolarGrid *gridP        - (in) grid on ground to
 *                                                  interpolate to
 *        struct Radar_PolarData *datP                - (in/out) structure to
 *                                                  receive interpolated values.
 *                                                  Should have been initialized
 *                                                  with Radar_PolarDataInit.
 *
 * Results:
 *         Return value is true if function succeeds, otherwise false.
 *
 * Side effects:
 *         The input Radar_PolarData structure receives a copy of the grid values
 *         and interpolated data.  Memory is allocated in the structure.  It
 *         should eventually be freed with a call to Radar_PolarDataFree.
 *
 *------------------------------------------------------------------------
 */

int Dorade_SweepToGrnd(struct Dorade_Sweep *swpP, unsigned p,
        const struct Radar_PolarGrid *gridP, struct Radar_PolarData *datP)
{
    float **swpRngP = NULL;                        /* Distance along ground
                                                 * to below a bin in sweep */
    Angle *swpAzP;                                /* Array of azimuths from
                                                 * sweep */
    float *f, *f1, *fe;                                /* Loop parameters */
    Angle *a, *a1, *ae;
    int i, j;                                        /* Azimuth, cell index */
    int il, ir;                                        /* Left and right azimuth
                                                 * indeces */
    int jl, jr;                                        /* Left and right gate
                                                 * indeces */
    int i0;                                        /* Index in datP->az */
    int j0;                                        /* Index in swpRngP */
    size_t sz;

    if ( !swpP ) {
        ErrMsg_Append("Dorade_SweepToGrnd called with bogus sweep.\n");
        return 0;
    }
    if ( !gridP ) {
        ErrMsg_Append("Dorade_SweepToGrnd called with bogus grid.\n");
        return 0;
    }
    if ( !datP ) {
        ErrMsg_Append("Dorade_SweepToGrnd called without place to "
                "send data.\n");
        return 0;
    }
    if (p >= swpP->n_parms) {
        ErrMsg_Append("Dorade_InterpSweep parameter index out of "
                "range.\n");
        return 0;
    }

    /*
     * Allocate and put values into swpRngP.  It will be dimensioned
     * [swpP->n_rays][swpP->n_cells].
     */

    swpRngP = (float **)MALLOC(swpP->n_rays * sizeof(float *));
    sz = swpP->n_rays * swpP->n_cells * sizeof(float);
    swpRngP[0] = (float *)MALLOC(sz);
    for (f = swpRngP[0], fe = f + swpP->n_rays * swpP->n_cells;
            f < fe; f++) {
        *f = Radar_NoData();
    }
    for (i = 1; i < swpP->n_rays; i++) {
        swpRngP[i] = swpRngP[i - 1] + swpP->n_cells;
    }
    swpAzP = (Angle *)MALLOC(swpP->n_rays * sizeof(Angle));
    for (i = 0; i < swpP->n_rays; i++) {
        if ( swpP->rayHdrP[i].good ) {
            for (j = 0; j < swpP->n_cells; j++) {
                swpRngP[i][j] = Radar_GrndRng(gridP->rng[j],
                        swpP->rayHdrP[i].altitude_msl,
                        swpP->rayHdrP[i].elevation);
            }
            swpAzP[i] = swpP->rayHdrP[i].azimuth;
        }
    }

    /*
     * Allocate the polar gridded data structure and copy grid information
     * into it.
     */

    if ( !Radar_PolarDataSetAlloc(datP, swpP->radar_location,
                gridP->n_az, gridP->n_rng) ) {
        ErrMsg_Append("Could not allocate memory to receive polar data.");
        goto error;
    }
    datP->n_az = gridP->n_az;
    for (a = gridP->az, a1 = datP->az, ae = a + gridP->n_az;
            a < ae; a++, a1++) {
        *a1 = *a;
    }
    datP->n_rng = gridP->n_rng;
    for (f = gridP->rng, f1 = datP->rng, fe = f + gridP->n_rng;
            f < fe; f++, f1++) {
        *f1 = *f;
    }

    /*
     * Find points in destination grid that are surrounded by sweep points
     * and interpolate.
     */

    for (i0 = 0; i0 < datP->n_az; i0++) {
        for (il = 0, ir = 1 ; ir < swpP->n_rays; il++, ir++) {
            if ( !swpP->rayHdrP[il].good
                    || !swpP->rayHdrP[ir].good ) {
                continue;
            }
            if (LonBtwn1(datP->az[i0], swpP->rayHdrP[il].azimuth,
                        swpP->rayHdrP[ir].azimuth)) {

                /*
                 * Ray i0 in the destination grid is surrounded
                 * by rays il and ir in the sweep.  Now step along
                 * azimuth datP->az[i0] and find the gates that
                 * surround each step.
                 */

                for (j0 = 0; j0 < datP->n_rng; j0++) {
                    jl = Radar_ArrSearch(swpRngP[il], swpP->n_cells,
                            datP->rng[j0]);
                    jr = Radar_ArrSearch(swpRngP[ir], swpP->n_cells,
                            datP->rng[j0]);
                    if (jl >= 0 && jr >= 0) {
                        /*
                         * Point at i0,j0 has points around it with known
                         * values.
                         */

                        double z00 = swpP->dat[p][il][jl];
                        double z01 = swpP->dat[p][il][jl + 1];
                        double z10 = swpP->dat[p][ir][jr];
                        double z11 = swpP->dat[p][ir][jr + 1];
                        Angle az0 = swpP->rayHdrP[il].azimuth;
                        Angle az1 = swpP->rayHdrP[ir].azimuth;

                        double r00 = swpRngP[il][jl];
                        double r01 = swpRngP[il][jl + 1];
                        double r10 = swpRngP[ir][jl];
                        double r11 = swpRngP[ir][jl + 1];

                        Angle az = datP->az[i0];
                        double r = datP->rng[j0];

                        double t, dr0_inv, dr1_inv;

                        az = DomainLon(az, az0);
                        az1 = DomainLon(az1, az0);
                        if (       Radar_ValIsNoData(z00)
                                || Radar_ValIsNoData(z01)
                                || Radar_ValIsNoData(z10)
                                || Radar_ValIsNoData(z11)
                                || fabs(r01 - r00) < FLT_EPSILON
                                || fabs(r11 - r10) < FLT_EPSILON) {
                            continue;
                        }
                        t = AngleToRad(az - az0) / AngleToRad(az1 - az0);
                        dr0_inv = 1.0 / (r01 - r00);
                        dr1_inv = 1.0 / (r11 - r10);
                        datP->dat[i0][j0]
                            = (1 - t) * (r01 - r) * dr0_inv * z00
                            + t * (r11 - r) * dr1_inv * z10
                            + t * (r - r10) * dr1_inv * z11
                            + (1 - t) * (r - r00) * dr0_inv * z01;

                    }
                }
            }
        }
    }

    /*
     * Clean up and go.
     */

    FREE(swpAzP);
    return 1;
    ErrMsg_Append("Interpolation failed.\n");
    goto error;
    FREE(swpRngP[0]);
    FREE(swpRngP);
    return 1;

error:
    if (swpRngP) {
        if (swpRngP[0]) {
            FREE(swpRngP[0]);
        }
        FREE(swpRngP);
    }
    Radar_PolarDataFree(datP);
    return 0;
}

/*
 *----------------------------------------------------------------------
 *
 * Dorade_DestroyLib --
 *
 *         This function cleans up.  It should be called when the application exits.
 *
 * Results:
 *         Memory is freed.
 *
 *----------------------------------------------------------------------
 */

void Dorade_DestroyLib(void)
{
    if (info) {
        FREE(info);
    }
    info = NULL;
    if (inBufP) {
        FREE(inBufP);
    }
    in_sz = 0;
}
