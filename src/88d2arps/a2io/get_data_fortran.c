/*
 *	Get_data_fortran.c
 *
 *	This is a Fortran interface version.
 *
 *	Do *NOT* include the "len" argument in any Fortran calls!
 *	It is supplied automatically.
 */
 
#pragma ident	"@(#)get_data_fortran.c	5.1	03/01/05	CAPS"

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <strings.h>

get_data_field_f( num, data, n, irc, len )
int *num;
float data[];
int *n;
int *irc;
{
	*irc = get_data_field( *num, data, *n );
	return;
}

get_data_field_f_( num, data, n, irc, len )
int *num;
float data[];
int *n;
int *irc;
{
	get_data_field_f( num, data, n, irc, len );
	return;
}

get_data_field_raw_f( num, data, n, irc, len )
int *num;
unsigned short int data[];
int *n;
int *irc;
{
	*irc = get_data_field_raw( *num, data, *n );
	return;
}

get_data_field_raw_f_( num, data, n, irc, len )
int *num;
unsigned short int data[];
int *n;
int *irc;
{
	get_data_field_raw_f( num, data, n, irc, len );
	return;
}

get_data_field_uchar_f( num, data, n, irc, len )
int *num;
unsigned char data[];
int *n;
int *irc;
{
	*irc = get_data_field_uchar( *num, data, *n );
	return;
}

get_data_field_uchar_f_( num, data, n, irc, len )
int *num;
unsigned char data[];
int *n;
int *irc;
{
	get_data_field_uchar_f( num, data, n, irc, len );
	return;
}

