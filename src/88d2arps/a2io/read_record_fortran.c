/*
 *	Read_record.c
 *
 *	This is a Fortran interface version.
 *
 */
 
#pragma ident	"@(#)read_record_fortran.c	5.1	03/01/05	CAPS"

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <strings.h>

read_record_f( irc )
int *irc;
{
	*irc = read_record();
	return;
}

read_record_f_( irc )
int *irc;
{
	read_record_f( irc );
	return;
}

set_console_messages_f( n )
int *n;
{
	set_console_messages( *n );
	return;
}

set_console_messages_f_( n )
int *n;
{
	set_console_messages_f( n );
	return;
}
