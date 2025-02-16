/*
 *	Read_radial.c
 *
 *	This is a Fortran interface version.
 *
 */
 
#pragma ident	"@(#)read_radial_fortran.c	5.1	03/01/05	CAPS"

read_radial_f( irc )
int *irc;
{
	*irc = read_radial();
	return;
}

read_radial_f_( irc )
int *irc;
{
	read_radial_f( irc );
	return;
}
