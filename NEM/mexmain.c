/*\

    mexmain.c

    A mexFunction() interface to compile a matlab routine which is
    executable from matlab.

    Jan 1998

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

Vers-mod  Date         Who  Description

1.05-a    25-JAN-1998  MD   Create to make an interface to matlab
1.05-b    30-JAN-1998  MD   Prototype of called mainfunc() in mainfunc.h
\*/


#include "mex.h"
#include "mainfunc.h"   /* mainfunc() */

#include <stdio.h>      /* printf, ... */
#include <string.h>     /* strlen, ... */


#define  SEPAR_S    " \t"
#define  CMDNAME_S  "mexprog"


/* Local functions */

static void make_args
(
 const char* cmd,      /* I : command line to analyse */
 int*        argc_p,   /* O : number of arguments */
 char***     argv_p    /* O : array of argc arguments (allocated) */
) ;

static int fetch_args  /* returns number of counted strings */
( 
 const char* cmd,  /* I : command line to analyze */
 int         nmax, /* I : number of allocated args in argv */
 char**      argv  /* O : table [0..nmax-1] strings (each string alloc here)*/
) ;



/* 
 *   Function called from matlab.
 *   A string parameter is expected = the command line string
 *   to invoke the executable program.
 */

void mexFunction(
                 int nlhs, Matrix *plhs[],
                 int nrhs, Matrix *prhs[]
		 )
{
  int    len ;   /* length of given argument string */
  char*  arg_s ; /* given argument string - allocated locally */
  char*  cmd ;   /* cmd + given argument string - allocated locally */

  int    argc ;  /* number of arguments + 1 in string */
  char** argv ;  /* array of the arguments (+ "command") - allocated locally */
  int    iarg ;  /* counter from 0 to argc-1 to free elts of argv */

  int    sts ;   /* status returned from mainfunc */

  /* If args given in a string */
  if ( nrhs >= 1 )
    {
      /* Allocate and get string */
      len = mxGetN( prhs[0] ) + strlen( CMDNAME_S ) + 1 ;
      arg_s = mxCalloc( len, sizeof( char ) ) ;
      if ( arg_s == NULL )
	mexErrMsgTxt( "Could not allocate string of args" ) ;
      mxGetString( prhs[0], arg_s, len ) ;

      cmd = mxCalloc( strlen( arg_s ) + strlen( CMDNAME_S ) + 1 + 1, 
		      sizeof( char ) ) ;
      if ( cmd == NULL )
	mexErrMsgTxt( "Could not allocate 2nd string of args" ) ;
      strcpy( cmd, CMDNAME_S ) ;
      strcat( cmd, " " ) ;
      strcat( cmd, arg_s ) ;

      /* Translate string into mainfunc args */
      make_args( cmd, & argc, & argv ) ;

      /* Free strings */
      mxFree( cmd ) ;
      mxFree( arg_s ) ;
    }
  else  /* no arg */
    {
      /* Translate no arg into mainfunc args */
      make_args( CMDNAME_S, & argc, & argv ) ;
    }

  /* Call function */
  printf( "=== Call mainfunc with %d args :\n===  " , argc ) ;
  for ( iarg = 0 ; iarg < argc ; iarg ++ )
    {
      printf( "[%s] ", argv[ iarg ] ) ;
    }
  printf( "\n" ) ;

  sts = mainfunc( argc, (const char* *) argv ) ;

  printf( "=== Returned from mainfunc (status = %d) \n", sts ) ;

  /* Free each element of array of args and the array of pointers itself */
/*  for ( iarg = 0 ; iarg < argc ; iarg ++ )
/*    {
/*	mxFree( argv[ iarg ] ) ;
/*    }
/*  mxFree( argv ) ;
 */

}



/* ------------------------------------------------------------------------ */
static void make_args
(
 const char* cmd,      /* I : command line to analyse */
 int*        argc_p,   /* O : number of arguments */
 char***     argv_p    /* O : array of argc arguments (allocated) */
)
/*\

    A command line CMD is taken as input.  Each distinct string within
    CMD is copied into a slot of (*ARGV_P).  After the call, the table
    (*ARGV_P) contain (*ARGC_P) allocated strings.

\*/
/* ------------------------------------------------------------------------ */
{
  char** tabstr ; /* table of n strings - allocate and copy in (*argv_p) */  

  /*
   *  1 : Count number of strings 
   */
  (*argc_p) = fetch_args( cmd, 0, NULL ) ;
  
  /* 
   *  2 : Allocate and fill array of these strings 
   */
  (*argv_p) = mxCalloc( (*argc_p), sizeof( char* ) ) ;
  if ( tabstr == NULL )
    mexErrMsgTxt( "Could not allocate array of args" ) ;

  fetch_args( cmd, (*argc_p), (*argv_p) ) ;
}


/* ------------------------------------------------------------------------ */
static int fetch_args  /* returns number of counted strings */
( 
 const char* cmd,  /* I : command line to analyze */
 int         nmax, /* I : number of allocated args in argv */
 char**      argv  /* O : table [0..nmax-1] strings (each string alloc here)*/
)
/*\ 

    A command line string CMD is taken as input.  If ARGV is not NULL,
    a string is allocated in each of the NMAX slots of ARGV, and the
    distinct strings of CMD are copied in those slots.  The function
    returns the number of distinct strings in CMD.

    Typically, first call fetch_args() to count in N the number of
    strings of CMD, then allocate a table of N pointers, then call
    again fetch_args() to copy the strings of CMD into the table.

\*/
/* ------------------------------------------------------------------------ */
{
  char*  copy_s ; /* copy of input string - allocated locally */
  char*  this_s ; /* pointer to current arg in copy_s */

  int    n ;      /* number of given strings - to copy in (*argc_p)*/

  copy_s = mxCalloc( strlen( cmd ) + 1 , sizeof( char ) ) ;
  if ( copy_s == NULL )
    mexErrMsgTxt( "Could not allocate copy of string of args" ) ;

  strncpy( copy_s, cmd, strlen( cmd ) + 1 ) ;
  n = 0 ;
  for ( this_s = strtok( copy_s , SEPAR_S ) ;
	this_s != NULL ;
	this_s = strtok( NULL , SEPAR_S ) )
    {
      if ( ( argv != NULL ) && ( n < nmax ) )
	{
	  argv[ n ] = mxCalloc( strlen( this_s ) + 1 , sizeof( char ) ) ;
	  if ( argv[ n ] == NULL )
	    mexErrMsgTxt( "Could not allocate an elt of array of args" ) ;

	  strncpy( argv[ n ] , this_s , strlen( this_s ) + 1 ) ;
	}
      
      n ++ ;
    }

  mxFree( copy_s ) ;

  return n ;
}

