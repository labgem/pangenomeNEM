/*\
    mexmemo.c

    Specific routines for MATLAB memory allocation

    Jan 1998

    Van Mo DANG       
    Universite de Technologie de Compiegne
    URA CNRS 817 

    Vers-mod  Date         Description

    1.05-a    26-JAN-1998  Creation
\*/


#include "mex.h"    /* mxCalloc, mxFree */
#include <stdio.h>  /* stderr */



/* ------------------------------------------------------------------- */
void* GenAlloc
(
 size_t       nelem,        /* I : number of elements to allocate */ 
 size_t       elsize,       /* I : size in bytes of each element */
 int          doexit,       /* I : 1 if failure exits, 0 if only return NULL */
 const char*  where,        /* I : name of calling function */
 const char*  what          /* I : name of allocated object */
)
{
  void *result = mxCalloc( nelem , elsize );
  if ( result != NULL )
    {
      return result ;
    }
  else
  {
    fprintf(stderr, "Fatal: in %s, no memory for %s (%d elements size %d)\n",
	    where, what, nelem, elsize);
    if ( doexit )
      mexErrMsgTxt( " " );
    else
      return result ;
  }
}
/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
void GenFree( void* ptr )
{
  if ( ptr != NULL )
    mxFree( ptr ) ;
}
/* ------------------------------------------------------------------- */
