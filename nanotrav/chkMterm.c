/**CFile***********************************************************************
 
   FileName    [chkMterm.c]

   PackageName [ntr]

   Synopsis    [Functions to check that the minterm counts have not
   changed.]

   Description [Functions to check that the minterm counts have not
   changed during reordering.<p>
   Internal procedures included in this module:
		<ul>
		<li> check_minterms()
		</ul>
   Static procedures included in this module:
		<ul>
		<li> stFree()
		</ul>]

   SeeAlso     []

   Author      [Fabio Somenzi]
		   
   Copyright [ This file was created at the University of Colorado at
   Boulder.  The University of Colorado at Boulder makes no warranty
   about the suitability of this software for any purpose.  It is
   presented on an AS IS basis.]

******************************************************************************/

#include "ntr.h"

/*---------------------------------------------------------------------------*/
/* Constant declarations                                                     */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Stucture declarations                                                     */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Type declarations                                                         */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Variable declarations                                                     */
/*---------------------------------------------------------------------------*/

#ifndef lint
static char rcsid[] UTIL_UNUSED = "$Id: chkMterm.c,v 1.6 2004/01/01 07:06:06 fabio Exp fabio $";
#endif

/*---------------------------------------------------------------------------*/
/* Macro declarations                                                        */
/*---------------------------------------------------------------------------*/

/**AutomaticStart*************************************************************/

/*---------------------------------------------------------------------------*/
/* Static function prototypes                                                */
/*---------------------------------------------------------------------------*/

static enum st_retval stFree (char *key, char *value, char *arg);

/**AutomaticEnd***************************************************************/

/*---------------------------------------------------------------------------*/
/* Definition of exported functions                                          */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Definition of internal functions                                          */
/*---------------------------------------------------------------------------*/
    

/**Function********************************************************************

  Synopsis    [Check that minterm counts have not changed.]

  Description [Counts the minterms in the global functions of the
  primary outputs of the network passed as argument.
  When it is calld with the second argument set to NULL, it allocates
  a symbol table and stores, for each output, the minterm count. If
  an output does not have a BDD, it stores a NULL pointer for it.
  If it is called with a non-null second argument, it assumes that
  the symbol table contains the minterm counts measured previously
  and it compares the new counts to the old ones. Finally, it frees
  the symbol table.
  check_minterms is designed so that it can be called twice: once before
  reordering, and once after reordering.
  Returns a pointer to the symbol table on the first invocation and NULL
  on the second invocation.]

  SideEffects [None]

  SeeAlso     []

******************************************************************************/
st_table *
checkMinterms(
  BnetNetwork * net,
  DdManager * dd,
  st_table * previous)
{
    BnetNode *po;
    int numPi;
    char *name;
    double *count, newcount, *oldcount;
    int flag,err,i;

    numPi = net->ninputs;

    if (previous == NULL) {
	previous = st_init_table(strcmp,st_strhash);
	if (previous == NULL) {
	    (void) printf("checkMinterms out-of-memory\n");
	    return(NULL);
	}
	for (i = 0; i < net->noutputs; i++) {
	    if (!st_lookup(net->hash,net->outputs[i],&po)) {
		exit(2);
	    }
	    name = net->outputs[i];
	    if (po->dd != NULL) {
		count = ALLOC(double,1);
		*count = Cudd_CountMinterm(dd,po->dd,numPi);
		err = st_insert(previous, name, (char *) count);
	    } else {
		err = st_insert(previous, name, NULL);
	    }
	    if (err) {
		(void) printf("Duplicate input name (%s)\n",name);
		return(NULL);
	    }
	}
	return(previous);
    } else {
	flag = 0;
	if (st_count(previous) != net->noutputs) {
	    (void) printf("Number of outputs has changed from %d to %d\n",
	    st_count(previous), net->noutputs);
	    flag = 1;
	}
	for (i = 0; i < net->noutputs; i++) {
	    if (!st_lookup(net->hash,net->outputs[i],&po)) {
		exit(2);
	    }
	    name = net->outputs[i];
	    if (st_lookup(previous,name,&oldcount)) {
		if (po->dd != NULL) {
		    newcount = Cudd_CountMinterm(dd,po->dd,numPi);
		    if (newcount != *oldcount) {
			(void) printf("Number of minterms of %s has changed from %g to %g\n",name,*oldcount,newcount);
			flag = 1;
		    }
		} else {
		    if (oldcount != NULL) {
			(void) printf("Output %s lost its BDD!\n",name);
			flag = 1;
		    }
		}
	    } else {
		(void) printf("Output %s is new!\n",name);
		flag = 1;
	    }
	}
	/*st_foreach(previous,(enum st_retval)stFree,NULL);*/
	st_foreach(previous,(ST_PFSR)stFree,NULL);
	st_free_table(previous);
	if (flag) {
	    return((st_table *) 1);
	} else {
	    return(NULL);
	}
    }

} /* end of checkMinterms */


/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/


/**Function********************************************************************

  Synopsis    [Frees the data of the symbol table.]

  Description []

  SideEffects [None]

  SeeAlso     []

*****************************************************************************/
static enum st_retval
stFree(
  char *key,
  char *value,
  char *arg)
{
    if (value != NULL) {
	FREE(value);
    }
    return(ST_CONTINUE);

} /* end of stFree */


