/**CHeaderFile*****************************************************************

  FileName    [cstringstream.h]

  PackageName [cstringstream]

  Synopsis    [Package for simple stringstreams in C.]

  Description [Package for simple stringstreams in C.]

  SeeAlso     []

  Author      [Fabio Somenzi]

  Copyright   [Copyright (c) 2014, Regents of the University of Colorado

  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  Neither the name of the University of Colorado nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.]

  Revision    [$Id$]

******************************************************************************/
#ifndef CSTRINGSTREAM_H_
#define CSTRINGSTREAM_H_

/*---------------------------------------------------------------------------*/
/* Type declarations                                                         */
/*---------------------------------------------------------------------------*/

typedef struct _cstringstream * cstringstream;
typedef struct _cstringstream const * const_cstringstream;

/*---------------------------------------------------------------------------*/
/* Function prototypes                                                       */
/*---------------------------------------------------------------------------*/

/* Return a new cstringstream with an empty string.
 * Return NULL if creation fails. */
cstringstream newStringStream(void);
/* Free cstringstream ss. */
void deleteStringStream(cstringstream ss);
/* Clear the contents of cstringstream ss.
 * Return 0 if succesful and -1 if ss is an invalid pointer. */
int clearStringStream(cstringstream ss);
/* Copies cstringstream src to a new cstringstream.  Return 0 if succesful
 * and -1 if src is an invalid pointer or memory allocation fails. */
cstringstream copyStringStream(const_cstringstream src);
/* Change the size of cstringstream ss.
 * Return 0 if successful and -1 if resizing fails. */
int resizeStringStream(cstringstream ss, size_t newSize);
/* Write the size of cstringstream ss to the location pointed by num.
 * Return 0 if succesful and -1 if ss is an invalid pointer. */
int sizeStringStream(const_cstringstream ss, size_t * num);
/* Write the i-th element of cstringstream ss to the location
 * pointed by c.  Return 0 if successful and -1 otherwise. */
int getStringStream(const_cstringstream ss, size_t i, char * c);
/* Add char c at the end of cstringstream ss.  Return 0 if
 * successful and -1 otherwise. */
int appendCharStringStream(cstringstream ss, char c);
/* Add string s at the end of cstringstream ss.  Return 0 if
 * successful and -1 otherwise. */
int appendStringStringStream(cstringstream ss, char const * s);
/* Add int/unsigned/long/unsigned long/double at the end of cstringstream ss.
 * Return 0 if successful and -1 otherwise. */
int appendIntStringStream(cstringstream ss, int d);
int appendUnsignedStringStream(cstringstream ss, unsigned u);
int appendLongStringStream(cstringstream ss, long ld);
int appendUnsignedLongStringStream(cstringstream ss, unsigned long lu);
int appendDoubleStringStream(cstringstream ss, double g);
/* Set the i-th element of cstringstream ss to c.  Return 0 if
 * successful and -1 otherwise.  The i-th element of ss
 * must already exist.  */
int putStringStream(cstringstream ss, size_t index, char c);
/* Return a NULL terminated string from the contents of
 * cstringstream ss.  In case of failure, return NULL.
 * The returned string must be freed by the caller. */
char * stringFromStringStream(const_cstringstream ss);
#endif
