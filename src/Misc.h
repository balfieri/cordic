// Copyright (c) 2014-2019 Robert A. Alfieri
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
// 

// Misc.h - miscellaneous functions
//
#ifndef _Misc_h
#define _Misc_h

//------------------------------------------
// Standard C/C++ Functions assumed to be available on all platforms
//------------------------------------------
#include "math.h"
#include <stdio.h>	
#include <stdlib.h>
#include <string.h>

#include <cmath>
#include <iostream>
#include <iomanip>

//------------------------------------------
// Debugging
//------------------------------------------
#define dprintf if ( 0 ) printf
#define dassert(expr) if ( !(expr) ) { printf( "ERROR: assertion failure: '%s' at %s:%d\n", #expr, __FILE__, __LINE__ ); my_exit( 1 ); }
#define error(msg) { printf( "ERROR: %s (%s:%d)\n", msg, __FILE__, __LINE__ ); my_exit( 1 ); }

extern void my_exit( int status );

//------------------------------------------
// Math
//------------------------------------------
extern int  int_log2( int n );
extern bool int_is_pow2( int n );
extern int  int_sqrt( int n );
extern int  rand_n( int n );  // returns 0 .. n-1
extern int  heads( void );    // returns 0 or 1

#endif
