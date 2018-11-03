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
#include "Misc.h"

void my_exit( int status )  // can set a breakpoint here
{
    exit( 1 );
}

int int_log2( int n )
{
    dassert( n > 0 );

    n += n-1;

    int result = 0;
    while( n != 1 ) 
    {
        n >>= 1;
        result++;
    }
    return result;
}

bool int_is_pow2( int n )
{
    return (1 << int_log2( n )) == n;
}

int int_sqrt( int n )
{
    return int( sqrt( double( n ) ) );
}

int rand_n( int n )
{
    int r = rand() & 0x7fffffff;
    return r % n;
}

int heads( void )
{
    return rand_n( 2 );
}
