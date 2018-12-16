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
// test_mpint.cpp - test multi-precision integer class (which has limited semantics on purpose)
//
#include "mpint.h"

int main( int argc, const char * argv[] )
{
    //---------------------------------------------------------------------------
    // Process command-line arguments after applying defaults.
    //---------------------------------------------------------------------------
    int int_w = 63;                            

    for( int i = 1; i < argc; i++ )
    {
        if ( strcmp( argv[i], "-int_w" ) == 0 ) {
            int_w = std::atoi( argv[++i] );
        } else {
            std::cout << "ERROR: unknown option " << argv[i] << "\n";
            exit( 1 );
        }
    }
    mpint::implicit_int_w_set( int_w );
    std::cout << "int_w=" << int_w << "\n";
    
    //---------------------------------------------------------------------------
    // Simple stuff.
    //---------------------------------------------------------------------------
    mpint x = 123456789;
    std::cout << "x should be 123456789: " << x << "\n";
    x = mpint::to_mpint( "1125899906842623" );
    std::cout << "x should be 1125899906842623 (max 63-bit signed int): " << x << "\n";
    return 0;
}
