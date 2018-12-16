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
    (void)argc;
    (void)argv;

    //---------------------------------------------------------------------------
    // Some small numbers that fit in 64 bits.
    //---------------------------------------------------------------------------
    {
        int int_w = 64;                            
        mpint::implicit_int_w_set( int_w );
        std::cout << "\nint_w=" << int_w << "\n";
        mpint y0 = 123456789;
        mpint y1 = 937183901;
        std::cout << "y0 should be 123456789: " << y0 << "\n";
        std::cout << "y1 should be 937183901: " << y1 << "\n";
        mpint y = y0 + y1;
        std::cout << "sum y="  << y << "\n";
        std::cout << "subtract y1=" << y1 << "\n";
        mpint z = y - y1;
        std::cout << "should get y0=" << z << "\n";
        iassert( z == y0, "y - y1 != y0" );

        y0 = mpint::to_mpint( "9223372036854775807" );
        std::cout << "y0 should be 9223372036854775807 (may positive 64-bit int): " << y0 << "\n";
        y0 = mpint(1) << 63;
        std::cout << "y0 should be 1 << 63 (may negative 64-bit int): " << y0 << "\n";
    }

    //---------------------------------------------------------------------------
    // Use medium-sized numbers (> 64 bits).
    //---------------------------------------------------------------------------
    {
        int int_w = 256;
        mpint::implicit_int_w_set( int_w );
        std::cout << "\nint_w=" << int_w << "\n";
        mpint y0 = mpint::to_mpint( "748310705147189301347207" );
        mpint y1 = mpint::to_mpint( "07084217431764890316589383" );
        std::cout << "         y0=" << y0 << "\n" << y0.to_string(2) << "\n";
        std::cout << "         y1=" << y1 << "\n" << y1.to_string(2) << "\n";
        mpint y  = y0 + y1;
        std::cout << "          y="  << y  << "\n" << y.to_string(2) << "\n";
        std::cout << "subtract y1=" << y1 << "\n";
        y1 = -y1;
        std::cout << "add    -y1=" << y1 << "\n" << y1.to_string(2) << "\n";
        mpint y11 = -y1;
        std::cout << "note:  --y1=" << y11 << "\n" << y11.to_string(2) << "\n";
        mpint z = y + y1;
        std::cout << "expect   y0=" << z << "\n" << z.to_string(2) << "\n";
        iassert( z == y0, "y - y1 != y0" );
    }

    //---------------------------------------------------------------------------
    // Use some really big numbers.
    //---------------------------------------------------------------------------
    {
        int int_w = 1024;
        mpint::implicit_int_w_set( int_w );
        std::cout << "\nint_w=" << int_w << "\n";
        mpint y0 = mpint::to_mpint( "7483107051471893013472076086841320964319066409318609463216094321" );
        mpint y1 = mpint::to_mpint( "07084217431764890316589361985618902365849612894613098648906312098460983216" );
        std::cout << "y0=" << y0 << "\n";
        std::cout << "y1=" << y1 << "\n";
        mpint y  = y0 + y1;
        std::cout << "y="  << y << "\n";
        std::cout << "subtract y1=" << y1 << "\n";
        mpint z = y - y1;
        std::cout << "should get y0=" << z << "\n";
        iassert( z == y0, "y - y1 != y0" );
    }
    return 0;
}
