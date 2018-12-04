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
// test_basic.cpp - basic black-box test of Cordic.h math functions
//
#include "Cordic.h"

using FLT = double;                                     // later, use a more precise float type
using T   = int64_t;                                    // integer container to hold encoded numbers

#include "test_helpers.h"                               // must be included after FLT is defined

int main( int argc, const char * argv[] )
{
    //---------------------------------------------------------------------------
    // Process command-line arguments after applying defaults.
    //---------------------------------------------------------------------------
    int int_w  = 7;                             // fixed-point for now, max integer is 127
    int frac_w = 52;                            // same as double
    FLT TOL = 1.0 / FLT( 1LL << (frac_w-13) );  // would like this to be much smaller
    bool     new_bugs = false;                  // by default, don't run new bugs
    uint32_t loop_cnt = 2;                      // by default, run first 2 iterations

    for( int i = 1; i < argc; i++ )
    {
        if ( strcmp( argv[i], "-int_w" ) == 0 ) {
            int_w = std::atoi( argv[++i] );
        } else if ( strcmp( argv[i], "-frac_w" ) == 0 ) {
            frac_w = std::atoi( argv[++i] );
        } else if ( strcmp( argv[i], "-tol" ) == 0 ) {
            TOL = std::atof( argv[++i] );
        } else if ( strcmp( argv[i], "-new_bugs" ) == 0 ) {
            new_bugs = true;
        } else if ( strcmp( argv[i], "-loop_cnt" ) == 0 ) {
            loop_cnt = std::atoi( argv[++i] );
        } else {
            std::cout << "ERROR: unknown option " << argv[i] << "\n";
            exit( 1 );
        }
    }
    std::cout << "int_w=" << int_w << " frac_w=" << frac_w << " tol=" << TOL << "\n\n";

    //---------------------------------------------------------------------------
    // Allocate do_reduce=true and do_reduce=false Cordic objects.
    //---------------------------------------------------------------------------
    Cordic<T, FLT> * cordicr  = new Cordic( int_w, frac_w, true );     // with arg reduction
    Cordic<T, FLT> * cordicnr = new Cordic( int_w, frac_w, false );    // without arg reduction

    //---------------------------------------------------------------------------
    // New and fixed bugs.
    //---------------------------------------------------------------------------
    bool do_reduce = true;      // let routines handle the general case
    if ( new_bugs ) {
        //---------------------------------------------------------------------------
        // Put new bugs here.
        // Once fixed, they will be moved to below.
        //---------------------------------------------------------------------------
        FLT x = 0.810431798013170871;
        FLT y = 0.681807431807431031;
        do_op2(  "3) normh(x,y)",        normh,  normh,         x+0.0, y, do_reduce );          // y/x has to be less than a certain amount
        do_op1(  "   sqrt(x)",           sqrt,   std::sqrt,     x+3.0   , do_reduce );
    }
    //---------------------------------------------------------------------------
    // Put fixed bugs here so they get regressed.
    //---------------------------------------------------------------------------
    do_op2(  "2) mul",                   mul,    mul,           1.45, 0.4782,  do_reduce );
    do_op1(  "1) log",                   log,    std::log,      1.53,          do_reduce );

    //---------------------------------------------------------------------------
    // Run through all operations quickly with do_reduce=false and do_reduce=true.
    //---------------------------------------------------------------------------
    for( uint32_t i = 0; i < loop_cnt; i++ )
    {
        do_reduce = i != 0;

        FLT x  = 0.681807431807431031 + 3*i;
        FLT xs = 0.681807431807431031 + 1.23*i;
        FLT y  = 0.810431798013170871 + 3*i;
        FLT w  = 0.103301038084310970 + 3*i;
        FLT b  = M_E * 1.1;

        //                          cordic   reference
        do_op3(  "x*y + w",          mad,     mad,            x, y, w, do_reduce );
        do_op2(  "x*y",              mul,     mul,            x, y, do_reduce );
        if ( x > 0.0 ) {
            do_op3(  "y/x + w",          dad,     dad,            y, x, w, do_reduce );
            do_op2(  "y/x",              div,     div,            y, x, do_reduce );
            do_op1(  "1/x",              one_over,one_over,       x   , do_reduce );
            if ( !do_reduce ) {
            do_op1(  "sqrt(x)",          sqrt,    std::sqrt,      x   , do_reduce );
            do_op1(  "one_over_sqrt(x)", one_over_sqrt, one_over_sqrt, x, do_reduce );
            }
        }
        do_op1(  "exp(x)",           exp,     std::exp,       x   , do_reduce );
        do_op2(  "pow(x,y)",         pow,     std::pow,       b, y, do_reduce );
        do_op1(  "pow2(x)",          pow2,    pow2,           x   , do_reduce );
        do_op1(  "pow10(x)",         pow10,   pow10,          xs  , true      );
        if ( x > 0.0 ) {
            do_op1(  "log(x)",       log,     std::log,       x   , true      );
            do_op2(  "logb(x,b)",    logb,    logb,           x, b, true      );
            do_op1(  "log2(x)",      log2,    log2,           x   , true      );
            do_op1(  "log10(x)",     log10,   log10,          x   , true      );
        }
        do_op1(  "sin(x)",           sin,     std::sin,       x   , do_reduce );
        do_op1(  "cos(x)",           cos,     std::cos,       x   , do_reduce );
        do_op12( "sin_cos(x)",       sin_cos, sin_cos,        x   , do_reduce );
        if ( std::cos(x) != 0.0 ) {
            do_op1(  "tan(x)",       tan,     std::tan,       x   , do_reduce );
        }
        if ( x >= -1.0 && x <= 1.0 ) {
            do_op1(  "asin(x)",      asin,    std::asin,      x   , do_reduce );
            do_op1(  "acos(x)",      acos,    std::acos,      x   , do_reduce );
            do_op1(  "atan(x)",      atan,    std::atan,      x   , do_reduce );
        }
        do_op2(  "atan2(y,x)",       atan2,   std::atan2,     y, x, true      );
        do_op1(  "sinh(x)",          sinh,    std::sinh,      x   , do_reduce );
        do_op1(  "cosh(x)",          cosh,    std::cosh,      x   , do_reduce );
        do_op12( "sinh_cosh(x)",     sinh_cosh,sinh_cosh,     x   , do_reduce );
        do_op1(  "tanh(x)",          tanh,    std::tanh,      x   , do_reduce );
        do_op1(  "asinh(x)",         asinh,   std::asinh,     x   , true      );
        if ( x >= 1.0 ) {
            do_op1(  "acosh(x)",     acosh,   std::acosh,     x   , true      );
        }
        if ( x >= -1.0 && x <= 1.0 ) {
            do_op1(  "atanh(x)",     atanh,   std::atanh,     x   , true      );
        }
        if ( y/x >= -1.0 && y/x <= 1.0 ) {
            do_op2(  "atanh2(y,x)",  atanh2,  atanh2,         y, x, true      );
        }
        do_op2(  "norm(x,y)",        norm,    norm,           x, y, do_reduce );
        if ( 0 && y >= x ) {
            do_op2(  "normh(x,y)",   normh,   normh,          y, x, do_reduce );
        }
        do_op22( "rect_to_polar(x,y)", rect_to_polar, rect_to_polar, x, y, do_reduce );
        do_op22( "polar_to_rect(x,y)", polar_to_rect, polar_to_rect, x, y, do_reduce );
    }

    return 0;
}
