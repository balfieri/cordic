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
#include "freal.h"                                      // not used yet, just here to test build
#include "mpint.h"

using FLT  = double;                                    // later, use a more precise float type
using T    = __int128_t;                                // integer container to hold encoded numbers (int64_t also works)
using real = freal<>;                                   // real number (not used yet)

#include "test_helpers.h"                               // must be included after FLT is defined

int main( int argc, const char * argv[] )
{
    //---------------------------------------------------------------------------
    // Process command-line arguments after applying defaults.
    //---------------------------------------------------------------------------
    int int_w  = 7;                             // fixed-point for now, max integer is 127
    int frac_w = 52;                            // same as double
    FLT TOL = 1.0 / FLT( 1LL << (frac_w-7) );   // we'd like this to be 1/(1 << (frac_w-1)) 
                                                // not clear the CPU is doing the ops correctly either
    bool     new_bugs = false;                  // by default, don't run new bugs
    uint32_t loop_cnt = 2;                      

    real unused = real::make_fixed( int_w, frac_w, 1.78302 ); // smoke test on freal type; freal is not used yet

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
        } else if ( strcmp( argv[i], "-log" ) == 0 ) {
            Cordic<T,FLT>::logger_set( new Logger<T,FLT>( Cordic<T,FLT>::op_to_str, "" ) );
        } else if ( strcmp( argv[i], "-log_file" ) == 0 ) {
            Cordic<T,FLT>::logger_set( new Logger<T,FLT>( Cordic<T,FLT>::op_to_str, argv[++i] ) );
        } else {
            std::cout << "ERROR: unknown option " << argv[i] << "\n";
            exit( 1 );
        }
    }
    std::cout << "int_w=" << int_w << " frac_w=" << frac_w << " tol=" << TOL << "\n\n";

    //---------------------------------------------------------------------------
    // Allocate do_reduce=true and do_reduce=false Cordic objects.
    //---------------------------------------------------------------------------
    Cordic<T, FLT> * cordicr  = new Cordic<T,FLT>( int_w, frac_w, true );     // with arg reduction
    Cordic<T, FLT> * cordicnr = new Cordic<T,FLT>( int_w, frac_w, false );    // without arg reduction

    //---------------------------------------------------------------------------
    // New and fixed bugs.
    //---------------------------------------------------------------------------
    bool do_reduce = true;      // let routines handle the general case
    if ( new_bugs ) {
        //---------------------------------------------------------------------------
        // Put new bugs here.
        // Once fixed, they will be moved to below.
        //---------------------------------------------------------------------------
    }
    //---------------------------------------------------------------------------
    // Put fixed bugs here so they get regressed.
    //---------------------------------------------------------------------------
    do_op2(  "7) pow",         pow,    std::pow,      0.004999995231628418, 0.45454543828964233, do_reduce );
    do_op2(  "6) norm",        norm,   norm,          0.70710676908493042, 0.70710664987564087, do_reduce );
    do_op1(  "5) cos",         cos,    std::cos,      1.6214,           do_reduce );
    do_op1(  "4) sqrt",        sqrt,   std::sqrt,     3.8104,           do_reduce );
    do_op2(  "3) normh",       normh,  normh,         0.8104, 0.6818,   do_reduce );          
    do_op2(  "2) mul",         mul,    mul,           1.45, 0.4782,     do_reduce );
    do_op1(  "1) log",         log,    std::log,      1.53,             do_reduce );

    //---------------------------------------------------------------------------
    // Run through all operations quickly with do_reduce=true.
    // Trying to run with do_reduce=false is painful to deal with.
    //---------------------------------------------------------------------------
    for( uint32_t i = 0; i < loop_cnt; i++ )
    {
        FLT x  = 0.681807431807431031 + 3*i;
        FLT xs = 0.681807431807431031 + 1.23*i;
        FLT y  = 0.810431798013170871 + 3*i;
        FLT w  = 0.103301038084310970 + 3*i;
        FLT b  = M_E * 1.1;
        switch( i ) 
        {
            case 0:
            case 1:
            {
                x  = 0.681807431807431031 + 3*i;
                xs = 0.681807431807431031 + 1.23*i;
                y  = 0.810431798013170871 + 3*i;
                w  = 0.103301038084310970 + 3*i;
                b  = M_E * 1.1;
            }
            break;
        }

        //                          cordic   reference
        do_op3(  "x*y + w",          mad,     mad,            x, y, w, do_reduce );
        do_op2(  "x*y",              mul,     mul,            x, y, do_reduce );
        if ( x > 0.0 ) {
            do_op3(  "y/x + w",          dad,     dad,            y, x, w, do_reduce );
            do_op2(  "y/x",              div,     div,            y, x, do_reduce );
            do_op1(  "1/x",              rcp,     rcp,            x   , do_reduce );
            do_op1(  "sqrt(x)",          sqrt,    std::sqrt,      x   , do_reduce );
            do_op1(  "rsqrt(x)",         rsqrt,   rsqrt,          x+1.0, do_reduce );
            do_op1(  "cbrt(x)",          cbrt,    std::cbrt,      x   , do_reduce );
            do_op1(  "rcbrt(x)",         rcbrt,   rcbrt,          x+1.0, do_reduce );
        }
        do_op2(  "fdim(x,y)",        fdim,    std::fdim,      x, y, do_reduce );
        do_op2(  "fmax(x,y)",        fmax,    std::fmax,      x, y, do_reduce );
        do_op2(  "fmin(x,y)",        fmin,    std::fmin,      x, y, do_reduce );

        do_op1(  "exp(x)",           exp,     std::exp,       x   , do_reduce );
        do_op1(  "expm1(x)",         expm1,   std::expm1,     x   , do_reduce );
        do_op1(  "exp2(x)",          exp2,    std::exp2,      x   , do_reduce );
        do_op1(  "exp10(x)",         exp10,   exp10,          xs  , do_reduce );
        do_op2(  "pow(x,y)",         pow,     std::pow,       b, y, do_reduce );
        if ( x > 0.0 ) {
            do_op1(  "log(x)",       log,     std::log,       x   , true      );
            do_op2(  "logb(x,b)",    logb,    logb,           x, b, true      );
            do_op1(  "log2(x)",      log2,    log2,           x   , true      );
            do_op1(  "log10(x)",     log10,   log10,          x   , true      );
        }
        do_op1(  "sin(x)",           sin,     std::sin,       x   , do_reduce );
        do_op1(  "cos(x)",           cos,     std::cos,       x   , do_reduce );
        do_op12( "sincos(x)",        sincos,  sincos,         x   , do_reduce );
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
        do_op12( "sinhcosh(x)",      sinhcosh,sinhcosh,       x   , do_reduce );
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
        do_op2(  "norm(x,y)",        norm,    norm,           y, x, do_reduce );
        if ( y >= x ) {
            do_op2(  "normh(x,y)",   normh,   normh,          y, x, true      );
        }
        do_op22( "rect_to_polar(x,y)", rect_to_polar, rect_to_polar, x, y, do_reduce );
        do_op22( "polar_to_rect(x,y)", polar_to_rect, polar_to_rect, x, y, do_reduce );
    }

    delete cordicr;
    delete cordicnr;

    return 0;
}
