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
// test_basic.cpp - basic black-box test of freal.h math functions
//
#include "freal.h"                                      // not used yet, just here to test build
#include "Analysis.h"
#include "AnalysisLight.h"
#include "mpint.h"

#include "test_helpers.h"                               // must be included after FLT is defined

int main( int argc, const char * argv[] )
{
    //---------------------------------------------------------------------------
    // Process command-line arguments after applying defaults.
    //---------------------------------------------------------------------------
    bool is_float = true;                       // floating-point or fixed-point?
    int  exp_or_int_w  = 8;                     // exponent width for floating point; for fixed-point max integer width
    int  frac_w = 23;                           // same as float
    FLT TOL = 1.0 / FLT( 1LL << (frac_w-3) );   // we'd like this to be 1/(1 << (frac_w-1))  (in most cases, it is)
                                                // not clear the CPU is doing the ops correctly either
    uint32_t loop_cnt = 2;                      

    for( int i = 1; i < argc; i++ )
    {
        if ( strcmp( argv[i], "-exp_or_int_w" ) == 0 ) {
            exp_or_int_w = std::atoi( argv[++i] );
        } else if ( strcmp( argv[i], "-is_float" ) == 0 ) {
            is_float = atoi( argv[++i] );
        } else if ( strcmp( argv[i], "-frac_w" ) == 0 ) {
            frac_w = std::atoi( argv[++i] );
        } else if ( strcmp( argv[i], "-tol" ) == 0 ) {
            TOL = std::atof( argv[++i] );
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
    std::cout << "exp_or_int_w=" << exp_or_int_w << " frac_w=" << frac_w << " tol=" << TOL << "\n\n";

    //---------------------------------------------------------------------------
    // Set up default freal type to use for implicit conversions.
    //---------------------------------------------------------------------------
    freal::implicit_to_set( exp_or_int_w, frac_w, is_float );
    freal::implicit_from_set( true );

    //---------------------------------------------------------------------------
    // Run through all operations.
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
                x  = 0.681807431807431031 + 3*i;
                xs = 0.681807431807431031 + 1.23*i;
                y  = 0.810431798013170871 + 3*i;
                w  = 0.103301038084310970 + 3*i;
                b  = M_E * 1.1;
                break;

            default:
                break;
        }

        //                          freal   reference
        do_op3(  "x*y + w",          fma,     fma,            x, y, w );
        do_op2(  "x*y",              mul,     mul,            x, y );
        if ( x > 0.0 ) {
            do_op3(  "y/x + w",          fda,     fda,            y, x, w );
            do_op2(  "y/x",              div,     div,            y, x );
            do_op1(  "1/x",              rcp,     rcp,            x    );
            do_op1(  "sqrt(x)",          sqrt,    std::sqrt,      x    );
            do_op1(  "rsqrt(x)",         rsqrt,   rsqrt,          x+1.0 );
            do_op1(  "cbrt(x)",          cbrt,    std::cbrt,      x    );
            do_op1(  "rcbrt(x)",         rcbrt,   rcbrt,          x+1.0 );
        }
        do_op2(  "fdim(x,y)",        fdim,    std::fdim,      x, y );
        do_op2(  "fmax(x,y)",        fmax,    std::fmax,      x, y );
        do_op2(  "fmin(x,y)",        fmin,    std::fmin,      x, y );

        do_op1(  "exp(x)",           exp,     std::exp,       x    );
        do_op1(  "expm1(x)",         expm1,   std::expm1,     x    );
        do_op1(  "exp2(x)",          exp2,    std::exp2,      x    );
        do_op1(  "exp10(x)",         exp10,   exp10,          xs   );
        do_op2(  "pow(x,y)",         pow,     std::pow,       b, y );
        if ( x > 0.0 ) {
            do_op1(  "log(x)",       log,     std::log,       x );
            do_op1(  "log1p(x)",     log1p,   std::log1p,     x );
            do_op2(  "log(x,b)",     log,     log,            x, b );
            do_op1(  "log2(x)",      log2,    log2,           x );
            do_op1(  "log10(x)",     log10,   log10,          x );
        }
        do_op1(  "sin(x)",           sin,     std::sin,       x    );
        do_op1(  "cos(x)",           cos,     std::cos,       x    );
        do_op12( "sincos(x)",        sincos,  sincos,         x    );
        if ( std::cos(x) != 0.0 ) {
            do_op1(  "tan(x)",       tan,     std::tan,       x    );
        }
        if ( x >= -1.0 && x <= 1.0 ) {
            do_op1(  "asin(x)",      asin,    std::asin,      x    );
            do_op1(  "acos(x)",      acos,    std::acos,      x    );
            do_op1(  "atan(x)",      atan,    std::atan,      x    );
        }
        do_op2(  "atan2(y,x)",       atan2,   std::atan2,     y, x );
        do_op1(  "sinh(x)",          sinh,    std::sinh,      x    );
        do_op1(  "cosh(x)",          cosh,    std::cosh,      x    );
        do_op12( "sinhcosh(x)",      sinhcosh,sinhcosh,       x    );
        do_op1(  "tanh(x)",          tanh,    std::tanh,      x    );
        do_op1(  "asinh(x)",         asinh,   std::asinh,     x    );
        if ( x >= 1.0 ) {
            do_op1(  "acosh(x)",     acosh,   std::acosh,     x    );
        }
        if ( x >= -1.0 && x <= 1.0 ) {
            do_op1(  "atanh(x)",     atanh,   std::atanh,     x    );
        }
        if ( y/x >= -1.0 && y/x <= 1.0 ) {
            do_op2(  "atanh2(y,x)",  atanh2,  atanh2,         y, x );
        }
        do_op2(  "hypot(x,y)",       hypot,   hypot,          y, x );
        if ( y >= x ) {
            do_op2(  "hypoth(x,y)",  hypoth,  hypoth,         y, x );
        }
        do_op22( "rect_to_polar(x,y)", rect_to_polar, rect_to_polar, x, y );
        do_op22( "polar_to_rect(x,y)", polar_to_rect, polar_to_rect, x, y );
    }

    //---------------------------------------------------------------------------
    // Put fixed bugs here so they get regressed.
    //---------------------------------------------------------------------------
    std::cout << "\nPAST BUGS:\n";

    do_op1(     "1) log",               log,    std::log,      1.53 );
    do_op2(     "2) mul",               mul,    mul,           1.45, 0.4782 );
    do_op2(     "3) hypoth",            hypoth, hypoth,        0.8104, 0.6818 );
    do_op1(     "4) sqrt",              sqrt,   std::sqrt,     3.8104 );
    do_op1(     "5) cos",               cos,    std::cos,      1.6214 );
    do_op2(     "6) hypot",             hypot,  hypot,         0.70710676908493042, 0.70710664987564087 );
    do_op2(     "7) pow",               pow,    std::pow,      0.004999995231628418, 0.45454543828964233 );
    do_op2(     "8) pow",               pow,    std::pow,      0.0, 0.4 );
    do_op12(    "9) sincos",            sincos, sincos,        0.887265 );
    do_op12(    "10) sincos",           sincos, sincos,        1.5707963891327381 );
    do_op22sc(  "11) sincos",           sincos, sincos,        0.7853982001543045, -0.64947649836540222 );
    do_op2(     "12) atan2(y,x)",       atan2,  std::atan2,    0.00018440188432577997, -0.62388521432876587 );
    do_op2(     "13) atan2(y,x)",       atan2,  std::atan2,    -0.000000059604651880817983, -0.84004300832748413 );
    do_op2(     "14) atan2(y,x)",       atan2,  std::atan2,    0.85865706205368042, 0 );
    do_op2(     "15) x*y",              mul,    mul,           -0.00000000000000000010980813584353493, 0.00000000000000000010980813584353493 );
    do_op2(     "16) x*y",              mul,    mul,           -2.02836e-24, -2.02836e-24 );
    do_op1(     "17) sqrt",             sqrt,   std::sqrt,     1.000000204890966415405273437500 );
    if ( 0 )  // not ready
    do_op1(     "18) atanh",            atanh,  std::atanh,    1.000000204890966415405273437500 );  
    do_op1(     "19) atanh",            atanh,  std::atanh,    1.000000204890966415405273437500/2.0 );  
    do_op2(     "20) atanh2",           atanh2, atanh2,        1.000229761004447937011718750000/2.0, 1.000000204890966415405273437500 );
    do_op2(     "21) hypot",            hypot,  hypot,         1.000000204890966415405273437500, 1.000229761004447937011718750000 );
    do_op2(     "22) hypoth",           hypoth, hypoth,        1.000000204890966415405273437500, 1.000229761004447937011718750000/2.0 );
    do_op1(     "23) rsqrt",            rsqrt,  rsqrt,         1.000058617442846298217773437500 );
    do_op1(     "24) exp",              exp,    std::exp,      1.000000204890966415405273437500 );
    do_op2(     "25) atan2",            atan2,  std::atan2,    1.000229761004447937011718750000, -1.000000204890966415405273437500 );
    do_op1(     "26) log",              log,    std::log,      1.000000204890966415405273437500 );
    do_op1(     "27) asin",             asin,   std::asin,     1.000000204890966415405273437500 );
    do_op1(     "28) acos",             acos,   std::acos,     1.000000204890966415405273437500 );
    do_op2(     "29) hypoth",           hypoth, hypoth,        1.0, 1.000000204890966415405273437500 );
    do_op2(     "30) add",              add,    add,           1.0, 1.000000204890966415405273437500 );
    do_op2(     "31) sub",              sub,    sub,           1.0, 1.000000204890966415405273437500 );
    do_op1(     "32) asin",             asin,   std::asin,     1.000000204890966415405273437500/2 );
    do_op1(     "33) acos",             acos,   std::acos,     1.000000204890966415405273437500/2 );
    do_op2(     "34) pow",              pow,    std::pow,      1.000000204890966415405273437500, 1.000229761004447937011718750000*8.0 );
    do_op12(    "35) sincos",           sincos, sincos,        1.000000204890966415405273437500 );

    std::cout << "PASSED\n";
    return 0;
}
