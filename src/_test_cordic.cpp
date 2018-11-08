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
// test_cordic.cpp - test Cordic math
//
#include "Cordic.h"

using FLT = double;
using FP  = int64_t;
constexpr int64_t  FW  = 56;
constexpr FLT      TOL = 1.0 / FLT( 1LL << (FW/4) ); 

FP to_fp( FLT x )
{
    return FP( x * FLT(FP(1) << FP(FW)) );
}

FLT to_flt( FP x )
{
    return FLT( x ) / FLT(FP(1) << FP(FW));
}

#define do_op1( str, cordic_fn, exp_fn, fltx )                          \
{                                                                       \
    FP  fpx  = to_fp( fltx );			                        \
    FP  fpz  = cordic_fn( fpx );			                \
    FLT fltz = to_flt( fpz );			                        \
    FLT flte = exp_fn( fltx );			                        \
    FLT flterr = std::abs( flte-fltz );			                \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << "\n" << #str << "\n";			                \
    std::cout << "Input:    " << std::setw(30) << fltx << "(fltx)\n";	\
    std::cout << "Expected: " << std::setw(30) << flte << "\n";		\
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";		\
    std::cout << "Diff:     " << std::setw(30) << flterr << "\n";	\
    dassert( flterr <= TOL );			                        \
}    

#define do_op12( str, cordic_fn, exp_fn, fltx )                         \
{                                                                       \
    FP  fpx  = to_fp( fltx );			                        \
    FP  fpz1, fpz2;                                                     \
    cordic_fn( fpx, fpz1, fpz2 );			                \
    FLT fltz1 = to_flt( fpz1 );			                        \
    FLT fltz2 = to_flt( fpz2 );			                        \
    FLT flte1, flte2;                                                   \
    exp_fn( fltx, flte1, flte2 );			                \
    FLT flterr1 = std::abs( flte1-fltz1 );			        \
    FLT flterr2 = std::abs( flte2-fltz2 );			        \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << "\n" << #str << "\n";			                \
    std::cout << "Input:    " << std::setw(30) << fltx << "(fltx)\n";	\
    std::cout << "Expected: " << std::setw(30) << flte1 << "," << flte2 << "\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz1 << "," << fltz2 << "\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr1 << "," << flterr2 << "\n"; \
    dassert( flterr1 <= TOL );			                        \
    dassert( flterr2 <= TOL );			                        \
}    

#define do_op2( str, cordic_fn, exp_fn, fltx, flty )                    \
{                                                                       \
    FP  fpx  = to_fp( fltx );			                        \
    FP  fpy  = to_fp( flty );			                        \
    FP  fpz  = cordic_fn( fpx, fpy );			                \
    FLT fltz = to_flt( fpz );			                        \
    FLT flte = exp_fn( fltx, flty );			                \
    FLT flterr = std::abs( flte-fltz );			                \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << "\n" << #str << "\n";			                \
    std::cout << "Input:    " << std::setw(30) << fltx << "(x) " << flty << "(y)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte << "\n";		\
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";		\
    std::cout << "Diff:     " << std::setw(30) << flterr << "\n";	\
    dassert( flterr <= TOL );			                        \
}    

FLT mul( FLT x, FLT y ) { return x*y; }
FLT div( FLT x, FLT y ) { return y/x; }
FLT pow2( FLT x )       { return std::pow( 2.0, x ); }
FLT pow10( FLT x )      { return std::pow( 10.0, x ); }
FLT logb( FLT x, FLT y ){ return std::log( x ) / std::log( y ); }
FLT log2( FLT x )       { return std::log( x ) / std::log( 2.0 ); }
FLT log10( FLT x )      { return std::log( x ) / std::log( 10.0 ); }
void sin_cos( FLT x, FLT& si, FLT& co ) { si = std::sin( x ); co = std::cos( x ); }

int main( int argc, const char * argv[] )
{
    Cordic<FP, FW> cordic;
    std::cout << "tol: " << TOL << "\n";

    FLT x = 0.681807431807431031;
    FLT y = 0.810431798013170871;

    do_op2( "x*y",              cordic.mul,     mul,            x, y );
    do_op2( "y/x",              cordic.div,     div,            x, y );
    do_op1( "sqrt(x)",          cordic.sqrt,    std::sqrt,      x    );
    do_op1( "exp(x)",           cordic.exp,     std::exp,       x    );
    do_op2( "pow(x,y)",         cordic.pow,     std::pow,       x, y );
    do_op1( "pow2(x)",          cordic.pow2,    pow2,           x    );
    do_op1( "pow10(x)",         cordic.pow10,   pow10,          x    );
    do_op1( "log(x)",           cordic.log,     std::log,       x    );
    do_op2( "logb(x,b)",        cordic.logb,    logb,           x, y );
    do_op1( "log2(x)",          cordic.log2,    log2,           x    );
    do_op1( "log10(x)",         cordic.log10,   log10,          x    );
    do_op1( "sin(x)",           cordic.sin,     std::sin,       x    );
    do_op1( "cos(x)",           cordic.cos,     std::cos,       x    );
    do_op12("sin_cos(x)",       cordic.sin_cos, sin_cos,        x    );
    do_op1( "tan(x)",           cordic.tan,     std::tan,       x    );
    do_op1( "asin(x)",          cordic.asin,    std::sin,       x    );
    do_op1( "acos(x)",          cordic.acos,    std::cos,       x    );
    do_op1( "atan(x)",          cordic.atan,    std::atan,      x    );
    do_op2( "atan2(y,x)",       cordic.atan2,   std::atan2,     y, x );

#if 0
    // norm( x, y ) = sqrt( x^2 + y^2 )
    fltx = 0.956728943106177373;
    flty = 0.708473170947310947;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    fpz  = cordic.norm( fpx, fpy );
    fltz = to_flt( fpz );
    flte = std::sqrt(fltx*fltx + flty*flty);
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nnorm(x, y)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(x), " << flty << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // normh( x, y ) = sqrt( x^2 - y^2 )
    fltx = 0.708473170947310947;
    flty = 0.556728943106177373;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    fpz  = cordic.normh( fpx, fpy );
    fltz = to_flt( fpz );
    flte = std::sqrt(fltx*fltx - flty*flty);
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nnormh(x, y)\n";
    std::cout << "Input:    " << std::setw(30) << fltx  << "(x), " << flty  << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // rect_to_polar
    fltx = 0.456728943106177373;
    flty = 0.708473170947310947;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    cordic.rect_to_polar( fpx, fpy, fpz0, fpz1 );
    fltz0 = to_flt( fpz0 );
    fltz1 = to_flt( fpz1 );
    flte0 = sqrt(fltx*fltx + flty*flty);
    flte1 = atan2( flty, fltx );
    flterr0 = std::abs( flte0-fltz0 );
    flterr1 = std::abs( flte1-fltz1 );

    std::cout.precision(24);
    std::cout << "\nrect_to_polar(x, y)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(x), " << flty << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte0 << ", " << flte1 << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz0 << ", " << fltz1 << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr0 << ", " << flterr1 << "\n";
    dassert( flterr0 <= TOL );
    dassert( flterr1 <= TOL );

    // polar_to_rect
    fltx = 0.456728943106177373;
    flty = 0.708473170947310947;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    cordic.polar_to_rect( fpx, fpy, fpz0, fpz1 );
    fltz0 = to_flt( fpz0 );
    fltz1 = to_flt( fpz1 );
    flte0 = fltx*cos(flty);
    flte1 = fltx*sin(flty);
    flterr0 = std::abs( flte0-fltz0 );
    flterr1 = std::abs( flte1-fltz1 );

    std::cout.precision(24);
    std::cout << "\npolar_to_rect(r, a)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(r), " << flty << "(a)\n";
    std::cout << "Expected: " << std::setw(30) << flte0 << ", " << flte1 << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz0 << ", " << fltz1 << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr0 << ", " << flterr1 << "\n";
    dassert( flterr0 <= TOL );
    dassert( flterr1 <= TOL );

    // sinh() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.sinh( fpx );
    fltz = to_flt( fpz );
    flte = std::sinh( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nsinh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // cosh() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.cosh( fpx );
    fltz = to_flt( fpz );
    flte = std::cosh( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\ncosh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // sinh( x ), cosh( x )
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    cordic.sinh_cosh( fpx, fpz0, fpz1 );
    fltz0 = to_flt( fpz0 );
    fltz1 = to_flt( fpz1 );
    flte0 = std::sinh( fltx );
    flte1 = std::cosh( fltx );
    flterr0 = std::abs( flte0-fltz0 );
    flterr1 = std::abs( flte1-fltz1 );

    std::cout.precision(24);
    std::cout << "\nsinh_cosh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx  << "\n";
    std::cout << "Expected: " << std::setw(30) << flte0 << ", " << flte1 << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz0 << ", " << fltz1 << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr0 << ", " << flterr1 << "\n";
    dassert( flterr0 <= TOL );
    dassert( flterr1 <= TOL );

    // tanh() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.tanh( fpx );
    fltz = to_flt( fpz );
    flte = std::tanh( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\ntanh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // asinh() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.asinh( fpx );
    fltz = to_flt( fpz );
    flte = std::asinh( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nasinh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // acosh() 
    //
    fltx = 1.591370341781322;
    fpx  = to_fp( fltx );
    fpz  = cordic.acosh( fpx );
    fltz = to_flt( fpz );
    flte = std::acosh( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nacosh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // atanh() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.atanh( fpx );
    fltz = to_flt( fpz );
    flte = std::atanh( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\natanh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // atanh2() 
    //
    fltx = 0.709831990704326039;
    flty = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    fpz  = cordic.atanh2( fpy, fpx );
    fltz = to_flt( fpz );
    flte = std::atanh( flty/fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\natanh2(y, x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(x), " << flty << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );
#endif
    return 0;
}
