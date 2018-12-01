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
using FP  = int64_t;
constexpr int int_w = 7;                                // fixed-point for now
constexpr int frac_w = 56;                              // same as double
constexpr FLT TOL = 1.0 / FLT( 1LL << (frac_w-12) );    // would like this to be much smaller

// some useful macros to avoid redundant typing
//
#define do_op1( str, c_fn, exp_fn, fltx, do_reduce )                    \
{                                                                       \
    auto c = do_reduce ? cordicr : cordicnr;                            \
    FP  fpx  = c->to_fp( fltx );			                \
    FP  fpz  = c->c_fn( fpx );		                                \
    FLT fltz = c->to_flt( fpz );	                                \
    FLT flte = exp_fn( fltx );			                        \
    FLT flterr = std::abs( flte-fltz );			                \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Input:    " << std::setw(30) << fltx << "(fltx)\n";	\
    std::cout << "Expected: " << std::setw(30) << flte << "\n";		\
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";		\
    std::cout << "Diff:     " << std::setw(30) << flterr << "\n\n";	\
    cassert( flterr <= TOL );			                        \
}    

#define do_op12( str, c_fn, exp_fn, fltx, do_reduce )                   \
{                                                                       \
    auto c = do_reduce ? cordicr : cordicnr;                            \
    FP  fpx  = c->to_fp( fltx );			                \
    FP  fpz1, fpz2;                                                     \
    c->c_fn( fpx, fpz1, fpz2 );		                                \
    FLT fltz1 = c->to_flt( fpz1 );		                        \
    FLT fltz2 = c->to_flt( fpz2 );		                        \
    FLT flte1, flte2;                                                   \
    exp_fn( fltx, flte1, flte2 );			                \
    FLT flterr1 = std::abs( flte1-fltz1 );			        \
    FLT flterr2 = std::abs( flte2-fltz2 );			        \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Reduce:   " << do_reduce << "\n";                     \
    std::cout << "Input:    " << std::setw(30) << fltx << "(fltx)\n";	\
    std::cout << "Expected: " << std::setw(30) << flte1 << ", " << flte2 << "\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz1 << ", " << fltz2 << "\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr1 << ", " << flterr2 << "\n\n"; \
    cassert( flterr1 <= TOL );			                        \
    cassert( flterr2 <= TOL );			                        \
}    

#define do_op2( str, c_fn, exp_fn, fltx, flty, do_reduce )              \
{                                                                       \
    auto c = do_reduce ? cordicr : cordicnr;                            \
    FP  fpx  = c->to_fp( fltx );			                \
    FP  fpy  = c->to_fp( flty );			                \
    FP  fpz  = c->c_fn( fpx, fpy );	                                \
    FLT fltz = c->to_flt( fpz );			                \
    FLT flte = exp_fn( fltx, flty );			                \
    FLT flterr = std::abs( flte-fltz );			                \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Reduce:   " << do_reduce << "\n";                     \
    std::cout << "Input:    " << std::setw(30) << fltx << "(x) " << flty << "(y)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte << "\n";		\
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";		\
    std::cout << "Diff:     " << std::setw(30) << flterr << "\n\n";	\
    cassert( flterr <= TOL );			                        \
}    

#define do_op22( str, c_fn, exp_fn, fltx, flty, do_reduce )             \
{                                                                       \
    auto c = do_reduce ? cordicr : cordicnr;                            \
    FP  fpx  = c->to_fp( fltx );			                \
    FP  fpy  = c->to_fp( flty );			                \
    FP  fpz1, fpz2;                                                     \
    c->c_fn( fpx, fpy, fpz1, fpz2 );	                                \
    FLT fltz1 = c->to_flt( fpz1 );		                        \
    FLT fltz2 = c->to_flt( fpz2 );		                        \
    FLT flte1, flte2;                                                   \
    exp_fn( fltx, flty, flte1, flte2 );			                \
    FLT flterr1 = std::abs( flte1-fltz1 );			        \
    FLT flterr2 = std::abs( flte2-fltz2 );			        \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Reduce:   " << do_reduce << "\n";                     \
    std::cout << "Input:    " << std::setw(30) << fltx << "(fltx) " << flty << "(flty)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte1 << ", " << flte2 << "\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz1 << ", " << fltz2 << "\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr1 << ", " << flterr2 << "\n\n"; \
    cassert( flterr1 <= TOL );			                        \
    cassert( flterr2 <= TOL );			                        \
}    

#define do_op3( str, c_fn, exp_fn, fltx, flty, fltw, do_reduce )        \
{                                                                       \
    auto c = do_reduce ? cordicr : cordicnr;                            \
    FP  fpx  = c->to_fp( fltx );			                \
    FP  fpy  = c->to_fp( flty );			                \
    FP  fpw  = c->to_fp( fltw );			                \
    FP  fpz  = c->c_fn( fpx, fpy, fpw );	                        \
    FLT fltz = c->to_flt( fpz );			                \
    FLT flte = exp_fn( fltx, flty, fltw );			        \
    FLT flterr = std::abs( flte-fltz );			                \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Reduce:   " << do_reduce << "\n";                     \
    std::cout << "Input:    " << std::setw(30) << fltx << "(x) " << flty << "(y) " << fltw << "(w)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte << "\n";		\
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";		\
    std::cout << "Diff:     " << std::setw(30) << flterr << "\n\n";	\
    cassert( flterr <= TOL );			                        \
}    

// FLT wrapper routines for those that are not in std::
//
FLT  mad( FLT x, FLT y, FLT w ) { return x*y + w; }
FLT  mul( FLT x, FLT y ) { return x*y; }
FLT  dad( FLT x, FLT y, FLT w ) { return x/y + w; }
FLT  div( FLT x, FLT y ) { return x/y; }
FLT  one_over( FLT x )   { return 1.0/x; }
FLT  one_over_sqrt( FLT x ) { return 1.0 / std::sqrt( x ); }
FLT  pow2( FLT x )       { return std::pow( 2.0, x ); }
FLT  pow10( FLT x )      { return std::pow( 10.0, x ); }
FLT  logb( FLT x, FLT y ){ return std::log( x ) / std::log( y ); }
FLT  log2( FLT x )       { return std::log( x ) / std::log( 2.0 ); }
FLT  log10( FLT x )      { return std::log( x ) / std::log( 10.0 ); }
void sin_cos( FLT x, FLT& si, FLT& co ) { si = std::sin( x ); co = std::cos( x ); }
void sinh_cosh( FLT x, FLT& si, FLT& co ) { si = std::sinh( x ); co = std::cosh( x ); }
FLT  atanh2( FLT y, FLT x ){ return std::atanh( y/x); }
FLT  norm( FLT x, FLT y ){ return std::sqrt( x*x + y*y ); }
FLT  normh( FLT x, FLT y ){ return std::sqrt( x*x - y*y ); }
void rect_to_polar( FLT x, FLT y, FLT& r, FLT& a ) { r = std::sqrt( x*x + y*y ); a = atan2( y, x ); }
void polar_to_rect( FLT r, FLT a, FLT& x, FLT& y ) { x = r*std::cos( a ); y = r*std::sin( a ); }

int main( int argc, const char * argv[] )
{
    Cordic<FP, FLT> * cordicr  = new Cordic( int_w, frac_w, true );     // with arg reduction
    Cordic<FP, FLT> * cordicnr = new Cordic( int_w, frac_w, false );    // without arg reduction
    std::cout << "tol: " << TOL << "\n";

    //---------------------------------------------------------------------------
    // Put new bugs here, numbered, most recent first so that 
    // fixed bugs get added to this basic regression.
    //---------------------------------------------------------------------------
    bool do_reduce = true;      // let routines handle the general case
    do_op2(  "2) mul",                   mul,    mul,           1.45, 0.4782,  do_reduce );
    do_op1(  "1) log",                   log,    std::log,      1.53,          do_reduce );

    //---------------------------------------------------------------------------
    // Run through all operations quickly with do_reduce=false and do_reduce=true.
    // Not thorough at all.
    //---------------------------------------------------------------------------
    for( uint32_t i = 0; i < 1; i++ )
    {
        do_reduce = i;

        FLT x = 0.681807431807431031 + 3*i;
        FLT y = 0.810431798013170871 + 3*i;
        FLT w = 0.103301038084310970 + 3*i;
        FLT b = M_E * 1.1;

        //                          cordic   reference
        do_op3(  "x*y + w",          mad,     mad,            x, y, w, do_reduce );
        do_op2(  "x*y",              mul,     mul,            x, y, do_reduce );
        do_op3(  "y/x + w",          dad,     dad,            y, x, w, do_reduce );
        do_op2(  "y/x",              div,     div,            y, x, do_reduce );
        do_op1(  "1/x",              one_over,one_over,       x   , do_reduce );
        if ( !do_reduce ) { // reduce not working yet
        do_op1(  "sqrt(x)",          sqrt,    std::sqrt,      x   , do_reduce );
        do_op1(  "one_over_sqrt(x)", one_over_sqrt, one_over_sqrt, x, do_reduce );
        }
        
        do_op1(  "exp(x)",           exp,     std::exp,       x   , do_reduce );
        do_op2(  "pow(x,y)",         pow,     std::pow,       b, y, do_reduce );
        do_op1(  "pow2(x)",          pow2,    pow2,           x   , do_reduce );
        do_op1(  "pow10(x)",         pow10,   pow10,          x   , true      );
        do_op1(  "log(x)",           log,     std::log,       x   , true      );
        do_op2(  "logb(x,b)",        logb,    logb,           1.76380274379013, 1.439028043178590, true );
        do_op1(  "log2(x)",          log2,    log2,           x   , true      );
        do_op1(  "log10(x)",         log10,   log10,          x   , true      );

        do_op1(  "sin(x)",           sin,     std::sin,       x   , do_reduce );
        do_op1(  "cos(x)",           cos,     std::cos,       x   , do_reduce );
        do_op12( "sin_cos(x)",       sin_cos, sin_cos,        x   , do_reduce );
        do_op1(  "tan(x)",           tan,     std::tan,       x   , do_reduce );
        if ( !do_reduce ) { // reduce not working yet
        do_op1(  "asin(x)",          asin,    std::asin,      x   , do_reduce );
        do_op1(  "acos(x)",          acos,    std::acos,      x   , do_reduce );
        do_op1(  "atan(x)",          atan,    std::atan,      x   , do_reduce );
        do_op2(  "atan2(y,x)",       atan2,   std::atan2,     y, x, do_reduce );
        do_op1(  "sinh(x)",          sinh,    std::sinh,      x   , do_reduce );
        do_op1(  "cosh(x)",          cosh,    std::cosh,      x   , do_reduce );
        do_op12( "sinh_cosh(x)",     sinh_cosh,sinh_cosh,     x   , do_reduce );
        do_op1(  "tanh(x)",          tanh,    std::tanh,      x   , do_reduce );
        do_op1(  "asinh(x)",         asinh,   std::asinh,     x   , do_reduce );
        do_op1(  "acosh(x)",         acosh,   std::acosh,     1.591370341781322, do_reduce );
        do_op1(  "atanh(x)",         atanh,   std::atanh,     x   , do_reduce );
        do_op2(  "atanh2(y,x)",      atanh2,  atanh2,         0.456728943106177373, 0.709831990704326039, do_reduce );
        }
        do_op2(  "norm(x,y)",        norm,    norm,           x, y, do_reduce );
        do_op2(  "normh(x,y)",       normh,   normh,          0.708473170947310947, 0.556728943106177373, do_reduce );

        do_op22( "rect_to_polar(x,y)", rect_to_polar, rect_to_polar, x, y, do_reduce );
        do_op22( "polar_to_rect(x,y)", polar_to_rect, polar_to_rect, x, y, false );
    }

    return 0;
}
