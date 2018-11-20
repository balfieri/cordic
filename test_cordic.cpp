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
constexpr int INT_W  = 7;
constexpr int FRAC_W = 56;
constexpr FLT TOL = 1.0 / FLT( 1LL << (FRAC_W/2) ); 

#define do_op1( str, cordic_fn, exp_fn, fltx, do_reduce )               \
{                                                                       \
    FP  fpx  = cordic.to_fp( fltx );			                \
    FP  fpz  = cordic_fn( fpx, do_reduce );		                \
    FLT fltz = cordic.to_flt( fpz );	                                \
    FLT flte = exp_fn( fltx );			                        \
    FLT flterr = std::abs( flte-fltz );			                \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Input:    " << std::setw(30) << fltx << "(fltx)\n";	\
    std::cout << "Expected: " << std::setw(30) << flte << "\n";		\
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";		\
    std::cout << "Diff:     " << std::setw(30) << flterr << "\n\n";	\
    dassert( flterr <= TOL );			                        \
}    

#define do_op12( str, cordic_fn, exp_fn, fltx, do_reduce )              \
{                                                                       \
    FP  fpx  = cordic.to_fp( fltx );			                \
    FP  fpz1, fpz2;                                                     \
    cordic_fn( fpx, fpz1, fpz2, do_reduce );		                \
    FLT fltz1 = cordic.to_flt( fpz1 );		                        \
    FLT fltz2 = cordic.to_flt( fpz2 );		                        \
    FLT flte1, flte2;                                                   \
    exp_fn( fltx, flte1, flte2 );			                \
    FLT flterr1 = std::abs( flte1-fltz1 );			        \
    FLT flterr2 = std::abs( flte2-fltz2 );			        \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Input:    " << std::setw(30) << fltx << "(fltx)\n";	\
    std::cout << "Expected: " << std::setw(30) << flte1 << "," << flte2 << "\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz1 << "," << fltz2 << "\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr1 << "," << flterr2 << "\n\n"; \
    dassert( flterr1 <= TOL );			                        \
    dassert( flterr2 <= TOL );			                        \
}    

#define do_op2( str, cordic_fn, exp_fn, fltx, flty, do_reduce )         \
{                                                                       \
    FP  fpx  = cordic.to_fp( fltx );			                \
    FP  fpy  = cordic.to_fp( flty );			                \
    FP  fpz  = cordic_fn( fpx, fpy, do_reduce );	                \
    FLT fltz = cordic.to_flt( fpz );			                \
    FLT flte = exp_fn( fltx, flty );			                \
    FLT flterr = std::abs( flte-fltz );			                \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Input:    " << std::setw(30) << fltx << "(x) " << flty << "(y)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte << "\n";		\
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";		\
    std::cout << "Diff:     " << std::setw(30) << flterr << "\n\n";	\
    dassert( flterr <= TOL );			                        \
}    

#define do_op22( str, cordic_fn, exp_fn, fltx, flty, do_reduce )        \
{                                                                       \
    FP  fpx  = cordic.to_fp( fltx );			                \
    FP  fpy  = cordic.to_fp( flty );			                \
    FP  fpz1, fpz2;                                                     \
    cordic_fn( fpx, fpy, fpz1, fpz2, do_reduce );	                \
    FLT fltz1 = cordic.to_flt( fpz1 );		                        \
    FLT fltz2 = cordic.to_flt( fpz2 );		                        \
    FLT flte1, flte2;                                                   \
    exp_fn( fltx, flty, flte1, flte2 );			                \
    FLT flterr1 = std::abs( flte1-fltz1 );			        \
    FLT flterr2 = std::abs( flte2-fltz2 );			        \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Input:    " << std::setw(30) << fltx << "(fltx) " << flty << "(flty)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte1 << "," << flte2 << "\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz1 << "," << fltz2 << "\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr1 << "," << flterr2 << "\n\n"; \
    dassert( flterr1 <= TOL );			                        \
    dassert( flterr2 <= TOL );			                        \
}    

#define do_op3( str, cordic_fn, exp_fn, fltx, flty, fltw, do_reduce )   \
{                                                                       \
    FP  fpx  = cordic.to_fp( fltx );			                \
    FP  fpy  = cordic.to_fp( flty );			                \
    FP  fpw  = cordic.to_fp( fltw );			                \
    FP  fpz  = cordic_fn( fpx, fpy, fpw, do_reduce );	                \
    FLT fltz = cordic.to_flt( fpz );			                \
    FLT flte = exp_fn( fltx, flty, fltw );			        \
    FLT flterr = std::abs( flte-fltz );			                \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Input:    " << std::setw(30) << fltx << "(x) " << flty << "(y) " << fltw << "(w)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte << "\n";		\
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";		\
    std::cout << "Diff:     " << std::setw(30) << flterr << "\n\n";	\
    dassert( flterr <= TOL );			                        \
}    

FLT mad( FLT x, FLT y, FLT w ) { return x*y + w; }
FLT mul( FLT x, FLT y ) { return x*y; }
FLT dad( FLT x, FLT y, FLT w ) { return x/y + w; }
FLT div( FLT x, FLT y ) { return x/y; }
FLT one_over( FLT x )   { return 1.0/x; }
FLT one_over_sqrt( FLT x ) { return 1.0 / std::sqrt( x ); }
FLT pow2( FLT x )       { return std::pow( 2.0, x ); }
FLT pow10( FLT x )      { return std::pow( 10.0, x ); }
FLT logb( FLT x, FLT y ){ return std::log( x ) / std::log( y ); }
FLT log2( FLT x )       { return std::log( x ) / std::log( 2.0 ); }
FLT log10( FLT x )      { return std::log( x ) / std::log( 10.0 ); }
void sin_cos( FLT x, FLT& si, FLT& co ) { si = std::sin( x ); co = std::cos( x ); }
void sinh_cosh( FLT x, FLT& si, FLT& co ) { si = std::sinh( x ); co = std::cosh( x ); }
FLT atanh2( FLT y, FLT x ){ return std::atanh( y/x); }
FLT norm( FLT x, FLT y ){ return std::sqrt( x*x + y*y ); }
FLT normh( FLT x, FLT y ){ return std::sqrt( x*x - y*y ); }
void rect_to_polar( FLT x, FLT y, FLT& r, FLT& a ) { r = std::sqrt( x*x + y*y ); a = atan2( y, x ); }
void polar_to_rect( FLT r, FLT a, FLT& x, FLT& y ) { x = r*std::cos( a ); y = r*std::sin( a ); }

int main( int argc, const char * argv[] )
{
    Cordic<FP, INT_W, FRAC_W> cordic;
    std::cout << "tol: " << TOL << "\n";

    for( uint32_t i = 0; i < 2; i++ )
    {
        bool do_reduce = i;

        FLT x = 0.681807431807431031;
        FLT y = 0.810431798013170871;
        FLT w = 0.103301038084310970;
        FLT b = M_E * 1.1;

        do_op3( "x*y + w",          cordic.mad,     mad,            x, y, w, do_reduce );
        do_op2( "x*y",              cordic.mul,     mul,            x, y, do_reduce );
        do_op3( "y/x + w",          cordic.dad,     dad,            y, x, w, do_reduce );
        do_op2( "y/x",              cordic.div,     div,            y, x, do_reduce );
        do_op1( "1/x",              cordic.one_over,one_over,       x   , do_reduce );
        do_op1( "sqrt(x)",          cordic.sqrt,    std::sqrt,      x   , do_reduce );
        do_op1( "one_over_sqrt(x)", cordic.one_over_sqrt, one_over_sqrt, x, do_reduce );
        do_op1( "exp(x)",           cordic.exp,     std::exp,       x   , do_reduce );
        do_op2( "pow(x,y)",         cordic.pow,     std::pow,       b, y, do_reduce );
        do_op1( "pow2(x)",          cordic.pow2,    pow2,           x   , do_reduce );
        do_op1( "pow10(x)",         cordic.pow10,   pow10,          x   , true      );
        do_op1( "log(x)",           cordic.log,     std::log,       x   , do_reduce );
        do_op2( "logb(x,b)",        cordic.logb,    logb,           b, y, do_reduce );
        do_op1( "log2(x)",          cordic.log2,    log2,           x   , do_reduce );
        do_op1( "log10(x)",         cordic.log10,   log10,          x   , do_reduce );
        do_op1( "sin(x)",           cordic.sin,     std::sin,       x   , do_reduce );
        do_op1( "cos(x)",           cordic.cos,     std::cos,       x   , do_reduce );
        do_op12("sin_cos(x)",       cordic.sin_cos, sin_cos,        x   , do_reduce );
        do_op1( "tan(x)",           cordic.tan,     std::tan,       x   , do_reduce );
        do_op1( "asin(x)",          cordic.asin,    std::asin,      x   , do_reduce );
        do_op1( "acos(x)",          cordic.acos,    std::acos,      x   , do_reduce );
        do_op1( "atan(x)",          cordic.atan,    std::atan,      x   , false );
        do_op2( "atan2(y,x)",       cordic.atan2,   std::atan2,     y, x, false );
        do_op1( "sinh(x)",          cordic.sinh,    std::sinh,      x   , do_reduce );
        do_op1( "cosh(x)",          cordic.cosh,    std::cosh,      x   , do_reduce );
        do_op12("sinh_cosh(x)",     cordic.sinh_cosh,sinh_cosh,     x   , do_reduce );
        do_op1( "tanh(x)",          cordic.tanh,    std::tanh,      x   , do_reduce );
        do_op1( "asinh(x)",         cordic.asinh,   std::asinh,     x   , do_reduce );
        do_op1( "acosh(x)",         cordic.acosh,   std::acosh,     1.591370341781322, do_reduce );
        do_op1( "atanh(x)",         cordic.atanh,   std::atanh,     x   , false );
        do_op2( "atanh2(y,x)",      cordic.atanh2,  atanh2,         0.456728943106177373, 0.709831990704326039, false );
        do_op2( "norm(x,y)",        cordic.norm,    norm,           x, y, do_reduce );
        do_op2( "normh(x,y)",       cordic.normh,   normh,          0.708473170947310947, 0.556728943106177373, do_reduce );

        do_op22("rect_to_polar(x,y)", cordic.rect_to_polar,rect_to_polar, x, y, do_reduce );
        do_op22("polar_to_rect(x,y)", cordic.polar_to_rect,polar_to_rect, x, y, false );

        do_op2( "x*y",              cordic.mul,     mul,            0.0001, 1.999999, do_reduce );
        do_op2( "x/y",              cordic.div,     div,            0.0003, 1.999999, do_reduce );
        do_op2( "x/y",              cordic.div,     div,            0.0003, 0.000555, do_reduce );
        do_op1( "sqrt(1.99999)",    cordic.sqrt,    sqrt,           1.99999, do_reduce );
        do_op1( "exp(0.5)",         cordic.exp,     exp,            0.5, do_reduce );
        do_op1( "exp(1)",           cordic.exp,     exp,            1.0, do_reduce );
        do_op1( "log(2.71)",        cordic.log,     log,            2.71, do_reduce );
        do_op1( "log(1.00)",        cordic.log,     log,            1.00, do_reduce );
        do_op1( "log(0.50)",        cordic.log,     log,            0.50, do_reduce );
    }

    return 0;
}
