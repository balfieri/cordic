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
// test_helpers.h - helper macros and wrapper functions used by all tests
//
#ifndef _test_helpers_h
#define _test_helpers_h

#include "freal.h"
#include <cstdint>

// some useful macros to avoid redundant typing
//

static inline FLT tolerance( uint32_t frac_w, FLT expected, FLT tol, int32_t& tol_lg2 )
{
    //------------------------------------------------------------
    // Fixed-point numbers > 1 are going to have less precision.
    // tol is the tolerance if the |result| is <= 1.
    //------------------------------------------------------------
    (void)frac_w; // unused
    FLT exp_abs = std::fabs( expected );
    tol_lg2 = int32_t( std::log2( tol ) - 0.5 );
    if ( exp_abs > 1.0 ) tol_lg2 += int32_t( std::log2( exp_abs ) ) + 1;
    return std::pow( 2.0, tol_lg2 );
}

#define do_op1( str, c_fn, exp_fn, fltx )                               \
{                                                                       \
    std::cout << "\n" << #str << "\n";			                \
    freal tx   = fltx;			                                \
    freal tz   = std::c_fn( tx );		                        \
    FLT   fltz = tz;	                                                \
    FLT   flte = exp_fn( fltx );			                \
    freal te   = std::rfrac( freal( flte ) );                           \
    freal terr = (tz >= te) ? (tz-te) : (te-tz);                        \
    FLT flterr = terr;                                                  \
    int32_t tol_lg2;                                                    \
    FLT tol    = tolerance( te.frac_w(), flte, TOL, tol_lg2 );          \
    int32_t tol_bits = te.frac_w() + tol_lg2;                           \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Input:    " << std::setw(30) << fltx   << " (x)\n";	\
    std::cout << "Tol:      " << std::setw(30) << tol    << " (" << tol_bits                << " bits)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte   << " (" << std::to_bstring( te )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz   << " (" << std::to_bstring( tz )   << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr << " (" << std::to_bstring( terr ) << ")\n\n"; \
    cassert( (std::isnan(fltz) && std::isnan(flte)) || flterr <= tol, "outside tolerance" );			\
}    

#define do_op12( str, c_fn, exp_fn, fltx )                              \
{                                                                       \
    std::cout << "\n" << #str << "\n";			                \
    freal tx  = fltx;			                                \
    freal tz1, tz2;                                                     \
    std::c_fn( tx, tz1, tz2 );		                                \
    FLT fltz1 = tz1;		                                        \
    FLT fltz2 = tz2;		                                        \
    FLT flte1, flte2;                                                   \
    exp_fn( fltx, flte1, flte2 );			                \
    freal te1  = std::rfrac( freal(flte1) );                            \
    freal te2  = std::rfrac( freal(flte2) );                            \
    freal terr1 = (tz1 >= te1) ? (tz1-te1) : (te1-tz1);                 \
    freal terr2 = (tz2 >= te2) ? (tz2-te2) : (te2-tz2);                 \
    FLT flterr1 = terr1;                                                \
    FLT flterr2 = terr2;                                                \
    int32_t tol1_lg2;                                                   \
    int32_t tol2_lg2;                                                   \
    FLT tol1  = tolerance( te1.frac_w(), flte1, TOL, tol1_lg2 );        \
    FLT tol2  = tolerance( te2.frac_w(), flte2, TOL, tol2_lg2 );        \
    int32_t tol1_bits = te1.frac_w() + tol1_lg2;                        \
    int32_t tol2_bits = te2.frac_w() + tol2_lg2;                        \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Input:    " << std::setw(30) << fltx << " (x)\n";	\
    std::cout << "Tol:      " << std::setw(30) << tol1    << " (" << tol1_bits                << " bits)\n"; \
    std::cout << "Tol:      " << std::setw(30) << tol2    << " (" << tol2_bits                << " bits)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte1   << " (" << std::to_bstring( te1 )   << ")\n"; \
    std::cout << "Expected: " << std::setw(30) << flte2   << " (" << std::to_bstring( te2 )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz1   << " (" << std::to_bstring( tz1 )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz2   << " (" << std::to_bstring( tz2 )   << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr1 << " (" << std::to_bstring( terr1 ) << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr2 << " (" << std::to_bstring( terr2 ) << ")\n\n"; \
    cassert( flterr1 <= tol1, "outside tolerance" );			\
    cassert( flterr2 <= tol2, "outside tolerance" );			\
}    

#define do_op2( str, c_fn, exp_fn, fltx, flty )                         \
{                                                                       \
    std::cout << "\n" << #str << "\n";			                \
    freal tx  = fltx;			                                \
    freal ty  = flty;			                                \
    freal tz  = std::c_fn( tx, ty );	                                \
    FLT fltz = tz;			                                \
    FLT flte = exp_fn( fltx, flty );                                    \
    freal te   = std::rfrac( freal(flte) );                             \
    freal terr = (tz >= te) ? (tz-te) : (te-tz);                        \
    FLT flterr = terr;                                                  \
    int32_t tol_lg2;                                                    \
    FLT tol  = tolerance( te.frac_w(), flte, TOL, tol_lg2 );            \
    int32_t tol_bits = te.frac_w() + tol_lg2;                           \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Input:    " << std::setw(30) << fltx << " (x)\n";     \
    std::cout << "Input:    " << std::setw(30) << flty << " (y)\n";     \
    std::cout << "Tol:      " << std::setw(30) << tol    << " (" << tol_bits                << " bits)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte   << " (" << std::to_bstring( te )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz   << " (" << std::to_bstring( tz )   << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr << " (" << std::to_bstring( terr ) << ")\n\n"; \
    cassert( (std::isnan(fltz) && std::isnan(flte)) || flterr <= tol, "outside tolerance" );			\
}    

#define do_op22( str, c_fn, exp_fn, fltx, flty )                        \
{                                                                       \
    std::cout << "\n" << #str << "\n";			                \
    freal tx  = fltx;			                                \
    freal ty  = flty;			                                \
    freal tz1, tz2;                                                     \
    std::c_fn( tx, ty, tz1, tz2 );	                                \
    FLT fltz1 = tz1;		                                        \
    FLT fltz2 = tz2;		                                        \
    FLT flte1, flte2;                                                   \
    exp_fn( fltx, flty, flte1, flte2 );			                \
    freal te1  = std::rfrac( freal(flte1) );                            \
    freal te2  = std::rfrac( freal(flte2) );                            \
    freal terr1 = (tz1 >= te1) ? (tz1-te1) : (te1-tz1);                 \
    freal terr2 = (tz2 >= te2) ? (tz2-te2) : (te2-tz2);                 \
    FLT flterr1 = terr1;                                                \
    FLT flterr2 = terr2;                                                \
    int32_t tol1_lg2;                                                   \
    int32_t tol2_lg2;                                                   \
    FLT tol1  = tolerance( te1.frac_w(), flte1, TOL, tol1_lg2 );        \
    FLT tol2  = tolerance( te2.frac_w(), flte2, TOL, tol2_lg2 );        \
    int32_t tol1_bits = te1.frac_w() + tol1_lg2;                        \
    int32_t tol2_bits = te2.frac_w() + tol2_lg2;                        \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Input:    " << std::setw(30) << fltx << " (x)\n";     \
    std::cout << "Input:    " << std::setw(30) << flty << " (y)\n";     \
    std::cout << "Tol:      " << std::setw(30) << tol1    << " (" << tol1_bits                << " bits)\n"; \
    std::cout << "Tol:      " << std::setw(30) << tol2    << " (" << tol2_bits                << " bits)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte1   << " (" << std::to_bstring( te1 )   << ")\n"; \
    std::cout << "Expected: " << std::setw(30) << flte2   << " (" << std::to_bstring( te2 )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz1   << " (" << std::to_bstring( tz1 )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz2   << " (" << std::to_bstring( tz2 )   << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr1 << " (" << std::to_bstring( terr1 ) << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr2 << " (" << std::to_bstring( terr2 ) << ")\n\n"; \
    cassert( flterr1 <= tol1, "outside tolerance" );			\
    cassert( flterr2 <= tol2, "outside tolerance" );			\
}    

#define do_op22sc( str, c_fn, exp_fn, fltx, flty )                      \
{                                                                       \
    std::cout << "\n" << #str << "\n";			                \
    freal tx  = fltx;			                                \
    freal ty  = flty;			                                \
    freal tz1, tz2;                                                     \
    std::c_fn( tx, tz1, tz2, ty );	                                \
    FLT fltz1 = tz1;		                                        \
    FLT fltz2 = tz2;		                                        \
    FLT flte1, flte2;                                                   \
    exp_fn( fltx, flte1, flte2, flty );			                \
    freal te1  = std::rfrac( freal(flte1) );                            \
    freal te2  = std::rfrac( freal(flte2) );                            \
    freal terr1 = (tz1 >= te1) ? (tz1-te1) : (te1-tz1);                 \
    freal terr2 = (tz2 >= te2) ? (tz2-te2) : (te2-tz2);                 \
    FLT flterr1 = terr1;                                                \
    FLT flterr2 = terr2;                                                \
    int32_t tol1_lg2;                                                   \
    int32_t tol2_lg2;                                                   \
    FLT tol1  = tolerance( te1.frac_w(), flte1, TOL, tol1_lg2 );        \
    FLT tol2  = tolerance( te2.frac_w(), flte2, TOL, tol2_lg2 );        \
    int32_t tol1_bits = te1.frac_w() + tol1_lg2;                        \
    int32_t tol2_bits = te2.frac_w() + tol2_lg2;                        \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Input:    " << std::setw(30) << fltx << " (x)\n";     \
    std::cout << "Input:    " << std::setw(30) << flty << " (y)\n";     \
    std::cout << "Tol:      " << std::setw(30) << tol1    << " (" << tol1_bits                << " bits)\n"; \
    std::cout << "Tol:      " << std::setw(30) << tol2    << " (" << tol2_bits                << " bits)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte1   << " (" << std::to_bstring( te1 )   << ")\n"; \
    std::cout << "Expected: " << std::setw(30) << flte2   << " (" << std::to_bstring( te2 )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz1   << " (" << std::to_bstring( tz1 )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz2   << " (" << std::to_bstring( tz2 )   << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr1 << " (" << std::to_bstring( terr1 ) << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr2 << " (" << std::to_bstring( terr2 ) << ")\n\n"; \
    cassert( flterr1 <= tol1, "outside tolerance" );			\
    cassert( flterr2 <= tol2, "outside tolerance" );			\
}    

#define do_op3( str, c_fn, exp_fn, fltx, flty, fltw )                   \
{                                                                       \
    std::cout << "\n" << #str << "\n";			                \
    freal tx   = fltx;			                                \
    freal ty   = flty;			                                \
    freal tw   = fltw;			                                \
    freal tz   = std::c_fn( tx, ty, tw );	                        \
    FLT fltz   = tz;			                                \
    FLT flte   = exp_fn( fltx, flty, fltw );			        \
    freal te   = std::rfrac( freal(flte) );                             \
    freal terr = (tz >= te) ? (tz-te) : (te-tz);                        \
    FLT flterr = terr;                                                  \
    int32_t tol_lg2;                                                    \
    FLT tol  = tolerance( te.frac_w(), flte, TOL, tol_lg2 );            \
    int32_t tol_bits = te.frac_w() + tol_lg2;                           \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Input:    " << std::setw(30) << fltx << " (x)\n";     \
    std::cout << "Input:    " << std::setw(30) << flty << " (y)\n";     \
    std::cout << "Input:    " << std::setw(30) << fltz << " (z)\n";     \
    std::cout << "Tol:      " << std::setw(30) << tol    << " (" << tol_bits                << " bits)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte   << " (" << std::to_bstring( te )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz   << " (" << std::to_bstring( tz )   << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr << " (" << std::to_bstring( terr ) << ")\n\n"; \
    cassert( (std::isnan(fltz) && std::isnan(flte) || flterr <= tol), "outside tolerance" );			\
}    

// FLT wrapper routines for those that are not in std::
//
FLT  add( FLT x, FLT y ) { return x+y; }
FLT  sub( FLT x, FLT y ) { return x-y; }
FLT  fma( FLT x, FLT y, FLT w ) { return x*y + w; }
FLT  mul( FLT x, FLT y ) { return x*y; }
FLT  fda( FLT x, FLT y, FLT w ) { return x/y + w; }
FLT  div( FLT x, FLT y ) { return x/y; }
FLT  rcp( FLT x )        { return 1.0/x; }
FLT  rsqrt( FLT x )      { return 1.0 / std::sqrt( x ); }
FLT  rcbrt( FLT x )      { return 1.0 / std::cbrt( x ); }
FLT  exp10( FLT x )      { return std::pow( 10.0, x ); }
FLT  log( FLT x, FLT y ) { return std::log( x ) / std::log( y ); }
FLT  log2( FLT x )       { return std::log( x ) / std::log( 2.0 ); }
FLT  log10( FLT x )      { return std::log( x ) / std::log( 10.0 ); }
void sincos( FLT x, FLT& si, FLT& co, FLT r=1.0 )   { si = r*std::sin( x );  co = r*std::cos( x ); }
void sinhcosh( FLT x, FLT& si, FLT& co, FLT r=1.0 ) { si = r*std::sinh( x ); co = r*std::cosh( x ); }
FLT  atanh2( FLT y, FLT x ){ return std::atanh( y/x); }
FLT  hypot( FLT x, FLT y ){ return std::sqrt( x*x + y*y ); }
FLT  hypoth( FLT x, FLT y ){ return std::sqrt( x*x - y*y ); }
void rect_to_polar( FLT x, FLT y, FLT& r, FLT& a ) { r = std::sqrt( x*x + y*y ); a = atan2( y, x ); }
void polar_to_rect( FLT r, FLT a, FLT& x, FLT& y ) { x = r*std::cos( a ); y = r*std::sin( a ); }

#endif
