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

#include "Cordic.h"

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

#define do_op1( str, c_fn, exp_fn, fltx, do_reduce )                    \
{                                                                       \
    auto c = do_reduce ? cordicr : cordicnr;                            \
    T   tx  = c->to_t( fltx );			                        \
    T   tz  = c->c_fn( tx );		                                \
    FLT fltz = c->to_flt( tz );	                                        \
    FLT flte = exp_fn( fltx );			                        \
    T   te   = c->to_t( flte );                                         \
    T   terr = (tz >= te) ? (tz-te) : (te-tz);                          \
    FLT flterr = c->to_flt( terr );                                     \
    int32_t tol_lg2;                                                    \
    FLT tol  = tolerance( c->frac_w(), flte, TOL, tol_lg2 );            \
    int32_t tol_bits = c->frac_w() + tol_lg2;                           \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Reduce:   " << do_reduce << "\n";                     \
    std::cout << "Input:    " << std::setw(30) << fltx   << " (x)\n";	\
    std::cout << "Tol:      " << std::setw(30) << tol    << " (" << tol_bits              << " bits)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte   << " (" << c->to_bstring( te )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz   << " (" << c->to_bstring( tz )   << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr << " (" << c->to_bstring( terr ) << ")\n\n"; \
    cassert( flterr <= tol, "outside tolerance" );			\
}    

#define do_op12( str, c_fn, exp_fn, fltx, do_reduce )                   \
{                                                                       \
    auto c = do_reduce ? cordicr : cordicnr;                            \
    T   tx  = c->to_t( fltx );			                        \
    T   tz1, tz2;                                                       \
    c->c_fn( tx, tz1, tz2 );		                                \
    FLT fltz1 = c->to_flt( tz1 );		                        \
    FLT fltz2 = c->to_flt( tz2 );		                        \
    FLT flte1, flte2;                                                   \
    exp_fn( fltx, flte1, flte2 );			                \
    T   te1  = c->to_t( flte1 );                                        \
    T   te2  = c->to_t( flte2 );                                        \
    T   terr1 = (tz1 >= te1) ? (tz1-te1) : (te1-tz1);                   \
    T   terr2 = (tz2 >= te2) ? (tz2-te2) : (te2-tz2);                   \
    FLT flterr1 = c->to_flt( terr1 );                                   \
    FLT flterr2 = c->to_flt( terr2 );                                   \
    int32_t tol1_lg2;                                                   \
    int32_t tol2_lg2;                                                   \
    FLT tol1  = tolerance( c->frac_w(), flte1, TOL, tol1_lg2 );         \
    FLT tol2  = tolerance( c->frac_w(), flte2, TOL, tol2_lg2 );         \
    int32_t tol1_bits = c->frac_w() + tol1_lg2;                         \
    int32_t tol2_bits = c->frac_w() + tol2_lg2;                         \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Reduce:   " << do_reduce << "\n";                     \
    std::cout << "Input:    " << std::setw(30) << fltx << " (x)\n";	\
    std::cout << "Tol:      " << std::setw(30) << tol1    << " (" << tol1_bits              << " bits)\n"; \
    std::cout << "Tol:      " << std::setw(30) << tol2    << " (" << tol2_bits              << " bits)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte1   << " (" << c->to_bstring( te1 )   << ")\n"; \
    std::cout << "Expected: " << std::setw(30) << flte2   << " (" << c->to_bstring( te2 )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz1   << " (" << c->to_bstring( tz1 )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz2   << " (" << c->to_bstring( tz2 )   << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr1 << " (" << c->to_bstring( terr1 ) << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr2 << " (" << c->to_bstring( terr2 ) << ")\n\n"; \
    cassert( flterr1 <= tol1, "outside tolerance" );			\
    cassert( flterr2 <= tol2, "outside tolerance" );			\
}    

#define do_op2( str, c_fn, exp_fn, fltx, flty, do_reduce )              \
{                                                                       \
    auto c = do_reduce ? cordicr : cordicnr;                            \
    T   tx  = c->to_t( fltx );			                        \
    T   ty  = c->to_t( flty );			                        \
    T   tz  = c->c_fn( tx, ty );	                                \
    FLT fltz = c->to_flt( tz );			                        \
    FLT flte = exp_fn( fltx, flty );                                    \
    T   te   = c->to_t( flte );                                         \
    T   terr = (tz >= te) ? (tz-te) : (te-tz);                          \
    FLT flterr = c->to_flt( terr );                                     \
    int32_t tol_lg2;                                                    \
    FLT tol  = tolerance( c->frac_w(), flte, TOL, tol_lg2 );            \
    int32_t tol_bits = c->frac_w() + tol_lg2;                           \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Reduce:   " << do_reduce << "\n";                     \
    std::cout << "Input:    " << std::setw(30) << fltx << " (x)\n";     \
    std::cout << "Input:    " << std::setw(30) << flty << " (y)\n";     \
    std::cout << "Tol:      " << std::setw(30) << tol    << " (" << tol_bits              << " bits)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte   << " (" << c->to_bstring( te )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz   << " (" << c->to_bstring( tz )   << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr << " (" << c->to_bstring( terr ) << ")\n\n"; \
    cassert( flterr <= tol, "outside tolerance" );			\
}    

#define do_op22( str, c_fn, exp_fn, fltx, flty, do_reduce )             \
{                                                                       \
    auto c = do_reduce ? cordicr : cordicnr;                            \
    T   tx  = c->to_t( fltx );			                        \
    T   ty  = c->to_t( flty );			                        \
    T   tz1, tz2;                                                       \
    c->c_fn( tx, ty, tz1, tz2 );	                                \
    FLT fltz1 = c->to_flt( tz1 );		                        \
    FLT fltz2 = c->to_flt( tz2 );		                        \
    FLT flte1, flte2;                                                   \
    exp_fn( fltx, flty, flte1, flte2 );			                \
    T   te1  = c->to_t( flte1 );                                        \
    T   te2  = c->to_t( flte2 );                                        \
    T   terr1 = (tz1 >= te1) ? (tz1-te1) : (te1-tz1);                   \
    T   terr2 = (tz2 >= te2) ? (tz2-te2) : (te2-tz2);                   \
    FLT flterr1 = c->to_flt( terr1 );                                   \
    FLT flterr2 = c->to_flt( terr2 );                                   \
    int32_t tol1_lg2;                                                   \
    int32_t tol2_lg2;                                                   \
    FLT tol1  = tolerance( c->frac_w(), flte1, TOL, tol1_lg2 );         \
    FLT tol2  = tolerance( c->frac_w(), flte2, TOL, tol2_lg2 );         \
    int32_t tol1_bits = c->frac_w() + tol1_lg2;                         \
    int32_t tol2_bits = c->frac_w() + tol2_lg2;                         \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Reduce:   " << do_reduce << "\n";                     \
    std::cout << "Input:    " << std::setw(30) << fltx << " (x)\n";     \
    std::cout << "Input:    " << std::setw(30) << flty << " (y)\n";     \
    std::cout << "Tol:      " << std::setw(30) << tol1    << " (" << tol1_bits              << " bits)\n"; \
    std::cout << "Tol:      " << std::setw(30) << tol2    << " (" << tol2_bits              << " bits)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte1   << " (" << c->to_bstring( te1 )   << ")\n"; \
    std::cout << "Expected: " << std::setw(30) << flte2   << " (" << c->to_bstring( te2 )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz1   << " (" << c->to_bstring( tz1 )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz2   << " (" << c->to_bstring( tz2 )   << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr1 << " (" << c->to_bstring( terr1 ) << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr2 << " (" << c->to_bstring( terr2 ) << ")\n\n"; \
    cassert( flterr1 <= tol1, "outside tolerance" );			\
    cassert( flterr2 <= tol2, "outside tolerance" );			\
}    

#define do_op3( str, c_fn, exp_fn, fltx, flty, fltw, do_reduce )        \
{                                                                       \
    auto c = do_reduce ? cordicr : cordicnr;                            \
    T   tx  = c->to_t( fltx );			                        \
    T   ty  = c->to_t( flty );			                        \
    T   tw  = c->to_t( fltw );			                        \
    T   tz  = c->c_fn( tx, ty, tw );	                                \
    FLT fltz = c->to_flt( tz );			                        \
    FLT flte = exp_fn( fltx, flty, fltw );			        \
    T   te   = c->to_t( flte );                                         \
    T   terr = (tz >= te) ? (tz-te) : (te-tz);                          \
    FLT flterr = c->to_flt( terr );                                     \
    int32_t tol_lg2;                                                    \
    FLT tol  = tolerance( c->frac_w(), flte, TOL, tol_lg2 );            \
    int32_t tol_bits = c->frac_w() + tol_lg2;                           \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Reduce:   " << do_reduce << "\n";                     \
    std::cout << "Input:    " << std::setw(30) << fltx << " (x)\n";     \
    std::cout << "Input:    " << std::setw(30) << flty << " (y)\n";     \
    std::cout << "Input:    " << std::setw(30) << fltz << " (z)\n";     \
    std::cout << "Tol:      " << std::setw(30) << tol    << " (" << tol_bits              << " bits)\n"; \
    std::cout << "Expected: " << std::setw(30) << flte   << " (" << c->to_bstring( te )   << ")\n"; \
    std::cout << "Actual:   " << std::setw(30) << fltz   << " (" << c->to_bstring( tz )   << ")\n"; \
    std::cout << "Diff:     " << std::setw(30) << flterr << " (" << c->to_bstring( terr ) << ")\n\n"; \
    cassert( flterr <= tol, "outside tolerance" );			\
}    

// FLT wrapper routines for those that are not in std::
//
FLT  mad( FLT x, FLT y, FLT w ) { return x*y + w; }
FLT  mul( FLT x, FLT y ) { return x*y; }
FLT  dad( FLT x, FLT y, FLT w ) { return x/y + w; }
FLT  div( FLT x, FLT y ) { return x/y; }
FLT  rcp( FLT x )        { return 1.0/x; }
FLT  rsqrt( FLT x )      { return 1.0 / std::sqrt( x ); }
FLT  pow2( FLT x )       { return std::pow( 2.0, x ); }
FLT  exp10( FLT x )      { return std::pow( 10.0, x ); }
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

#endif
