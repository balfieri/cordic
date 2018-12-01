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
#define do_op1( str, c_fn, exp_fn, fltx, do_reduce )                    \
{                                                                       \
    auto c = do_reduce ? cordicr : cordicnr;                            \
    T   tx  = c->to_t( fltx );			                        \
    T   tz  = c->c_fn( tx );		                                \
    FLT fltz = c->to_flt( tz );	                                        \
    FLT flte = exp_fn( fltx );			                        \
    FLT flterr = std::abs( flte-fltz );			                \
			                                                \
    std::cout.precision(24);			                        \
    std::cout << #str << "\n";			                        \
    std::cout << "Reduce:   " << do_reduce << "\n";                     \
    std::cout << "Input:    " << std::setw(30) << fltx << "(fltx)\n";	\
    std::cout << "Expected: " << std::setw(30) << flte << "\n";		\
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";		\
    std::cout << "Diff:     " << std::setw(30) << flterr << "\n\n";	\
    cassert( flterr <= TOL );			                        \
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
    T   tx  = c->to_t( fltx );			                        \
    T   ty  = c->to_t( flty );			                        \
    T   tz  = c->c_fn( tx, ty );	                                \
    FLT fltz = c->to_flt( tz );			                        \
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
    T   tx  = c->to_t( fltx );			                        \
    T   ty  = c->to_t( flty );			                        \
    T   tz1, tz2;                                                       \
    c->c_fn( tx, ty, tz1, tz2 );	                                \
    FLT fltz1 = c->to_flt( tz1 );		                        \
    FLT fltz2 = c->to_flt( tz2 );		                        \
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
    T   tx  = c->to_t( fltx );			                        \
    T   ty  = c->to_t( flty );			                        \
    T   tw  = c->to_t( fltw );			                        \
    T   tz  = c->c_fn( tx, ty, tw );	                                \
    FLT fltz = c->to_flt( tz );			                        \
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

#endif
