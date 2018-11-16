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
#ifndef _Cordic_h
#define _Cordic_h

#include "Misc.h"

// T      = some signed integer type that can hold fixed-point values (e.g., int64_t)
// INT_W  = integer width to left of fixed decimal point (not including sign)
// FRAC_W = fraction width to right of fixed decimal point
// FLT    = some floating-point type that can hold constants of the desired precision (default is double)
//
template< typename T, int INT_W, int FRAC_W, typename FLT=double >              
class Cordic
{
public:
    //-----------------------------------------------------
    // Constructor
    //
    // nc == number of iterations for circular
    // nh == number of iterations for hyperbolic
    // nl == number of iterations for linear
    //-----------------------------------------------------
    Cordic( uint32_t nc=FRAC_W, uint32_t nh=FRAC_W, uint32_t nl=FRAC_W );
    ~Cordic();

    //-----------------------------------------------------
    // Queries
    //-----------------------------------------------------
    const T ZERO    = 0;
    const T ONE     = T(1) << T(FRAC_W);
    const T QUARTER = T(1) << T(FRAC_W-2);
    const uint32_t MAX_INT = (1 << INT_W)-1;

    T       to_fp( FLT x ) const;
    FLT     to_flt( const T& x ) const;

    uint32_t n_circular( void ) const;
    uint32_t n_linear( void ) const;
    uint32_t n_hyperbolic( void ) const;

    T gain( void ) const;               // circular
    T gainh( void ) const;              // hyperbolic
    T one_over_gain( void ) const;      // circular
    T one_over_gainh( void ) const;     // hyperbolic

    //-----------------------------------------------------
    // The CORDIC functions
    //-----------------------------------------------------

    // circular rotation mode after step n:
    //      x = gain*(x0*cos(z0) - y0*sin(z0))
    //      y = gain*(y0*cos(z0) + x0*sin(z0))
    //      z = 0
    //
    void circular_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const;

    // circular vectoring mode after step n:
    //      x = gain*sqrt(x0^2 + y0^2)
    //      y = 0
    //      z = z0 + atan( y0/x0 )
    //
    void circular_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const;

    // hyperbolic rotation mode after step n:
    //      x = gain*(x0*cosh(z0) + y0*sinh(z0))
    //      y = gain*(y0*cosh(z0) + x0*sinh(z0))
    //      z = 0
    //
    void hyperbolic_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const;

    // hyperbolic vectoring mode after step n:
    //      x = gain*sqrt(x0^2 - y0^2)
    //      y = 0
    //      z = z0 + atanh( y0/x0 )
    //
    void hyperbolic_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const;

    // linear rotation mode after step n:
    //      x = x0
    //      y = y0 + x0*z0
    //      z = 0
    //
    void linear_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const;

    // linear vectoring mode after step n:
    //      x = x0
    //      y = 0
    //      z = z0 + y0/x0
    //
    void linear_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const;

    //-----------------------------------------------------
    // Well-Known Math Functions Implemented Using CORDIC
    //
    // (2) means requires 2 applications of a CORDIC algorithm.
    //
    // do_reduce=true means that arguments need to be range-reduced.
    //-----------------------------------------------------

    T    mul( const T& x, const T& y, const T addend = T(0), bool do_reduce=false ) const;      // x*y + addend
    T    div( const T& y, const T& x, const T addend = T(0), bool do_reduce=false ) const;      // y/x + addend
    T    sqrt( const T& x, bool do_reduce=false ) const;                                        // normh( x+0.25, x-0.25 )
    T    exp( const T& x, bool do_reduce=false ) const;                                         // e^x
    T    pow( const T& b, const T& x, bool do_reduce=false ) const;                             // exp(x * log(b))              (3)
    T    powc( const FLT& b, const T& x, bool do_reduce=false ) const;                          // exp(x * log(b))  b=const     (2)
    T    pow2( const T& x, bool do_reduce=false ) const;                                        // exp(x * log(2))              (2)
    T    pow10( const T& x, bool do_reduce=false ) const;                                       // exp(x * log(10))             (2)
    T    log( const T& x, bool do_reduce=false ) const;                                         // 2*atan2(x-1, x+1)    
    T    logb( const T& x, const T& b, bool do_reduce=false ) const;                            // log(x)/log(b)                (3)
    T    logc( const T& x, const FLT& b, bool do_reduce=false ) const;                          // log(x)/log(b)    b=const     (2)
    T    log2( const T& x, bool do_reduce=false ) const;                                        // log(x)/log(2)                (2)
    T    log10( const T& x, bool do_reduce=false ) const;                                       // log(x)/log(10)               (2)

    T    sin( const T& x, bool do_reduce=false ) const;                                         // sin(x)
    T    cos( const T& x, bool do_reduce=false ) const;                                         // cos(x)
    void sin_cos( const T& x, T& si, T& co, bool do_reduce=false ) const;                       // si=sin(x), co=cos(x)
    T    tan( const T& x, bool do_reduce=false ) const;                                         // sin(x) / cos(x)              (2)
    T    asin( const T& x, bool do_reduce=false ) const;                                        // atan2(x, sqrt(1 - x^2))      (2)
    T    acos( const T& x, bool do_reduce=false ) const;                                        // atan2(sqrt(1 - x^2), x)      (2)
    T    atan( const T& x, bool do_reduce=false ) const;                                        // atan(x)
    T    atan2( const T& y, const T& x, bool do_reduce=false ) const;                           // atan2(y, x)                  

    void polar_to_rect( const T& r, const T& a, T& x, T& y, bool do_reduce=false ) const;       // x=r*cos(a), y=r*sin(a)
    void rect_to_polar( const T& x, const T& y, T& r, T& a, bool do_reduce=false ) const;       // r=sqrt(x^2 + y^2), a=atan2(y, x)
    T    norm( const T& x, const T& y, bool do_reduce=false ) const;                            // sqrt(x^2 + y^2)
    T    normh( const T& x, const T& y, bool do_reduce=false ) const;                           // sqrt(x^2 - y^2)

    T    sinh( const T& x, bool do_reduce=false ) const;                                        // sinh(x)
    T    cosh( const T& x, bool do_reduce=false ) const;                                        // cosh(x)
    void sinh_cosh( const T& x, T& sih, T& coh, bool do_reduce=false ) const;                   // sih=sinh(x), coh=cosh(x)
    T    tanh( const T& x, bool do_reduce=false ) const;                                        // sinh(x) / cosh(x)            (2)
    T    asinh( const T& x, bool do_reduce=false ) const;                                       // log(x + sqrt(1 + x^2))       (2)
    T    acosh( const T& x, bool do_reduce=false ) const;                                       // log(x + sqrt(x^2 - 1))       (2)
    T    atanh( const T& x, bool do_reduce=false ) const;                                       // atanh(x)
    T    atanh2( const T& y, const T& x, bool do_reduce=false ) const;                          // atanh2(y, x)

    //-----------------------------------------------------
    // Argument Range Reduction Routines
    //
    // If you take the default above for do_reduce=true, then
    // you need not call these yourself.
    //
    // All inputs must be non-negative.
    //-----------------------------------------------------
    void reduce_arg( T& x, int32_t& x_lshift, bool shift_x=true, bool normalize=false ) const;  // 0.0 .. <2.0
    void reduce_mul_args( T& x, T& y, int32_t& x_lshift, int32_t& y_lshift ) const;             // reduce_arg x and y
    void reduce_div_args( T& x, T& y, int32_t& x_lshift, int32_t& y_lshift ) const;             // reduce_arg x, normalize_arg y
    void reduce_sqrt_arg( T& x, int32_t& x_lshift ) const;                                      // reduce_arg but lshift must pow-of-2
    void reduce_exp_arg( FLT b, T& x, T& factor ) const;                                        // b=const_base, multiply exp(x) by factor
    void reduce_log_arg( T& x, T& addend ) const;                                               // reduce_arg x, addend to log(x) 
    void reduce_norm_args( T& x, T& y, int32_t& lshift ) const;                                 // reduce_arg x and y with same lshift
    void reduce_angle( T& a, uint32_t& quadrant ) const;                                        // to 0 .. pi/2

private:
    class Impl;
    std::unique_ptr<Impl> impl;
};

#endif
