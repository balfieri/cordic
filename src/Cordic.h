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

// T  = some signed integer type that can hold fixed-point values (e.g., int64_t)
// FW = fraction width of fixed point values                    
//
template< typename T, int FW >              
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
    Cordic( uint32_t nc=FW, uint32_t nh=FW, uint32_t nl=FW );
    ~Cordic();

    //-----------------------------------------------------
    // Queries
    //-----------------------------------------------------
    const T ZERO    = 0;
    const T ONE     = T(1) << T(FW);
    const T QUARTER = T(1) << T(FW-2);

    uint32_t n_circular( void ) const;
    uint32_t n_linear( void ) const;
    uint32_t n_hyperbolic( void ) const;

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
    //-----------------------------------------------------

    T    mul( const T& x, const T& y, const T addend = T(0) ) const;            // x*y + addend
    T    div( const T& y, const T& x, const T addend = T(0) ) const;            // y/x + addend
    T    sqrt( const T& x ) const;                                              // sqrt(x)
    T    exp( const T& x ) const;                                               // e^x
    T    pow( const T& b, const T& x ) const;                                   // exp(x * log(b))              (3)
    T    pow2( const T& x ) const;                                              // exp(x * log(2))              (2)
    T    pow10( const T& x ) const;                                             // exp(x * log(10))             (2)
    T    log( const T& x ) const;                                               // 2*atan2(x-1, x+1)    
    T    logb( const T& x, const T& b ) const;                                  // log(x)/log(b)                (3)
    T    log2( const T& x ) const;                                              // log(x)/log(2)                (2)
    T    log10( const T& x ) const;                                             // log(x)/log(10)               (2)

    T    sin( const T& x ) const;                                               // sin(x)
    T    cos( const T& x ) const;                                               // cos(x)
    void sin_cos( const T& x, T& si, T& co ) const;                             // si=sin(x), co=cos(x)
    T    tan( const T& x ) const;                                               // sin(x) / cos(x)              (2)
    T    asin( const T& x ) const;                                              // atan2(x, sqrt(1 - x^2))      (2)
    T    acos( const T& x ) const;                                              // atan2(sqrt(1 - x^2), x)      (2)
    T    atan( const T& x ) const;                                              // atan(x)
    T    atan2( const T& y, const T& x ) const;                                 // atan2(y, x)                  

    void polar_to_rect( const T& r, const T& a, T& x, T& y ) const;             // x=r*cos(a), y=r*sin(a)
    void rect_to_polar( const T& x, const T& y, T& r, T& a ) const;             // r=sqrt(x^2 + y^2), a=atan2(y, x)
    T    norm( const T& x, const T& y ) const;                                  // sqrt(x^2 + y^2)
    T    normh( const T& x, const T& y ) const;                                 // sqrt(x^2 - y^2)

    T    sinh( const T& x ) const;                                              // sinh(x)
    T    cosh( const T& x ) const;                                              // cosh(x)
    void sinh_cosh( const T& x, T& sih, T& coh ) const;                         // sih=sinh(x), coh=cosh(x)
    T    tanh( const T& x ) const;                                              // sinh(x) / cosh(x)            (2)
    T    asinh( const T& x ) const;                                             // log(x + sqrt(1 + x^2))       (2)
    T    acosh( const T& x ) const;                                             // log(x + sqrt(x^2 - 1))       (2)
    T    atanh( const T& x ) const;                                             // atanh(x)
    T    atanh2( const T& y, const T& x ) const;                                // atanh2(y, x)

private:
    class Impl;
    std::unique_ptr<Impl> impl;
};

#endif
