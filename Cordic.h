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

// T      = some signed integer type that can hold fixed-point values (default is int64_t)
// FLT    = some floating-point type that can hold constants of the desired precision (default is double)
//
template< typename T=int64_t, typename FLT=double >              
class Cordic
{
public:
    //-----------------------------------------------------
    // Constructor
    //
    // int_w  = fixed-point integer width
    // frac_w = fixed-point fraction width
    // 1+int_w+frac_w must fit in T
    //
    // So the most significant bit is the sign, followed by int_w bits of integer, followed by frac_w bits of fraction.
    //
    // do_reduce indicates whether Cordic routines should first do argument reduction, which is usually the case.
    // If you need a mix of reduce and no-reduce, then please allocate two different Cordic objects with different do_reduce
    // settings but the same other parameters.
    //
    // nc == number of iterations for circular   (0 == use FRAC_W)
    // nh == number of iterations for hyperbolic (0 == use FRAC_W)
    // nl == number of iterations for linear     (0 == use FRAC_W)
    //-----------------------------------------------------
    Cordic( uint32_t int_w, uint32_t frac_w, bool do_reduce=true, uint32_t nc=0, uint32_t nh=0, uint32_t nl=0 );
    ~Cordic();

    //-----------------------------------------------------
    // Conversions
    //-----------------------------------------------------
    T       to_fp( FLT x ) const;               // float to fixed-point
    FLT     to_flt( const T& x ) const;         // fixed-point to float

    T       make_fp( bool sign, T i, T f );     // encode a fixed-point value using sign, integer part i, and fractional part f

    //-----------------------------------------------------
    // Well-Known Math Functions Implemented Using CORDIC
    //
    // (2) means requires 2 applications of a CORDIC algorithm.              functionality
    //-----------------------------------------------------               ---------------------------
    T    mad( const T& x, const T& y, const T addend ) const;             // x*y + addend
    T    mul( const T& x, const T& y ) const;                             // x*y 
    T    dad( const T& y, const T& x, const T addend = T(0) ) const;      // y/x + addend
    T    div( const T& y, const T& x ) const;                             // y/x
    T    one_over( const T& x ) const;                                    // 1/x
    T    sqrt( const T& x ) const;                                        // normh( x+0.25, x-0.25 )
    T    one_over_sqrt( const T& x ) const;                               // 1/sqrt 

    T    exp( const T& x ) const;                                         // e^x
    T    pow( const T& b, const T& x ) const;                             // exp(x * log(b))              (3)
    T    powc( const FLT& b, const T& x ) const;                          // exp(x * log(b))  b=const     (2)
    T    pow2( const T& x ) const;                                        // exp(x * log(2))              (2)
    T    pow10( const T& x ) const;                                       // exp(x * log(10))             (2)
    T    log( const T& x ) const;                                         // 2*atan2(x-1, x+1)    
    T    logb( const T& x, const T& b ) const;                            // log(x)/log(b)                (3)
    T    logc( const T& x, const FLT& b ) const;                          // log(x)/log(b)    b=const     (2)
    T    log2( const T& x ) const;                                        // log(x)/log(2)                (2)
    T    log10( const T& x ) const;                                       // log(x)/log(10)               (2)

    T    sin( const T& x, const T * r=nullptr  ) const;                   // r*sin(x)                   (default r is 1)
    T    cos( const T& x, const T * r=nullptr ) const;                    // r*cos(x)                   (default r is 1)
    void sin_cos( const T& x, T& si, T& co, const T * r=nullptr ) const;  // si=r*sin(x), co=r*cos(x)   (default r is 1)
    T    tan( const T& x ) const;                                         // sin(x) / cos(x)              (2)
    T    asin( const T& x ) const;                                        // atan2(x, sqrt(1 - x^2))      (2)
    T    acos( const T& x ) const;                                        // atan2(sqrt(1 - x^2), x)      (2)
    T    atan( const T& x ) const;                                        // atan(x)
    T    atan2( const T& y, const T& x ) const;                           // atan2(y, x)                  

    void polar_to_rect( const T& r, const T& a, T& x, T& y ) const;       // sin_cos(a, x, y, &r)  
    void rect_to_polar( const T& x, const T& y, T& r, T& a ) const;       // r=sqrt(x^2 + y^2), a=atan2(y, x)
    T    norm( const T& x, const T& y ) const;                            // sqrt(x^2 + y^2)
    T    normh( const T& x, const T& y ) const;                           // sqrt(x^2 - y^2)

    T    sinh( const T& x, const T * r=nullptr ) const;                   // r*sinh(x), also r*(e^x - e^-x)/2  (default r is 1)
    T    cosh( const T& x, const T * r=nullptr ) const;                   // r*cosh(x), also r*(e^x + e^-x)/2  (default r is 1)
    void sinh_cosh( const T& x, T& sih, T& coh, const T * r=nullptr ) const;// sih=r*sinh(x), coh=r*cosh(x)    (default r is 1)
    T    tanh( const T& x ) const;                                        // sinh(x) / cosh(x)            (2)
    T    asinh( const T& x ) const;                                       // log(x + sqrt(x^2 + 1))       (2)
    T    acosh( const T& x ) const;                                       // log(x + sqrt(x^2 - 1))       (2)
    T    atanh( const T& x ) const;                                       // atanh(x)
    T    atanh2( const T& y, const T& x ) const;                          // atanh2(y, x)

    //-----------------------------------------------------
    // Useful Identities
    //
    // exp(x+y)         = exp(x) * exp(y)
    // log(x)           = 2*atan2(x-1, x+1)              
    // log(x)           = atanh2(x^2 - 1, x^2 + 1) 
    // log(x*y)         = log(x) + log(y)
    // log(x/y)         = log(x) - log(y)
    // sqrt(x)          = sqrt( (x+0.25)^2 - (x-0.25)^2 )               (allows use of normh)
    //
    // sin(-x)          = -sin(x)
    // sin(x+y)         = sin(x)*cos(y) + cos(x)*sin(y)
    // cos(x+y)         = cos(x)*sin(y) - sin(x)*cos(y)
    // cos(-x)          = cos(x)
    // tan(-x)          = -tan(x)
    // tan(x+y)         = (tan(x) + tan(y)) / (1 - tan(x)*tan(y))
    //
    // asin(-x)         = -asin(x)
    // asin(x)          = atan2(x, sqrt(1 - x^2))
    // asin(x+y)        = 
    // acos(-x)         = acos(-x)
    // acos(x)          = atan2(sqrt(1 - x^2), x)
    // acos(x+y)        = PI/2 - asin(x+y)
    // atan(-x)         = -atan(x)
    // atan(x+y)        = atan((x+y) / (1 - x*y))
    // atan(x*y)        = asin(x + y)
    //
    // sinh(-x)         = -sinh(x)
    // sinh(x)          = (e^x - e^-x)/2
    // sinh(x+y)        = sinh(x)*cosh(y) + cosh(x)*sinh(y)
    // cosh(-x)         = cosh(x)
    // cosh(x)          = (e^x + e^-x)/2
    // cosh(x+y)        = cosh(x)*cosh(y) - sinh(x)*sinh(y)
    // tanh(-x)         = -tanh(x)
    // tanh(x)          = (e^x - e^-x) / (e^x + e^-x)
    //
    // asinh(-x)        = -asinh(x)
    // asinh(x)         = log(x + sqrt(x^2 + 1))
    // asinh(x)         = atanh(x / sqrt(1 + x^2))
    // acosh(-x)        = [illegal, x must be >= 1]
    // acosh(x)         = log(x + sqrt(x^2 - 1))
    // acosh(x)         = abs( asinh(sqrt(x^2 - 1)) )
    // atanh(-x)        = -atanh(x)
    // atanh(x)         = log((1+x)/(1-x))/2 = log(1+x)/2 - log(1-x)/2    (note: x must be between -1 and 1)
    // atanh(x)         = asinh(x / sqrt(1 - x^2)) 
    // atanh(x)         = +/- acosh(1 / sqrt(1 - x^2))
    //-----------------------------------------------------

    //-----------------------------------------------------
    //-----------------------------------------------------
    //-----------------------------------------------------
    //-----------------------------------------------------
    //-----------------------------------------------------
    //
    // WHAT FOLLOWS ARE ADVANCED ROUTINES THAT MOST PEOPLE WON'T EVER NEED
    // 
    //-----------------------------------------------------
    //-----------------------------------------------------
    //-----------------------------------------------------
    //-----------------------------------------------------
    //-----------------------------------------------------

    //-----------------------------------------------------
    // Constants 
    //-----------------------------------------------------
    uint32_t int_w( void ) const;
    uint32_t frac_w( void ) const;
    T maxint( void ) const;             // largest positive integer
    T zero( void ) const;               // 0.0
    T one( void ) const;                // 1.0
    T quarter( void ) const;            // 0.25
    T gain( void ) const;               // circular
    T gainh( void ) const;              // hyperbolic
    T one_over_gain( void ) const;      // circular
    T one_over_gainh( void ) const;     // hyperbolic

    uint32_t n_circular( void ) const;  // nc
    uint32_t n_linear( void ) const;    // nh
    uint32_t n_hyperbolic( void ) const;// nl

    //-----------------------------------------------------
    // The basic CORDIC functions that all the above math functions ultimately use.
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
    // These version are used internally, but making them available publically.
    // In general, you should only call the earlier routines.
    //-----------------------------------------------------
    T    mad( const T& x, const T& y, const T addend, bool do_reduce ) const; // same but override do_reduce
    T    mul( const T& x, const T& y, bool do_reduce ) const;             // same but override do_reduce
    T    dad( const T& y, const T& x, const T addend, bool do_reduce ) const; // same but override do_reduce
    T    div( const T& y, const T& x, bool do_reduce ) const;             // same but override do_reduce
    T    log( const T& x, bool do_reduce ) const;                         // 2*atan2(x-1, x+1, do_reduce)    
    T    atan2(  const T& y, const T& x, bool do_reduce, bool x_is_one=false, T * r=nullptr ) const; // same but override do_reduce 
    T    atanh2( const T& y, const T& x, bool do_reduce, bool x_is_one=false ) const; // same but override do_reduce

    //-----------------------------------------------------
    // Argument Range Reduction Routines
    //
    // If your Cordic object was allocated with do_reduce==false, then
    // you are still free to call these routines to perform reduction where you want.
    // However, it's easier to just allocate another Cordic object with do_reduce=true
    // and use the appropriate one depending on whether you want reduction or not.
    //-----------------------------------------------------
    void reduce_arg( T& x, int32_t& x_lshift, bool& sign, bool shift_x=true, bool normalize=false ) const; // 0.0 .. <2.0
    void reduce_mul_args( T& x, T& y, int32_t& x_lshift, int32_t& y_lshift, bool& sign ) const; // reduce_arg x and y
    void reduce_div_args( T& x, T& y, int32_t& x_lshift, int32_t& y_lshift, bool& sign ) const; // reduce_arg x, normalize_arg y
    void reduce_sqrt_arg( T& x, int32_t& x_lshift ) const;                                      // reduce_arg but lshift must pow-of-2
    void reduce_exp_arg( FLT b, T& x, T& factor, bool& sign ) const;                            // b=const_base, mul/div exp(x) by factor
    void reduce_log_arg( T& x, T& addend ) const;                                               // reduce_arg x, addend to log(x) 
    void reduce_norm_args( T& x, T& y, int32_t& lshift ) const;                                 // reduce_arg x and y with same lshift
    void reduce_angle_arg( T& a, uint32_t& quadrant, bool& sign ) const;                        // to 0 .. pi/2

private:
    class Impl;
    std::unique_ptr<Impl> impl;
};

#endif
