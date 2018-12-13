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

#include <cmath>
#include <iostream>
#include <iomanip>

#include "Logger.h"

#ifdef DEBUG_LEVEL
static constexpr uint32_t debug = DEBUG_LEVEL;
#else
static constexpr uint32_t debug = 0;
#endif

#define cassert(expr, msg) if ( !(expr) ) \
                { std::cout << "ERROR: assertion failure: " << (msg) << " at " << __FILE__ << ":" << __LINE__ << "\n"; exit( 1 ); }


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
    //-----------------------------------------------------
    Cordic( uint32_t int_w,                     // fixed-point integer width
            uint32_t frac_w,                    // fixed-point fraction width
            bool     do_reduce=true,            // whether to do range reduction by default
            uint32_t n=0 );                     // number of iterations used for CORDIC proper (0 == default == frac_w)
    ~Cordic();

    //-----------------------------------------------------
    // Explicit Conversions
    //-----------------------------------------------------
    T       to_t( FLT x ) const;                // FLT to T encoded value
    FLT     to_flt( const T& x, bool can_log=false ) const; // T encoded value to FLT
    std::string to_string( const T& x ) const;  // T to std::string in decimal floating-point format
    std::string to_rstring( const T& _x ) const;// T to std::string in raw decimal integer format 
    std::string to_bstring( const T& x ) const; // T to std::string in binary format, like "1 001 101101011010"

    T       make_fixed( bool sign, T i, T f );  // encode a fixed-point    value using sign, integer  part i, and fractional part f
    T       make_float( bool sign, T e, T f );  // encode a floating-point value using sign, exponent part 3, and fractional part f

    //-----------------------------------------------------
    // Constants 
    //-----------------------------------------------------
    uint32_t int_w( void ) const;                       // int_w  from above
    uint32_t frac_w( void ) const;                      // frac_w from above
    uint32_t n( void ) const;                           // n      from above

    T maxval( void ) const;                             // maximum positive encoded value 
    T minval( void ) const;                             // minimum positive encoded value
    T maxint( void ) const;                             // largest positive integer (just integer part, does not include fraction)
    T zero( void ) const;                               // encoded 0.0
    T one( void ) const;                                // encoded 1.0
    T half( void ) const;                               // encoded 0.5
    T quarter( void ) const;                            // encoded 0.25
    T sqrt2( void ) const;                              // encoded sqrt(2)
    T sqrt2_div_2( void ) const;                        // encoded sqrt(2)/2
    T pi( void ) const;                                 // encoded PI
    T pi_div_2( void ) const;                           // encoded PI/2
    T pi_div_4( void ) const;                           // encoded PI/4
    T two_div_pi( void ) const;                         // encoded 2/PI
    T four_div_pi( void ) const;                        // encoded PI/4
    T e( void ) const;                                  // encoded natural exponent

    //-----------------------------------------------------
    // Well-Known Math Functions Implemented Using CORDIC
    //
    // (2) means requires 2 applications of a CORDIC algorithm.              functionality
    //-----------------------------------------------------               ---------------------------
    void constructed( const T& x ) const;                                 // so we can log creation of x
    void destructed( const T& x ) const;                                  // so we can log destruction of x
    T&   assign( T& x, const T& y ) const;                                // x = y  (this exists so we can log assignments)

    bool isgreater( const T& x, const T& y ) const;                       // x > y
    bool isgreaterequal( const T& x, const T& y ) const;                  // x >= y
    bool isless( const T& x, const T& y ) const;                          // x > y
    bool islessequal( const T& x, const T& y ) const;                     // x >= y
    bool islessgreater( const T& x, const T& y ) const;                   // x != y (but returns false if at least one is NaN)
    bool isunordered( const T& x, const T& y ) const;                     // returns true if x is a NaN OR y is a NaN
    bool isunequal( const T& x, const T& y ) const;                      // x != y (returns true if only one is NaN)
    bool isequal( const T& x, const T& y ) const;                         // x == y (returns true if both are NaN)

    T    abs( const T& x ) const;                                         // |x|
    T    neg( const T& x ) const;                                         // -x
    T    floor( const T& x ) const;                                       // largest  integral value <= x
    T    ceil( const T& x ) const;                                        // smallest integral value >= x

    T    add( const T& x, const T& y ) const;                             // x+y 
    T    sub( const T& x, const T& y ) const;                             // x-y 
    T    mad( const T& x, const T& y, const T& addend ) const;            // x*y + addend
    T    fma( const T& x, const T& y, const T& addend ) const;            // mad( x, y, addend )        (same thing)
    T    mul( const T& x, const T& y ) const;                             // x*y 
    T    lshift( const T& x, int y ) const;                               // x << y
    T    rshift( const T& x, int y ) const;                               // x >> y
    T    dad( const T& y, const T& x, const T addend = T(0) ) const;      // y/x + addend
    T    div( const T& y, const T& x ) const;                             // y/x
    T    rcp( const T& x ) const;                                         // 1/x
    T    sqrt( const T& x ) const;                                        // normh( x+1, x-1 ) / 2
    T    rsqrt( const T& x ) const;                                       // 1/sqrt 
    T    cbrt( const T& x ) const;                                        // x^(1/3)  (but x can be negative)
    T    rcbrt( const T& x ) const;                                       // 1/cbrt

    T    fdim( const T& x, const T& y ) const;                            // (x >= y) ? (x-y) : 0
    T    fmax( const T& x, const T& y ) const;                            // max(x, y)
    T    fmin( const T& x, const T& y ) const;                            // min(x, y)

    T    exp( const T& x ) const;                                         // e^x
    T    expm1( const T& x ) const;                                       // e^x - 1
    T    exp2( const T& x ) const;                                        // 2^x
    T    exp10( const T& x ) const;                                       // 10^x
    T    pow( const T& b, const T& x ) const;                             // b^x  = exp(x * log(b))              (3)
    T    powc( const FLT& b, const T& x ) const;                          // b^x  = exp(x * log(b))  b=const     (2)
    T    log( const T& x ) const;                                         // 2*atan2(x-1, x+1)    
    T    logb( const T& x, const T& b ) const;                            // log(x)/log(b)                (3)
    T    logc( const T& x, const FLT& b ) const;                          // log(x)/log(b)    b=const     (2)
    T    log2( const T& x ) const;                                        // log(x)/log(2)                (2)
    T    log10( const T& x ) const;                                       // log(x)/log(10)               (2)

    T    sin( const T& x, const T * r=nullptr  ) const;                   // r*sin(x)                   (default r is 1)
    T    cos( const T& x, const T * r=nullptr ) const;                    // r*cos(x)                   (default r is 1)
    void sincos( const T& x, T& si, T& co, const T * r=nullptr ) const;   // si=r*sin(x), co=r*cos(x)   (default r is 1)
    T    tan( const T& x ) const;                                         // sin(x) / cos(x)              (2)
    T    asin( const T& x ) const;                                        // atan2(x, sqrt(1 - x^2))      (2)
    T    acos( const T& x ) const;                                        // atan2(sqrt(1 - x^2), x)      (2)
    T    atan( const T& x ) const;                                        // atan(x)
    T    atan2( const T& y, const T& x ) const;                           // atan2(y, x)                  

    void polar_to_rect( const T& r, const T& a, T& x, T& y ) const;       // sincos(a, x, y, &r)  
    void rect_to_polar( const T& x, const T& y, T& r, T& a ) const;       // r=sqrt(x^2 + y^2), a=atan2(y, x)
    T    norm( const T& x, const T& y ) const;                            // sqrt(x^2 + y^2)
    T    hypot( const T& x, const T& y ) const;                           // norm(x, y);  (same thing)
    T    normh( const T& x, const T& y ) const;                           // sqrt(x^2 - y^2)

    T    sinh( const T& x, const T * r=nullptr ) const;                   // r*sinh(x), also r*(e^x - e^-x)/2  (default r is 1)
    T    cosh( const T& x, const T * r=nullptr ) const;                   // r*cosh(x), also r*(e^x + e^-x)/2  (default r is 1)
    void sinhcosh( const T& x, T& sih, T& coh, const T * r=nullptr ) const;// sih=r*sinh(x), coh=r*cosh(x)    (default r is 1)
    T    tanh( const T& x ) const;                                        // sinh(x) / cosh(x)            (2)
    T    asinh( const T& x ) const;                                       // log(x + sqrt(x^2 + 1))       (2)
    T    acosh( const T& x ) const;                                       // log(x + sqrt(x^2 - 1))       (2)
    T    atanh( const T& x ) const;                                       // atanh(x)
    T    atanh2( const T& y, const T& x ) const;                          // atanh(y/x)

    //-----------------------------------------------------
    // Bob's Collection of Math Identities (some are used in the implementation, most are not)
    //
    // sqrt(x)          = sqrt((x+0.25)^2 - (x-0.25)^2)         allows use of hyperbolic_vectoring
    // sqrt(x)          = sqrt((x+1)^2 - (x-1)^2) / 2           ditto, but easier to make sure |atanh((x-1)/(x+1))| is small enough
    // sqrt(x*y)        = sqrt((x+y)^2 - (x-y)^2) / 2           ditto, y doesn't matter
    // sqrt(x*y)        = sqrt(|x|) * sqrt(|y|)                 assuming x*y >= 0
    // sqrt(x^2 - y^2)  = sqrt((m+d)^2 - (m-d)^2) = 2*sqrt(m*d) where m=(x+y)/2, d=(x-y)/2 
    //                                                          factor m*d=p*s so that p is a power-of-2 and s is within -1 .. 1
    //                                                          2*sqrt(m*d) = 2*sqrt(p) * sqrt(s) = sqrt(p) * sqrt((s+1)^2 - (s-1)^2)
    //                                                          if log2(p) is even and >= 2, then sqrt(p) = 2^(log2(p)/2) = some integer
    // sqrt(1+x) - 1    = x / (sqrt(1+x)+1)
    // 1-sqrt(1-x)      = x / (sqrt(1-x)+1)
    // 1-(1-x)^2        = x * (2-x)
    //
    // pow(b,x)         = exp(log(b) * x)
    // exp(x)           = sinh(x) + cosh(x)
    // exp(-x)          = 1/exp(x)
    // exp(x+y)         = exp(x) * exp(y)
    // exp(ix)          = cos(x) + i*sin(x)                     Euler's Formula, i = sqrt(-1)
    // exp(i*pi) + 1    = 0                                     Euler's Identity
    // exp(x)-1         = tanh(x/2)*(exp(x)+1)                  expm1(x)
    // exp(x)-1         = u = exp(x); (u == 1) ? x : ((u-1) == -1) ? -1 : (u-1)*x/log(u)
    //
    // log(x)           = 2*atanh2(x-1, x+1)              
    // log(x)           = atanh2((x^2 - 1, x^2 + 1) 
    // log(x*y)         = log(x) + log(y)
    // log(x/4)         = atanh2(x-0.25, x+0.25) 
    // log(x/y)         = log(x) - log(y)
    // log(x/y)         = 2*atanh2(x-y, x+y)  
    // log(1+x)         = 2*atanh(x/(x+2))                      log1p()
    //
    // sin(-x)          = -sin(x)
    // sin(x)           = sin(i*pi/2 + f)                       where i is an integer and |f| <= pi/4
    // sin(x+y)         = sin(x)*cos(y) + cos(x)*sin(y)         
    // sin(x+pi/4)      = sqrt(2)/2 * (sin(x) + cos(x))
    // sin(i*pi/2 + f)  = +/-sin(f) or +/-cos(f)
    // sin(x*pi)        = sin(2x * pi/2) = sin((i+f)*pi/2) = +/-sin(f*pi/2) or +/-cos(f*pi/2)
    // sin(x)           = Im(e^(ix)) = (e^(ix) + e^(-ix)) / 2
    // sin(ix)          = i*sinh(x)
    // sin(x)/x         = (1 + x*x*(1/6)) == 1) ? 1 : sin(x)/x
    //
    // cos(-x)          = cos(x)
    // cos(x+y)         = cos(x)*sin(y) - sin(x)*cos(y)
    // cos(x+pi/4)      = sqrt(2)/2 * (cos(x) - sin(x))
    // cos(i*pi/2 + f)  = +/-cos(f) or +/-sin(f)
    // cos(x*pi)        = cos(2x * pi/2) = cos((i+f)*pi/2) = +/-cos(f*pi/2) or +/-sin(f*pi/2)
    // cos(x)           = Re(e^(ix)) = (e^(ix) - e^(-ix)) / 2i
    // cos(ix)          = cosh(x)
    // 1-cos(x)         = 2*sin(x/2)^2
    // (1-cos(x))/x     = (1 + x*x) == 1) ? 0.5*x : cosf1(x)/x
    //
    // tan(-x)          = -tan(x)
    // tan(x+y)         = (tan(x) + tan(y)) / (1 - tan(x)*tan(y))
    //
    // asin(-x)         = -asin(x)
    // asin(x)          = atan2(x, sqrt(1 - x^2))
    //
    // acos(-x)         = acos(-x)
    // acos(x)          = atan2(sqrt(1 - x^2), x)
    // acos(x+y)        = pi/2 - asin(x+y)
    // acos(1-x)        = 2*asin(x/2)
    // acos(dot(u, v))  = 2*atan2( |u-v|, |u+v| )
    // acos(dot(u, v))  = (dot(u,v) < 0) ? (pi - 2*asin(|-v-u|/2)) : 2*asin(|v-u|/2)
    //
    // atan(x)          = atan2(x, 1)                           this is true only when second argument is > 0 (see below)
    // atan(-x)         = -atan(x)
    // atan(1/x)        = pi/2 - atan(x)                        if x > 0
    // atan(x)          = asin(x / sqrt(1 + x^2))
    // atan(x)          = 2*atan(x / (1 + sqrt(1 + x^2)))
    // atan2(y,x)       = 2*atan(y / (sqrt(x^2 + y^2) + x))     if x >  0    (x < 0 can have inflated rounding errors, so...)
    // atan2(y,x)       = 2*atan((sqrt(x^2 + y^2) - x) / y)     if x <= 0 && y != 0
    // atan2(y,x)       = pi                                    if x <  0 && y == 0
    // atan2(y,x)       = undefined                             if x == 0 && y == 0
    //
    // pi               = 4*atan(1)                             but low 2 bits will end up as 0 for a fixed-point number
    // pi               = acos(-1)                              but that just uses atan
    // pi               = pi()                                  function in this class to return a precomputed high-precision version
    //
    // sinh(-x)         = -sinh(x)
    // sinh(x)          = (e^x - e^-x)/2
    // sinh(x+y)        = sinh(x)*cosh(y) + cosh(x)*sinh(y)
    // sinh(x)          = expm1(x) * (exp1(x)+2)/(expm1(x)+1) / 2
    //
    // cosh(-x)         = cosh(x)
    // cosh(x)          = (e^x + e^-x)/2
    // cosh(x+y)        = cosh(x)*cosh(y) + sinh(x)*sinh(y)
    //
    // tanh(-x)         = -tanh(x)
    // tanh(x)          = (e^x - e^-x) / (e^x + e^-x)
    // tanh(x)          = expm1(2*x) / (expm1(2*x) + 2)
    //
    // asinh(-x)        = -asinh(x)
    // asinh(x)         = log(x + sqrt(x^2 + 1))
    // asinh(x)         = atanh(x / sqrt(1 + x^2))
    // asinh(x)         = log1p(x * (1 + x/(sqrt(x^2 + 1)+1))
    //
    // acosh(-x)        = [illegal, x must be >= 1]
    // acosh(x)         = log(x + sqrt(x^2 - 1))
    // acosh(x)         = abs( asinh(sqrt(x^2 - 1)) )
    //
    // atanh(-x)        = -atanh(x)
    // atanh(x)         = log((1+x)/(1-x))/2 = log(1+x)/2 - log(1-x)/2    note: x must be between -1 and 1
    // atanh(x)         = asinh(x / sqrt(1 - x^2)) 
    // atanh(x)         = acosh(1 / sqrt(1 - x^2))  (+/-)
    // atanh(x)         = log1p(2*x/(1-x)) / 2
    //
    // slerp(v0,v1,t)   = sin((1-t)*Ang)/sin(Ang)*v0 + sin(t*Ang)/sin(Ang)*v1           (Ang = angle between v0 and v1
    //                    [see other equations when Ang ~= 0 or Ang ~= 180]
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
    // The basic CORDIC functions that all the above math functions ultimately use.
    //-----------------------------------------------------

    // circular rotation mode results after step n:
    //      x = gain*(x0*cos(z0) - y0*sin(z0))          gain=1.64676...
    //      y = gain*(y0*cos(z0) + x0*sin(z0))
    //      z = 0
    //
    // input ranges allowed:
    //      -1    <= x0 <= 1
    //      -1    <= y0 <= 1
    //      -PI/4 <= z0 <= PI/4
    //
    // output ranges:
    //      -sqrt(2) <= x <= sqrt(2)
    //      -sqrt(2) <= y <= sqrt(2)
    //      0        <= z <= 0
    //
    void circular_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const;

    // circular vectoring mode results after step n:
    //      x = gain*sqrt(x0^2 + y0^2)                  gain=1.64676...
    //      y = 0
    //      z = z0 + atan( y0/x0 )
    //
    // input ranges allowed:
    //      -3    <= x0 <= 3
    //      -1    <= y0 <= 1
    //      -PI/4 <= z0 <= PI/4
    //      |atan(y/x)| <= PI/4
    //
    // output ranges:
    //      0    <= x <= sqrt(2)
    //      0    <= y <= 0
    //      -PI  <= z <= PI     (if z0 == 0)
    //
    void circular_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const;
    void circular_vectoring_xy( const T& x0, const T& y0, T& x, T& y ) const;  // if z not needed

    // hyperbolic rotation mode results after step n:
    //      x = gain*(x0*cosh(z0) + y0*sinh(z0))        gain=0.828159...
    //      y = gain*(y0*cosh(z0) + x0*sinh(z0))
    //      z = 0
    //
    // input ranges allowed:
    //      -2    <= x0 <= 2
    //      -2    <= y0 <= 2
    //      |z0|  <= 1.1182...
    //
    // output ranges:
    //      1  <= x <= 2
    //      -2 <= y <= 2
    //      0  <= z <= 0
    //
    void hyperbolic_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const;

    // hyperbolic vectoring mode results after step n:
    //      x = gain*sqrt(x0^2 - y0^2)                  gain=0.828159...
    //      y = 0
    //      z = z0 + atanh( y0/x0 )
    //
    // input ranges allowed:
    //      -2  <= x0 <= 2
    //      -2  <= y0 <= 2
    //      -PI <= z0 <= PI
    //      |atanh(y0/x0)| <= 1.1182...
    //
    // output ranges:
    //      -2    <= x <= 2
    //      0     <= y <= 0
    //      -PI/2 <= z <= PI/2  (if z0 == 0)
    //
    void hyperbolic_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const;
    void hyperbolic_vectoring_xy( const T& x0, const T& y0, T& x, T& y ) const;  // if z not needed

    // linear rotation mode results after step n:
    //      x = x0
    //      y = y0 + x0*z0
    //      z = 0
    //
    // input ranges allowed:
    //      -2    <= x0 <= 2
    //      -2    <= y0 <= 2
    //      |z0|  <= 1
    //
    // output ranges:
    //      -2    <= x <= 2
    //      -2    <= y <= 2
    //      0     <= z <= 0
    //
    void linear_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const;

    // linear vectoring mode results after step n:
    //      x = x0
    //      y = 0
    //      z = z0 + y0/x0
    //
    // input ranges allowed:
    //      -2      <= x0 <= 2
    //      -2      <= y0 <= 2
    //      |y0/x0| <= 1
    //
    // output ranges:
    //      -2    <= x <= 2
    //      0     <= y <= 0
    //      -PI/2 <= z <= PI/2  (if z0 == 0)
    //
    void linear_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const;
    void linear_vectoring_xy( const T& x0, const T& y0, T& x, T& y ) const;   // if z not needed

    //-----------------------------------------------------
    // More Constants 
    //-----------------------------------------------------
    T circular_rotation_gain( void ) const;              // circular_rotation()     x result with x0=1, y0=0, z0=0 (1.646760258121066...)
    T circular_vectoring_gain( void ) const;             // circular_vectoring()    x result with x0=1, y0=0, z0=0 (1.646760258121066...)
    T hyperbolic_rotation_gain( void ) const;            // hyperbolic_rotation()   x result with x0=1, y0=0, z0=0 (0.828159360960214...)
    T hyperbolic_vectoring_gain( void ) const;           // hyperbolic_vectoring()  x result with x0=1, y0=0, z0=0 (0.828159360960214...)
    T circular_rotation_one_over_gain( void ) const;     // 1.0/circular_rotation_gain()                           (0.607252935008881...)
    T circular_vectoring_one_over_gain( void ) const;    // 1.0/circular_vectoring_gain()                          (0.607252935008881...)
    T hyperbolic_rotation_one_over_gain( void ) const;   // 1.0/hyperbolic_rotation_gain()                         (1.207497067763074...)
    T hyperbolic_vectoring_one_over_gain( void ) const;  // 1.0/hyperbolic_vectoring_gain()                        (1.207497067763074...)

    T circular_angle_max( void ) const;                  // circular_vectoring      z result with x0=1, y0=1, z0=0 (0.785398163397447... == PI/4)
    T hyperbolic_angle_max( void ) const;                // hyperbolic_vectoring    z result with x0=1, y0=1, z0=0 (1.118173015526502... == [need to figure out])

    //-----------------------------------------------------
    // These version are used internally, but making them available publically.
    // In general, you should only call the earlier routines.
    //-----------------------------------------------------
    T    mad( const T& x, const T& y, const T addend, bool do_reduce ) const; // same but override do_reduce
    T    mul( const T& x, const T& y, bool do_reduce ) const;             // same but override do_reduce
    T    dad( const T& y, const T& x, const T addend, bool do_reduce ) const; // same but override do_reduce
    T    div( const T& y, const T& x, bool do_reduce ) const;             // same but override do_reduce
    T    log( const T& x, bool do_reduce ) const;                         // 2*atan2(x-1, x+1, do_reduce)    
    T    norm( const T& _x, const T& _y, bool do_reduce ) const;          
    T    atan2(  const T& y, const T& x, bool do_reduce, bool x_is_one=false, T * r=nullptr ) const; // same but override do_reduce 
    T    atanh2( const T& y, const T& x, bool do_reduce, bool x_is_one=false ) const; // same but override do_reduce
    void sincos( const T& x, T& si, T& co, bool do_reduce, bool need_si=true, bool need_co=true, const T * r=nullptr ) const;
    void sinhcosh( const T& x, T& sih, T& coh, bool do_reduce, bool need_sih=true, bool need_coh=true, const T * r=nullptr ) const;

    //-----------------------------------------------------
    // Argument Range Reduction Routines
    //
    // If your Cordic object was allocated with do_reduce==false, then
    // you are still free to call these routines to perform reduction where you want.
    // However, it's easier to just allocate another Cordic object with do_reduce=true
    // and use the appropriate one depending on whether you want reduction or not.
    //
    // Look at the implementation below to understand what these routines are doing.  
    //-----------------------------------------------------
    void reduce_arg( T& x, int32_t& x_lshift, bool& sign, bool shift_x=true, bool normalize=false, bool for_sqrt=false ) const;
    void reduce_mul_args( T& x, T& y, int32_t& x_lshift, int32_t& y_lshift, bool& sign ) const; 
    void reduce_div_args( T& y, T& x, int32_t& y_lshift, int32_t& x_lshift, bool& sign ) const;
    void reduce_sqrt_arg( T& x, int32_t& lshift ) const;
    void reduce_exp_arg( FLT b, T& x, T& factor ) const;
    void reduce_log_arg( T& x, T& addend ) const;                                            
    void reduce_atan2_args( T& y, T& x, bool& y_sign, bool& x_sign, bool& swapped, bool& is_pi ) const;     
    void reduce_norm_args( T& x, T& y, int32_t& lshift, bool& swapped ) const;
    void reduce_sincos_arg( T& a, uint32_t& quadrant, bool& sign, bool& did_minus_pi_div_4 ) const;
    void reduce_sinhcosh_arg( T& x, T& sinh_i, T& cosh_i, bool& sign ) const;                  

    //-----------------------------------------------------
    // Logging Support
    //
    // There are many reasons you might want to log Cordic operations.
    // You can use the default Logger (new Logger<T,FLT>( ... )) or supply your 
    // own Logger subclass.  
    //
    // NOTE: Logging is done globally, not just for one Cordic.  
    //       Thus the static methods here.
    //-----------------------------------------------------
    static void            logger_set( Logger<T,FLT> * logger );  // null means use the default logger
    static Logger<T,FLT> * logger_get( void );                    // returns current logger
    static std::string     op_to_str( uint16_t op );              // supply this to Logger constructor

    enum class OP
    {
        make_constant,
        use_constant,
        to_flt,
        assign,

        isgreater,
        isgreaterequal,
        isless,
        islessequal,
        islessgreater,
        isunordered,
        isunequal,
        isequal,

        abs,
        neg,
        floor,
        ceil,

        add,
        sub,
        mad,
        fma,
        mul,
        lshift,
        rshift,
        dad,
        div,
        rcp,
        sqrt,
        rsqrt,
        cbrt,
        rcbrt,
        fdim,
        fmax,
        fmin,

        exp,
        expm1,
        exp2,
        exp10,
        pow,
        powc,
        log,
        logb,
        logc,
        log2,
        log10,

        sin,
        cos,
        sincos,
        tan,
        asin,
        acos,
        atan,
        atan2,

        polar_to_rect,
        rect_to_polar,
        norm,
        hypot,
        normh,

        sinh,
        cosh,
        sinhcosh,
        tanh,
        asinh,
        acosh,
        atanh,
        atanh2,
    };

    static constexpr uint32_t OP_cnt = uint32_t(OP::atanh2) + 1;

private:
    class Impl;
    Impl * impl;

    static Logger<T,FLT> * logger;
};

//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
// 
// IMPLEMENTATION  IMPLEMENTATION  IMPLEMENTATION
//
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
template< typename T, typename FLT >
struct Cordic<T,FLT>::Impl
{
    uint32_t                    int_w;
    uint32_t                    frac_w;
    bool                        do_reduce;
    uint32_t                    n;

    T                           maxint;
    T                           maxval;
    T                           minval;
    T                           zero;
    T                           one;
    T                           half;
    T                           quarter;
    T                           sqrt2;
    T                           sqrt2_div_2;
    T                           pi;
    T                           pi_div_2;
    T                           pi_div_4;
    T                           two_div_pi;
    T                           four_div_pi;
    T                           e;
    T                           log2;                                   
    T                           log10;                                 

    T *                         circular_atan;                          // circular atan values
    T                           circular_rotation_gain;                 // circular rotation gain
    T                           circular_rotation_one_over_gain;        // circular rotation 1/gain
    T                           circular_vectoring_gain;                // circular vectoring gain
    T                           circular_vectoring_one_over_gain;       // circular vectoring 1/gain
    T                           circular_angle_max;                     // circular vectoring |z0| max value

    T *                         hyperbolic_atanh;                       // hyperbolic atanh values
    T                           hyperbolic_rotation_gain;               // hyperbolic rotation gain
    T                           hyperbolic_rotation_one_over_gain;      // hyperbolic rotation 1/gain
    T                           hyperbolic_vectoring_gain;              // hyperbolic vectoring gain
    T                           hyperbolic_vectoring_one_over_gain;     // hyperbolic vectoring 1/gain
    T                           hyperbolic_angle_max;                   // hyperbolic vectoring |z0| max value

    T *                         linear_pow2;                            // linear 2^(-i) values
    T *                         reduce_sinhcosh_sinh_i;                // for each possible integer value, sinh(i)
    T *                         reduce_sinhcosh_cosh_i;                // for each possible integer value, cosh(i)
    bool *                      reduce_sinhcosh_sinh_i_oflow;          // boolean indicating if this index creates too big of a number
    bool *                      reduce_sinhcosh_cosh_i_oflow;          // boolean indicating if this index creates too big of a number
    FLT *                       reduce_exp_factor;                      // for each possible integer value, std::exp(i)
    T *                         reduce_log_addend;                      // for each possible lshift value, log( 1 << lshift )
};

//-----------------------------------------------------
// Constructor
//-----------------------------------------------------
template< typename T, typename FLT >
Cordic<T,FLT>::Cordic( uint32_t int_w, uint32_t frac_w, bool do_reduce, uint32_t n )
{
    if ( n == 0 ) n = frac_w;
    if ( logger != nullptr ) logger->cordic_constructed( this, int_w, frac_w, n );

    cassert( (1+int_w+frac_w) <= (sizeof( T ) * 8), "1+int_w+frac_w does not fit in T container" );
    cassert( int_w  != 0, "int_w must be > 0 currently" );
    cassert( frac_w != 0, "frac_w must be > 0 currently" );

    impl = new Impl;

    impl->int_w   = int_w;
    impl->frac_w  = frac_w;
    impl->do_reduce = do_reduce;
    impl->n       = n;
    impl->maxint        = (T(1) << int_w) - 1;
    impl->one           = T(1) << frac_w;                       // required before calling to_t()

    constructed( impl->maxval        );                         // to keep logging analysis happy
    constructed( impl->minval        );
    constructed( impl->zero          );
    constructed( impl->one           );
    constructed( impl->half          );
    constructed( impl->quarter       );
    constructed( impl->sqrt2         );
    constructed( impl->sqrt2_div_2   );
    constructed( impl->pi            );
    constructed( impl->pi_div_2      );
    constructed( impl->pi_div_4      );
    constructed( impl->two_div_pi    );
    constructed( impl->four_div_pi   );
    constructed( impl->e             );
    constructed( impl->log2          );
    constructed( impl->log10         );

    assign( impl->maxval        , to_t( std::pow( 2.0, int_w ) - 1.0 ) );
    assign( impl->minval        , to_t( 1.0 / std::pow( 2.0, frac_w ) ) );
    assign( impl->zero          , to_t( 0.0 ) );
    assign( impl->one           , to_t( 1.0 ) );
    assign( impl->half          , to_t( 0.5 ) );
    assign( impl->quarter       , to_t( 0.25 ) );
    assign( impl->sqrt2         , to_t( std::sqrt( 2.0 ) ) );
    assign( impl->sqrt2_div_2   , to_t( std::sqrt( 2.0 ) / FLT(2.0) ) );
    assign( impl->pi            , to_t( std::acos( FLT(-1.0) ) ) );
    assign( impl->pi_div_2      , to_t( std::acos( FLT(-1.0) ) / FLT(2.0) ) );
    assign( impl->pi_div_4      , to_t( std::acos( FLT(-1.0) ) / FLT(4.0) ) );
    assign( impl->two_div_pi    , to_t( FLT(2.0) / std::acos( FLT(-1.0) ) ) );
    assign( impl->four_div_pi   , to_t( FLT(4.0) / std::acos( FLT(-1.0) ) ) );
    assign( impl->e             , to_t( std::exp( FLT(  1 ) ) ) );
    assign( impl->log2          , to_t( std::log( FLT(  2 ) ) ) );
    assign( impl->log10         , to_t( std::log( FLT( 10 ) ) ) );

    impl->circular_atan    = new T[n+1];
    impl->hyperbolic_atanh = new T[n+1];
    impl->linear_pow2      = new T[n+1];

    // compute atan/atanh table in high-resolution floating point
    //
    for( uint32_t i = 0; i <= n; i++ )
    {
        assign( impl->linear_pow2[i], T(1) << (frac_w-i) );
        FLT pow2                  = to_flt( impl->linear_pow2[i] );
        FLT a                     = std::atan( pow2 );
        FLT ah                    = std::atanh( pow2 );
        assign( impl->circular_atan[i], to_t( a ) );
        assign( impl->hyperbolic_atanh[i], (i == 0) ? to_t(-1) : to_t( ah ) );

        if ( debug ) printf( "i=%2d a=%30.27g ah=%30.27g y=%30.27g\n", i, double(a), double(ah), double(pow2) );
    }

    // calculate max |z0| angle allowed
    T xx, yy, zz;
    assign( impl->circular_angle_max, impl->one );     // to avoid triggering assert
    assign( impl->hyperbolic_angle_max, impl->zero );  // to disable assert
    circular_vectoring(   impl->one,     impl->one, impl->zero, xx, yy, impl->circular_angle_max );
    hyperbolic_vectoring( to_t(0.5), impl->one, impl->zero, xx, yy, impl->hyperbolic_angle_max );
    if ( debug ) std::cout << "circular_angle_max="             << std::setw(30) << to_flt(impl->circular_angle_max) << "\n";
    if ( debug ) std::cout << "hyperbolic_angle_max="           << std::setw(30) << to_flt(impl->hyperbolic_angle_max) << "\n";
    
    // calculate gain by plugging in x=1,y=0,z=0 into CORDICs
    circular_rotation(    impl->one, impl->zero, impl->zero, impl->circular_rotation_gain,    yy, zz );
    circular_vectoring(   impl->one, impl->zero, impl->zero, impl->circular_vectoring_gain,   yy, zz );
    hyperbolic_rotation(  impl->one, impl->zero, impl->zero, impl->hyperbolic_rotation_gain,  yy, zz );
    hyperbolic_vectoring( impl->one, impl->zero, impl->zero, impl->hyperbolic_vectoring_gain, yy, zz );

    // calculate 1/gain which are the multiplication factors
    assign( impl->circular_rotation_one_over_gain,    to_t( FLT(1) / to_flt(impl->circular_rotation_gain) ) );
    assign( impl->circular_vectoring_one_over_gain,   to_t( FLT(1) / to_flt(impl->circular_vectoring_gain) ) );
    assign( impl->hyperbolic_rotation_one_over_gain,  to_t( FLT(1) / to_flt(impl->hyperbolic_rotation_gain) ) );
    assign( impl->hyperbolic_vectoring_one_over_gain, to_t( FLT(1) / to_flt(impl->hyperbolic_vectoring_gain) ) );
    if ( debug ) std::cout << "circular_rotation_gain="             << std::setw(30) << to_flt(impl->circular_rotation_gain) << "\n";
    if ( debug ) std::cout << "circular_vectoring_gain="            << std::setw(30) << to_flt(impl->circular_vectoring_gain) << "\n";
    if ( debug ) std::cout << "hyperbolic_rotation_gain="           << std::setw(30) << to_flt(impl->hyperbolic_rotation_gain) << "\n";
    if ( debug ) std::cout << "hyperbolic_vectoring_gain="          << std::setw(30) << to_flt(impl->hyperbolic_vectoring_gain) << "\n";
    if ( debug ) std::cout << "circular_rotation_one_over_gain="    << std::setw(30) << to_flt(impl->circular_rotation_one_over_gain) << "\n";
    if ( debug ) std::cout << "circular_vectoring_one_over_gain="   << std::setw(30) << to_flt(impl->circular_vectoring_one_over_gain) << "\n";
    if ( debug ) std::cout << "hyperbolic_rotation_one_over_gain="  << std::setw(30) << to_flt(impl->hyperbolic_rotation_one_over_gain) << "\n";
    if ( debug ) std::cout << "hyperbolic_vectoring_one_over_gain=" << std::setw(30) << to_flt(impl->hyperbolic_vectoring_one_over_gain) << "\n";

    // construct LUT used by reduce_sinhcosh_arg();
    // use integer part plus 0.25 bit of fraction
    cassert( int_w <= 24, "too many cases to worry about" );
    uint32_t N = 1 << (2+int_w);
    T *        addend       = new T[N];
    uint32_t * quadrant     = new uint32_t[N];
    bool *     odd_pi_div_4 = new bool[N];
    T *        sinh_i       = new T[N];
    T *        cosh_i       = new T[N];
    bool *     sinh_i_oflow = new bool[N];
    bool *     cosh_i_oflow = new bool[N];
    impl->reduce_sinhcosh_sinh_i       = sinh_i;
    impl->reduce_sinhcosh_cosh_i       = cosh_i;
    impl->reduce_sinhcosh_sinh_i_oflow = sinh_i_oflow;
    impl->reduce_sinhcosh_cosh_i_oflow = cosh_i_oflow;
    const FLT PI       = M_PI;
    const FLT PI_DIV_2 = PI / 2.0;
    const FLT PI_DIV_4 = PI / 4.0;
    const T   MASK     = (T(1) << (int_w+T(1)))-T(1);  // include 0.5 bit of fraction
    const T   MAX      = (T(1) << (int_w+frac_w))-T(1);
    const FLT MAX_F    = to_flt( MAX );
    for( T i = 0; i <= MASK; i++ )
    {
        FLT i_f = FLT(i) / 4.0;

        FLT sinh_i_f    = std::sinh( i_f );
        FLT cosh_i_f    = std::cosh( i_f );
        sinh_i_oflow[i] = sinh_i_f > MAX_F;
        cosh_i_oflow[i] = cosh_i_f > MAX_F;
        assign( sinh_i[i], sinh_i_oflow[i] ? MAX : to_t( std::sinh( i_f ) ) );
        assign( cosh_i[i], cosh_i_oflow[i] ? MAX : to_t( std::cosh( i_f ) ) );
        if ( debug ) std::cout << "reduce_sinhcosh_arg LUT: i_f=" << i_f << " sinh_i=" << to_flt(sinh_i[i]) << " cosh_i=" << to_flt(cosh_i[i]) << 
                                  " sinh_i_oflow=" << sinh_i_oflow[i] << " cosh_i_oflow=" << cosh_i_oflow[i] << "\n";
    }

    // construct LUT used by reduce_exp_arg()
    // values for negative integers come first.
    FLT * factor = new FLT[2*N];
    impl->reduce_exp_factor = factor;
    T MIN_INT = -impl->maxint - 1;
    for( T i = MIN_INT; i <= impl->maxint; i++ )
    {
        T index = i - MIN_INT;
        factor[index] = std::exp(FLT(i));
        if ( debug ) std::cout << "reduce_exp_arg LUT: factor[" << to_rstring(i) << "]=" << factor[index] << 
                                  " index=" << to_rstring(index) << "\n";
    }

    // construct LUT used by reduce_log_arg()
    addend = new T[frac_w+int_w];
    impl->reduce_log_addend = addend;
    for( int32_t i = -frac_w; i <= int32_t(int_w); i++ )
    {
        double addend_f = std::log( std::pow( 2.0, double( i ) ) );
        assign( addend[frac_w+i], to_t( addend_f ) );
        if ( debug ) std::cout << "addend[]=0x" << std::hex << to_rstring(addend[frac_w+i]) << "\n" << std::dec;
        if ( debug ) std::cout << "reduce_log_arg LUT: addend[" << to_rstring(i) << "]=" << to_flt(addend[frac_w+i]) << 
                                  " addend_f=" << addend_f << "\n";
    }
}

template< typename T, typename FLT >
Cordic<T,FLT>::~Cordic( void )
{
    if ( logger != nullptr ) logger->cordic_destructed( this );

    delete impl->circular_atan;
    delete impl->hyperbolic_atanh;
    delete impl->linear_pow2;
    delete impl->reduce_sinhcosh_sinh_i;
    delete impl->reduce_sinhcosh_cosh_i;
    delete impl->reduce_sinhcosh_sinh_i_oflow;
    delete impl->reduce_sinhcosh_cosh_i_oflow;
    delete impl->reduce_exp_factor;
    delete impl->reduce_log_addend;
    delete impl;
    impl = nullptr;
}

//-----------------------------------------------------
// Logging
//-----------------------------------------------------
template< typename T, typename FLT >
Logger<T,FLT> * Cordic<T,FLT>::logger = nullptr;

template< typename T, typename FLT >
void Cordic<T,FLT>::logger_set( Logger<T,FLT> * _logger )
{
    logger = _logger;
}    

template< typename T, typename FLT >
Logger<T,FLT> * Cordic<T,FLT>::logger_get( void )
{
    return logger;
}    

template< typename T, typename FLT >
std::string Cordic<T,FLT>::op_to_str( uint16_t op )
{
    #define _ocase( op ) case OP::op: return #op;
    
    switch( OP( op ) )
    {
        _ocase( make_constant )
        _ocase( use_constant )
        _ocase( to_flt )
        _ocase( assign )

        _ocase( isgreater )
        _ocase( isgreaterequal )
        _ocase( isless )
        _ocase( islessequal )
        _ocase( islessgreater )
        _ocase( isunordered )
        _ocase( isunequal )
        _ocase( isequal )

        _ocase( abs )
        _ocase( neg )
        _ocase( floor )
        _ocase( ceil )

        _ocase( add )
        _ocase( sub )
        _ocase( mad )
        _ocase( fma )
        _ocase( mul )
        _ocase( lshift )
        _ocase( rshift )
        _ocase( dad )
        _ocase( div )
        _ocase( rcp )
        _ocase( sqrt )
        _ocase( rsqrt )

        _ocase( exp )
        _ocase( exp2 )
        _ocase( exp10 )
        _ocase( pow )
        _ocase( powc )
        _ocase( log )
        _ocase( logb )
        _ocase( logc )
        _ocase( log2 )
        _ocase( log10 )

        _ocase( sin )
        _ocase( cos )
        _ocase( sincos )
        _ocase( tan )
        _ocase( asin )
        _ocase( acos )
        _ocase( atan )
        _ocase( atan2 )

        _ocase( polar_to_rect )
        _ocase( rect_to_polar )
        _ocase( norm )
        _ocase( hypot )
        _ocase( normh )

        _ocase( sinh )
        _ocase( cosh )
        _ocase( sinhcosh )
        _ocase( tanh )
        _ocase( asinh )
        _ocase( acosh )
        _ocase( atanh )
        _ocase( atanh2 )
        default: return "<unknown OP>";
    }
}

#define _log1( op, opnd1 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op1( uint16_t(Cordic<T,FLT>::OP::op), &opnd1 )
#define _log2( op, opnd1, opnd2 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op2( uint16_t(Cordic<T,FLT>::OP::op), &opnd1, &opnd2 )
#define _log2i( op, opnd1, opnd2 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op2( uint16_t(Cordic<T,FLT>::OP::op), &opnd1, opnd2 )
#define _log2f( op, opnd1, opnd2 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op2( uint16_t(Cordic<T,FLT>::OP::op), &opnd1, opnd2 )
#define _log3( op, opnd1, opnd2, opnd3 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op3( uint16_t(Cordic<T,FLT>::OP::op), &opnd1, &opnd2, &opnd3 )
#define _log4( op, opnd1, opnd2, opnd3, opnd4 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op4( uint16_t(Cordic<T,FLT>::OP::op), &opnd1, &opnd2, &opnd3, &opnd4 )

//-----------------------------------------------------
// Constants
//-----------------------------------------------------
template< typename T, typename FLT >
inline uint32_t Cordic<T,FLT>::int_w( void ) const
{
    return impl->int_w;
}

template< typename T, typename FLT >
inline uint32_t Cordic<T,FLT>::frac_w( void ) const
{
    return impl->frac_w;
}

template< typename T, typename FLT >
inline uint32_t Cordic<T,FLT>::n( void ) const
{
    return impl->n;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::maxint( void ) const
{
    _log1( use_constant, impl->maxint );
    return impl->maxint;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::maxval( void ) const
{
    _log1( use_constant, impl->maxval );
    return impl->maxval;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::minval( void ) const
{
    _log1( use_constant, impl->minval );
    return impl->minval;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::zero( void ) const
{
    _log1( use_constant, impl->zero );
    return impl->zero;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::one( void ) const
{
    _log1( use_constant, impl->one ); 
    return impl->one;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::half( void ) const
{
    _log1( use_constant, impl->half ); 
    return impl->half;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::quarter( void ) const
{
    _log1( use_constant, impl->quarter ); 
    return impl->quarter;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sqrt2( void ) const
{
    _log1( use_constant, impl->sqrt2 ); 
    return impl->sqrt2;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sqrt2_div_2( void ) const
{
    _log1( use_constant, impl->sqrt2_div_2 ); 
    return impl->sqrt2_div_2;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pi( void ) const
{
    _log1( use_constant, impl->pi ); 
    return impl->pi;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pi_div_2( void ) const
{
    _log1( use_constant, impl->pi_div_2 ); 
    return impl->pi_div_2;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pi_div_4( void ) const
{
    _log1( use_constant, impl->pi_div_4 ); 
    return impl->pi_div_4;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::two_div_pi( void ) const
{
    _log1( use_constant, impl->two_div_pi ); 
    return impl->two_div_pi;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::four_div_pi( void ) const
{
    _log1( use_constant, impl->four_div_pi ); 
    return impl->four_div_pi;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::e( void ) const
{
    _log1( use_constant, impl->e ); 
    return impl->e;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::circular_rotation_gain( void ) const
{
    _log1( use_constant, impl->circular_rotation_gain ); 
    return impl->circular_rotation_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::circular_vectoring_gain( void ) const
{
    _log1( use_constant, impl->circular_vectoring_gain ); 
    return impl->circular_vectoring_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hyperbolic_rotation_gain( void ) const
{
    _log1( use_constant, impl->hyperbolic_rotation_gain ); 
    return impl->hyperbolic_rotation_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hyperbolic_vectoring_gain( void ) const
{
    _log1( use_constant, impl->hyperbolic_vectoring_gain ); 
    return impl->hyperbolic_vectoring_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::circular_rotation_one_over_gain( void ) const
{
    return impl->circular_rotation_one_over_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::circular_vectoring_one_over_gain( void ) const
{
    _log1( use_constant, impl->circular_vectoring_one_over_gain ); 
    return impl->circular_vectoring_one_over_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hyperbolic_rotation_one_over_gain( void ) const
{
    _log1( use_constant, impl->hyperbolic_rotation_one_over_gain ); 
    return impl->hyperbolic_rotation_one_over_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hyperbolic_vectoring_one_over_gain( void ) const
{
    _log1( use_constant, impl->hyperbolic_vectoring_one_over_gain ); 
    return impl->hyperbolic_vectoring_one_over_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::circular_angle_max( void ) const
{
    _log1( use_constant, impl->circular_angle_max ); 
    return impl->circular_angle_max;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hyperbolic_angle_max( void ) const
{
    _log1( use_constant, impl->hyperbolic_angle_max ); 
    return impl->hyperbolic_angle_max;
}

//-----------------------------------------------------
// Conversion
//-----------------------------------------------------
template< typename T, typename FLT >
inline T Cordic<T,FLT>::to_t( FLT _x ) const
{
    FLT x = _x;
    bool is_neg = x < 0.0;
    if ( is_neg ) x = -x;
    cassert( T(x) < (T(1) << impl->int_w), "to_t: integer part of |x| " + std::to_string(x) + " does not fit in int_w bits" ); 
    T x_t = x * FLT( impl->one ); // need to round
    if ( is_neg ) x_t = -x_t;
    _log2f( make_constant, x_t, _x );
    return x_t;
}

template< typename T, typename FLT >
inline FLT Cordic<T,FLT>::to_flt( const T& _x, bool can_log ) const
{
    T x = _x;
    bool is_neg = x < 0;
    if ( is_neg ) x = -x;
    FLT x_f = FLT( x ) / FLT( impl->one );
    if ( is_neg ) x_f = -x_f;
    if ( can_log ) _log2f( to_flt, _x, x_f );
    return x_f;
}

template< typename T, typename FLT >
inline std::string Cordic<T,FLT>::to_string( const T& x ) const
{
    // floating-point representation
    return std::to_string( to_flt( x ) );  
}

template< typename T, typename FLT >
inline std::string Cordic<T,FLT>::to_rstring( const T& _x ) const
{
    // raw integer representation
    T x = _x;
    bool sign = x < 0;
    if ( sign ) x = -x;
    int32_t twidth = 8 * sizeof( T );
    int32_t bwidth = impl->int_w + impl->frac_w;
    int32_t dwidth = FLT(bwidth) / std::ceil( std::log2( 10 ) );
    char s[1024];
    cassert( (dwidth+2) < int32_t(sizeof(s)), "to_rstring: need to make s buffer bigger" );
    memset( s, 0, dwidth+1 );
    int32_t m = dwidth;
    int32_t i;
    for( i = bwidth-1; i >= 0; i-- )
    {
        // avoiding divides 
        bool carry = (x >> i) & 1; 
        int32_t j;
        for( j = dwidth-1; j >= (m + 1) || carry; j-- ) {
            int32_t d = 2 * s[j] + carry;
            carry = d > 9;
            s[j] = carry ? (d - 10) : d;
        }
        m = j;
    }

    // skip leading zeros
    for( i = 0; i < (dwidth-1); i++ ) 
    {
        if ( s[i] != 0 ) break;  
    }

    // convert to characters
    for( int32_t j = i; j < dwidth; j++ ) 
    {
        s[j] += '0';
    }

    // add optional sign
    if ( sign ) {
        i--;
        s[i] = '-';
    }        

    return s+i;
}

template< typename T, typename FLT >
inline std::string Cordic<T,FLT>::to_bstring( const T& _x ) const
{
    // binary representation
    T x = _x;
    uint32_t width = 1 + impl->int_w + impl->frac_w;
    std::string bs = "";
    for( uint32_t i = 0; i < width; i++ )
    {
        if ( i == (impl->int_w+impl->frac_w) || i == impl->frac_w ) bs = " " + bs;
        const char * b = (x & 1) ? "1" : "0";
        bs = b + bs;
        x >>= 1;
    }
    return bs;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::make_fixed( bool sign, T i, T f )
{
    cassert( i >= 0 && i <= impl->maxint              , "make_fixed integer part must be in range 0 .. impl->maxint" );
    cassert( f >= 0 && f <= ((T(1) << impl->frac_w)-1), "make_fixed fractional part must be in range 0 .. (1 << frac_w)-1" );

    return (T(sign) << (impl->int_w + impl->frac_w)) |
           (T(i)    << impl->frac_w)             |
           (T(f)    << 0);
}

//-----------------------------------------------------
// The CORDIC Functions
//-----------------------------------------------------
template< typename T, typename FLT >
void Cordic<T,FLT>::circular_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
{
    //-----------------------------------------------------
    // input ranges allowed:
    //      -1  <= x0 <= 1
    //      -1  <= y0 <= 1
    //      |z0| <= 0.7854...
    //-----------------------------------------------------
    const T ONE = impl->one;
    const T ANGLE_MAX = circular_angle_max();
    if ( debug ) std::cout << "circular_rotation begin: x0,y0,z0=[ " << to_flt(x0) << ", " << to_flt(y0) << ", " << to_flt(z0) << "]\n";
    cassert( x0 >= -ONE       && x0 <= ONE,       "circular_rotation x0 must be in the range -1 .. 1" );
    cassert( y0 >= -ONE       && y0 <= ONE,       "circular_rotation y0 must be in the range -1 .. 1" );
    cassert( z0 >= -ANGLE_MAX && z0 <= ANGLE_MAX, "circular_rotation |z0| must be <= circular_angle_max()" );

    //-----------------------------------------------------
    // d = (z >= 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctan(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->n;
    for( uint32_t i = 0; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "circular_rotation: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int(z >= impl->zero) );
        if ( z >= T(0) ) {
            xi = x - (y >> i);
            yi = y + (x >> i);
            zi = z - impl->circular_atan[i];
        } else {
            xi = x + (y >> i);
            yi = y - (x >> i);
            zi = z + impl->circular_atan[i];
        }
        x = xi;
        y = yi;
        z = zi;
    }
}

template< typename T, typename FLT >
void Cordic<T,FLT>::circular_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
{
    //-----------------------------------------------------
    // input ranges allowed:
    //      -3  <= x0 <= 3
    //      -1  <= y0 <= 1
    //      -PI <= z0 <= PI
    //      |atan(y0/x0)| <= 0.7854...
    //-----------------------------------------------------
    const T ONE = impl->one;
    const T THREE = 3*ONE;
    const T PI  = impl->pi;
    const T ANGLE_MAX = circular_angle_max();
    if ( debug ) std::cout << "circular_vectoring begin: x0,y0,z0=[ " << to_flt(x0) << ", " << to_flt(y0) << ", " << to_flt(z0) << "]\n";
    cassert( x0 >= -THREE && x0 <= THREE, "circular_vectoring x0 must be in the range -3 .. 3" );
    cassert( y0 >= -ONE   && y0 <= ONE  , "circular_vectoring y0 must be in the range -1 .. 1" );
    cassert( z0 >= -PI    && z0 <= PI   , "circular_vectoring z0 must be in the range -PI .. PI" );
    cassert( std::abs( std::atan( to_flt(y0) / to_flt(x0) ) ) <= to_flt(ANGLE_MAX),
                                        "circular_vectoring |atan(y0/x0)| must be <= circular_angle_max()" );

    //-----------------------------------------------------
    // d = (y < 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctan(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->n;
    for( uint32_t i = 0; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "circular_vectoring: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int((x < impl->zero) != (y < impl->zero)) );
        if ( y < T(0) ) {
            xi = x - (y >> i);
            yi = y + (x >> i);
            zi = z - impl->circular_atan[i];
        } else {
            xi = x + (y >> i);
            yi = y - (x >> i);
            zi = z + impl->circular_atan[i];
        }
        x = xi;
        y = yi;
        z = zi;
    }
}

template< typename T, typename FLT >
void Cordic<T,FLT>::circular_vectoring_xy( const T& x0, const T& y0, T& x, T& y ) const
{
    //-----------------------------------------------------
    // input ranges allowed:
    //      -3  <= x0 <= 3
    //      -1  <= y0 <= 1
    //-----------------------------------------------------
    const T ONE = impl->one;
    const T THREE = 3*ONE;
    if ( debug ) std::cout << "circular_vectoring_xy begin: x0,y0=[ " << to_flt(x0) << ", " << to_flt(y0) << "\n";
    cassert( x0 >= -THREE && x0 <= THREE, "circular_vectoring_xy x0 must be in the range -3 .. 3" );
    cassert( y0 >= -ONE   && y0 <= ONE  , "circular_vectoring_xy y0 must be in the range -1 .. 1" );

    //-----------------------------------------------------
    // d = (y < 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    //-----------------------------------------------------
    x = x0;
    y = y0;
    uint32_t n = impl->n;
    for( uint32_t i = 0; i <= n; i++ )
    {
        T xi;
        T yi;
        if ( debug ) printf( "circular_vectoring_xy: i=%d xy=[%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), int(y < impl->zero) );
        if ( y < T(0) ) {
            xi = x - (y >> i);
            yi = y + (x >> i);
        } else {
            xi = x + (y >> i);
            yi = y - (x >> i);
        }
        x = xi;
        y = yi;
    }
}

template< typename T, typename FLT >
void Cordic<T,FLT>::hyperbolic_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
{
    //-----------------------------------------------------
    // input ranges allowed:
    //      -1  <= x0 <= 1
    //      -1  <= y0 <= 1
    //      |z0| <= 1.1182...
    //-----------------------------------------------------
    const T TWO = impl->one << 1;
    const T ANGLE_MAX = hyperbolic_angle_max();
    if ( debug ) std::cout << "hyperbolic_rotation begin: x0,y0,z0=[ " << to_flt(x0) << ", " << to_flt(y0) << ", " << to_flt(z0) << "]\n";
    cassert( x0 >= -TWO       && x0 <= TWO,       "hyperbolic_rotation x0 must be in the range -2 .. 2" );
    cassert( y0 >= -TWO       && y0 <= TWO,       "hyperbolic_rotation y0 must be in the range -2 .. 2" );
    cassert( z0 >= -ANGLE_MAX && z0 <= ANGLE_MAX, "hyperbolic_rotation |z0| must be <= hyperbolic_angle_max()" );

    //-----------------------------------------------------
    // d = (z >= 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctanh(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->n;
    uint32_t next_dup_i = 4;     
    for( uint32_t i = 1; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "hyperbolic_rotation: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int(z >= impl->zero) );
        if ( z >= T(0) ) {
            xi = x + (y >> i);
            yi = y + (x >> i);
            zi = z - impl->hyperbolic_atanh[i];
        } else {
            xi = x - (y >> i);
            yi = y - (x >> i);
            zi = z + impl->hyperbolic_atanh[i];
        }
        x = xi;
        y = yi;
        z = zi;

        if ( i == next_dup_i ) {
            // for hyperbolic, we must duplicate iterations 4, 13, 40, 121, ..., 3*i+1
            next_dup_i = 3*i + 1;
            i--;
        }
    }
}

template< typename T, typename FLT >
void Cordic<T,FLT>::hyperbolic_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
{
    //-----------------------------------------------------
    // input ranges allowed:
    //      -2  <= x0 <= 2
    //      -2  <= y0 <= 2
    //      -PI <= z0 <= PI
    //      |atanh(y0/x0)| <= 1.1182...
    //-----------------------------------------------------
    const T TWO = impl->one << 1;
    const T PI  = impl->pi;
    const T ANGLE_MAX = hyperbolic_angle_max();
    if ( debug ) std::cout << "hyperbolic_vectoring begin: x0,y0,z0=[ " << to_flt(x0) << ", " << to_flt(y0) << ", " << to_flt(z0) << "]\n";
    cassert( x0 >= -TWO && x0 <= TWO, "hyperbolic_vectoring x0 must be in the range -2 .. 2" );
    cassert( y0 >= -TWO && y0 <= TWO, "hyperbolic_vectoring y0 must be in the range -2 .. 2" );
    cassert( z0 >= -PI  && z0 <= PI , "hyperbolic_vectoring z0 must be in the range -PI .. PI" );
    cassert( (ANGLE_MAX == 0 || std::abs( std::atanh( to_flt(y0) / to_flt(x0) ) ) <= to_flt(ANGLE_MAX)),
                                        "hyperbolic_vectoring |atanh(y0/x0)| must be <= hyperbolic_angle_max()" );

    //-----------------------------------------------------
    // d = (y < 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctanh(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->n;
    uint32_t next_dup_i = 4;     
    for( uint32_t i = 1; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "hyperbolic_vectoring: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, 
                             to_flt(x), to_flt(y), to_flt(z), int((x < impl->zero) != (y < impl->zero)) );
        if ( y < T(0) ) {
            xi = x + (y >> i);
            yi = y + (x >> i);
            zi = z - impl->hyperbolic_atanh[i];
        } else {
            xi = x - (y >> i);
            yi = y - (x >> i);
            zi = z + impl->hyperbolic_atanh[i];
        }
        x = xi;
        y = yi;
        z = zi;

        if ( i == next_dup_i ) {
            // for hyperbolic, we must duplicate iterations 4, 13, 40, 121, ..., 3*i+1
            next_dup_i = 3*i + 1;
            i--;
        }
    }
}

template< typename T, typename FLT >
void Cordic<T,FLT>::hyperbolic_vectoring_xy( const T& x0, const T& y0, T& x, T& y ) const
{
    //-----------------------------------------------------
    // input ranges allowed:
    //      -2  <= x0 <= 2
    //      -2  <= y0 <= 2
    //-----------------------------------------------------
    const T TWO = impl->one << 1;
    const T PI  = impl->pi;
    if ( debug ) std::cout << "hyperbolic_vectoring_xy begin: x0,y0=[ " << to_flt(x0) << ", " << to_flt(y0) << "\n";
    cassert( x0 >= -TWO && x0 <= TWO, "hyperbolic_vectoring_xy x0 must be in the range -2 .. 2" );
    cassert( y0 >= -TWO && y0 <= TWO, "hyperbolic_vectoring_xy y0 must be in the range -2 .. 2" );

    //-----------------------------------------------------
    // d = (y < 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    //-----------------------------------------------------
    x = x0;
    y = y0;
    uint32_t n = impl->n;
    uint32_t next_dup_i = 4;     
    for( uint32_t i = 1; i <= n; i++ )
    {
        T xi;
        T yi;
        if ( debug ) printf( "hyperbolic_vectoring_xy: i=%d xy=[%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), int(y < impl->zero) );
        if ( y < T(0) ) {
            xi = x + (y >> i);
            yi = y + (x >> i);
        } else {
            xi = x - (y >> i);
            yi = y - (x >> i);
        }
        x = xi;
        y = yi;

        if ( i == next_dup_i ) {
            // for hyperbolic, we must duplicate iterations 4, 13, 40, 121, ..., 3*i+1
            next_dup_i = 3*i + 1;
            i--;
        }
    }
}

template< typename T, typename FLT >
void Cordic<T,FLT>::linear_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
{
    //-----------------------------------------------------
    // input ranges allowed:
    //      -2    <= x0 <= 2
    //      -2    <= y0 <= 2
    //      |z0|  <= 1
    //-----------------------------------------------------
    const T ONE = impl->one;
    const T TWO = ONE << 1;
    if ( debug ) std::cout << "linear_rotation begin: x0,y0,z0=[ " << to_flt(x0) << ", " << to_flt(y0) << ", " << to_flt(z0) << "]\n";
    cassert( x0 >= -TWO && x0 <= TWO, "linear_rotation x0 must be in the range -2 .. 2" );
    cassert( y0 >= -TWO && y0 <= TWO, "linear_rotation y0 must be in the range -2 .. 2" );
    //cassert( z0 >= -ONE && z0 <= ONE, "linear_rotation z0 must be in the range -1 .. 1" );
    
    //-----------------------------------------------------
    // d = (z >= 0) ? 1 : -1
    // xi = x 
    // yi = y + d*(x >> i)
    // zi = z - d*2^(-i)
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->n;
    for( uint32_t i = 0; i <= n; i++ )
    {
        if ( debug ) printf( "linear_rotation: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int(z >= impl->zero) );
        T yi;
        T zi;
        if ( z >= T(0) ) {
            yi = y + (x >> i);
            zi = z - impl->linear_pow2[i];
        } else {
            yi = y - (x >> i);
            zi = z + impl->linear_pow2[i];
        }
        y = yi;
        z = zi;
    }
}

template< typename T, typename FLT >
void Cordic<T,FLT>::linear_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
{
    //-----------------------------------------------------
    // input ranges allowed:
    //      -2      <= x0 <= 2
    //      -2      <= y0 <= 2
    //      |y0/x0| <= 1
    //-----------------------------------------------------
    const T ONE = impl->one;
    const T TWO = ONE << 1;
    if ( debug ) std::cout << "linear_vectoring begin: x0,y0,z0=[ " << to_flt(x0) << ", " << to_flt(y0) << ", " << to_flt(z0) << "]\n";
    cassert( x0 >= -TWO && x0 <= TWO, "linear_vectoring x0 must be in the range -2 .. 2, got " + to_string(x0) );
    cassert( y0 >= -TWO && y0 <= TWO, "linear_vectoring y0 must be in the range -2 .. 2, got " + to_string(y0) );
    //cassert( std::abs( to_flt(y0) / to_flt(x0) ) <= FLT(1.0) &&
    //                                    "linear_vectoring y0/x0 must be in the range -1 .. 1" );
    
    //-----------------------------------------------------
    // d = (y < 0) ? 1 : -1
    // xi = x
    // yi = y + d*(x >> i)
    // zi = z - d*2^(-i)
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->n;
    for( uint32_t i = 0; i <= n; i++ )
    {
        if ( debug ) printf( "linear_vectoring: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", 
                             i, to_flt(x), to_flt(y), to_flt(z), int(y < impl->zero) );
        T yi;
        T zi;
        if ( y < T(0) ) {
            yi = y + (x >> i);
            zi = z - impl->linear_pow2[i];
        } else {
            yi = y - (x >> i);
            zi = z + impl->linear_pow2[i];
        }
        y = yi;
        z = zi;
    }
}

template< typename T, typename FLT >
void Cordic<T,FLT>::linear_vectoring_xy( const T& x0, const T& y0, T& x, T& y ) const
{
    //-----------------------------------------------------
    // input ranges allowed:
    //      -2      <= x0 <= 2
    //      -2      <= y0 <= 2
    //-----------------------------------------------------
    const T TWO = impl->one << 1;
    if ( debug ) std::cout << "linear_vectoring_xy begin: x0,y0=[ " << to_flt(x0) << ", " << to_flt(y0) << "\n";
    cassert( x0 >= -TWO && x0 <= TWO, "linear_vectoring_xy x0 must be in the range -2 .. 2" );
    cassert( y0 >= -TWO && y0 <= TWO, "linear_vectoring_xy y0 must be in the range -2 .. 2" );
    
    //-----------------------------------------------------
    // d = (y < 0) ? 1 : -1
    // xi = x
    // yi = y + d*(x >> i)
    //-----------------------------------------------------
    x = x0;
    y = y0;
    uint32_t n = impl->n;
    for( uint32_t i = 0; i <= n; i++ )
    {
        if ( debug ) printf( "linear_vectoring_xy: i=%d xy=[%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), int(y < impl->zero) );
        T yi;
        if ( y < T(0) ) {
            yi = y + (x >> i);
        } else {
            yi = y - (x >> i);
        }
        y = yi;
    }
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::constructed( const T& x ) const
{
    if ( logger != nullptr ) logger->constructed( &x, this );
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::destructed( const T& x ) const
{
    if ( logger != nullptr ) logger->destructed( &x, this );
}

template< typename T, typename FLT >
inline T&   Cordic<T,FLT>::assign( T& x, const T& y ) const
{
    _log2( assign, x, y );
    x = y;
    return x;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isgreater( const T& x, const T& y ) const
{
    _log2( isgreater, x, y );
    return x > y;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isgreaterequal( const T& x, const T& y ) const
{
    _log2( isgreaterequal, x, y );
    return x >= y;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isless( const T& x, const T& y ) const
{
    _log2( isless, x, y );
    return x < y;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::islessequal( const T& x, const T& y ) const
{
    _log2( islessequal, x, y );
    return x <= y;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::islessgreater( const T& x, const T& y ) const
{
    _log2( islessgreater, x, y );
    return x < y || x > y;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isunordered( const T& x, const T& y ) const
{
    _log2( isunordered, x, y );
    return !(x == y || x != y);
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isunequal( const T& x, const T& y ) const
{
    _log2( isunequal, x, y );
    return x != y;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isequal( const T& x, const T& y ) const
{
    _log2( isequal, x, y );
    return x == y;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::abs( const T& x ) const
{
    _log1( abs, x );
    T    x_abs  = x;
    bool x_sign = x_abs < T(0);
    if ( x_sign ) x_abs = -x;
    T    sign_mask = x_abs >> (impl->int_w + impl->frac_w);
    cassert( (sign_mask == T(0) || sign_mask == T(-1)), "abs caused overflow" ); 
    return x_abs;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::neg( const T& x ) const
{
    _log1( neg, x );
    bool x_sign = x < 0;
    T    x_neg  = -x;
    T    sign_mask = x_neg >> (impl->int_w + impl->frac_w);
    cassert( (x == 0 || sign_mask == (x_sign ? T(0) : T(-1))), "neg caused overflow" ); 
    return x_neg;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::floor( const T& x ) const
{
    _log1( floor, x );
    T frac_mask = (T(1) << impl->frac_w) - 1;
    if ( (x & frac_mask) == 0 ) {
        return x;
    } else if ( x > 0 ) {
        return x & ~frac_mask;
    } else {
        return (x & ~frac_mask) - impl->one;
    }
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::ceil( const T& x ) const
{
    _log1( ceil, x );
    T frac_mask = (T(1) << impl->frac_w) - 1;
    if ( (x & frac_mask) == 0 ) {
        return x;
    } else if ( x > 0 ) {
        return (x & ~frac_mask) + impl->one;
    } else {
        return x & ~frac_mask;
    }
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::add( const T& x, const T& y ) const
{
    _log2( add, x, y );
    bool x_sign = x < T(0);
    bool y_sign = y < T(0);
    T    sum    = x + y;
    T    sign_mask = sum >> (impl->int_w + impl->frac_w);
    cassert( sign_mask == T(0) || sign_mask == T(-1), "add caused overflow" );
    return sum;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sub( const T& x, const T& y ) const
{
    _log2( sub, x, y );
    bool x_sign = x < T(0);
    bool y_sign = y < T(0);
    T    sum    = x - y;
    T    sign_mask = sum >> (impl->int_w + impl->frac_w);
    cassert( sign_mask == 0 || sign_mask == T(-1), "sub caused overflow" );
    return sum;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::mad( const T& _x, const T& _y, const T addend, bool do_reduce ) const
{
    _log3( mad, _x, _y, addend );
    T x = _x;
    T y = _y;
    if ( debug ) std::cout << "mad begin: x_orig=" << to_flt(x) << " y_orig=" << to_flt(y) << " addend=" << to_flt(addend) << " do_reduce=" << do_reduce << "\n";
    cassert( do_reduce || addend >= 0, "mad addend must be non-negative" );
    int32_t x_lshift;
    int32_t y_lshift;
    bool    sign;
    if ( do_reduce ) reduce_mul_args( x, y, x_lshift, y_lshift, sign );

    T xx, yy, zz;
    linear_rotation( x, do_reduce ? impl->zero : addend, y, xx, yy, zz );
    if ( do_reduce ) {
        yy = lshift( yy, x_lshift + y_lshift );
        yy += addend;
        if ( sign ) yy = -yy;
    }
    if ( debug ) std::cout << "mad end: x_orig=" << to_flt(_x) << " y_orig=" << to_flt(_y) << 
                              " addend=" << to_flt(addend) << " do_reduce=" << do_reduce << 
                              " x_reduced=" << to_flt(x) << " y_reduced=" << to_flt(y) << 
                              " yy=" << to_flt(yy) << " x_lshift=" << x_lshift << " y_lshift=" << y_lshift << 
                              " sign=" << sign << "\n";
    return yy;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::mad( const T& x, const T& y, const T& addend ) const
{
    return mad( x, y, addend, impl->do_reduce );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::fma( const T& x, const T& y, const T& addend ) const
{
    return mad( x, y, addend );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::mul( const T& x, const T& y ) const
{
    return mad( x, y, impl->zero );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::mul( const T& x, const T& y, bool do_reduce ) const
{
    return mad( x, y, impl->zero, do_reduce );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::lshift( const T& x, int ls ) const
{
    _log2i( lshift, x, ls );
    cassert( x >= 0, "lshift x should be non-negative" );
    if ( ls > 0 ) {
        //-----------------------------------------------------
        // For now, crap out if we overflow.
        // At some point, we'll have options to saturate or set a flag in the container.
        //-----------------------------------------------------
        int32_t ls_max = impl->int_w;
        uint32_t i = x >> impl->frac_w;
        cassert( i <= impl->maxint, "lshift x integer part should be <= impl->maxint"  );
        while( i != 0 ) 
        {
            ls_max--;
            i >>= 1;
        }
        if ( ls > ls_max ) {
            std::cout << "lshift x << " << ls << " will overflow x\n";
            exit( 1 );
        }
        return x << ls;
    } else if ( ls < 0 ) {
        return x >> -ls;
    } else {
        return x;
    }
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::rshift( const T& x, int rs ) const
{
    return lshift( x, -rs );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::dad( const T& _y, const T& _x, const T addend, bool do_reduce ) const
{
    _log3( dad, _y, _x, addend );
    T x = _x;
    T y = _y;
    if ( debug ) std::cout << "dad begin: x_orig=" << to_flt(x) << " y_orig=" << to_flt(y) << " addend=" << to_flt(addend) << " do_reduce=" << do_reduce << "\n";
    cassert( x != 0 , "dad x (denominator) must be non-zero" );
    cassert( do_reduce || addend >= 0, "dad addend must be non-negative (need to fix this soon)" );
    int32_t x_lshift;
    int32_t y_lshift;
    bool    sign;
    if ( do_reduce ) reduce_div_args( y, x, y_lshift, x_lshift, sign );

    T xx, yy, zz;
    linear_vectoring( x, y, do_reduce ? impl->zero : addend, xx, yy, zz );
    if ( do_reduce ) {
        zz = lshift( zz, y_lshift-x_lshift );
        zz += addend;
        if ( sign ) zz = -zz;
    }
    if ( debug ) std::cout << "dad end: x_orig=" << to_flt(_x) << " y_orig=" << to_flt(_y) << 
                              " addend=" << to_flt(addend) << " do_reduce=" << do_reduce <<
                              " x_reduced=" << to_flt(x) << " y_reduced=" << to_flt(y) << 
                              " zz=" << to_flt(zz) << " x_lshift=" << x_lshift << " y_lshift=" << y_lshift << 
                              " sign=" << sign << "\n";
    return zz;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::dad( const T& _y, const T& _x, const T addend ) const
{
    return dad( _y, _x, addend, impl->do_reduce );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::div( const T& y, const T& x ) const
{
    return dad( y, x, impl->zero );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::div( const T& y, const T& x, bool do_reduce ) const
{
    return dad( y, x, impl->zero, do_reduce );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::rcp( const T& x ) const
{
    _log1( rcp, x );
    return div( impl->one, x );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::sqrt( const T& _x ) const
{ 
    //-----------------------------------------------------
    // Identities:
    //     sqrt(x) = sqrt((x+1)^2 - (x-1)^2) / 2          
    // Strategy:
    //     reduce_sqrt_arg() will factor x=p*s where p is a power-of-2 
    //     where log2(p) is even >= 2 and s is between 0.5..1.
    //     sqrt(x) = sqrt(p) * sqrt(s) = 2^(log2(p)/2) * sqrt((s+1)^2 - (s-1)^2)/2
    //     Use hyperbolic_vectoring() for sqrt((s+1)^2 - (s-1)^2).
    //     Then lshift that by log2(p)/2 - 1.
    //-----------------------------------------------------
    _log1( sqrt, _x );
    T x = _x;
    if ( debug ) std::cout << "sqrt begin: x_orig=" << to_flt(_x) << " do_reduce=" << impl->do_reduce << "\n";
    int32_t ls;
    if ( impl->do_reduce ) reduce_sqrt_arg( x, ls );

    T xx, yy;
    hyperbolic_vectoring_xy( x+impl->one, x-impl->one, xx, yy );  // gain*sqrt((s+1)^2 - (s-1)^2)
    xx = mul( xx, hyperbolic_vectoring_one_over_gain(), false );   // sucks that we have to do this
    if ( impl->do_reduce ) xx = lshift( xx, ls );                  // log2(p)/2 - 1

    if ( debug ) std::cout << "sqrt end: x_orig=" << to_flt(_x) << " x_reduced=s=" << to_flt(x) << " do_reduce=" << impl->do_reduce << " xx=" << to_flt(xx) << "\n";
    return xx;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::rsqrt( const T& x ) const
{ 
    _log1( rsqrt, x );
    if ( debug ) std::cout << "rsqrt begin: x_orig=" << to_flt(x) << " do_reduce=" << impl->do_reduce << "\n";
    cassert( x != 0, "rsqrt x must not be 0" );

    // There might be a better way, but exp(-log(x)/2) is probably not it
    return div( impl->one, sqrt( x ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::cbrt( const T& _x ) const
{ 
    _log1( cbrt, _x );
    T x = _x;
    bool sign = x < 0;
    if ( sign ) x = -x;
    const T THREE = (impl->one << 1) | impl->one;
    T r = exp( div( log( x ), THREE ) );
    if ( sign ) r = -r;
    return r;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::rcbrt( const T& x ) const
{ 
    _log1( rcbrt, x );
    if ( debug ) std::cout << "rcbrt begin: x_orig=" << to_flt(x) << " do_reduce=" << impl->do_reduce << "\n";

    // There might be a better way, but exp(-log(x)/3) is probably not it
    return div( impl->one, cbrt( x ) );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::fdim( const T& x, const T& y ) const
{ 
    _log2( fdim, x, y );
    return (x >= y) ? (x - y) : 0;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::fmax( const T& x, const T& y ) const
{ 
    _log2( fmax, x, y );
    return (x >= y) ? x : y;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::fmin( const T& x, const T& y ) const
{ 
    _log2( fmin, x, y );
    return (x < y) ? x : y;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::exp( const T& _x ) const
{ 
    //-----------------------------------------------------
    // Identities:
    //     Assume: x = i + f  (integer plus fraction)
    //     exp(i+f) = exp(i) * exp(f)
    //     pow(b,x) = log(b) * exp(x) = [log(b)*exp(i)] * exp(f)
    //
    // Strategy:
    //     Find i such that x-i is in -1 .. 1.
    //     Because x can be negative, so can i.
    //     exp(i) comes from a pre-built LUT kept in FLT
    //     so we can multiply it by log(e)==1 before converting to type T and
    //     then multiplying by exp(f) here.
    //-----------------------------------------------------
    _log1( exp, _x );
    T x = _x;
    T factor;
    if ( impl->do_reduce ) reduce_exp_arg( M_E, x, factor );  // x=log(f) factor=log(e)*exp(i)

    T xx, yy, zz;
    hyperbolic_rotation( hyperbolic_rotation_one_over_gain(), hyperbolic_rotation_one_over_gain(), x, xx, yy, zz );
    if ( impl->do_reduce ) {
        if ( debug ) std::cout << "exp mid: b=" << M_E << " x_orig=" << to_flt(_x) << " f=reduced_x=" << to_flt(x) << 
                                  " exp(f)=" << to_flt(xx) << " log(b)*log(i)=" << to_flt(factor) << "\n";
        xx = mul( xx, factor, true );
    }
    if ( debug ) std::cout << "exp: x_orig=" << to_flt(_x) << " reduced_x=" << to_flt(x) << " exp=" << to_flt(xx) << "\n";
    return xx;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::expm1( const T& x ) const
{ 
    return exp( x ) - one();            // see identities for more accurate way to do this
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::exp2( const T& x ) const
{ 
    return powc( 2.0, x );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::exp10( const T& x ) const
{ 
    return powc( 10.0, x );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pow( const T& b, const T& x ) const
{ 
    _log2( pow, b, x );
    if ( b == impl->zero ) return impl->zero;
    cassert( b >= 0, "pow base b must be non-negative" );
    return exp( mul( x, log( b, true ) ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::powc( const FLT& b, const T& x ) const
{ 
    _log2f( powc, x, b );
    cassert( b > 0, "powc b must be positive" );
    const FLT log_b_f = std::log( b );
    cassert( log_b_f >= 0.0, "powc log(b) must be non-negative" );
    const T   log_b   = to_t( log_b_f );
    return exp( mul( x, log_b ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::log( const T& _x, bool do_reduce ) const
{ 
    _log1( log, _x );
    T x = _x;
    cassert( x > 0, "log: x must be positive" );
    T addend;
    if ( do_reduce ) reduce_log_arg( x, addend );
    T dv = div( x-impl->one, x+impl->one, false );
    T lg = atanh( dv ) << 1;
    if ( do_reduce ) lg += addend;
    if ( debug ) std::cout << "log: x_orig=" << to_flt(_x) << " reduced_x=" << to_flt(x) << " log=" << to_flt(lg) << "\n";
    return lg;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::log( const T& _x ) const
{ 
    return log( _x, impl->do_reduce );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::logb( const T& x, const T& b ) const
{ 
    _log2( logb, x, b );
    cassert( b > 0, "logb b must be positive" );
    return div( log(x), log(b) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::logc( const T& x, const FLT& b ) const
{ 
    _log2f( logc, x, b );
    cassert( b > 0.0, "logc b must be positive" );
    const FLT  one_over_log_b_f = FLT(1) / std::log( b );
    const T    one_over_log_b   = to_t( one_over_log_b_f );
          T    log_x            = log( x );
    const bool log_x_sign       = log_x < 0;
    if ( log_x_sign ) log_x = -log_x;
    T z = mul( log_x, one_over_log_b );
    if ( log_x_sign ) z = -z;
    return z;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::log2( const T& x ) const
{ 
    return logc( x, 2.0 );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::log10( const T& x ) const
{ 
    return logc( x, 10.0 );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sin( const T& x, const T * r ) const
{ 
    T si;
    T co;
    sincos( x, si, co, impl->do_reduce, true, false, r );
    return si;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::cos( const T& x, const T * r ) const
{ 
    T si;
    T co;
    sincos( x, si, co, impl->do_reduce, false, true, r );
    return co;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::sincos( const T& x, T& si, T& co, const T * r ) const             
{
    sincos( x, si, co, impl->do_reduce, true, true, r );
}

template< typename T, typename FLT >
void Cordic<T,FLT>::sincos( const T& _x, T& si, T& co, bool do_reduce, bool need_si, bool need_co, const T * _r ) const             
{ 
    if ( _r != nullptr ) {
        _log4( sincos, _x, si, co, *_r );
    } else {
        _log3( sincos, _x, si, co );
    }
    T x = _x;
    uint32_t quadrant;
    bool x_sign;
    bool did_minus_pi_div_4;
    if ( do_reduce ) {
        //-----------------------------------------------------
        // reduce_sincos_arg() will get x in the range 0 .. PI/2 and tell us the quadrant.
        // It will then check if x is still > PI/4 and, if so, subtract PI/4 and
        // set did_minus_pi_div_4.  If did_minus_pi_div_4 is true, then we need
        // to do some adjustments below after the cordic routine completes.
        //-----------------------------------------------------
        reduce_sincos_arg( x, quadrant, x_sign, did_minus_pi_div_4 );
    }

    T r = circular_rotation_one_over_gain();
    int32_t r_lshift;
    bool r_sign = false;
    if ( _r != nullptr ) {
        r = *_r;
        if ( do_reduce ) reduce_arg( r, r_lshift, r_sign );
        r = mul( r, circular_rotation_one_over_gain(), true );  
    }

    T zz;
    circular_rotation( r, impl->zero, x, co, si, zz );
    if ( do_reduce ) {
        //-----------------------------------------------------
        // If did_minus_pi_div_4 is true, then we need to perform this
        // modification for sin and cos:
        //
        // sin(x+PI/4) = sqrt(2)/2 * ( sin(x) + cos(x) )
        // cos(x+PI/4) = sqrt(2)/2 * ( cos(x) - sin(x) )
        //-----------------------------------------------------
        if ( did_minus_pi_div_4 ) {
            T si_new = mul( impl->sqrt2_div_2, si+co, true );
            T co_new = mul( impl->sqrt2_div_2, co-si, true );
            si = si_new;
            co = co_new;
        }

        //-----------------------------------------------------
        // Next, make adjustments for the quadrant and r multiply.
        //-----------------------------------------------------
        if ( quadrant&1 ) {
            T tmp = co;
            co = si;
            si = tmp;
        }
        if ( _r != nullptr ) {
            if ( need_si ) si <<= r_lshift;
            if ( need_co ) co <<= r_lshift;
        }
        if ( need_si && (r_sign ^ x_sign ^ (quadrant >= 2)) )                  si = -si;
        if ( need_co && (r_sign ^          (quadrant == 1) || quadrant == 2) ) co = -co;
    }
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::tan( const T& x ) const
{ 
    _log1( tan, x );
    T si, co;
    sincos( x, si, co );
    return div( si, co );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::asin( const T& x ) const
{ 
    _log1( asin, x );
    cassert( x >= -impl->one && x <= impl->one, "asin x must be between -1 and 1" );
    return atan2( x, normh( impl->one, x ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::acos( const T& x ) const
{ 
    _log1( acos, x );
    cassert( x >= -impl->one && x <= impl->one, "acos x must be between -1 and 1" );
    return atan2( normh( impl->one, x ), x, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::atan( const T& x ) const
{ 
    cassert( x >= -impl->one && x <= impl->one, "atan x must be between -1 and 1" );
    return atan2( x, impl->one, impl->do_reduce, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::atan2( const T& y, const T& x ) const
{ 
    return atan2( y, x, impl->do_reduce );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::atan2( const T& _y, const T& _x, bool do_reduce, bool x_is_one, T * r ) const
{ 
    _log2( atan2, _y, _x );
    T y = _y;
    T x = _x;

    //-----------------------------------------------------
    // Identities:
    //     atan2(y,x)       = undefined                             if x == 0 && y == 0
    //     atan2(y,x)       = PI                                    if x <  0 && y == 0
    //     atan2(y,x)       = 2*atan(y / (sqrt(x^2 + y^2) + x))     if x >  0    
    //     atan2(y,x)       = 2*atan((sqrt(x^2 + y^2) + |x|) / y)   if x <= 0 && y != 0
    //     atan(1/x)        = PI/2 - atan(x)                        if x >  0
    // Strategy:
    //     Use reduce_atan2_args() to reduce y and x and get y_sign and x_sign.
    //     Return PI if we're done.
    //     Do 2*atan( y / (norm(x, y) + x) )    if x > 0
    //     Do 2*atan( (norm(x, y) + x) / y )    if x <= 0
    //     When using atan2 for the latter, if the numerator is larger than
    //     the denominator, then use PI/2 - atan(x/y)
    //-----------------------------------------------------
    if ( debug ) std::cout << "atan2 begin: y=" << to_flt(y) << " x=" << to_flt(x) << " do_reduce=" << do_reduce << " x_is_one=" << x_is_one << "\n";
    cassert( (x != 0 || y != 0), "atan2: x or y needs to be non-zero for result to be defined" );
    T xx, yy, zz;
    if ( r != nullptr ) *r = norm( _x, _y, do_reduce );  // optimize this later with below norm() of reduced x,y
    if ( do_reduce ) {
        bool y_sign;
        bool x_sign;
        bool is_pi;
        bool swapped;
        reduce_atan2_args( y, x, y_sign, x_sign, swapped, is_pi );
        if ( is_pi ) {
            if ( debug ) std::cout << "atan2 end: y=" << to_flt(_y) << " x=" << to_flt(_x) << " do_reduce=" << do_reduce << 
                                      " x_is_one=" << x_is_one << 
                                      " zz=PI" << " r=" << ((r != nullptr) ? to_flt(*r) : to_flt(impl->zero)) << "\n";
            return impl->pi;
        }

        const T norm_plus_x = norm( x, y, true ) + x;
        if ( debug ) std::cout << "atan2 cordic begin: y=y=" << to_flt(y) << " x=norm_plus_x=" << to_flt(norm_plus_x) << " swapped=" << swapped << "\n";
        circular_vectoring( norm_plus_x, y, impl->zero, xx, yy, zz );
        zz <<= 1;
        if ( swapped ) zz = impl->pi_div_2 - zz;
        if ( y_sign ) zz = -zz;
        if ( debug ) std::cout << "atan2 cordic end: zz=" << to_flt(zz) << "\n";
    } else {
        circular_vectoring( x, y, impl->zero, xx, yy, zz );
    }
    if ( debug ) std::cout << "atan2 end: y=" << to_flt(_y) << " x=" << to_flt(_x) << " do_reduce=" << do_reduce << " x_is_one=" << x_is_one << 
                              " zz=" << to_flt(zz) << " r=" << ((r != nullptr) ? to_flt(*r) : to_flt(impl->zero)) << "\n";
    return zz;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::polar_to_rect( const T& r, const T& a, T& x, T& y ) const
{
    _log4( polar_to_rect, r, a, x, y );
    if ( debug ) std::cout << "polar_to_rect begin: r=" << to_flt(r) << " a=" << to_flt(a) << " do_reduce=" << impl->do_reduce << "\n";
    sincos( a, y, x, true, true, true, &r );
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::rect_to_polar( const T& x, const T& y, T& r, T& a ) const
{
    _log4( rect_to_polar, x, y, r, a );
    if ( debug ) std::cout << "rect_to_polar begin: x=" << to_flt(x) << " y=" << to_flt(y) << " do_reduce=" << impl->do_reduce << "\n";
    a = atan2( y, x, true, false, &r );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::norm( const T& x, const T& y ) const
{
    return norm( x, y, impl->do_reduce );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hypot( const T& x, const T& y ) const
{
    return norm( x, y, impl->do_reduce );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::norm( const T& _x, const T& _y, bool do_reduce ) const
{
    _log2( norm, _x, _y );
    T x = _x;
    T y = _y;
    if ( debug ) std::cout << "norm begin: x=" << to_flt(x) << " y=" << to_flt(y) << " do_reduce=" << do_reduce << "\n";
    int32_t ls;
    bool    swapped;  // unused
    if ( do_reduce ) reduce_norm_args( x, y, ls, swapped );

    T xx, yy, zz;
    circular_vectoring_xy( x, y, xx, yy );
    xx = mul( xx, circular_vectoring_one_over_gain(), true );
    if ( do_reduce ) xx = lshift( xx, ls );
    if ( debug ) std::cout << "norm end: x=" << to_flt(x) << " y=" << to_flt(y) << " do_reduce=" << do_reduce << " xx=" << to_flt(xx) << "\n";
    return xx;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::normh( const T& x, const T& y ) const
{
    //-----------------------------------------------------
    // Identities:
    //     sqrt(x^2 - y^2) = sqrt((x+y)(x-y))
    // Strategy:
    //     Try this easy way, though I suspect there will be issues.
    //-----------------------------------------------------
    _log2( normh, x, y );
    if ( debug ) std::cout << "normh begin: x=" << to_flt(x) << " y=" << to_flt(y) << " do_reduce=" << impl->do_reduce << "\n";
    cassert( x >= y, "normh x must be >= y" );
    return sqrt( mul( x+y, x-y ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sinh( const T& x, const T * r ) const
{ 
    T sih;
    T coh;
    sinhcosh( x, sih, coh, impl->do_reduce, true, false, r );
    return sih;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::cosh( const T& x, const T * r ) const
{ 
    T sih;
    T coh;
    sinhcosh( x, sih, coh, impl->do_reduce, false, true, r );
    return coh;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::sinhcosh( const T& x, T& sih, T& coh, const T * r ) const
{ 
    sinhcosh( x, sih, coh, impl->do_reduce, true, true, r );
}

template< typename T, typename FLT >
void Cordic<T,FLT>::sinhcosh( const T& _x, T& sih, T& coh, bool do_reduce, bool need_sih, bool need_coh, const T * _r ) const
{ 
    if ( _r != nullptr ) {
        _log4( sinhcosh, _x, sih, coh, *_r );
    } else {
        _log3( sinhcosh, _x, sih, coh );
    }

    //-----------------------------------------------------
    // Identities:
    //     sinh(-x)         = -sinh(x)
    //     sinh(x+y)        = sinh(x)*cosh(y) + cosh(x)*sinh(y)    
    //     cosh(-x)         = cosh(x)
    //     cosh(x+y)        = cosh(x)*cosh(y) - sinh(x)*sinh(y)     
    // Strategy:
    //     split abs(x) into i + f
    //     use LUT for sinh(i) and cosh(i)
    //     run cordic on f
    //     do the multiplications
    //     fix sign of sih
    //-----------------------------------------------------
    T x = _x;
    T sinh_i; 
    T cosh_i;
    bool sign;
    if ( do_reduce ) reduce_sinhcosh_arg( x, sinh_i, cosh_i, sign );  

    T r = hyperbolic_rotation_one_over_gain();
    int32_t r_lshift;
    bool r_sign = false;
    if ( _r != nullptr ) {
        r = *_r;
        if ( do_reduce ) reduce_arg( r, r_lshift, r_sign );
        r = mul( r, hyperbolic_rotation_one_over_gain(), false );  // should not need to reduce (I think)
    }

    T sinh_f;
    T cosh_f;
    T zz;
    hyperbolic_rotation( r, impl->zero, x, cosh_f, sinh_f, zz );
    if ( do_reduce ) {
        if ( need_sih ) {
            sih = mul( sinh_f, cosh_i, true ) + mul( cosh_f, sinh_i, true );
            if ( sign ) sih = -sih;
        }
        if ( need_coh ) {
            coh = mul( cosh_f, cosh_i, true ) + mul( sinh_f, sinh_i, true );
        }
    } else {
        sih = sinh_f;
        coh = cosh_f;
    }
}

template< typename T, typename FLT >
T Cordic<T,FLT>::tanh( const T& x ) const
{ 
    _log1( tanh, x );
    T sih, coh;
    sinhcosh( x, sih, coh );
    return div( sih, coh );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::asinh( const T& x ) const
{ 
    _log1( asinh, x );
    return log( x + norm( x, impl->one ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::acosh( const T& x ) const
{ 
    _log1( acosh, x );
    cassert( x >= impl->one, "acosh x must be >= 1" );
    return log( x + normh( x, impl->one ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::atanh( const T& x ) const
{ 
    return atanh2( x, impl->one, impl->do_reduce, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::atanh2( const T& y, const T& x ) const             
{ 
    return atanh2( y, x, impl->do_reduce );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::atanh2( const T& _y, const T& _x, bool do_reduce, bool x_is_one ) const             
{ 
    _log2( atanh2, _y, _x );
    T y = _y;
    T x = _x;

    //-----------------------------------------------------
    // Identities:
    //     atan(-x) = -atan(x)
    //     abs(y/x) must be between 0 and 1
    // Strategy:
    //     reduce y and x for division
    //     sum of their lshifts should be 0
    //-----------------------------------------------------
    int32_t y_lshift = 0;
    int32_t x_lshift = 0;
    bool    sign = false;
    if ( do_reduce ) {
        bool y_sign;
        bool x_sign = false;
        reduce_arg( y, y_lshift, y_sign );
        if ( !x_is_one ) reduce_arg( x, x_lshift, x_sign, true ); 
        sign = y_sign != x_sign;
    }
    int32_t lshift = x_lshift + y_lshift;
    cassert( lshift == 0, "atanh2: abs(y/x) must be between 0 and 1" );

    T xx, yy, zz;
    hyperbolic_vectoring( x, y, impl->zero, xx, yy, zz );
    if ( sign ) zz = -zz;
    return zz;
}

template< typename T, typename FLT >
void Cordic<T,FLT>::reduce_arg( T& x, int32_t& x_lshift, bool& sign, bool shift_x, bool normalize, bool for_sqrt ) const
{
    T x_orig = x;
    sign = x < 0;
    if ( sign ) x = -x;
    x_lshift = 0;
    T other = for_sqrt ? impl->half : impl->one;
    while( x > other ) 
    {
        x_lshift++;
        if ( shift_x ) {
            x >>= 1;
        } else {
            other <<= 1;
        }
    }
    while( normalize && x > 0 && x < other )
    {
        x_lshift--;
        if ( shift_x ) {
            x <<= 1;
        } else {
            other >>= 1;
        }
    }
    if ( debug ) std::cout << "reduce_arg: x_orig=" << to_flt(x_orig) << " x_reduced=" << to_flt(x) << " x_lshift=" << x_lshift << "\n"; 
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_mul_args( T& x, T& y, int32_t& x_lshift, int32_t& y_lshift, bool& sign ) const
{
    if ( debug ) std::cout << "reduce_mul_args: x_orig=" << to_flt(x) << " y_orig=" << to_flt(y) << "\n";
    bool x_sign;
    bool y_sign;
    reduce_arg( x, x_lshift, x_sign );
    reduce_arg( y, y_lshift, y_sign );
    sign = x_sign ^ y_sign;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_div_args( T& y, T& x, int32_t& y_lshift, int32_t& x_lshift, bool& sign ) const
{
    if ( debug ) std::cout << "reduce_div_args: x_orig=" << to_flt(x) << " y_orig=" << to_flt(y) << "\n";
    bool x_sign;
    bool y_sign;
    reduce_arg( y, y_lshift, y_sign );
    reduce_arg( x, x_lshift, x_sign, true, true );
    sign = x_sign ^ y_sign;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_sqrt_arg( T& x, int32_t& lshift ) const
{
    //-----------------------------------------------------
    // Identities:
    //     sqrt(x) = sqrt((x+1)^2 - (x-1)^2) / 2          
    // Strategy:
    //     Factor x=p*s where p is a power-of-2 
    //     where log2(p) is even >= 2 and s is between 0.5..1.
    //     sqrt(x) = sqrt(p) * sqrt(s) = 2^(log2(p)/2) * sqrt((s+1)^2 - (s-1)^2)/2
    //     Use hyperbolic_vectoring() for sqrt((s+1)^2 - (s-1)^2).
    //     Then lshift that by log2(p)/2 - 1.
    //-----------------------------------------------------
    cassert( x >= 0, "reduce_sqrt_arg x must be non-negative" );
    T x_orig = x;
    bool sign;
    reduce_arg( x, lshift, sign, true, true, true );
    if ( lshift & 1 ) {
        // make it even
        lshift++;
        x >>= 1;
    }
    lshift = lshift/2 - 1;
    if ( debug ) std::cout << "reduce_sqrt_arg: x_orig=" << to_flt(x_orig) << " x_reduced=s=" << to_flt(x) << " lshift=" << to_flt(lshift) << "\n";
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_exp_arg( FLT b, T& x, T& factor ) const
{
    //-----------------------------------------------------
    // Identities:
    //     Assume: x = i + f  (integer plus fraction)
    //     exp(i+f) = exp(i) * exp(f)
    //     pow(b,x) = log(b) * exp(x) = [log(b)*exp(i)] * exp(f)
    //
    // Strategy:
    //     Find i such that x-i is in -1 .. 1.
    //     Because x can be negative, so can i.
    //     exp(i) comes from a pre-built LUT kept in FLT
    //     so we can multiply it by log(b) before converting to type T and
    //     then multiplying by exp(f) in the caller.
    //-----------------------------------------------------
    const T TWO = impl->one << 1;
    const T MININT = -impl->maxint - 1;
    T x_orig = x;
    if ( debug ) std::cout << "reduce_exp_arg: b=" << b << " x_orig=" << to_flt(x_orig) << "\n";
    const FLT * factors_f = impl->reduce_exp_factor;
    T   i         = x >> impl->frac_w;  // can be + or -
    T   index     = -MININT + i;
    FLT factor_f  = std::log(b) * factors_f[index];   // could build per-b factors_f[] LUT with multiply already done
    factor        = to_t( factor_f );
    x            -= i << impl->frac_w;
    if ( debug ) std::cout << "reduce_exp_arg: b=" << b << " x_orig=" << to_flt(x_orig) << 
                              " i=" << to_rstring(i) << " index=" << to_rstring(index) << " log(b)*exp(i)=" << to_flt(factor) << 
                              " f=x_reduced=" << to_flt(x) << "\n"; 
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_log_arg( T& x, T& addend ) const 
{
    //-----------------------------------------------------
    // log(ab) = log(a) + log(b)
    // 
    // Normalize x so that it's in 1.00 .. 2.00.
    // Then addend = log(1 << lshift).
    //-----------------------------------------------------
    cassert( x >= 0, "log argument may not be negative for fixed-point numbers" );
    T x_orig = x;
    int32_t x_lshift;
    bool x_sign;
    reduce_arg( x, x_lshift, x_sign, true, true, true );
    const T * addends = impl->reduce_log_addend;
    addend = addends[impl->frac_w+x_lshift];
    if ( debug ) std::cout << "reduce_log_arg: x_orig=" << to_flt(x_orig) << " x_reduced=" << to_flt(x) << " addend=" << to_flt(addend) << "\n"; 
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_atan2_args( T& y, T& x, bool& y_sign, bool& x_sign, bool& swapped, bool& is_pi ) const
{
    //-----------------------------------------------------
    // Identities:
    //     atan2(y,x)       = undefined                             if x == 0 && y == 0
    //     atan2(y,x)       = PI                                    if x <  0 && y == 0
    //     atan2(y,x)       = 2*atan(y / (sqrt(x^2 + y^2) + x))     if x >  0    
    //     atan2(y,x)       = 2*atan((sqrt(x^2 + y^2) + |x|) / y)   if x <= 0 && y != 0
    // Strategy:
    //     Use reduce_norm_arg().
    //-----------------------------------------------------
    const T y_orig = y;
    const T x_orig = x;
    cassert( (x != 0 || y != 0), "atan2: x or y needs to be non-zero for result to be defined" );

    x_sign = x < 0;
    y_sign = y < 0;
    if ( x_sign && y == 0 ) {
        // PI
        is_pi = true;
        if ( debug ) std::cout << "reduce_atan2_args: xy_orig=[" << to_flt(x_orig) << "," << to_flt(y_orig) << "]" << " PI returned\n";
    } else {
        is_pi = false;
        int32_t lshift;
        reduce_norm_args( x, y, lshift, swapped );
        if ( debug ) std::cout << "reduce_atan2_args: xy_orig=[" << to_flt(x_orig) << "," << to_flt(y_orig) << "]" << 
                                  " xy_reduced=[" << to_flt(x) << "," << to_flt(y) << "] y_sign=" << y_sign << " x_sign=" << x_sign << 
                                  " lshift=" << lshift << " swapped=" << swapped << "\n";
    }

}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_norm_args( T& x, T& y, int32_t& lshift, bool& swapped ) const
{
    //-----------------------------------------------------
    // Must shift both x and y by max( x_lshift, y_lshift ).
    // If x or y is negative, it's fine to negate them
    // because squaring them anyway.
    //-----------------------------------------------------
    if ( x < 0 ) x = -x;
    if ( y < 0 ) y = -y;
    T x_orig = x;
    T y_orig = y;
    int32_t x_lshift;
    int32_t y_lshift;
    bool x_sign;
    bool y_sign;
    reduce_arg( x, x_lshift, x_sign, false );
    reduce_arg( y, y_lshift, y_sign, false );   
    lshift = (x_lshift > y_lshift) ? x_lshift : y_lshift;
    x >>= lshift;
    y >>= lshift;
    swapped = x < y;
    if ( swapped ) {
        T tmp = x;
        x = y;
        y = tmp;
    }
    if ( debug ) std::cout << "reduce_norm_args: xy_orig=[" << to_flt(x_orig) << "," << to_flt(y_orig) << "]" << 
                                               " xy_reduced=[" << to_flt(x) << "," << to_flt(y) << "] lshift=" << lshift << " swapped=" << swapped << "\n"; 
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_sincos_arg( T& a, uint32_t& quad, bool& sign, bool& did_minus_pi_div_4 ) const
{
    //-----------------------------------------------------
    // Compute a * 4/PI, then take integer part of that.
    // Subtract int_part*PI/4 from a.
    // If we ended up subtracting an odd PI/4, then set the flag.
    //-----------------------------------------------------
    const T a_orig = a;
    sign = a < 0;
    if ( sign ) a = -a;
    const T m = mul( a, impl->four_div_pi );
    const T i = m >> impl->frac_w;
    const T s = i * impl->pi_div_4;
    a        -= s;
    quad      = (i >> 1) & 3;
    did_minus_pi_div_4 = i & 1;
    if ( debug ) std::cout << "reduce_sincos_arg: a_orig=" << to_flt(a_orig) << " m=" << to_flt(m) << " i=" << to_rstring(i) << 
                              " subtract=" << to_flt(s) << " a_reduced=" << to_flt(a) << 
                              " quadrant=" << quad << " did_minus_pi_div_4=" << did_minus_pi_div_4 << "\n"; 

}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_sinhcosh_arg( T& x, T& sinh_i, T& cosh_i, bool& sign ) const
{
    //-----------------------------------------------------
    // Identities:
    //     sinh(-x)         = -sinh(x)
    //     sinh(x+y)        = sinh(x)*cosh(y) + cosh(x)*sinh(y)    
    //     cosh(-x)         = cosh(x)
    //     cosh(x+y)        = cosh(x)*cosh(y) - sinh(x)*sinh(y)     
    // Strategy:
    //     split abs(x) into i + f
    //     use LUT for sinh(i) and cosh(i)
    //     run cordic on f
    //     do the multiplications
    //     fix sign of sih
    //-----------------------------------------------------
    T x_orig = x;
    sign = x < 0;
    if ( sign ) x = -x;

    const T MASK = (T(1) << (impl->int_w+2)) - T(1);  // include 0.25 bit of fraction
    T i = (x >> (impl->frac_w-2)) & MASK;
    x   = x & ((T(1) << (impl->frac_w-2))-T(1));
    const T *    sinh_i_vals   = impl->reduce_sinhcosh_sinh_i;
    const T *    cosh_i_vals   = impl->reduce_sinhcosh_cosh_i;
    const bool * sinh_i_oflows = impl->reduce_sinhcosh_sinh_i_oflow;
    const bool * cosh_i_oflows = impl->reduce_sinhcosh_cosh_i_oflow;
    cassert( !sinh_i_oflows[i], "reduce_sinhcosh_arg x will cause an overflow for sinh" );
    cassert( !cosh_i_oflows[i], "reduce_sinhcosh_arg x will cause an overflow for cosh" );
    sinh_i = sinh_i_vals[i];
    cosh_i = cosh_i_vals[i];
    if ( debug ) std::cout << "reduce_sinhcosh_arg: x_orig=" << to_flt(x_orig) << " sinh_i[" << to_rstring(i) << 
                              "]=" << to_flt(sinh_i) << " coshh_i[" << to_rstring(i) << "]=" << to_flt(cosh_i) << 
                              " x_reduced=" << to_flt(x) << "\n";
}

template class Cordic<int64_t, double>;

#endif
