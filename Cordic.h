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
#ifdef issubnormal
#undef issubnormal
#endif
#include <cfenv>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <cstring>

#include "Logger.h"

#ifdef DEBUG_LEVEL
static constexpr uint32_t debug = DEBUG_LEVEL;
#else
static constexpr uint32_t debug = 0;
#endif

#ifdef NO_ASSERT
static constexpr bool do_asserts = false;
#else
static constexpr bool do_asserts = true;
#endif

#define cassert(expr, msg) if ( do_asserts && !(expr) ) \
                { std::cout << "ERROR: assertion failure: " << (msg) << " at " << __FILE__ << ":" << __LINE__ << "\n"; exit( 1 ); }

// extended rounding modes beyond FE_TONEAREST etc.
//
static constexpr int FE_NOROUND      = 0x1000;  // perform no rounding at all; leave guard bits alone
static constexpr int FE_AWAYFROMZERO = 0x2000;  // signbit(x) ? trunc(x) : ceil(x);   and clear guard bits

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
    // int_exp_w      = fixed-point integer width  OR floating-point exponent width
    // frac_w         = fixed-point fraction width OR floating-point mantissa width
    // is_float       = true=floating-point, false=fixed-point
    // guard_w        = fixed-point guard bits width (default is calculated as log2(frac_w)
    // 1+int_exp_w+frac_w+guard_w must fit in T
    //
    // So the most significant bit is the sign, followed by int_exp_w bits of integer/exponent, followed by frac_w+guard_w bits of fraction.
    //-----------------------------------------------------
    Cordic( uint32_t int_exp_w,                 // fixed-point integer width  OR floating-point exponent width
            uint32_t frac_w,                    // fixed-point fraction width OR floating-point mantissa width
            bool     is_float=true,             // true=floating-point, false=fixed-point
            uint32_t guard_w=-1,                // number of guard bits used for CORDIC proper (-1 == default == log2(frac_w))
            uint32_t n=-1 );                    // number of iterations used for CORDIC proper (-1 == default == frac_w)
    ~Cordic();
    
    void log_constructed( void );               // for when Cordic constructed before logging enabled

    //-----------------------------------------------------
    // Construction
    //-----------------------------------------------------
    T    make_fixed( bool sign, const T& i, const T& f ) const;         // fixed-point    value with sign, integer  part i, and fractional part f
    T    make_float( bool sign, const T& e, const T& f ) const;         // floating-point value using sign, exponent part 3, and mantissa part f

    //-----------------------------------------------------
    // Explicit Conversions
    //-----------------------------------------------------
    T           to_t( FLT x, bool is_final=false, bool to_fixed=false ) const; // FLT to T encoded value; lsb is never rounded
    FLT         to_flt( const T& x ) const;                                    // T encoded value to FLT
    std::string to_string(  const T& x, bool from_fixed=false ) const;         // T to std::string in decimal floating-point format
    std::string to_rstring( const T& x, bool from_fixed=false ) const;         // T to std::string in raw decimal integer format 
    std::string to_bstring( const T& x, bool from_fixed=false ) const;         // T to std::string in binary format, like "1 001 101101011010"

    //-----------------------------------------------------
    // Constants (T ones are never rounded, so call rfrac() if you want them rounded)
    //-----------------------------------------------------
    bool     is_float( void ) const;                    // is_float from above
    uint32_t int_w( void ) const;                       // int_w   from above
    uint32_t exp_w( void ) const;                       // exp_w   from above
    uint32_t frac_w( void ) const;                      // frac_w  from above
    uint32_t guard_w( void ) const;                     // guard_w from above
    uint32_t w( void ) const;                           // 1 + int_exp_w + frac_w + guard_w (i.e., overall width)
    uint32_t n( void ) const;                           // n       from above
    T maxint( void ) const;                             // largest positive integer (just integer part)

    T max( void ) const;                                // encoded maximum positive value 
    T min( void ) const;                                // encoded minimum positive value
    T denorm_min( void ) const;                         // fixed-point: min()
    T lowest( void ) const;                             // encoded most negative value
    T epsilon( void ) const;                            // difference between 1 and least value greater than 1
    T round_error( void ) const;                        // maximum rounding error
    T zero( void ) const;                               // encoded 0.0
    T one( void ) const;                                // encoded 1.0
    T neg_one( void ) const;                            // encoded -1.0
    T two( void ) const;                                // encoded 2.0
    T half( void ) const;                               // encoded 0.5
    T quarter( void ) const;                            // encoded 0.25
    T sqrt2( void ) const;                              // encoded sqrt(2)
    T sqrt2_div_2( void ) const;                        // encoded sqrt(2)/2
    T pi( void ) const;                                 // encoded PI
    T tau( void ) const;                                // encoded 2*PI
    T pi_div_2( void ) const;                           // encoded PI/2
    T pi_div_4( void ) const;                           // encoded PI/4
    T one_div_pi( void ) const;                         // encoded 1/PI
    T two_div_pi( void ) const;                         // encoded 2/PI
    T four_div_pi( void ) const;                        // encoded PI/4
    T e( void ) const;                                  // encoded natural exponent
    T nan( const char * arg ) const;                    // to_t( std::strtod( arg, nullptr )
    T quiet_NaN( void ) const;                          // to_t( std::numeric_limits<FLT>::quiet_NaN() )
    T signaling_NaN( void ) const;                      // to_t( std::numeric_limits<FLT>::signaling_NaN() )
    T infinity( void ) const;                           // to_t( std::numeric_limits<FLT>::infinity() )
    T ninfinity( void ) const;                          // to_t( -std::numeric_limits<FLT>::infinity() )

    //-----------------------------------------------------
    // Well-Known Math Functions Implemented Using CORDIC
    //
    // (2) means requires 2 applications of a CORDIC algorithm.              functionality
    //-----------------------------------------------------               ---------------------------

    // deconstruction
    bool signbit( const T& x ) const;                                   // x < 0
    T    frexp( const T& x, int * e ) const;                            // return normalized fraction and exponent
    T    modf( const T& x, T * i ) const;                               // decompose into fraction and integer (both with sign of x), still encoded 
    int  ilogb( const T& x ) const;                                     // unbiased exponent as int
    T    logb( const T& x ) const;                                      // unbiased exponent, still encoded
    int  fpclassify( const T& x ) const;                                // FP_ZERO, FP_NORMAL, FP_SUBNORMAL, FP_INFINITE, FP_NAN
    bool iszero( const T& x ) const;                                    // fpclassify() is FP_ZERO
    bool isfinite( const T& x ) const;                                  // fpclassify() is FP_ZERO or FP_NORMAL or FP_SUBNORMAL
    bool isinf( const T& x ) const;                                     // fpclassify() is FP_INFINITE
    bool isnan( const T& x ) const;                                     // fpclassify() is FP_NAN
    bool isnormal( const T& x ) const;                                  // fpclassify() is FP_NORMAL
    bool issubnormal( const T& x ) const;                               // fpclassify() is FP_SUBNORMAL


    // rounding
    int  fesetround( int round );                                       // set rounding mode to FE_{DOWNWARD,UPWARD,TOWARDZDERO,TONEAREST}; return 0
    int  fegetround( void ) const;                                      // returns current rounding mode
    T    nextafter( const T& from, const T& to ) const;                 // (from == to) ?      to  : (from +/- min toward to)
    T    nexttoward( const T& from, long double to ) const;             // (from == to) ? to_t(to) : (from +/- min toward to)
    T    floor( const T& x ) const;                                     // largest  integral value <= x
    T    ceil( const T& x ) const;                                      // smallest integral value >= x
    T    trunc( const T& x ) const;                                     // nearest  integral value toward 0
    T    extend( const T& x ) const;                                    // nearest  integral value away from 0
    T    round( const T& x ) const;                                     // nearest  integral value; halfway cases away from 0 
    long lround( const T& x ) const;                                    // same as round() except returns just the raw integer part as long 
    long long llround( const T& x ) const;                              // same as round() except returns just the raw integer part as long long
    T    iround( const T& x ) const;                                    // same as round() except returns just the raw integer part as T
    T    rint( const T& x, int rmode=-1 ) const;                        // nearest integral value according to rounding mode (-1 means fegetround()):
                                                                        //    FE_NOROUND:      x                (extension)
                                                                        //    FE_DOWNWARD:     floor(x)
                                                                        //    FE_UPWARD:       ceil(x)
                                                                        //    FE_TOWWARDZERO:  trunc(x)
                                                                        //    FE_AWAYFROMZERO: extend(x)        (extension)
                                                                        //    FE_TONEAREST:    round(x)
    long lrint( const T& x ) const;                                     // same as rint() except returns just the raw integer part as long
    long long llrint( const T& x ) const;                               // same as rint() except returns just the raw integer part as long long
    T    irint( const T& x ) const;                                     // same as rint() except returns just the raw integer part as T
    T    nearbyint( const T& x ) const;                                 // same as rint() but never raises FE_INEXACT
    T    floorfrac( const T& x ) const;                                 // largest  value <= x                       (and clear guard bits)
    T    ceilfrac( const T& x ) const;                                  // smallest value >= x                       (and clear guard bits)
    T    truncfrac( const T& x ) const;                                 // nearest  value toward 0                   (and clear guard bits)
    T    extendfrac( const T& x ) const;                                // nearest  value away from 0                (and clear guard bits)
    T    roundfrac( const T& x ) const;                                 // nearest  value; halfway cases away from 0 (and clear guard bits)
    T    rfrac( const T& x, int rmode=-1 ) const;                       // use guard bits to round value according to rounding mode (-1 means fegetround()):
                                                                        //    FE_NOROUND:      x                (extension)
                                                                        //    FE_DOWNWARD:     floorfrac(x)
                                                                        //    FE_UPWARD:       ceilfrac(x)
                                                                        //    FE_TOWWARDZERO:  truncfrac(x)
                                                                        //    FE_AWAYFROMZERO: extendfrac(x)    (extension)
                                                                        //    FE_TONEAREST:    roundfrac(x)

    // basic arithmetic
    T    abs( const T& x ) const;                                       // |x|
    T    neg( const T& x ) const;                                       // -x
    T    copysign( const T& x, const T& y ) const;                      // |x| with sign of y
    T    add( const T& x, const T& y ) const;                           // x+y 
    T    sub( const T& x, const T& y ) const;                           // x-y 
    T    scalbn( const T& x, int y ) const;                             // x * 2^y    (i.e., left-shift)
    T    scalbnn( const T& x, int y ) const;                            // x * 2^(-y) (i.e., right-shift)
    T    ldexp( const T& x, int y ) const;                              // x * 2^y    (same as scalbn)
    T    fma( const T& x, const T& y, const T& addend ) const;          // x*y + addend
    T    mul( const T& x, const T& y ) const;                           // x*y 
    T    mulc( const T& x, const T& c ) const;                          // x*c where c is known to be a constant
    T    sqr( const T& x ) const;                                       // x*x
    T    fda( const T& y, const T& x, const T& addend ) const;          // y/x + addend
    T    div( const T& y, const T& x ) const;                           // y/x
    T    remainder( const T& y, const T& x ) const;                     // IEEE 754-style remainder: y - n*x (where n is nearest int)
    T    fmod( const T& y, const T& x ) const;                          // rem = remainder( |y|, |x| ); if (rem < 0) rem += x; rem = copysign(rem, x)
    T    remquo( const T& y, const T& x, int * quo ) const;             // remainder(), and *quo receives sign and at least 3 LSBs of y/x
    T    rcp( const T& x ) const;                                       // 1/x

    // comparisons
    int  compare( const T& _x, const T& _y ) const;                     // -2 if one NaN, -4 if two NaNs, -1 if x < y, 0 if x == y, 1 if x > y
    bool isgreater( const T& x, const T& y ) const;                     // x > y
    bool isgreaterequal( const T& x, const T& y ) const;                // x >= y
    bool isless( const T& x, const T& y ) const;                        // x > y
    bool islessequal( const T& x, const T& y ) const;                   // x >= y
    bool islessgreater( const T& x, const T& y ) const;                 // x != y (but returns false if at least one is NaN)
    bool isunordered( const T& x, const T& y ) const;                   // returns true if x is a NaN OR y is a NaN
    bool isunequal( const T& x, const T& y ) const;                     // x != y (returns true if only one is NaN)
    bool isequal( const T& x, const T& y ) const;                       // x == y (returns true if both are NaN)
    T    fdim( const T& x, const T& y ) const;                          // (x >= y) ? (x-y) : 0
    T    fmax( const T& x, const T& y ) const;                          // max(x, y)
    T    fmin( const T& x, const T& y ) const;                          // min(x, y)

    // elementary functions
    T    sqrt( const T& x ) const;                                      // hypoth( x+1, x-1 ) / 2
    T    rsqrt( const T& x ) const;                                     // 1.0 / sqrt( x )
    T    rsqrt_orig( const T& x ) const;                                // x^(-1/2) = exp(log(x)/-2)
    T    cbrt( const T& x ) const;                                      // x^(1/3)  = exp(log(x)/3)
    T    rcbrt( const T& x ) const;                                     // x^(-1/3) = exp(log(x)/-3)

    T    exp( const T& x ) const;                                       // e^x
    T    expm1( const T& x ) const;                                     // e^x - 1
    T    expc( const FLT& b, const T& x ) const;                        // b^x  = exp(x * log(b))  b=const     (2)
    T    exp2( const T& x ) const;                                      // 2^x
    T    exp10( const T& x ) const;                                     // 10^x
    T    pow( const T& b, const T& x ) const;                           // b^x  = exp(x * log(b))              (3)
    T    log( const T& x ) const;                                       // 2*atanh2(x-1, x+1)    
    T    log( const T& x, const T& b ) const;                           // log(x)/log(b)                (3)
    T    log1p( const T& x ) const;                                     // 2*atanh2(x, x+2) = log(x+1)
    T    logc( const T& x, const FLT& b ) const;                        // log(x)/log(b)    b=const     (2)
    T    log2( const T& x ) const;                                      // log(x)/log(2)                (2)
    T    log10( const T& x ) const;                                     // log(x)/log(10)               (2)

    T    deg2rad( const T& x ) const;                                   // x * PI / 180
    T    rad2deg( const T& x ) const;                                   // x * 180 / PI
    T    sin( const T& x, const T * r=nullptr  ) const;                 // r*sin(x)                   (default r is 1)
    T    sinpi( const T& x, const T * r=nullptr ) const;                // r*sin(x*PI)                (defualt r is 1
    T    cos( const T& x, const T * r=nullptr ) const;                  // r*cos(x)                   (default r is 1)
    T    cospi( const T& x, const T * r=nullptr ) const;                // r*cos(x*PI)                (defualt r is 1
    void sincos( const T& x, T& si, T& co, const T * r=nullptr ) const; // si=r*sin(x), co=r*cos(x)   (default r is 1)
    void sinpicospi( const T& x, T& si, T& co, const T * r=nullptr ) const;// si=r*sin(x*PI), co=r*cos(x*PI) (default r is 1)
    T    tan( const T& x ) const;                                       // sin(x) / cos(x)              (2)
    T    tanpi( const T& x ) const;                                     // sin(x*PI) / cos(x*PI)        (2)
    T    asin( const T& x ) const;                                      // atan2(x, sqrt(1 - x^2))      (2)
    T    acos( const T& x ) const;                                      // atan2(sqrt(1 - x^2), x)      (2)
    T    atan( const T& x ) const;                                      // atan(x)
    T    atan2( const T& y, const T& x ) const;                         // atan2(y, x)                  

    void polar_to_rect( const T& r, const T& a, T& x, T& y ) const;     // sincos(a, x, y, &r)  
    void rect_to_polar( const T& x, const T& y, T& r, T& a ) const;     // r=sqrt(x^2 + y^2), a=atan2(y, x)
    T    hypot( const T& x, const T& y ) const;                         // sqrt(x^2 + y^2)  (Euclidean hypot)
    T    hypoth( const T& x, const T& y ) const;                        // sqrt(x^2 - y^2)  ("hyperbolic hypot")

    T    sinh( const T& x, const T * r=nullptr ) const;                 // r*sinh(x), also r*(e^x - e^-x)/2  (default r is 1)
    T    cosh( const T& x, const T * r=nullptr ) const;                 // r*cosh(x), also r*(e^x + e^-x)/2  (default r is 1)
    void sinhcosh( const T& x, T& sih, T& coh, const T * r=nullptr ) const;// sih=r*sinh(x), coh=r*cosh(x)    (default r is 1)
    T    tanh( const T& x ) const;                                      // sinh(x) / cosh(x)            (2)
    T    asinh( const T& x ) const;                                     // log(x + sqrt(x^2 + 1))       (2)
    T    acosh( const T& x ) const;                                     // log(x + sqrt(x^2 - 1))       (2)
    T    atanh( const T& x ) const;                                     // atanh(x)
    T    atanh2( const T& y, const T& x ) const;                        // atanh(y/x)

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
    // exp2(x)          = 2^x = exp(log(2) * x) 
    // exp(x)           = exp2(log2(e) * x)
    // exp2(i+f)        = exp2(i) * exp2(f)                     i = integer part, f = fractional remainder
    //                  = exp2(f) << i    
    //                  = exp(log(2)*f) << i
    // exp(x)           = sinh(x) + cosh(x)                     if x is already reduced, use hyperbolic CORDIC directly to get this sum
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
    // log(b^i)         = i*log(b)
    // log(2^i + f)     = i*log(2) + log(f)                     i=integer f=remainder
    // log(x/4)         = atanh((x-0.25)/(x+0.25)) 
    // log(x/y)         = log(x) - log(y)
    // log(x/y)         = 2*atanh((x-y)/(x+y))  
    // log(1+x)         = 2*atanh(x/(x+2))                      log1p()
    //
    // sin(-x)          = -sin(x)
    // sin(x)           = cos(x - pi/2)
    // sin(x)           = sin(i*pi/2 + f)                       where i is an integer and |f| <= pi/4
    // sin(x+y)         = sin(x)*cos(y) + cos(x)*sin(y)         
    // sin(x-y)         = sin(x)*cos(y) - cos(x)*sin(y)
    // sin(x+y)         = 2*sin(x)*cos(y) - sin(x-y)
    // sin(2x)          = 2*sin(x)*cos(x)                       but it's cheaper to just compute sin(x+x)
    // sin(2x)          = sin(x) * (2*cos(x))
    // sin^2(x)         = (1 - cos(2x)) / 2
    // sin(3x)          = sin(x) * (4*cos^2(x) - 1)
    // sin(nx)          = 2*sin((n-1)*x)*cos(x) - sin((n-2)*x)
    // sin(x+pi/4)      = sqrt(2)/2 * (sin(x) + cos(x))
    // sin(i*pi/2 + f)  = +/-sin(f) or +/-cos(f)
    // sin(x*pi)        = sin(2x * pi/2) = sin((i+f)*pi/2) = +/-sin(f*pi/2) or +/-cos(f*pi/2)
    // sin(x)           = Im(e^(ix)) = (e^(ix) + e^(-ix)) / 2
    // sin(ix)          = i*sinh(x)
    // sin(x)/x         = (1 + x*x*(1/6)) == 1) ? 1 : sin(x)/x
    // sin(x)*sin(y)    = (cos(x-y) - cos(x+y)) / 2
    //
    // cos(-x)          = cos(x)
    // cos(x)           = sin(x + pi/2)
    // cos(x+y)         = cos(x)*cos(y) - sin(x)*sin(y)
    // cos(x+pi/4)      = sqrt(2)/2 * (cos(x) - sin(x))
    // cos(i*pi/2 + f)  = +/-cos(f) or +/-sin(f)                i=integer, f=frac
    // cos(x*pi)        = cos(2x * pi/2) = cos((i+f)*pi/2) = +/-cos(f*pi/2) or +/-sin(f*pi/2)
    // cos(2x)          = -sin^2(x) + cos^2(x)
    // cos(2x)          = 2*cos^2(x) - 1                        but it's cheaper to just compute cos(x+x)
    // cos^2(x)         = (cos(2x) + 1) / 2                     this could come in handy
    // cos(3x)          = -3*cos(x) + 4*cos^3(x)                but it's cheaper to just compute cos(x+x+x)
    // cos(nx)          = 2*cos((n-1)*x)*cos(x) - cos((n-2)x)   
    // cos(x)           = Re(e^(ix)) = (e^(ix) - e^(-ix)) / 2i  i=sqrt(-1)
    // cos(ix)          = cosh(x)                               i=sqrt(-1)
    // 1-cos(x)         = 2*sin(x/2)^2
    // (1-cos(x))/x     = (1 + x*x) == 1) ? 0.5*x : cosf1(x)/x
    // cos(x)*cos(y)    = (cos(x-y) + cos(x+y)) / 2
    // sin(x)*cos(v)    = (sin(x+y) + sin(x-y)) / 2
    // cos(x)*sin(v)    = (sin(x+y) - sin(x-y)) / 2
    //
    // tan(-x)          = -tan(x)
    // tan(x+y)         = (tan(x) + tan(y)) / (1 - tan(x)*tan(y))
    // tan(x/2)         = sin(x) / (1 + cos(x))
    // tan(x/2 + PI/4)  = cos(x) / (1 - sin(x))
    // tan^2(x)         = (1 - cos(2x)) / (1 + cos(2x))
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
    // -3*pi/4          = atan(-1,-1)
    //
    // sinh(-x)         = -sinh(x)
    // sinh(x)          = (e^x - e^-x)/2
    //                  = (e^x - 1/e^x)/2
    // sinh(x+y)        = sinh(x)*cosh(y) + cosh(x)*sinh(y)             
    // sinh(2^i + f)    = sinh(2^i)*cosh(f) + cosh(2^i)*sinh(f)
    // sinh(2^i)        = (e^(2^i) - e^(-2^i))/2
    //                  = (2^(log2(e) << i) - 2^(-log2(e) << i))/2      let j = integer part of (log2(e) << i), g = fractional part
    //                  = ((2^g << j) - (2^(-g) >> j))/2
    //                  = (exp(log(2)*g) << (j-1)) - (1/exp(log(2)*g) >> (j+1))   
    //                    note: exp(log(2)*g) == hyperbolic_rotation()  
    //                    note: this whole method doesn't look more efficient than just using (e^x - 1/e^x)/2
    //
    // cosh(-x)         = cosh(x)
    // cosh(x)          = (e^x + e^-x)/2
    //                  = (e^x + 1/e^x)/2
    // cosh(x+y)        = cosh(x)*cosh(y) + sinh(x)*sinh(y)
    // cosh(2^i + f)    = cosh(2^i)*cosh(f) + sinh(2^i)*sinh(f)
    // cosh(2^i)        = (e^(2^i) + e^(-2^i))/2
    //                  = (2^(log2(e) << i) + 2^(-log2(e) << i))/2      let j = integer part of (log2(e) << i), g = fractional part
    //                  = ((2^g << j) + (2^(-g) >> j))/2
    //                  = (exp(log(2)*g) << (j-1)) + (1/exp(log(2)*g) >> (j+1))   
    //                    note: exp(log(2)*g) == hyperbolic_rotation() 
    //                    note: this whole method doesn't look more efficient than just using (e^x + 1/e^x)/2
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

    //-----------------------------------------------------
    // These version are used internally, but making them available publically.
    // In general, you should only call the earlier routines.
    //-----------------------------------------------------
    T    neg( const T& x, bool is_final ) const;                                      
    T    add( const T& x, const T& y, bool is_final ) const;                 
    T    sub( const T& x, const T& y, bool is_final ) const;                 
    T    scalbn( const T& x, int y, bool is_final ) const;                             
    T    fma_fda( bool is_fma, const T& x, const T& y, const T& addend, bool is_final ) const;
    T    mul( const T& x, const T& y, bool is_final ) const;                 
    T    mulc( const T& x, const T& c, bool is_final ) const;
    T    sqr( const T& x, bool is_final ) const;
    T    div( const T& y, const T& x, bool is_final ) const;                  
    T    sqrt( const T& x, bool is_final ) const;                              
    T    exp( const T& x, bool is_final, FLT b=M_E ) const;                              
    T    log( const T& x, bool is_final ) const;                             
    T    log1p( const T& x, bool is_final ) const;                          
    T    hypot( const T& x, const T& y, bool is_final ) const;          
    T    hypoth( const T& x, const T& y, bool is_final ) const;          
    T    atan2(  const T& y, const T& x, bool is_final, bool x_is_one, T * r ) const; 
    T    atanh2( const T& y, const T& x, bool is_final, bool x_is_one ) const; 
    void sincos( bool times_pi, const T& x, T& si, T& co, bool is_final, bool need_si, bool need_co, const T * r ) const;
    void sinhcosh( const T& x, T& sih, T& coh, bool is_final, bool need_sih, bool need_coh, const T * r ) const;

    //-----------------------------------------------------
    // Argument Range Reduction Routines
    //
    // These are called to reduce arguments from normal fixed/float form to 
    // reduced fixed form appropriate for core CORDIC routins.
    //-----------------------------------------------------
    enum class EXP_CLASS
    {
        ZERO,                           // can be used to save power
        NORMAL,                         // for fixed-point CORDIC, this is always the class
        SUBNORMAL,                      // these others apply to floating-point (is_float=true)
        INFINITE,
        NOT_A_NUMBER,
    };

    std::string to_str( const EXP_CLASS& c ) const
    {
        switch( c )
        {
            case EXP_CLASS::ZERO:               return "ZERO";
            case EXP_CLASS::NORMAL:             return "NORMAL";
            case EXP_CLASS::SUBNORMAL:          return "SUBNORMAL";
            case EXP_CLASS::INFINITE:           return "INFINITY";
            case EXP_CLASS::NOT_A_NUMBER:       return "NAN";
            default:                            return "<unknown>";
        }
    }

    FLT  _to_flt( const T& x, bool is_final=false, bool from_fixed=false, bool allow_debug=false ) const; // T encoded value to FLT
    EXP_CLASS classify( const T& x ) const;                                              // returns exp class only
    void deconstruct( T& x, EXP_CLASS& x_exp_class, int32_t& x_exp, bool& sign, bool allow_debug=true ) const;  // x will end up as fixed-point for _is_float=true
    void reconstruct( T& x, EXP_CLASS  x_exp_class, int32_t  x_exp, bool  sign ) const;  // x will end up as float value for _is_float=true

    void reduce_add_args( T& x, T& y, EXP_CLASS& x_exp_class, int32_t& x_exp, bool& x_sign, EXP_CLASS& y_exp_class, int32_t& y_exp, bool& y_sign ) const; 
    void reduce_mul_div_args( bool is_fma, T& x, T& y, EXP_CLASS& x_exp_class, int32_t& x_exp, EXP_CLASS& y_exp_class, int32_t& y_exp, bool& sign ) const; 
    void reduce_sqrt_arg( T& x, EXP_CLASS& x_exp_class, int32_t& x_exp, bool& x_sign ) const;
    void reduce_exp_arg( FLT b, T& x, int32_t& i, EXP_CLASS& x_exp_class, bool& x_sign ) const;
    void reduce_log_arg( T& x, EXP_CLASS& x_exp_class, bool& x_sign, T& addend ) const;
    void reduce_hypot_args( T& x, T& y, EXP_CLASS& exp_class, int32_t& exp, bool& swapped ) const;
    void reduce_sincos_arg( bool times_pi, T& a, uint32_t& quadrant, EXP_CLASS& exp_class, bool& sign, bool& did_minus_pi_div_4 ) const;

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

    void constructed( const T& x ) const;                         // so we can log creation of x
    void destructed( const T& x ) const;                          // so we can log destruction of x
    T&   assign( T& x, const T& y ) const;                        // x = y  (this exists so we can log assignments)
    T&   pop_value( T& x, const T& y ) const;                     // x = y  (where y is top of stack for logging)
    bool pop_bool( bool ) const;                                  // so we can log consumption of bool

    enum class OP
    {
        push_constant,
        to_flt,
        assign,
        pop_value,
        pop_bool,

        frexp,
        modf,
        ilogb,
        logb,

        nextafter,
        nexttoward,
        floor,
        ceil,
        trunc,
        extend,
        round,
        lround,
        llround,
        iround,
        rint,
        lrint,
        llrint,
        irint,
        nearbyint,
        floorfrac,
        ceilfrac,
        truncfrac,
        extendfrac,
        roundfrac,
        rfrac,

        abs,
        neg,
        copysign,
        add,
        sub,
        fma,
        mul,
        mulc,
        sqr,
        scalbn,
        ldexp,
        fda,
        div,
        rcp,
        sqrt,
        rsqrt,
        cbrt,
        rcbrt,

        iszero,
        isgreater,
        isgreaterequal,
        isless,
        islessequal,
        islessgreater,
        isunordered,
        isunequal,
        isequal,

        fdim,
        fmax,
        fmin,

        exp,
        expm1,
        expc,
        exp2,
        exp10,
        pow,
        log,
        log1p,
        logn,
        logc,
        log2,
        log10,

        deg2rad,
        rad2deg,
        sin,
        sinpi,
        cos,
        cospi,
        sincos,
        sinpicospi,
        tan,
        tanpi,

        asin,
        acos,
        atan,
        atan2,

        polar_to_rect,
        rect_to_polar,
        hypot,
        hypoth,

        sinh,
        cosh,
        sinhcosh,
        tanh,
        asinh,
        acosh,
        atanh,
        atanh2,

        // these are here for convenience, but have nothing to do with CORDIC
        sram_rd,
        sram_wr,
        dram_rd,
        dram_wr,
    };

    static constexpr uint32_t OP_cnt = uint32_t(OP::dram_wr) + 1;







//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
// 
//
//  ___                 _                           _        _   _
// |_ _|_ __ ___  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _| |_(_) ___  _ __
//  | || '_ ` _ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | __| |/ _ \| '_ \
//  | || | | | | | |_) | |  __/ | | | | |  __/ | | | || (_| | |_| | (_) | | | |
// |___|_| |_| |_| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_|
//               |_|
//
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------

private:
    bool                        _is_float;
    uint32_t                    _int_w;
    uint32_t                    _exp_w;
    uint32_t                    _frac_w;
    uint32_t                    _guard_w;
    uint32_t                    _frac_guard_w;  // sum of the two
    T                           _frac_guard_mask;
    T                           _guard_mask;
    uint32_t                    _exp_mask;
    int32_t                     _exp_bias;
    int32_t                     _exp_unbiased_min;
    int32_t                     _exp_unbiased_max;
    uint32_t                    _w;
    uint32_t                    _n;
    int                         _rounding_mode;

    T                           _quiet_NaN_fxd;
    T                           _maxint;
    T                           _max;
    T                           _min;
    T                           _min_fxd;
    T                           _lowest;
    T                           _zero;
    T                           _zero_fxd;
    T                           _neg_zero;
    T                           _quarter;
    T                           _third;
    T                           _neg_third;
    T                           _half;
    T                           _half_fxd;
    T                           _one;
    T                           _one_fxd;
    T                           _neg_one;
    T                           _two;
    T                           _two_fxd;
    T                           _four;
    T                           _four_fxd;
    T                           _sqrt2;
    T                           _sqrt2_div_2;
    T                           _pi;
    T                           _neg_pi;
    T                           _pi_fxd;
    T                           _tau;
    T                           _pi_div_2;
    T                           _pi_div_4;
    T                           _one_div_pi;
    T                           _two_div_pi;
    T                           _four_div_pi;
    T                           _e;
    T                           _log2;                                   
    T                           _log10;                                 

    T *                         _circular_atan_fxd;                      // circular atan values
    T                           _circular_rotation_gain_fxd;             // circular rotation gain
    T                           _circular_rotation_one_over_gain_fxd;    // circular rotation 1/gain
    T                           _circular_rotation_one_over_gain;        // circular rotation 1/gain
    T                           _circular_vectoring_gain_fxd;            // circular vectoring gain
    T                           _circular_vectoring_one_over_gain_fxd;   // circular vectoring 1/gain
    T                           _circular_vectoring_one_over_gain;       // circular vectoring 1/gain
    T                           _circular_angle_max_fxd;                 // circular vectoring |z0| max value

    T *                         _hyperbolic_atanh_fxd;                   // hyperbolic atanh values
    T                           _hyperbolic_rotation_gain_fxd;           // hyperbolic rotation gain
    T                           _hyperbolic_rotation_one_over_gain_fxd;  // hyperbolic rotation 1/gain
    T                           _hyperbolic_rotation_one_over_gain;      // hyperbolic rotation 1/gain
    T                           _hyperbolic_vectoring_gain_fxd;          // hyperbolic vectoring gain
    T                           _hyperbolic_vectoring_one_over_gain_fxd; // hyperbolic vectoring 1/gain
    T                           _hyperbolic_vectoring_one_over_gain;     // hyperbolic vectoring 1/gain
    T                           _hyperbolic_angle_max_fxd;               // hyperbolic vectoring |z0| max value

    static Logger<T,FLT> * logger;
};

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
        _ocase( push_constant )
        _ocase( to_flt )
        _ocase( assign )
        _ocase( pop_value )
        _ocase( pop_bool )

        _ocase( frexp )
        _ocase( modf )
        _ocase( ilogb )
        _ocase( logb )

        _ocase( nextafter )
        _ocase( nexttoward )
        _ocase( floor )
        _ocase( ceil )
        _ocase( trunc )
        _ocase( extend )
        _ocase( round )
        _ocase( lround )
        _ocase( llround )
        _ocase( iround )
        _ocase( rint )
        _ocase( lrint )
        _ocase( llrint )
        _ocase( irint )
        _ocase( nearbyint )
        _ocase( floorfrac )
        _ocase( ceilfrac )
        _ocase( truncfrac )
        _ocase( extendfrac )
        _ocase( roundfrac )
        _ocase( rfrac )

        _ocase( abs )
        _ocase( neg )
        _ocase( copysign )
        _ocase( add )
        _ocase( sub )
        _ocase( fma )
        _ocase( mul )
        _ocase( mulc )
        _ocase( sqr )
        _ocase( scalbn )
        _ocase( ldexp )
        _ocase( fda )
        _ocase( div )
        _ocase( rcp )
        _ocase( sqrt )
        _ocase( rsqrt )

        _ocase( iszero )
        _ocase( isgreater )
        _ocase( isgreaterequal )
        _ocase( isless )
        _ocase( islessequal )
        _ocase( islessgreater )
        _ocase( isunordered )
        _ocase( isunequal )
        _ocase( isequal )

        _ocase( exp )
        _ocase( expm1 )
        _ocase( expc )
        _ocase( exp2 )
        _ocase( exp10 )
        _ocase( pow )
        _ocase( log )
        _ocase( logn )
        _ocase( log1p )
        _ocase( logc )
        _ocase( log2 )
        _ocase( log10 )

        _ocase( deg2rad )
        _ocase( rad2deg )
        _ocase( sin )
        _ocase( sinpi )
        _ocase( cos )
        _ocase( cospi )
        _ocase( sincos )
        _ocase( sinpicospi )
        _ocase( tanpi )

        _ocase( asin )
        _ocase( acos )
        _ocase( atan )
        _ocase( atan2 )

        _ocase( polar_to_rect )
        _ocase( rect_to_polar )
        _ocase( hypot )
        _ocase( hypoth )

        _ocase( sinh )
        _ocase( cosh )
        _ocase( sinhcosh )
        _ocase( tanh )
        _ocase( asinh )
        _ocase( acosh )
        _ocase( atanh )
        _ocase( atanh2 )

        _ocase( sram_rd )
        _ocase( sram_wr )
        _ocase( dram_rd )
        _ocase( dram_wr )

        default: return "<unknown OP>";
    }
}

#define _log_1( op, opnd1 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op1( uint16_t(Cordic<T,FLT>::OP::op), &opnd1 )
#define _log_1i( op, opnd1 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op1( uint16_t(Cordic<T,FLT>::OP::op), opnd1 )
#define _log_1b( op, opnd1 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op1( uint16_t(Cordic<T,FLT>::OP::op), opnd1 )
#define _log_1f( op, opnd1 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op1( uint16_t(Cordic<T,FLT>::OP::op), opnd1 )
#define _log_2( op, opnd1, opnd2 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op2( uint16_t(Cordic<T,FLT>::OP::op), &opnd1, &opnd2 )
#define _log_2i( op, opnd1, opnd2 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op2( uint16_t(Cordic<T,FLT>::OP::op), &opnd1, opnd2 )
#define _log_2f( op, opnd1, opnd2 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op2( uint16_t(Cordic<T,FLT>::OP::op), &opnd1, opnd2 )
#define _log_3( op, opnd1, opnd2, opnd3 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op3( uint16_t(Cordic<T,FLT>::OP::op), &opnd1, &opnd2, &opnd3 )
#define _log_4( op, opnd1, opnd2, opnd3, opnd4 ) \
            if ( Cordic<T,FLT>::logger != nullptr ) Cordic<T,FLT>::logger->op4( uint16_t(Cordic<T,FLT>::OP::op), &opnd1, &opnd2, &opnd3, &opnd4 )
#define _logconst( c ) \
            constructed( c ); \
            _log_1f( push_constant, _to_flt(c) ); \
            _log_2i( pop_value, c, c );

//-----------------------------------------------------
// Constructor
//-----------------------------------------------------
template< typename T, typename FLT >
Cordic<T,FLT>::Cordic( uint32_t int_exp_w, uint32_t frac_w, bool is_float, uint32_t guard_w, uint32_t n )
{
    if ( n == uint32_t(-1) ) n = 1 + frac_w;
    if ( guard_w == uint32_t(-1) ) guard_w = std::ceil(std::log2(frac_w));
    if ( logger != nullptr ) logger->cordic_constructed( this, int_exp_w, frac_w, is_float, guard_w, n );

    cassert( (1+int_exp_w+frac_w+guard_w) <= (sizeof( T ) * 8), "1 + int_exp_w + frac_w + guard_w does not fit in T container" );
    cassert( int_exp_w != 0, "int_exp_w must be > 0" );
    cassert( frac_w    != 0, "frac_w must be > 0" );

    _is_float        = is_float;
    _int_w           = is_float ? 0         : int_exp_w;                
    _exp_w           = is_float ? int_exp_w : 0;
    _frac_w          = frac_w;
    _guard_w         = guard_w;
    _frac_guard_w    = frac_w + guard_w;
    _frac_guard_mask = (T(1) << _frac_guard_w) - 1;
    _guard_mask      = (T(1) << _guard_w) - 1;
    _exp_mask        = is_float ? ((1 << int_exp_w)-1) : 0;
    _exp_bias        = is_float ? ((_exp_mask >> 1) - 1) : 0;
    _exp_unbiased_min= is_float ? (0 - _exp_bias) : 0;
    _exp_unbiased_max= is_float ? ((1 << (int_exp_w-1))-1) : 0;
    _w               = 1 + int_exp_w + frac_w + guard_w;
    _n               = n;
    _rounding_mode   = FE_TONEAREST;
    _maxint          = is_float ? T(0) : ((T(1) << int_exp_w) - 1);

    // these must be done first because to_t() depends on some of them
    _zero_fxd        = 0;                                               
    _min_fxd         = T(1) << _guard_w;
    _half_fxd        = T(1) << (_frac_guard_w-1);
    _one_fxd         = T(1) << _frac_guard_w;                        
    _two_fxd         = T(2) << _frac_guard_w;
    _four_fxd        = T(4) << _frac_guard_w;
    _quiet_NaN_fxd   = T(1) << (_frac_guard_w-1);

    _max             = is_float ? make_float( false, _exp_mask-1, _frac_guard_mask ) : to_t( std::pow( 2.0, int_exp_w ) - 1.0 );
    _min             = is_float ? make_float( false, 0,           1                ) : to_t( 1.0 / std::pow( 2.0, frac_w ) );
    _lowest          = is_float ? make_float( true,  _exp_mask-1, _frac_guard_mask ) : (T(-1) << (_w-1));
    _zero            = to_t( 0.0 );
    _neg_zero        = to_t( -0.0 );
    _quarter         = to_t( 0.25 );
    _third           = to_t( 1.0 / 3.0 );
    _neg_third       = to_t( -1.0 / 3.0 );
    _half            = to_t( 0.5 );
    _one             = to_t( 1.0 );
    _neg_one         = to_t( -1.0 );
    _two             = to_t( 2.0 );
    _four            = to_t( 4.0 );
    _sqrt2           = to_t( std::sqrt( 2.0 ) );
    _sqrt2_div_2     = to_t( std::sqrt( 2.0 ) / FLT(2.0) );
    _pi              = to_t( std::acos( FLT(-1.0) ) );
    _pi_fxd          = to_t( std::acos( FLT(-1.0) ), false, true );
    _neg_pi          = to_t( -std::acos( FLT(-1.0) ) );
    _tau             = to_t( 2.0 * std::acos( FLT(-1.0) ) );
    _pi_div_2        = to_t( std::acos( FLT(-1.0) ) / FLT(2.0) );
    _pi_div_4        = to_t( std::acos( FLT(-1.0) ) / FLT(4.0) );
    _one_div_pi      = to_t( FLT(1.0) / std::acos( FLT(-1.0) ) );
    _two_div_pi      = to_t( FLT(2.0) / std::acos( FLT(-1.0) ) );
    _four_div_pi     = to_t( FLT(4.0) / std::acos( FLT(-1.0) ) );
    _e               = to_t( std::exp( FLT(  1 ) ) );
    _log2            = to_t( std::log( FLT(  2 ) ) );
    _log10           = to_t( std::log( FLT( 10 ) ) );

    _logconst( _zero );
    _logconst( _one  );

    _circular_atan_fxd    = new T[n+1];
    _hyperbolic_atanh_fxd = new T[n+1];

    // compute atan/atanh table in high-resolution floating point
    //
    FLT pow2 = 1.0;
    for( uint32_t i = 0; i <= n; i++ )
    {
        FLT a  = std::atan( pow2 );
        FLT ah = std::atanh( pow2 );
        _circular_atan_fxd[i] =    to_t( a, false, true );
        _hyperbolic_atanh_fxd[i] = (i == 0) ? to_t( -1, false, true ) : to_t( ah, false, true );

        if ( debug ) printf( "i=%2d a=%30.27g ah=%30.27g y=%30.27g\n", i, double(a), double(ah), double(pow2) );
        pow2 /= 2.0;
    }

    // calculate max |z0| angle allowed
    T xx, yy, zz;
    _circular_angle_max_fxd   = _one_fxd;   // to avoid triggering assert
    _hyperbolic_angle_max_fxd = _zero_fxd;  // to disable assert
    circular_vectoring(   _one_fxd,  _one_fxd, _zero_fxd, xx, yy, _circular_angle_max_fxd );
    hyperbolic_vectoring( _half_fxd, _one_fxd, _zero_fxd, xx, yy, _hyperbolic_angle_max_fxd );
    if ( debug ) std::cout << "circular_angle_max_fxd="                 << std::setw(30) << _to_flt(_circular_angle_max_fxd, false, true) << "\n";
    if ( debug ) std::cout << "hyperbolic_angle_max_fxd="               << std::setw(30) << _to_flt(_hyperbolic_angle_max_fxd, false, true) << "\n";
    
    // calculate gain by plugging in x=1,y=0,z=0 into CORDICs
    circular_rotation(    _one_fxd, _zero_fxd, _zero_fxd, _circular_rotation_gain_fxd,    yy, zz );
    circular_vectoring(   _one_fxd, _zero_fxd, _zero_fxd, _circular_vectoring_gain_fxd,   yy, zz );
    hyperbolic_rotation(  _one_fxd, _zero_fxd, _zero_fxd, _hyperbolic_rotation_gain_fxd,  yy, zz );
    hyperbolic_vectoring( _one_fxd, _zero_fxd, _zero_fxd, _hyperbolic_vectoring_gain_fxd, yy, zz );

    // calculate 1/gain_fxd which are the multiplication factors
    _circular_rotation_one_over_gain_fxd    = to_t( FLT(1) / _to_flt(_circular_rotation_gain_fxd,    false, true ),     false, true  );
    _circular_rotation_one_over_gain        = to_t( FLT(1) / _to_flt(_circular_rotation_gain_fxd,    false, true ),     false, false );
    _circular_vectoring_one_over_gain_fxd   = to_t( FLT(1) / _to_flt(_circular_vectoring_gain_fxd,   false, true ),     false, true  );
    _circular_vectoring_one_over_gain       = to_t( FLT(1) / _to_flt(_circular_vectoring_gain_fxd,   false, true ),     false, false );
    _hyperbolic_rotation_one_over_gain_fxd  = to_t( FLT(1) / _to_flt(_hyperbolic_rotation_gain_fxd,  false, true ),     false, true  );
    _hyperbolic_rotation_one_over_gain      = to_t( FLT(1) / _to_flt(_hyperbolic_rotation_gain_fxd,  false, true ),     false, false );
    _hyperbolic_vectoring_one_over_gain_fxd = to_t( FLT(1) / _to_flt(_hyperbolic_vectoring_gain_fxd, false, true ),     false, true  );
    _hyperbolic_vectoring_one_over_gain     = to_t( FLT(1) / _to_flt(_hyperbolic_vectoring_gain_fxd, false, true ),     false, false );
    if ( debug ) printf( "circular_rotation_gain_fxd:                   %016" FMT_LLX "   %.30f\n",  _circular_rotation_gain_fxd, _to_flt(_circular_rotation_gain_fxd, false, true) );
    if ( debug ) printf( "circular_vectoring_gain_fxd:                  %016" FMT_LLX "   %.30f\n",  _circular_vectoring_gain_fxd, _to_flt(_circular_vectoring_gain_fxd, false, true) );
    if ( debug ) printf( "hyperbolic_rotation_gain_fxd:                 %016" FMT_LLX "   %.30f\n",  _hyperbolic_rotation_gain_fxd, _to_flt(_hyperbolic_rotation_gain_fxd, false, true) );
    if ( debug ) printf( "hyperbolic_vectoring_gain_fxd:                %016" FMT_LLX "   %.30f\n",  _hyperbolic_vectoring_gain_fxd, _to_flt(_hyperbolic_vectoring_gain_fxd, false, true) );
    if ( debug ) printf( "circular_rotation_one_over_gain_fxd:          %016" FMT_LLX "   %.30f\n",  _circular_rotation_one_over_gain_fxd, _to_flt(_circular_rotation_one_over_gain_fxd, false, true) );
    if ( debug ) printf( "circular_vectoring_one_over_gain_fxd:         %016" FMT_LLX "   %.30f\n",  _circular_vectoring_one_over_gain_fxd, _to_flt(_circular_vectoring_one_over_gain_fxd, false, true) );
    if ( debug ) printf( "hyperbolic_rotation_one_over_gain_fxd:        %016" FMT_LLX "   %.30f\n",  _hyperbolic_rotation_one_over_gain_fxd, _to_flt(_hyperbolic_rotation_one_over_gain_fxd, false, true) );
    if ( debug ) printf( "hyperbolic_vectoring_one_over_gain_fxd:       %016" FMT_LLX "   %.30f\n",  _hyperbolic_vectoring_one_over_gain_fxd, _to_flt(_hyperbolic_vectoring_one_over_gain_fxd, false, true) );
}

template< typename T, typename FLT >
Cordic<T,FLT>::~Cordic( void )
{
    if ( logger != nullptr ) logger->cordic_destructed( this );

    delete _circular_atan_fxd;
    delete _hyperbolic_atanh_fxd;
}

template< typename T, typename FLT >
void Cordic<T,FLT>::log_constructed( void )
{
    if ( logger != nullptr ) logger->cordic_constructed( this, _int_w|_exp_w, _frac_w, _is_float, _guard_w, _n );
}

//-----------------------------------------------------
// Constants
//-----------------------------------------------------
template< typename T, typename FLT >
inline bool Cordic<T,FLT>::is_float( void ) const
{
    return _is_float;
}

template< typename T, typename FLT >
inline uint32_t Cordic<T,FLT>::int_w( void ) const
{
    return _int_w;
}

template< typename T, typename FLT >
inline uint32_t Cordic<T,FLT>::exp_w( void ) const
{
    return _exp_w;
}

template< typename T, typename FLT >
inline uint32_t Cordic<T,FLT>::frac_w( void ) const
{
    return _frac_w;
}

template< typename T, typename FLT >
inline uint32_t Cordic<T,FLT>::guard_w( void ) const
{
    return _guard_w;
}

template< typename T, typename FLT >
inline uint32_t Cordic<T,FLT>::w( void ) const
{
    return _w;
}

template< typename T, typename FLT >
inline uint32_t Cordic<T,FLT>::n( void ) const
{
    return _n;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::maxint( void ) const
{
    return _maxint;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::max( void ) const
{
    _log_1f( push_constant, to_flt(_max) );
    return _max;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::min( void ) const
{
    _log_1f( push_constant, to_flt(_min) );
    return _min;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::denorm_min( void ) const
{
    return min();
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::lowest( void ) const
{
    _log_1f( push_constant, to_flt(_lowest) );
    return _lowest;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::epsilon( void ) const
{
    return min();
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::round_error( void ) const
{
    return min();
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::zero( void ) const
{
    _log_1f( push_constant, to_flt(_zero) );
    return _zero;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::one( void ) const
{
    _log_1f( push_constant, to_flt(_one) ); 
    return _one;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::neg_one( void ) const
{
    _log_1f( push_constant, to_flt(_neg_one) ); 
    return _neg_one;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::two( void ) const
{
    _log_1f( push_constant, to_flt(_two) ); 
    return _two;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::half( void ) const
{
    _log_1f( push_constant, to_flt(_half) ); 
    return _half;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::quarter( void ) const
{
    _log_1f( push_constant, to_flt(_quarter) ); 
    return _quarter;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sqrt2( void ) const
{
    _log_1f( push_constant, to_flt(_sqrt2) ); 
    return _sqrt2;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sqrt2_div_2( void ) const
{
    _log_1f( push_constant, to_flt(_sqrt2_div_2) ); 
    return _sqrt2_div_2;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pi( void ) const
{
    _log_1f( push_constant, to_flt(_pi) ); 
    return _pi;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::tau( void ) const
{
    _log_1f( push_constant, to_flt(_tau) ); 
    return _tau;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pi_div_2( void ) const
{
    _log_1f( push_constant, to_flt(_pi_div_2) ); 
    return _pi_div_2;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pi_div_4( void ) const
{
    _log_1f( push_constant, to_flt(_pi_div_4) ); 
    return _pi_div_4;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::one_div_pi( void ) const
{
    _log_1f( push_constant, to_flt(_one_div_pi) ); 
    return _one_div_pi;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::two_div_pi( void ) const
{
    _log_1f( push_constant, to_flt(_two_div_pi) ); 
    return _two_div_pi;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::four_div_pi( void ) const
{
    _log_1f( push_constant, to_flt(_four_div_pi) ); 
    return _four_div_pi;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::e( void ) const
{
    _log_1f( push_constant, to_flt(_e) ); 
    return _e;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::nan( const char * arg ) const
{
    return to_t( FLT( std::strtod( arg, nullptr ) ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::quiet_NaN( void ) const
{
    return to_t( std::numeric_limits<FLT>::quiet_NaN() );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::signaling_NaN( void ) const
{
    return to_t( std::numeric_limits<FLT>::signaling_NaN() );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::infinity( void ) const
{
    return to_t( std::numeric_limits<FLT>::infinity() );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::ninfinity( void ) const
{
    return to_t( -std::numeric_limits<FLT>::infinity() );
}

//-----------------------------------------------------
// Conversion
//-----------------------------------------------------
template< typename T, typename FLT >
inline T Cordic<T,FLT>::to_t( FLT _x, bool is_final, bool to_fixed ) const
{
    if ( is_final && debug ) std::cout << "to_t begin: x=" << _x << " is_final=" << is_final << " to_fixed=" << to_fixed << "\n";
    FLT x = _x;
    T   x_t;
    if ( _is_float && !to_fixed ) {
        // FLOAT
        //
        // Deconstruct FLT into sign, exp, mantissa,
        // then reconstruct into new encoding.
        //
        bool      x_sign  = std::signbit( x );
        int       x_class = std::fpclassify( x );
        int32_t   x_exp   = 0;
        EXP_CLASS x_exp_class;
        switch( x_class )
        {
            case FP_ZERO:       
            {
                x_exp_class = EXP_CLASS::ZERO; 
                x_t         = 0;
                break;
            }

            case FP_INFINITE:
            {
                x_exp_class = EXP_CLASS::INFINITE;
                x_t         = 0;
                break;
            }

            case FP_SUBNORMAL:
            case FP_NORMAL:
            case FP_NAN:
            {
                // cheat for now and use double-precision to extract exponent and fraction
                //
                union {
                    double      f;
                    uint64_t    u;
                } xx;
                         xx.f     = x;
                uint32_t f_frac_w = 52;
                uint32_t f_exp_w  = 11;
                uint64_t x_u      = xx.u;
                         x_t      = (x_u >> 0)        & ((1LL << f_frac_w)-1);
                         x_exp    = (x_u >> f_frac_w) & ((1LL << f_exp_w)-1);
                         x_exp   -= (1 << (f_exp_w-1))-1;          // exp bias
                int32_t lshift = int32_t(_frac_guard_w) - int32_t(f_frac_w);
                if ( false && debug ) std::cout << "to_t: x_f=" << xx.f << std::hex << " x_u=" << x_u << 
                                          " x_t=" << x_t << std::dec << " x_unbiased_exp=" << x_exp << " lshift=" << lshift << "\n";
                x_t = (lshift >= 0) ? (x_t << lshift) : (x_t >> -lshift);
                cassert( (x_t & _frac_guard_mask) == x_t, "to_t: unexpected x_t after lshift/rshift" );

                if ( x_class == FP_NORMAL ) {
                    x_exp_class = EXP_CLASS::NORMAL;
                    x_t |= _one_fxd;
                } else if ( x_class == FP_SUBNORMAL ) {
                    x_exp_class = EXP_CLASS::SUBNORMAL;
                    x_exp = 0;
                } else {
                    x_exp_class = EXP_CLASS::NOT_A_NUMBER;
                    x_exp = 0;
                }
                if ( false && debug ) std::cout << "to_t: x_t_new=" << std::hex << x_t << std::dec << " x_exp_class=" << to_str(x_exp_class) <<
                                          " x_exp=" << x_exp << "\n";
                break;
            }

            default:
            {
                x_exp_class = EXP_CLASS::NOT_A_NUMBER;
                x_t = 666;
                break;
            }
        }

        reconstruct( x_t, x_exp_class, x_exp, x_sign );

    } else {
        // FIXED: convert to integer 
        //
        bool x_sign = x < FLT(0.0);
        if ( x_sign ) x = -x;
        cassert( _is_float ? (T(x) < (1 << 2)) : (T(x) < (T(1) << _int_w)),
                 "to_t: integer part of |x| " + std::to_string(x) + " does not fit in fixed-point int_w bits" ); 
        
        FLT x_f = x * FLT( _one_fxd );  // treat it as an integer
        switch( _rounding_mode )
        {
            case FE_NOROUND:                                                        break;
            case FE_DOWNWARD:                       x_f = std::floor( x_f );        break;
            case FE_UPWARD:                         x_f = std::ceil( x_f );         break;
            case FE_TOWARDZERO:                     x_f = std::trunc( x_f );        break;
            case FE_AWAYFROMZERO:                   x_f = std::signbit( x_f ) ? std::floor( x_f ) : std::ceil( x_f ); break;
            case FE_TONEAREST:                      x_f = std::round( x_f );        break;
            default:                                                                break;
        }
        x_t = x_f;
        if ( x_sign ) x_t = -x_t;
    } 
    if ( is_final ) _log_1f( push_constant, _x );       // note: we purposely don't round constants
    return x_t;
}

template< typename T, typename FLT >
inline FLT Cordic<T,FLT>::to_flt( const T& x ) const
{
    return _to_flt( x, true, false, true );
}

template< typename T, typename FLT >
inline FLT Cordic<T,FLT>::_to_flt( const T& _x, bool is_final, bool from_fixed, bool allow_debug ) const
{
    T    x = _x;
    bool x_sign;
    FLT  x_f;
    if ( _is_float && !from_fixed ) {
        // FLOAT
        //
        // Deconstruct into sign, exp_class, exp, fraction.
        // Then do the right thing based on the exp_class.
        // 
        EXP_CLASS x_exp_class;
        int32_t   x_exp;
        deconstruct( x, x_exp_class, x_exp, x_sign, allow_debug );
        switch( x_exp_class )
        {
            case EXP_CLASS::ZERO:                   
                x_f = 0.0;
                break;

            case EXP_CLASS::NORMAL: 
                x_f = std::scalbn( FLT( x ), x_exp - _frac_guard_w );
                break;

            case EXP_CLASS::SUBNORMAL:
                x_f = std::scalbn( FLT( x ), -_exp_bias - _frac_guard_w );
                break;

            case EXP_CLASS::INFINITE:
                x_f = std::numeric_limits<FLT>::infinity();
                break;

            case EXP_CLASS::NOT_A_NUMBER:
            default:
                x_f = std::numeric_limits<FLT>::quiet_NaN();   // FIXIT: retain frac bits 
                break;
        }
    } else {
        // FIXED: return x / _one
        //
        x_sign = x < 0;
        if ( x_sign ) x = -x;
        x_f = FLT( x ) / FLT( _one_fxd );
    } 
    if ( x_sign ) x_f = -x_f;
    if ( is_final ) _log_2f( to_flt, _x, x_f );
    return x_f;
}

template< typename T, typename FLT >
inline std::string Cordic<T,FLT>::to_string( const T& x, bool from_fixed ) const
{
    // floating-point representation
    return std::to_string( _to_flt( x, false, from_fixed ) );  
}

template< typename T, typename FLT >
inline std::string Cordic<T,FLT>::to_rstring( const T& _x, bool from_fixed ) const
{
    // raw integer representation
    T x = _x;
    bool sign = signbit( x );
    if ( sign && (!_is_float || from_fixed) ) x = -x;
    int32_t twidth = 8 * sizeof( T );
    int32_t bwidth = _w - 1;
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
inline std::string Cordic<T,FLT>::to_bstring( const T& _x, bool from_fixed ) const
{
    (void)from_fixed;
    // binary representation
    T x = _x;
    std::string bs = "";
    for( uint32_t i = 0; i < _w; i++ )
    {
        if ( i == (_w - 1) || i == _frac_guard_w || i == _guard_w ) bs = " " + bs;
        const char * b = (x & 1) ? "1" : "0";
        bs = b + bs;
        x >>= 1;
    }
    return bs;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::make_fixed( bool sign, const T& i, const T& f ) const
{
    cassert( !_is_float, "make_fixed may be called only for is_float=false Cordics" );
    cassert( i >= 0 && i <= _maxint, "make_fixed integer part must be in range 0 .. _maxint" );
    cassert( f >= 0 && f <= _frac_guard_mask, "make_fixed fractional part must be in range 0 .. (1 << (frac_w+guard_w))-1" );

    return (T(sign) << (_w - 1))      |
           (T(i)    << _frac_guard_w) |
           (T(f)    << 0);
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::make_float( bool sign, const T& e, const T& f ) const
{
    cassert( _is_float, "make_float may be called only for is_float=true Cordics" );
    cassert( e >= 0 && e <= _exp_mask, "make_float biased exponent part must be in range 0 .. _exp_mask, got " + std::to_string(e) );
    cassert( f >= 0 && f <= _frac_guard_mask, "make_float mantissa part must be in range 0 .. (1 << (frac_w+guard_w))-1" );

    return (T(sign) << (_w - 1))      |
           (T(e)    << _frac_guard_w) |
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
    const T ONE = _one_fxd;
    const T ANGLE_MAX = _circular_angle_max_fxd + 2*_min_fxd;
    if ( debug ) printf( "circular_rotation begin: xyz_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX ",0x%016" FMT_LLX "] xyz=[%.30f,%.30f,%.30f]\n",
                         x0, y0, z0, _to_flt(x0, false, true), _to_flt(y0, false, true), _to_flt(z0, false, true) );
    cassert( x0 >= -ONE       && x0 <= ONE,       "circular_rotation x0 must be in the range -1 .. 1" );
    cassert( y0 >= -ONE       && y0 <= ONE,       "circular_rotation y0 must be in the range -1 .. 1" );
    cassert( z0 >= -ANGLE_MAX && z0 <= ANGLE_MAX, "circular_rotation |z0| must be <= circular_angle_max (" +
                                                  to_string(ANGLE_MAX, true) + "), got z0=" + to_string(z0, true) );

    //-----------------------------------------------------
    // d = (z >= 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctan(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = _n;
    for( uint32_t i = 0; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "circular_rotation: i=%2d xyz_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX ",0x%016" FMT_LLX "] xyz=[%.30f,%.30f,%.30f] test=%d\n", 
                             i, x, y, z, _to_flt(x, false, true), _to_flt(y, false, true), _to_flt(z, false, true), int(z >= _zero_fxd) );
        if ( z >= 0 ) {
            xi = x - (y >> i);
            yi = y + (x >> i);
            zi = z - _circular_atan_fxd[i];
        } else {
            xi = x + (y >> i);
            yi = y - (x >> i);
            zi = z + _circular_atan_fxd[i];
        }
        x = xi;
        y = yi;
        z = zi;
    }

    //-----------------------------------------------------
    // circular rotation mode results after step n:
    //      x = gain*(x0*cos(z0) - y0*sin(z0))          gain=1.64676...
    //      y = gain*(y0*cos(z0) + x0*sin(z0))
    //      z = 0
    //-----------------------------------------------------
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
    const T ONE = _one_fxd;
    const T THREE = 3*ONE;
    const T PI  = _pi_fxd;
    const T ANGLE_MAX = _circular_angle_max_fxd + 2*_min_fxd;
    if ( debug ) printf( "circular_vectoring begin: xyz_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX ",0x%016" FMT_LLX "] xyz=[%.30f,%.30f,%.30f]\n",
                         x0, y0, z0, _to_flt(x0, false, true), _to_flt(y0, false, true), _to_flt(z0, false, true) );
    cassert( x0 >= -THREE && x0 <= THREE, "circular_vectoring x0 must be in the range -3 .. 3" );
    cassert( y0 >= -ONE   && y0 <= ONE  , "circular_vectoring y0 must be in the range -1 .. 1" );
    cassert( z0 >= -PI    && z0 <= PI   , "circular_vectoring z0 must be in the range -PI .. PI" );
    cassert( std::abs( std::atan( _to_flt(y0, false, true) / _to_flt(x0, false, true) ) ) <= _to_flt(ANGLE_MAX, false, true),
                                        "circular_vectoring |atan(y0/x0)| must be <= circular_angle_max" );

    //-----------------------------------------------------
    // d = (y < 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctan(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = _n;
    for( uint32_t i = 0; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "circular_vectoring: i=%2d xyz_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX ",0x%016" FMT_LLX "] xyz=[%.30f,%.30f,%.30f] test=%d\n", 
                             i, x, y, z, _to_flt(x, false, true), _to_flt(y, false, true), _to_flt(z, false, true), int((x < 0) != (y < 0)) );
        if ( y < 0 ) {
            xi = x - (y >> i);
            yi = y + (x >> i);
            zi = z - _circular_atan_fxd[i];
        } else {
            xi = x + (y >> i);
            yi = y - (x >> i);
            zi = z + _circular_atan_fxd[i];
        }
        x = xi;
        y = yi;
        z = zi;
    }

    //-----------------------------------------------------
    // circular vectoring mode results after step n:
    //      x = gain*sqrt(x0^2 + y0^2)                  gain=1.64676...
    //      y = 0
    //      z = z0 + atan( y0/x0 )
    //-----------------------------------------------------
}

template< typename T, typename FLT >
void Cordic<T,FLT>::circular_vectoring_xy( const T& x0, const T& y0, T& x, T& y ) const
{
    //-----------------------------------------------------
    // input ranges allowed:
    //      -3  <= x0 <= 3
    //      -1  <= y0 <= 1
    //-----------------------------------------------------
    const T ONE = _one_fxd;
    const T THREE = 3*ONE;
    if ( debug ) printf( "circular_vectoring_xy begin: xy_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX "] xy=[%.30f,%.30f]\n",
                         x0, y0, _to_flt(x0, false, true), _to_flt(y0, false, true) );
    cassert( x0 >= -THREE && x0 <= THREE, "circular_vectoring_xy x0 must be in the range -3 .. 3" );
    cassert( y0 >= -ONE   && y0 <= ONE  , "circular_vectoring_xy y0 must be in the range -1 .. 1" );

    //-----------------------------------------------------
    // d = (y < 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    //-----------------------------------------------------
    x = x0;
    y = y0;
    uint32_t n = _n;
    for( uint32_t i = 0; i <= n; i++ )
    {
        T xi;
        T yi;
        if ( debug ) printf( "circular_vectoring_xy: i=%d xy_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX "] xy=[%.30f,%.30f] test=%d\n", 
                             i, x, y, _to_flt(x, false, true), _to_flt(y, false, true), int(y < 0) );
        if ( y < 0 ) {
            xi = x - (y >> i);
            yi = y + (x >> i);
        } else {
            xi = x + (y >> i);
            yi = y - (x >> i);
        }
        x = xi;
        y = yi;
    }

    //-----------------------------------------------------
    // circular vectoring mode results after step n:
    //      x = gain*sqrt(x0^2 + y0^2)                  gain=1.64676...
    //      y = 0
    //-----------------------------------------------------
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
    const T TWO = _two_fxd;
    const T ANGLE_MAX = _hyperbolic_angle_max_fxd + 2*_min_fxd;
    if ( debug ) printf( "hyperbolic_rotation begin: xyz_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX ",0x%016" FMT_LLX "] xyz=[%.30f,%.30f,%.30f]\n",
                         x0, y0, z0, _to_flt(x0, false, true), _to_flt(y0, false, true), _to_flt(z0, false, true) );
    cassert( x0 >= -TWO       && x0 <= TWO,       "hyperbolic_rotation x0 must be in the range -2 .. 2" );
    cassert( y0 >= -TWO       && y0 <= TWO,       "hyperbolic_rotation y0 must be in the range -2 .. 2" );
    cassert( z0 >= -ANGLE_MAX && z0 <= ANGLE_MAX, "hyperbolic_rotation |z0| must be <= hyperbolic_angle_max (" + 
                                                  to_string(ANGLE_MAX, true) + "), got z0=" + to_string(z0, true) );

    //-----------------------------------------------------
    // d = (z >= 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctanh(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = _n;
    uint32_t next_dup_i = 4;     
    for( uint32_t i = 1; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "hyperbolic_rotation: i=%2d xyz_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX ",0x%016" FMT_LLX "] xyz=[%.30f,%.30f,%.30f] test=%d\n", 
                             i, x, y, z, _to_flt(x, false, true), _to_flt(y, false, true), _to_flt(z, false, true), int(z >= 0) );
        if ( z >= 0 ) {
            xi = x + (y >> i);
            yi = y + (x >> i);
            zi = z - _hyperbolic_atanh_fxd[i];
        } else {
            xi = x - (y >> i);
            yi = y - (x >> i);
            zi = z + _hyperbolic_atanh_fxd[i];
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

    //-----------------------------------------------------
    // hyperbolic rotation mode results after step n:
    //      x = gain*(x0*cosh(z0) + y0*sinh(z0))        gain=0.828159...
    //      y = gain*(y0*cosh(z0) + x0*sinh(z0))
    //      z = 0
    //-----------------------------------------------------
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
    const T TWO = _two_fxd;
    const T PI  = _pi_fxd;
    const T ANGLE_MAX = _hyperbolic_angle_max_fxd;
    if ( debug ) printf( "hyperbolic_vectoring begin: xyz_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX ",0x%016" FMT_LLX "] xyz=[%.30f,%.30f,%.30f]\n",
                         x0, y0, z0, _to_flt(x0, false, true), _to_flt(y0, false, true), _to_flt(z0, false, true) );
    cassert( x0 >= -TWO && x0 <= TWO, "hyperbolic_vectoring x0 must be in the range -2 .. 2" );
    cassert( y0 >= -TWO && y0 <= TWO, "hyperbolic_vectoring y0 must be in the range -2 .. 2" );
    cassert( z0 >= -PI  && z0 <= PI , "hyperbolic_vectoring z0 must be in the range -PI .. PI" );
    cassert( (ANGLE_MAX == 0 || std::abs( std::atanh( _to_flt(y0, false, true) / _to_flt(x0, false, true) ) ) <= _to_flt(ANGLE_MAX, false, true)),
                                        "hyperbolic_vectoring |atanh(y0/x0)| must be <= hyperbolic_angle_max=" + 
                                        std::to_string(_to_flt(ANGLE_MAX, false, true)) );

    //-----------------------------------------------------
    // d = (y < 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctanh(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = _n;
    uint32_t next_dup_i = 4;     
    for( uint32_t i = 1; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "hyperbolic_vectoring: i=%2d xyz_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX ",0x%016" FMT_LLX "] xyz=[%.30f,%.30f,%.30f] test=%d\n", 
                             i, x, y, z, _to_flt(x, false, true), _to_flt(y, false, true), _to_flt(z, false, true), int((x < 0) != (y < 0)) );
        if ( y < 0 ) {
            xi = x + (y >> i);
            yi = y + (x >> i);
            zi = z - _hyperbolic_atanh_fxd[i];
        } else {
            xi = x - (y >> i);
            yi = y - (x >> i);
            zi = z + _hyperbolic_atanh_fxd[i];
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

    //-----------------------------------------------------
    // hyperbolic vectoring mode results after step n:
    //      x = gain*sqrt(x0^2 - y0^2)                  gain=0.828159...
    //      y = 0
    //      z = z0 + atanh( y0/x0 )
    //-----------------------------------------------------
}

template< typename T, typename FLT >
void Cordic<T,FLT>::hyperbolic_vectoring_xy( const T& x0, const T& y0, T& x, T& y ) const
{
    //-----------------------------------------------------
    // input ranges allowed:
    //      -2  <= x0 <= 2
    //      -2  <= y0 <= 2
    //-----------------------------------------------------
    const T TWO = _two_fxd;
    const T PI  = _pi_fxd;
    if ( debug ) printf( "hyperbolic_vectoring_xy begin: xy_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX "] xy=[%.30f,%.30f]\n",
                         x0, y0, _to_flt(x0, false, true), _to_flt(y0, false, true) );
    cassert( x0 >= -TWO && x0 <= TWO, "hyperbolic_vectoring_xy x0 must be in the range -2 .. 2" );
    cassert( y0 >= -TWO && y0 <= TWO, "hyperbolic_vectoring_xy y0 must be in the range -2 .. 2" );

    //-----------------------------------------------------
    // d = (y < 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    //-----------------------------------------------------
    x = x0;
    y = y0;
    uint32_t n = _n;
    uint32_t next_dup_i = 4;     
    for( uint32_t i = 1; i <= n; i++ )
    {
        T xi;
        T yi;
        if ( debug ) printf( "hyperbolic_vectoring_xy: i=%2d xy_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX "] xy=[%.30f,%.30f] test=%d\n", 
                             i, x, y, _to_flt(x, false, true), _to_flt(y, false, true), int(y < 0) );
        if ( y < 0 ) {
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

    //-----------------------------------------------------
    // hyperbolic vectoring mode results after step n:
    //      x = gain*sqrt(x0^2 - y0^2)                  gain=0.828159...
    //      y = 0
    //-----------------------------------------------------
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
    const T ONE = _one_fxd;
    const T TWO = _two_fxd;
    if ( debug ) printf( "linear_rotation begin: xyz_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX ",0x%016" FMT_LLX "] xyz=[%.30f,%.30f,%.30f]\n",
                         x0, y0, z0, _to_flt(x0, false, true), _to_flt(y0, false, true), _to_flt(z0, false, true) );
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
    uint32_t n = _n;
    T pow2 = ONE;
    for( uint32_t i = 0; i <= n; i++, pow2 >>= 1 )
    {
        if ( debug ) printf( "linear_rotation: i=%2d xyz_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX ",0x%016" FMT_LLX "] xyz=[%.30f,%.30f,%.30f] test=%d\n", 
                             i, x, y, z, _to_flt(x, false, true), _to_flt(y, false, true), _to_flt(z, false, true), int(z >= 0) );
        T yi;
        T zi;
        if ( z >= 0 ) {
            yi = y + (x >> i);
            zi = z - pow2;
        } else {
            yi = y - (x >> i);
            zi = z + pow2;
        }
        y = yi;
        z = zi;
    }

    //-----------------------------------------------------
    // linear rotation mode results after step n:
    //      x = x0
    //      y = y0 + x0*z0
    //      z = 0
    //-----------------------------------------------------
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
    const T ONE = _one_fxd;
    const T TWO = _two_fxd;
    if ( debug ) printf( "linear_vectoring begin: xyz_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX ",0x%016" FMT_LLX "] xyz=[%.30f,%.30f,%.30f]\n",
                         x0, y0, z0, _to_flt(x0, false, true), _to_flt(y0, false, true), _to_flt(z0, false, true) );
    cassert( x0 >= -TWO && x0 <= TWO, "linear_vectoring x0 must be in the range -2 .. 2, got " + to_string(x0, true) );
    cassert( y0 >= -TWO && y0 <= TWO, "linear_vectoring y0 must be in the range -2 .. 2, got " + to_string(y0, true) );
    //cassert( std::abs( _to_flt(y0, false, true) / _to_flt(x0, false, true) ) <= FLT(1.0) &&
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
    uint32_t n = _n;
    T pow2 = ONE;
    for( uint32_t i = 0; i <= n; i++, pow2 >>= 1 )
    {
        if ( debug ) printf( "linear_vectoring: i=%2d xyz_fxd=[0x%016" FMT_LLX ",0x%016" FMT_LLX ",0x%016" FMT_LLX "] xyz=[%.30f,%.30f,%.30f] test=%d\n", 
                             i, x, y, z, _to_flt(x, false, true), _to_flt(y, false, true), _to_flt(z, false, true), int(y < 0) );
        T yi;
        T zi;
        if ( y < 0 ) {
            yi = y + (x >> i);
            zi = z - pow2;
        } else {
            yi = y - (x >> i);
            zi = z + pow2;
        }
        y = yi;
        z = zi;
    }

    //-----------------------------------------------------
    // linear vectoring mode results after step n:
    //      x = x0
    //      y = 0
    //      z = z0 + y0/x0
    //-----------------------------------------------------
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
    _log_2( assign, x, y );
    x = y;
    return x;
}

template< typename T, typename FLT >
inline T&   Cordic<T,FLT>::pop_value( T& x, const T& y ) const
{
    _log_2i( pop_value, x, y );
    x = y;
    return x;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::pop_bool( bool b ) const
{
    _log_1i( pop_bool, b );
    return b;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::signbit( const T& x ) const                                     
{
    return (x >> (_w-1)) & 1;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::frexp( const T& _x, int * e ) const
{
    _log_1( frexp, _x );
    switch( fpclassify( _x ) )
    {
        case FP_SUBNORMAL:
        case FP_NORMAL:
        {
            T x = _x;
            EXP_CLASS x_exp_class;
            int32_t   x_exp;
            bool      x_sign;
            deconstruct( x, x_exp_class, x_exp, x_sign );
            *e = x_exp;
            if ( x_sign ) x = -x;
            _log_2i( frexp, _x, T(x_exp) );
            return x;
        }

        default:
        {
            _log_2i( frexp, _x, T(0) );
            *e = 0;
            return _x;
        }
    }
}

template< typename T, typename FLT >
T Cordic<T,FLT>::modf( const T& _x, T * i ) const
{
    if ( debug ) std::cout << "modf begin: x=" << _to_flt(_x) << "\n";
    _log_1( modf, _x );
    T  x = _x;
      *i = _x;
    switch( fpclassify( _x ) )
    {
        case FP_ZERO:
        {
            // return +/- zero 
            break;
        }

        case FP_NAN:
        {
            // return +/- quiet_NaN for both
            x |= _quiet_NaN_fxd;
            break;
        }

        case FP_INFINITE:
        {
            // return +/- 0 and leave *i as +/- inf
            x &= ~((T(1) << (_w-1)) - 1);          // keep sign
            break;
        }
            
        case FP_SUBNORMAL:
        {
            // return subnormal and set *i to +/- zero
            *i &= ~((T(1) << (_w-1)) - 1);         // keep sign
            break;
        }

        case FP_NORMAL:
        {
            if ( _is_float ) {
                // FLOAT: take the easy way out for now using FLT  (FIXIT: do right way later)
                //
                FLT x_f  = _to_flt( x );
                T   x_i  = x_f;
                    *i   = to_t( FLT(x_i), false );
                    x_f -= x_i;
                    x    = to_t( x_f, false );
            } else {
                // FIXED: simple, but be careful to extend the sign
                //
                T x_sign = signbit( x );
                *i = x & ~_frac_guard_mask;
                x &= _frac_guard_mask;
                if ( x_sign ) x |= T(-1) << (_frac_w + _guard_w);
            }
            break;
        }

        default:
        {
            x |= _quiet_NaN_fxd;
            *i = x;
            break;
        }
    }
    return x;
}

template< typename T, typename FLT >
int Cordic<T,FLT>::ilogb( const T& x ) const
{
    int exp;
    (void)frexp( x, &exp );
    return exp;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::logb( const T& x ) const
{
    int exp = ilogb( x );
    return scalbn( _one, exp );
}

template< typename T, typename FLT >
inline int Cordic<T,FLT>::fpclassify( const T& _x ) const                                     
{
    T x = _x;
    EXP_CLASS x_exp_class;
    int32_t   x_exp;
    bool      x_sign;
    deconstruct( x, x_exp_class, x_exp, x_sign );
    switch( x_exp_class )
    {
        case EXP_CLASS::ZERO:                   return FP_ZERO;
        case EXP_CLASS::NORMAL:                 return FP_NORMAL;
        case EXP_CLASS::SUBNORMAL:              return FP_SUBNORMAL;
        case EXP_CLASS::INFINITE:               return FP_INFINITE;
        case EXP_CLASS::NOT_A_NUMBER:           return FP_NAN;
        default:                                return FP_NAN;
    }
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::iszero( const T& x ) const                                     
{
    return fpclassify( x ) == FP_ZERO;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isfinite( const T& x ) const                                     
{
    int c = fpclassify( x );
    return c != FP_INFINITE && c != FP_NAN;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isinf( const T& x ) const                                     
{
    return fpclassify( x ) == FP_INFINITE;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isnan( const T& x ) const                                     
{
    return fpclassify( x ) == FP_NAN;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isnormal( const T& x ) const                                     
{
    return fpclassify( x ) == FP_NORMAL;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::issubnormal( const T& x ) const                                     
{
    return fpclassify( x ) == FP_SUBNORMAL;
}

template< typename T, typename FLT >
inline int Cordic<T,FLT>::fesetround( int round )
{
    switch( round )
    {
        case FE_NOROUND:        
        case FE_DOWNWARD:        
        case FE_UPWARD:        
        case FE_TOWARDZERO:
        case FE_AWAYFROMZERO:
        case FE_TONEAREST:
            _rounding_mode = round;
            return 0;

        default:
            return -1;
    }
}

template< typename T, typename FLT >
inline int Cordic<T,FLT>::fegetround( void ) const
{
    return _rounding_mode;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::nextafter( const T& from, const T& to ) const
{
    _log_2( nextafter, from, to );
    if ( isequal( from, to ) ) {
        return to;
    } else if ( isless( from, to ) ) {
        return add( from, min() );
    } else {
        return sub( from, min() );
    }
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::nexttoward( const T& from, long double to ) const
{
    return nextafter( from, to_t(to) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::floor( const T& x ) const
{
    _log_1( floor, x );
    if ( (x & _frac_guard_mask) == 0 ) {
        return x;
    } else {
        // get integer part, still encoded;
        // subtract one() if negative
        T i;
        modf( x, &i );
        if ( signbit( i ) ) i = add( i, neg_one() );
        return i;
    }
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::ceil( const T& x ) const
{
    _log_1( ceil, x );
    if ( (x & _frac_guard_mask) == 0 ) {
        return x;
    } else {
        // get integer part, still encoded;
        // add one() if positive
        T i;
        modf( x, &i );
        if ( !signbit( i ) ) i = add( i, one() );
        return i;
    }
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::trunc( const T& x ) const
{
    _log_1( trunc, x );
    if ( (x & _frac_guard_mask) == 0 ) {
        return x;
    } else {
        // get integer part, still encoded;
        // add one() if negative
        T i;
        modf( x, &i );
        if ( signbit( i ) ) i = add( i, one() );
        return i;
    }
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::extend( const T& x ) const
{
    _log_1( extend, x );
    if ( (x & _frac_guard_mask) == 0 ) {
        return x;
    } else {
        // get integer part, still encoded;
        // sub one() if negative
        // add one() if positive
        T i;
        modf( x, &i );
        i = add( i, signbit( i ) ? -neg_one() : one() );
        return i;
    }
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::round( const T& x ) const
{
    _log_1( round, x );
    T i;
    T f = modf( x, &i );
    if ( isgreaterequal( f, half() ) ) f = add( f, signbit(f) ? neg_one() : one() );
    return f;
}

template< typename T, typename FLT >
inline long Cordic<T,FLT>::lround( const T& x ) const
{
    // use round(), then convert to FLT and then integer
    return FLT( round( x ) );
}

template< typename T, typename FLT >
inline long long Cordic<T,FLT>::llround( const T& x ) const
{
    // use round(), then convert to FLT and then integer
    return FLT( round( x ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::iround( const T& x ) const
{
    // same as round(), then convert to FLT and then integer
    return FLT( round( x ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::rint( const T& x, int rmode ) const
{
    if ( rmode < 0 ) rmode = _rounding_mode;

    switch( rmode )
    {
        case FE_NOROUND:                        return x;
        case FE_DOWNWARD:                       return floor( x );
        case FE_UPWARD:                         return ceil( x );
        case FE_TOWARDZERO:                     return trunc( x );
        case FE_AWAYFROMZERO:                   return extend( x );
        case FE_TONEAREST:                      return round( x );
        default:                                return x;                         // shouldn't happen
    }
}

template< typename T, typename FLT >
inline long Cordic<T,FLT>::lrint( const T& x ) const
{
    // use rint() then convert to FLT and then integer
    return FLT( rint( x ) );
}

template< typename T, typename FLT >
inline long long Cordic<T,FLT>::llrint( const T& x ) const
{
    // use rint() then convert to FLT and then integer
    return FLT( rint( x ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::irint( const T& x ) const
{
    // use rint() then convert to FLT and then integer
    return FLT( rint( x ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::nearbyint( const T& x ) const
{
    return rint( x );                   // needs to make sure FE_INEXACT doesn't get raised
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::rfrac( const T& _x, int rmode ) const
{
    if ( rmode < 0 ) rmode = _rounding_mode;

    T x = _x;
    T guard = x & _guard_mask;
    if ( rmode == FE_NOROUND || guard == 0 ) return x; // no need to round

    EXP_CLASS x_exp_class;
    int32_t   x_exp;
    bool      x_sign;
    deconstruct( x, x_exp_class, x_exp, x_sign );  // note that x will be non-negative after this

    x &= ~_guard_mask;

    switch( rmode )
    {
        case FE_DOWNWARD:       if ( x_sign )  x += _min_fxd;                           break;
        case FE_UPWARD:         if ( !x_sign ) x += _min_fxd;                           break;
        case FE_TOWARDZERO:                                                             break;
        case FE_AWAYFROMZERO:   x += x_sign ? -_min_fxd : _min_fxd;                     break;
        case FE_TONEAREST:      if ( guard >= (1 << (_guard_w-1)) ) x += _min_fxd;      break;
        default:                                                                        break;
    }

    if ( x_exp_class == EXP_CLASS::SUBNORMAL && x == 0 ) x_exp_class = EXP_CLASS::ZERO;
    reconstruct( x, x_exp_class, x_exp, x_sign );  
    return x;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::floorfrac( const T& x ) const    { return rfrac( x, FE_DOWNWARD ); }

template< typename T, typename FLT >
inline T Cordic<T,FLT>::ceilfrac( const T& x ) const     { return rfrac( x, FE_UPWARD ); }

template< typename T, typename FLT >
inline T Cordic<T,FLT>::truncfrac( const T& x ) const    { return rfrac( x, FE_TOWARDZERO ); }

template< typename T, typename FLT >
inline T Cordic<T,FLT>::extendfrac( const T& x ) const   { return rfrac( x, FE_AWAYFROMZERO ); }

template< typename T, typename FLT >
inline T Cordic<T,FLT>::roundfrac( const T& x ) const    { return rfrac( x, FE_TONEAREST ); }

template< typename T, typename FLT >
inline T Cordic<T,FLT>::neg( const T& x, bool is_final ) const
{
    if ( is_final ) _log_1( neg, x );
    T x_neg;
    if ( _is_float ) {
        // FLOAT: toggle sign bit
        x_neg = x ^ (1LL << (_w-1));
    } else {
        // FIXED: 2's complement
        bool x_sign = x < 0;
             x_neg  = -x;
        T    sign_mask = x_neg >> (_w - 1);
        cassert( (x == 0 || sign_mask == (x_sign ? T(0) : T(-1))), "neg caused overflow" ); 
    }
    return x_neg;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::neg( const T& x ) const
{
    return neg( x, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::abs( const T& x ) const
{
    _log_1( abs, x );
    T x_abs = x;
    if ( signbit( x_abs ) ) x_abs = neg( x_abs, false );
    return x_abs;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::copysign( const T& x, const T& y ) const
{
    _log_2( copysign, x, y );
    bool x_sign = signbit( x );
    bool y_sign = signbit( y );
    return (x_sign != y_sign) ? neg( x, false ) : x;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::add( const T& x, const T& y ) const
{

    return add( x, y, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sub( const T& x, const T& y, bool is_final ) const
{
    return add( x, neg( y, false ), is_final );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sub( const T& x, const T& y ) const
{
    return sub( x, y, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::add( const T& _x, const T& _y, bool is_final ) const
{
    if ( is_final ) _log_2( add, _x, _y );
    T x = _x;  // will also contain the result
    T y = _y;
    if ( _is_float ) {
        // FLOAT: normalize to largest of two exponents
        //
        EXP_CLASS x_exp_class;
        EXP_CLASS y_exp_class;
        int32_t   x_exp;
        int32_t   y_exp;
        bool      x_sign;
        bool      y_sign;
        reduce_add_args( x, y, x_exp_class, x_exp, x_sign, y_exp_class, y_exp, y_sign );

        // check for special cases
        //
        if ( x_exp_class == EXP_CLASS::NOT_A_NUMBER || y_exp_class == EXP_CLASS::ZERO ) {
            // x is the answer

        } else if ( x_exp_class == EXP_CLASS::ZERO || y_exp_class == EXP_CLASS::NOT_A_NUMBER ) {
            // y is the answer
            x_exp_class = y_exp_class;
            x_exp       = y_exp;
            x           = y;
            x_sign      = y_sign;

        } else if ( x_exp_class == EXP_CLASS::INFINITE ) {
            // x is the answer, but could change to quiet_NaN 
            if ( x_sign != y_sign && y_exp_class == EXP_CLASS::INFINITE ) {
                x_exp_class = EXP_CLASS::NOT_A_NUMBER;
                x |= _quiet_NaN_fxd;
            }

        } else if ( y_exp_class == EXP_CLASS::INFINITE ) {
            // y is the answer, but could change to quiet_NaN 
            if ( x_sign != y_sign && x_exp_class == EXP_CLASS::INFINITE ) {
                y_exp_class = EXP_CLASS::NOT_A_NUMBER;
                y |= _quiet_NaN_fxd;
            }
            x_exp_class = y_exp_class;
            x_exp       = y_exp;
            x           = y;
            x_sign      = y_sign;

        } else {
            //
            // add or sub aligned mantissas
            //
            if ( x_sign != y_sign ) y = -y;
            x += y;
            if ( x < 0 ) {
                x = -x;
                x_sign = !x_sign;
            }

            if ( x == 0 ) {
                x_exp_class = EXP_CLASS::ZERO;
                x_exp = 0;
            } else if ( x >= _two_fxd ) {
                cassert( x < _four_fxd, "add() sum of mantissas should have been less than 4" );
                cassert( x_exp_class == EXP_CLASS::NORMAL || x_exp_class == EXP_CLASS::SUBNORMAL, "add() expected NORMAL or SUBNORMAL class x value at this point" );
                x_exp++;
                x = (x >> 1) | (x & 1); // perhaps set sticky bit
                if ( x_exp >= _exp_unbiased_max ) {
                    x_exp_class = EXP_CLASS::INFINITE;
                    x_exp = _exp_mask;
                    x = 0;
                }
            } else if ( x >= _one_fxd && x_exp_class == EXP_CLASS::SUBNORMAL ) {
                x_exp_class = EXP_CLASS::NORMAL;
                x_exp = 0;
            } else if ( x < _one_fxd ) {
                while( x < _one_fxd && x_exp > _exp_unbiased_min )
                {
                    x_exp--;
                    x <<= 1;
                }
                x_exp_class = (x_exp == _exp_unbiased_min) ? EXP_CLASS::SUBNORMAL : EXP_CLASS::NORMAL;
            }
        }

        reconstruct( x, x_exp_class, x_exp, x_sign );

    } else {
        // FIXED: do simple integer add and check for overflow (this is the one place where fixed is cheaper than float)
        //
        bool x_sign = x < 0;
        bool y_sign = y < 0;
             x     += y;
        T    sign_mask = x >> (_w - 1);
        cassert( sign_mask == T(0) || sign_mask == T(-1), "add caused overflow" );
    }
    return x;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::scalbn( const T& _x, int ls, bool is_final ) const
{
    if ( is_final ) _log_2i( scalbn, _x, T(ls) );
    T         x = _x;
    EXP_CLASS x_exp_class;
    int32_t   x_exp;
    bool      x_sign;

    //-----------------------------------------------------
    // Deconstruct
    // x_exp += ls
    // Reconstruct
    //-----------------------------------------------------
    deconstruct( x, x_exp_class, x_exp, x_sign );
    x_exp += ls;
    reconstruct( x, x_exp_class, x_exp, x_sign );
    if ( is_final ) x = rfrac( x );
    return x;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::scalbn( const T& x, int ls ) const
{
    return scalbn( x, ls, true );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::scalbnn( const T& x, int rs ) const
{
    return scalbn( x, -rs, true );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::ldexp( const T& x, int y ) const
{
    return scalbn( x, y );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::fma_fda( bool is_fma, const T& _x, const T& _y, const T& addend, bool is_final ) const
{
    bool have_addend = !iszero( addend );
    if ( is_final ) {
        if (  have_addend &&  is_fma ) _log_3( fma, _x, _y, addend );
        if (  have_addend && !is_fma ) _log_3( fda, _y, _x, addend );   // note: div convention is y/x, not x/y
        if ( !have_addend &&  is_fma ) _log_2( mul, _x, _y );
        if ( !have_addend && !is_fma ) _log_2( div, _y, _x );
    }
    T x = _x;
    T y = _y;
    std::string kind = is_fma ? "fma" : "fda";
    if ( debug ) std::cout << kind << " begin: x_orig=" << _to_flt(x, is_final) << 
                                        " y_orig=" << _to_flt(y, is_final) << 
                                        " addend=" << _to_flt(addend, is_final) << "\n";
    EXP_CLASS x_exp_class;
    EXP_CLASS y_exp_class;
    int32_t   x_exp;
    int32_t   y_exp;
    EXP_CLASS rr_exp_class;
    int32_t   rr_exp = 0;
    bool      rr_sign;
    reduce_mul_div_args( is_fma, x, y, x_exp_class, x_exp, y_exp_class, y_exp, rr_sign );

    // check for special cases
    //
    bool do_rest = false;
    T rr;
    if ( x_exp_class == EXP_CLASS::NOT_A_NUMBER ) {
        // x is the answer
        rr_exp_class = EXP_CLASS::NOT_A_NUMBER;
        rr           = x;

    } else if ( y_exp_class == EXP_CLASS::NOT_A_NUMBER ) {
        // y is the answer
        rr_exp_class = EXP_CLASS::NOT_A_NUMBER;
        rr           = y;

    } else if ( x_exp_class == EXP_CLASS::INFINITE ) {
        if ( is_fma ) {
            if ( y_exp_class == EXP_CLASS::ZERO ) {
                rr_exp_class = EXP_CLASS::NOT_A_NUMBER;
            } else {
                rr_exp_class = EXP_CLASS::INFINITE;
            }
        } else {
            if ( y_exp_class == EXP_CLASS::INFINITE ) {
                rr_exp_class = EXP_CLASS::NOT_A_NUMBER;
            } else {
                rr_exp_class = EXP_CLASS::ZERO;
            }
        }

    } else if ( y_exp_class == EXP_CLASS::INFINITE ) {
        if ( is_fma ) {
            if ( x_exp_class == EXP_CLASS::ZERO ) {
                rr_exp_class = EXP_CLASS::NOT_A_NUMBER;
            } else { 
                rr_exp_class = EXP_CLASS::INFINITE;
            }
        } else {
            if ( y_exp_class == EXP_CLASS::INFINITE ) {
                rr_exp_class = EXP_CLASS::NOT_A_NUMBER;
            } else {
                rr_exp_class = EXP_CLASS::INFINITE;
            }
        }

    } else if ( x_exp_class == EXP_CLASS::ZERO ) {
        if ( is_fma ) {
            rr_exp_class = EXP_CLASS::ZERO;
        } else {
            rr_exp_class = EXP_CLASS::INFINITE;
        }

    } else if ( y_exp_class == EXP_CLASS::ZERO ) {
        rr_exp_class = EXP_CLASS::ZERO;

    } else {
        T xx, yy, zz;
        if ( is_fma ) {
            linear_rotation( x, _zero, y, xx, rr, zz );
        } else {
            linear_vectoring( x, y, _zero, xx, yy, rr );
        }
        rr_exp_class = x_exp_class;
        rr_exp = y_exp + (is_fma ? x_exp : -x_exp);
        if ( rr == 0 ) rr_exp_class = EXP_CLASS::ZERO;
        do_rest = true;
    }

    reconstruct( rr, rr_exp_class, rr_exp, rr_sign );
    if ( debug ) std::cout << kind << " mid: rr=" << _to_flt(rr, is_final) << "\n";
    
    if ( do_rest ) {
        if ( have_addend ) rr = add( rr, addend, false );
        if ( is_final ) rr = rfrac( rr );
    }

    if ( debug ) std::cout << kind << " end: x_orig=" << _to_flt(_x, is_final) << 
                              " y_orig=" << _to_flt(_y, is_final) << 
                              " have_addend=" << have_addend << " addend=" << _to_flt(addend, is_final) << 
                              " x_reduced=" << _to_flt(x, false, true) << " y_reduced=" << _to_flt(y, false, true) << 
                              " " << kind << "=" << _to_flt(rr, is_final) << " (0x" << std::hex << rr << std::dec << 
                              ") x_exp_class=" << to_str(x_exp_class) << " x_exp=" << x_exp << 
                              " y_exp_class=" << to_str(y_exp_class) << " y_exp=" << y_exp << "\n";
    return rr;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::fma( const T& x, const T& y, const T& addend ) const
{
    return fma_fda( true, x, y, addend, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::mul( const T& x, const T& y ) const
{
    return fma( x, y, _zero );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::mul( const T& x, const T& y, bool is_final ) const
{
    return fma_fda( true, x, y, _zero, is_final );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::mulc( const T& x, const T& c, bool is_final ) const
{
    if ( is_final ) _log_2i( mulc, x, c );
    T r = mul( x, c, false );     // later, change this to minimum shifts and adds
    if ( is_final ) r = rfrac( r );
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::mulc( const T& x, const T& c ) const
{
    return mulc( x, c, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sqr( const T& x, bool is_final ) const
{
    if ( is_final ) _log_1( sqr, x );
    T r = mul( x, x, false );
    if ( is_final ) r = rfrac( r );
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sqr( const T& x ) const
{
    return sqr( x, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::fda( const T& _y, const T& _x, const T& addend ) const
{
    return fma_fda( false, _x, _y, addend, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::div( const T& y, const T& x ) const
{
    return fda( y, x, _zero );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::div( const T& y, const T& x, bool is_final ) const
{
    return fma_fda( false, x, y, _zero, is_final );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::remainder( const T& y, const T& x ) const
{
    return div( y, x );                 // FIXIT: placeholder
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::fmod( const T& y, const T& x ) const
{
    return div( y, x );                 // FIXIT: placeholder
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::remquo( const T& y, const T& x, int * quo ) const
{
    return div( y, x );                 // FIXIT: placeholder
    *quo = 0;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::rcp( const T& x ) const
{
    return div( _one, x );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::sqrt( const T& _x, bool is_final ) const
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
    if ( is_final ) _log_1( sqrt, _x );
    T x = _x;
    if ( debug ) std::cout << "sqrt begin: x_orig=" << _to_flt(_x, is_final) << "\n";
    EXP_CLASS x_exp_class;
    int32_t   x_exp;
    bool      x_sign;
    reduce_sqrt_arg( x, x_exp_class, x_exp, x_sign );

    // check for special cases
    //
    bool do_rest = false;
    if ( x_sign || x_exp_class == EXP_CLASS::NOT_A_NUMBER ) {
        // quiet NaN
        x_exp_class = EXP_CLASS::NOT_A_NUMBER;
        x |= _quiet_NaN_fxd;
        reconstruct( x, x_exp_class, x_exp, x_sign );

    } else if ( x_exp_class == EXP_CLASS::ZERO || x_exp_class == EXP_CLASS::INFINITE ) {
        // x is the answer

    } else {
        T xx, yy;
        hyperbolic_vectoring_xy( x+_one_fxd, x-_one_fxd, xx, yy );        // gain*sqrt((s+1)^2 - (s-1)^2)
        reconstruct( xx, x_exp_class, 0, x_sign );
        x = mulc( xx, _hyperbolic_vectoring_one_over_gain, false );
        do_rest = true;
    }

    if ( debug ) std::cout << "sqrt mid: x=" << _to_flt(x, is_final) << " x_exp=" << x_exp << "\n";

    if ( do_rest ) {
        x = scalbn( x, x_exp, false );                             // log2(p)/2 - 1
        if ( debug ) std::cout << "sqrt mid2: x=" << _to_flt(x, is_final) << "\n";
        if ( is_final ) x = rfrac( x );
    }

    if ( debug ) std::cout << "sqrt end: x_orig=" << _to_flt(_x, is_final) << " sqrt=" << _to_flt(x, is_final) << "\n";
    return x;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sqrt( const T& x ) const
{ 
    return sqrt( x, true );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::rsqrt( const T& x ) const
{ 
    //-----------------------------------------------------
    // 1.0 / sqrt( x )
    //-----------------------------------------------------
    _log_1( rsqrt, x );
    T sq = sqrt( x, false );
    T r = div( _one, sq );
    if ( debug ) std::cout << "rsqrt end: x_orig=" << _to_flt(x) << " sqrt=" << _to_flt(sq) << " r=" << _to_flt(r, false) << "\n";
    return r;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::rsqrt_orig( const T& x ) const
{ 
    //-----------------------------------------------------
    // x^(-1/2) = exp( log(x) / -2 );
    //-----------------------------------------------------
    _log_1( rsqrt, x );
    T lg = log( x, false );
    T lg_div_m2 = scalbn( neg( lg, false ), -1, false );
    T r = exp( lg_div_m2, false );
    r = rfrac( r );
    if ( debug ) std::cout << "rsqrt end: x_orig=" << _to_flt(x) << " lg=" << _to_flt(lg) <<
                              " lg_div_m2=" << _to_flt(lg_div_m2) << " r=" << _to_flt(r, false) << "\n";
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::cbrt( const T& _x ) const
{ 
    // x^(1/3) = exp( log(x) / 3 );
    _log_1( cbrt, _x );
    T x = _x;
    bool sign = signbit( x );
    if ( sign ) x = neg( x, false );
    T r = exp( mulc( log( x, false ), _third, false ) );
    if ( sign ) r = neg( x, false );
    r = rfrac( r );
    if ( debug ) std::cout << "cbrt end: x_orig=" << _to_flt(_x, false) << 
                              " x_reduced=" << _to_flt(x, false, true, false) << " r=" << _to_flt(r, false) << "\n";
    return r;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::rcbrt( const T& _x ) const
{ 
    // x^(-1/3) = exp( log(x) / -3 );
    _log_1( rcbrt, _x );
    T x = _x;
    bool sign = signbit( x );
    if ( sign ) x = neg( x, false );
    T r = exp( mulc( log( x, false ), _neg_third, false ) ); 
    if ( sign ) r = neg( x, false );
    r = rfrac( r );
    if ( debug ) std::cout << "rcbrt end: x_orig=" << _to_flt(_x, false) << 
                              " x_reduced=" << _to_flt(x, false, true, false) << " r=" << _to_flt(r, false) << "\n";
    return r;
}

template< typename T, typename FLT >
inline int Cordic<T,FLT>::compare( const T& _x, const T& _y ) const
{
    T         x = _x;
    T         y = _y;
    EXP_CLASS x_exp_class;
    EXP_CLASS y_exp_class;
    int32_t   x_exp;
    int32_t   y_exp;
    bool      x_sign;
    bool      y_sign;
    deconstruct( x, x_exp_class, x_exp, x_sign );
    deconstruct( y, y_exp_class, y_exp, y_sign );

    if ( x_exp_class == EXP_CLASS::ZERO && y_exp_class == EXP_CLASS::ZERO ) {
        return 0;
    }
    if ( x_exp_class == EXP_CLASS::NOT_A_NUMBER || y_exp_class == EXP_CLASS::NOT_A_NUMBER) {
        return ((x_exp_class == EXP_CLASS::NOT_A_NUMBER) ? -2 : 0) + ((y_exp_class == EXP_CLASS::NOT_A_NUMBER) ? -2 : 0); // -2 or -4
    }
    if ( x_exp_class == EXP_CLASS::INFINITE && y_exp_class == EXP_CLASS::INFINITE ) {
        return (x_sign == y_sign) ? -2 : x_sign ? -1 : 1;
    }
    if ( x_exp_class == EXP_CLASS::INFINITE ) {
        return x_sign ? -1 : 1;
    }
    if ( y_exp_class == EXP_CLASS::INFINITE ) {
        return y_sign ? 1 : -1;
    }
    if ( x_exp_class == EXP_CLASS::ZERO ) {
        return y_sign ? 1 : -1;
    }
    if ( y_exp_class == EXP_CLASS::ZERO ) {
        return x_sign ? -1 : 1;
    }         
    if ( x_sign ^ y_sign ) {
        return x_sign ? -1 : 1;                                            
    }

    cassert( x_exp_class == EXP_CLASS::NORMAL || x_exp_class == EXP_CLASS::SUBNORMAL, "compare() unexpected x_exp_class at this point" );
    cassert( y_exp_class == EXP_CLASS::NORMAL || y_exp_class == EXP_CLASS::SUBNORMAL, "compare() unexpected y_exp_class at this point" );

    if ( x_exp == y_exp ) {
        if ( x == y ) return 0;
        if ( x >  y ) return x_sign ? -1 : 1;
        return x_sign ? 1 : -1;
    } else if ( x_exp > y_exp ) {
        return x_sign ? -1 : 1;
    } else {
        return x_sign ? 1 : -1;
    }
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isgreater( const T& x, const T& y ) const
{
    _log_2( isgreater, x, y );
    bool b = compare( x, y ) == 1;
    if ( debug ) std::cout << "isgreater end: x_orig=" << _to_flt(x) << " y_orig=" << _to_flt(y) << " isgreater=" << int(b) << "\n";
    return b;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isgreaterequal( const T& x, const T& y ) const
{
    _log_2( isgreaterequal, x, y );
    bool b = compare( x, y ) >= 0;
    if ( debug ) std::cout << "isgreaterequal end: x_orig=" << _to_flt(x) << " y_orig=" << _to_flt(y) << " isgreaterequal=" << int(b) << "\n";
    return b;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isless( const T& x, const T& y ) const
{
    _log_2( isless, x, y );
    bool b = compare( x, y ) == -1;
    if ( debug ) std::cout << "isless end: x_orig=" << _to_flt(x) << " y_orig=" << _to_flt(y) << " isless=" << int(b) << "\n";
    return b;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::islessequal( const T& x, const T& y ) const
{
    _log_2( islessequal, x, y );
    int cmp = compare( x, y );
    bool b = cmp == -1 || cmp == 0;
    if ( debug ) std::cout << "islessequal end: x_orig=" << _to_flt(x) << " y_orig=" << _to_flt(y) << " islessequal=" << int(b) << "\n";
    return b;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::islessgreater( const T& x, const T& y ) const
{
    _log_2( islessgreater, x, y );
    int cmp = compare( x, y );
    bool b = cmp == -1 || cmp == 1;
    if ( debug ) std::cout << "islessgreater end: x_orig=" << _to_flt(x) << " y_orig=" << _to_flt(y) << " islessgreater=" << int(b) << "\n";
    return b;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isunordered( const T& x, const T& y ) const
{
    _log_2( isunordered, x, y );
    bool b = compare( x, y ) <= -2;   // either is a NaN
    if ( debug ) std::cout << "isunordered end: x_orig=" << _to_flt(x) << " y_orig=" << _to_flt(y) << " isunordered=" << int(b) << "\n";
    return b;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isunequal( const T& x, const T& y ) const
{
    _log_2( isunequal, x, y );
    int cmp = compare( x, y );
    bool b = cmp != 0 && cmp != -4;   // two NaNs returns false
    if ( debug ) std::cout << "isunequal end: x_orig=" << _to_flt(x) << " y_orig=" << _to_flt(y) << " isunequal=" << int(b) << "\n";
    return b;
}

template< typename T, typename FLT >
inline bool Cordic<T,FLT>::isequal( const T& x, const T& y ) const
{
    _log_2( isequal, x, y );
    int cmp = compare( x, y );
    bool b = cmp == 0 || cmp == -4;   // two NaNs returns true
    if ( debug ) std::cout << "isequal end: x_orig=" << _to_flt(x) << " y_orig=" << _to_flt(y) << " isequal=" << int(b) << "\n";
    return b;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::fdim( const T& x, const T& y ) const
{ 
    _log_2( fdim, x, y );
    T r = isgreaterequal( x, y ) ? sub( x, y, true ) : _zero;
    return r;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::fmax( const T& x, const T& y ) const
{ 
    _log_2( fmax, x, y );
    T r = isgreaterequal( x, y ) ? x : y;
    return r;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::fmin( const T& x, const T& y ) const
{ 
    _log_2( fmin, x, y );
    T r = isless( x, y ) ? x : y;
    return r;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::exp( const T& _x, bool is_final, FLT b ) const
{ 
    //-----------------------------------------------------
    // Identities:
    //     exp(x)    = sinh(x) + cosh(x)
    //     pow(b,x)  = exp2(log2(b) * x)
    //     exp2(i+f) = exp2(i) * exp2(f)     i = integer part, f = fractional remainder
    //               = exp2(f) << i    
    //               = exp(log(2)*f) << i
    //
    // Strategy:
    //     Call reduce_exp_arg() to get i and x=log(b)*f.
    //     Call hyperbolic_rotation() to get sinh(x) + cosh(x) in one shot.
    //-----------------------------------------------------
    if ( debug ) std::cout << "exp begin: x_orig=" << _to_flt(_x, is_final) << " b=" << b << "\n";
    if ( is_final ) _log_1( exp, _x );
    T x = _x;
    int32_t i;
    EXP_CLASS x_exp_class;
    bool x_sign;
    reduce_exp_arg( b, x, i, x_exp_class, x_sign ); 

    bool do_rest = false;
    if ( x_exp_class == EXP_CLASS::ZERO ) {
        // one is the answer
        x = _one;

    } else if ( x_exp_class == EXP_CLASS::INFINITE ) {
        // 0 or +inf
        x = x_sign ? _zero : _x;

    } else if ( x_exp_class == EXP_CLASS::NOT_A_NUMBER ) {
        // x is the answer
        x = _x;
    
    } else {
        T xx, yy, zz;
        hyperbolic_rotation( _hyperbolic_rotation_one_over_gain_fxd, _hyperbolic_rotation_one_over_gain_fxd, x, xx, yy, zz );
        if ( debug ) std::cout << "exp mid: b=" << b << " x_orig=" << _to_flt(_x, is_final) << 
                                  " i=" << i << " exp(log(2)*f)=" << _to_flt(xx, false, true) << "\n";
        reconstruct( xx, x_exp_class, 0, false );
        x = scalbn( xx, i, false );
        if ( is_final ) x = rfrac( x );
    }

    if ( debug ) std::cout << "exp end: x_orig=" << _to_flt(_x, is_final) << " b=" << b << " exp=" << _to_flt(x, is_final) << "\n";
    return x;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::exp( const T& x ) const
{ 
    return exp( x, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::expm1( const T& x ) const
{ 
    //-----------------------------------------------------
    // Compute without rounding, then round.
    //-----------------------------------------------------
    return sub( exp( x, false ), one() );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::expc( const FLT& b, const T& x ) const
{ 
    return exp( x, true, b );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::exp2( const T& x ) const
{ 
    return expc( 2.0, x );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::exp10( const T& x ) const
{ 
    return expc( 10.0, x );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pow( const T& b, const T& x ) const
{ 
    if ( debug ) std::cout << "pow begin: b=" << _to_flt(b) << " x=" << _to_flt(x) << "\n";
    _log_2( pow, b, x );
    T lg_b = log( b, false );
    T m  = mul( x, lg_b, false );
    T r = exp( m, false );
    r = rfrac( r );
    if ( debug ) std::cout << "pow end: b=" << _to_flt(b) << " x=" << _to_flt(x) << " pow=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::log( const T& _x, bool is_final ) const
{ 
    //-----------------------------------------------------
    // log(x) = 2*atanh2(x-1, x+1);
    //-----------------------------------------------------
    if ( debug ) std::cout << "log begin: x_orig=" << _to_flt(_x) << " is_final=" << is_final << "\n";
    if ( is_final ) _log_1( log, _x );
    T x = _x;
    EXP_CLASS x_exp_class;
    bool x_sign;
    T addend;
    reduce_log_arg( x, x_exp_class, x_sign, addend );                        // does not deconstruct value
    T r;
    if ( x_exp_class == EXP_CLASS::ZERO ) {
        // -inf
        r = ninfinity();
        if ( debug ) std::cout << "log end: x_orig=" << _to_flt(_x) << " log=-inf" << "\n";

    } else if ( x_exp_class == EXP_CLASS::NOT_A_NUMBER ) {
        // _x is the answer
        r = _x;
        if ( debug ) std::cout << "log end: x_orig=" << _to_flt(_x) << " log=nan" << "\n";

    } else if ( x_sign ) {
        // < 0 => NaN
        r = quiet_NaN();
        if ( debug ) std::cout << "log end: x_orig=" << _to_flt(_x) << " log=nan" << "\n";

    } else if ( x_exp_class == EXP_CLASS::INFINITE ) {
        // x is already infinity
        r = _x;
        if ( debug ) std::cout << "log end: x_orig=" << _to_flt(_x) << " log=inf" << "\n";

    } else {
        T x_m1 = sub( x, _one, false );
        T x_p1 = add( x, _one, false );
        T dv   = div( x_m1, x_p1, false );
        T lg1  = atanh2( dv, _one, false, true );
        T lg2  = scalbn( lg1, 1, false );
          r    = add( lg2, addend, false );
        if ( is_final ) r = rfrac( r );
        if ( debug ) std::cout << "log end: x_orig=" << _to_flt(_x) << " x_m1=" << _to_flt(x_m1) << " x_p1=" << _to_flt(x_p1) <<
                                          " dv=" << _to_flt(dv) << " lg1=" << _to_flt(lg1) << " lg2=" << _to_flt(lg2) << " addend=" << _to_flt(addend) <<
                                          " log=" << _to_flt(r) << "\n";
    }
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::log( const T& _x ) const
{ 
    return log( _x, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::log1p( const T& _x, bool is_final ) const
{ 
    return log( add( _x, _one, is_final ), is_final );   
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::log1p( const T& _x ) const
{ 
    return log1p( _x, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::log( const T& x, const T& b ) const
{ 
    _log_2( logn, x, b );
    T lgx = log( x, false );
    T lgb = log( b, false );
    T r = div( lgx, lgb, false );
    r = rfrac( r );
    if ( debug ) std::cout << "log: b=" << _to_flt(b) << " x=" << _to_flt(x) << " log=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::logc( const T& x, const FLT& b ) const
{ 
    _log_2f( logc, x, b );
    const FLT  one_over_log_b_f = FLT(1) / std::log( b );
    const T    one_over_log_b   = to_t( one_over_log_b_f );
          T    log_x            = log( x, false );
    T r = mulc( log_x, one_over_log_b, false );
    r = rfrac( r );
    if ( debug ) std::cout << "logc: b=" << _to_flt(b) << " x=" << _to_flt(x) << " reduced_x=" << _to_flt(x, false, true) << " log=" << _to_flt(r) << "\n";
    return r;
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
inline T Cordic<T,FLT>::deg2rad( const T& x ) const
{
    _log_1( deg2rad, x );
    T _180 = to_t( FLT(180) );
    T r = mulc( x, _180, false );
      r = mulc( r, _one_div_pi, false );
      r = rfrac( r );
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::rad2deg( const T& x ) const
{
    _log_1( rad2deg, x );
    T ONE_DIV_180 = to_t( FLT(1) / FLT(180) );
    T r = mulc( x, _pi, false );
      r = mulc( r, ONE_DIV_180, false );
      r = rfrac( r );
    return r;
}

template< typename T, typename FLT >
void Cordic<T,FLT>::sincos( bool times_pi, const T& _x, T& si, T& co, bool is_final, bool need_si, bool need_co, const T * _r ) const             
{ 
    if ( is_final ) {
        if ( _r != nullptr ) {
            if ( times_pi ) {
                _log_4( sinpicospi, _x, si, co, *_r );
            } else {
                _log_4( sincos, _x, si, co, *_r );
            }
        } else {
            if ( times_pi ) {
                _log_3( sinpicospi, _x, si, co );
            } else {
                _log_3( sincos, _x, si, co );
            }
        }
    }
    T x = _x;

    //-----------------------------------------------------
    // reduce_sincos_arg() will get x in the range 0 .. PI/2 and tell us the quadrant.
    // It will then check if x is still > PI/4 and, if so, subtract PI/4 and
    // set did_minus_pi_div_4.  If did_minus_pi_div_4 is true, then we need
    // to do some adjustments below after the cordic routine completes.
    //-----------------------------------------------------
    uint32_t quadrant;
    EXP_CLASS x_exp_class;
    bool x_sign;
    bool did_minus_pi_div_4;
    reduce_sincos_arg( times_pi, x, quadrant, x_exp_class, x_sign, did_minus_pi_div_4 );

    // check for special cases
    //
    if ( x_exp_class == EXP_CLASS::ZERO ) {
        // si = 0
        // co = 1
        if ( need_si ) si = _x;
        if ( need_co ) co = signbit(_x) ? _neg_one : _one;

    } else if ( x_exp_class == EXP_CLASS::NOT_A_NUMBER || x_exp_class == EXP_CLASS::INFINITE ) {
        // NaN
        x_exp_class = EXP_CLASS::NOT_A_NUMBER;
        x |= _quiet_NaN_fxd;
        reconstruct( x, x_exp_class, 0, x_sign );
        if ( need_si ) si = x;
        if ( need_co ) co = x;

    } else {
        T zz;
        circular_rotation( _circular_rotation_one_over_gain_fxd, _zero, x, co, si, zz );

        if ( si < 0 ) si = 0;                   // FIXIT: temporary hack when x is tiny
        if ( co < 0 ) co = 0;
        reconstruct( si, (si == 0) ? EXP_CLASS::ZERO : EXP_CLASS::NORMAL, 0, false );
        reconstruct( co, (co == 0) ? EXP_CLASS::ZERO : EXP_CLASS::NORMAL, 0, false );

        //-----------------------------------------------------
        // If did_minus_pi_div_4 is true, then we need to perform this
        // modification for sin and cos:
        //
        // sin(x+PI/4) = sqrt(2)/2 * ( sin(x) + cos(x) )
        // cos(x+PI/4) = sqrt(2)/2 * ( cos(x) - sin(x) )
        //-----------------------------------------------------
        if ( did_minus_pi_div_4 ) {
            T si_new = mulc( add( si, co, false ), _sqrt2_div_2, false );
            T co_new = mulc( sub( co, si, false ), _sqrt2_div_2, false );
            si = si_new;
            co = co_new;
        }

        //-----------------------------------------------------
        // Next, make adjustments for the quadrant.
        //-----------------------------------------------------
        if ( quadrant&1 ) {
            T tmp = co;
            co = si;
            si = tmp;
        }
        if ( need_si && (x_sign ^ (quadrant >= 2)) )                  si = neg( si, false );
        if ( need_co && (         (quadrant == 1) || quadrant == 2) ) co = neg( co, false );

        if ( _r != nullptr ) {
            //-----------------------------------------------------
            // Keep it simple for now.
            // Don't attempt to include r in the circular_rotation().
            //-----------------------------------------------------
            if ( need_si ) si = mul( si, *_r, false );
            if ( need_co ) co = mul( co, *_r, false );
        }

        if ( is_final ) {
            if ( need_si ) si = rfrac( si );
            if ( need_co ) co = rfrac( co );
        }
    }
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sin( const T& x, const T * r ) const
{ 
    if ( r != nullptr ) {
        _log_2( sin, x, *r );
    } else {
        _log_1( sin, x );
    }
    T si;
    T co;
    sincos( false, x, si, co, false, true, false, r );
    si = rfrac( si );
    if ( debug ) std::cout << "sin end: x_orig=" << _to_flt(x) << " sin=" << _to_flt(si) << "\n";
    return si;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::cos( const T& x, const T * r ) const
{ 
    if ( r != nullptr ) {
        _log_2( cos, x, *r );
    } else {
        _log_1( cos, x );
    }
    T si;
    T co;
    sincos( false, x, si, co, false, false, true, r );
    co = rfrac( co );
    if ( debug ) std::cout << "cos end: x_orig=" << _to_flt(x) << " cos=" << _to_flt(co) << "\n";
    return co;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::sincos( const T& x, T& si, T& co, const T * r ) const             
{
    sincos( false, x, si, co, true, true, true, r );
    if ( debug ) std::cout << "sincos end: x_orig=" << _to_flt(x) << " sin=" << _to_flt(si) << " cos=" << _to_flt(co) << 
                              " r=" << ((r != nullptr) ? _to_flt(*r) : 1.0) << "\n";
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::tan( const T& x ) const
{ 
    _log_1( tan, x );
    T si, co;
    sincos( false, x, si, co, false, true, true, nullptr );
    T r = div( si, co, false );
    r = rfrac( r );
    if ( debug ) std::cout << "tan end: x_orig=" << _to_flt(x) << " tan=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sinpi( const T& x, const T * r ) const
{ 
    if ( r != nullptr ) {
        _log_2( sinpi, x, *r );
    } else {
        _log_1( sinpi, x );
    }
    T si;
    T co;
    sincos( true, x, si, co, false, true, false, r );
    si = rfrac( si );
    if ( debug ) std::cout << "sinpi end: x_orig=" << _to_flt(x) << " sinpi=" << _to_flt(si) << "\n";
    return si;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::cospi( const T& x, const T * r ) const
{ 
    if ( r != nullptr ) {
        _log_2( cospi, x, *r );
    } else {
        _log_1( cospi, x );
    }
    T si;
    T co;
    sincos( true, x, si, co, false, false, true, r );
    co = rfrac( co );
    if ( debug ) std::cout << "cospi end: x_orig=" << _to_flt(x) << " cospi=" << _to_flt(co) << "\n";
    return co;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::sinpicospi( const T& x, T& si, T& co, const T * r ) const             
{
    sincos( true, x, si, co, true, true, true, r );
    if ( debug ) std::cout << "sinpicospi end: x_orig=" << _to_flt(x) << " sinpi=" << _to_flt(si) << " cospi=" << _to_flt(co) << 
                              " r=" << ((r != nullptr) ? _to_flt(*r) : 1.0) << "\n";
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::tanpi( const T& x ) const
{ 
    _log_1( tan, x );
    T si, co;
    sincos( true, x, si, co, false, true, true, nullptr );
    T r = div( si, co, false );
    r = rfrac( r );
    if ( debug ) std::cout << "tanpi end: x_orig=" << _to_flt(x) << " tanpi=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::asin( const T& x ) const
{ 
    _log_1( asin, x );
    T nh = hypoth( _one, x, false );
    T r = atan2( x, nh, false, false, nullptr );
    r = rfrac( r );
    if ( debug ) std::cout << "asin end: x_orig=" << _to_flt(x) << " asin=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::acos( const T& x ) const
{ 
    _log_1( acos, x );
    T nh = hypoth( _one, x, false );
    T r = atan2( nh, x, false, false, nullptr );
    r = rfrac( r );
    if ( debug ) std::cout << "acos end: x_orig=" << _to_flt(x) << " acos=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::atan( const T& x ) const
{ 
    T r = atan2( x, _one, true, true, nullptr );
    if ( debug ) std::cout << "atan end: x_orig=" << _to_flt(x) << " atan=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::atan2( const T& y, const T& x ) const
{ 
    T r = atan2( y, x, true, false, nullptr );
    if ( debug ) std::cout << "atan2 end: y=" << _to_flt(y) << " x=" << _to_flt(x) << " atan2=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::atan2( const T& _y, const T& _x, bool is_final, bool x_is_one, T * r ) const
{ 
    if ( is_final ) _log_2( atan2, _y, _x );
    T y = _y;
    T x = _x;

    //-----------------------------------------------------
    // Identities:
    //     atan2(y,x)       = undefined                             if x == 0 && y == 0
    //     atan2(y,x)       = PI                                    if x <  0 && y == 0
    //     atan2(y,x)       = 2*atan(y / (sqrt(x^2 + y^2) + x))     if x >  0    
    //     atan2(y,x)       = 2*atan((sqrt(x^2 + y^2) - |x|) / y)   if x <= 0 && y != 0
    // Strategy:
    //     Look at signs and values.
    //     Return PI if we're done.
    //     Do 2*atan( y / (hypot(x, y) + x) )    if x > 0
    //     Do 2*atan( (hypot(x, y) - x) / y )    if x <= 0
    //     When using cordic atan2 for the latter, if the numerator is larger than
    //     the denominator, then use PI - atan(x/y)
    //-----------------------------------------------------
    if ( debug ) std::cout << "atan2 begin: y=" << _to_flt(y, is_final) << 
                                          " x=" << _to_flt(x, is_final) << 
                                          " x_is_one=" << x_is_one << "\n";
    if ( r != nullptr ) *r = hypot( _x, _y, false );  // FIXIT: optimize this later 

    bool      x_sign = signbit( x );
    bool      y_sign = signbit( y );
    EXP_CLASS x_exp_class = classify( x );
    EXP_CLASS y_exp_class = classify( y );
    bool      x_gt_0 = !x_sign && x_exp_class != EXP_CLASS::ZERO;
    
    // check for special cases
    if ( x_exp_class == EXP_CLASS::NOT_A_NUMBER ) {
        // x is the answer 
        return _x;

    } else if ( y_exp_class == EXP_CLASS::NOT_A_NUMBER ) {
        // y is the answer

    } else if ( (x_exp_class == EXP_CLASS::ZERO     && y_exp_class == EXP_CLASS::ZERO) ||
                (x_exp_class == EXP_CLASS::INFINITE && y_exp_class == EXP_CLASS::INFINITE) ) {
        // NaN
        return quiet_NaN();

    } else if ( x_sign && y_exp_class == EXP_CLASS::ZERO ) {
        // PI
        return _pi;

    } else if ( (x_gt_0  && x_exp_class == EXP_CLASS::INFINITE) ||
                (!x_gt_0 && y_exp_class == EXP_CLASS::INFINITE) ) {
        // atan(0) == 0
        return (x_sign ^ y_sign) ? _neg_zero : _zero;
    }

    // normal case
    //
    T hpx = hypot( x, y, false );
      hpx = add( hpx, x_gt_0 ? x : neg( x, false ), false );
    if ( debug ) std::cout << "atan2 mid: y=" << _to_flt(_y) << " x=" << _to_flt(_x) << " hpx=" << _to_flt(hpx) << "\n";

    // decide which is the numerator and which is the denominator
    //
    if ( x_gt_0 ) {
        x = hpx;
    } else {
        x = y;
        y = hpx;
    }

    bool      rr_sign = signbit(x) ^ signbit(y);
    EXP_CLASS exp_class;
    int32_t   exp;
    bool      swapped;
    reduce_hypot_args( x, y, exp_class, exp, swapped );

    // check for special cases
    T rr;
    bool did_neg = false;
    if ( exp_class == EXP_CLASS::ZERO || exp_class == EXP_CLASS::INFINITE ) {
        // 0 or inf
        // 
        rr = 0;
    } else if ( exp_class == EXP_CLASS::NOT_A_NUMBER ) {
        // NaN
        rr = x | _quiet_NaN_fxd;

    } else {
        cassert( exp_class == EXP_CLASS::NORMAL || exp_class == EXP_CLASS::SUBNORMAL,
                 "atan2: unexpected exp_class=" + to_str(exp_class) );
        exp_class = EXP_CLASS::NORMAL;

        // reduce
        T xx, yy;
        circular_vectoring( x, y, _zero, xx, yy, rr );

        did_neg = rr < 0;
        if ( did_neg ) {
            rr = -rr;
            rr_sign = !rr_sign;
        }
    }

    reconstruct( rr, exp_class, 0, rr_sign );
    rr = scalbn( rr, 1, false );
    if ( swapped ) rr = sub( did_neg ? _neg_pi : _pi, rr, false );
    if ( is_final ) rr = rfrac( rr );
    if ( debug ) std::cout << "atan2 end: y=" << _to_flt(_y) << " x=" << _to_flt(_x) << " x_is_one=" << x_is_one << 
                              " swapped=" << swapped << " did_neg=" << did_neg << " atan2=" << _to_flt(rr) << 
                              " r=" << ((r != nullptr) ? _to_flt(*r) : _to_flt(_zero)) << "\n";
    return rr;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::polar_to_rect( const T& r, const T& a, T& x, T& y ) const
{
    _log_4( polar_to_rect, r, a, x, y );
    if ( debug ) std::cout << "polar_to_rect begin: r=" << _to_flt(r) << " a=" << _to_flt(a) << "\n";
    sincos( false, a, y, x, false, true, true, &r );
    x = rfrac( x );
    y = rfrac( y );
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::rect_to_polar( const T& x, const T& y, T& r, T& a ) const
{
    _log_4( rect_to_polar, x, y, r, a );
    if ( debug ) std::cout << "rect_to_polar begin: x=" << _to_flt(x) << " y=" << _to_flt(y) << "\n";
    a = atan2( y, x, false, false, &r );
    r = rfrac( r );
    a = rfrac( a );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hypot( const T& _x, const T& _y, bool is_final ) const
{
    if ( is_final ) _log_2( hypot, _x, _y );
    T x = _x;
    T y = _y;
    if ( debug ) std::cout << "hypot begin: x=" << _to_flt(x, is_final) << " y=" << _to_flt(y, is_final) << "\n";
    EXP_CLASS exp_class;
    int32_t   exp;
    bool      swapped;  // unused
    reduce_hypot_args( x, y, exp_class, exp, swapped );

    T xx, yy, zz;
    if ( exp_class == EXP_CLASS::NORMAL || exp_class == EXP_CLASS::SUBNORMAL ) {
        circular_vectoring_xy( x, y, xx, yy );
    } 
    reconstruct( xx, exp_class, exp, false );
    xx = mulc( xx, _circular_vectoring_one_over_gain, false );
    if ( is_final ) xx = rfrac( xx );
    if ( debug ) std::cout << "hypot end: x_orig=" << _to_flt(_x, is_final) << 
                                        " y_orig=" << _to_flt(_y, is_final) << " hypot=" << _to_flt(xx, is_final) << "\n";
    return xx;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hypot( const T& x, const T& y ) const
{
    return hypot( x, y, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hypoth( const T& x, const T& y, bool is_final ) const
{
    //-----------------------------------------------------
    // Identities:
    //     sqrt(x^2 - y^2) = sqrt((x+y)(x-y))
    // Strategy:
    //     Try this easy way, though I suspect there will be issues.
    //-----------------------------------------------------
    if ( is_final ) _log_2( hypoth, x, y );
    if ( debug ) std::cout << "hypoth begin: x=" << _to_flt(x, is_final) << " y=" << _to_flt(y, is_final) << "\n";
    T x_p_y = add( x, y, false );
    T x_m_y = sub( x, y, false );
    T r = sqrt( mul( x_p_y, x_m_y, false ), false ); // FIXIT: go back to using CORDIC hyperbolic core
    if ( is_final ) r = rfrac( r );
    if ( debug ) std::cout << "hypoth end: x_orig=" << _to_flt(x, is_final) << 
                                         " y_orig=" << _to_flt(y, is_final) << " hypoth=" << _to_flt(r, is_final) << "\n";
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hypoth( const T& x, const T& y ) const
{
    return hypoth( x, y, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sinh( const T& x, const T * r ) const
{ 
    if ( r != nullptr ) {
        _log_2( sinh, x, *r );
    } else {
        _log_1( sinh, x );
    }
    T sih;
    T coh;
    sinhcosh( x, sih, coh, false, true, false, r );
    sih = rfrac( sih );
    if ( debug ) std::cout << "sinh end: x_orig=" << _to_flt(x) << " sinh=" << _to_flt(sih) << "\n";
    return sih;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::cosh( const T& x, const T * r ) const
{ 
    if ( r != nullptr ) {
        _log_2( cosh, x, *r );
    } else {
        _log_1( cosh, x );
    }
    T sih;
    T coh;
    sinhcosh( x, sih, coh, true, false, true, r );
    coh = rfrac( coh );
    if ( debug ) std::cout << "cosh end: x_orig=" << _to_flt(x) << " cosh=" << _to_flt(coh) << "\n";
    return coh;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::sinhcosh( const T& x, T& sih, T& coh, const T * r ) const
{ 
    sinhcosh( x, sih, coh, true, true, true, r );
    if ( debug ) std::cout << "sinhcosh end: x_orig=" << _to_flt(x) << " sinh=" << _to_flt(sih) << " cosh=" << _to_flt(coh) << 
                              " r=" << ((r != nullptr) ? _to_flt(*r) : 1.0) << "\n";
}

template< typename T, typename FLT >
void Cordic<T,FLT>::sinhcosh( const T& _x, T& sih, T& coh, bool is_final, bool need_sih, bool need_coh, const T * _r ) const
{ 
    if ( is_final ) {
        if ( _r != nullptr ) {
            _log_4( sinhcosh, _x, sih, coh, *_r );
        } else {
            _log_3( sinhcosh, _x, sih, coh );
        }
    }

    T x = _x;

    //-----------------------------------------------------
    // sinh(x) = (exp(x) - 1/exp(x))/2
    // cosh(x) = (exp(x) + 1/exp(x))/2
    //
    // Currently, this appears more efficient than other techniques that don't use a LUT.
    //-----------------------------------------------------
    T expx = exp( x, false );
    T one_over_expx = div( _one, expx, false );
    if ( _r != nullptr ) expx = mul( expx, *_r, false );
    if ( need_sih ) {
        sih = sub( expx, one_over_expx, false );
        sih = scalbn( sih, -1, false );
    }
    if ( need_coh ) {
        coh = add( expx, one_over_expx, false );
        coh = scalbn( coh, -1, false );
    }

    if ( is_final ) {
        if ( need_sih ) sih = rfrac( sih );
        if ( need_coh ) coh = rfrac( coh );
    }
}

template< typename T, typename FLT >
T Cordic<T,FLT>::tanh( const T& x ) const
{ 
    _log_1( tanh, x );
    T sih, coh;
    sinhcosh( x, sih, coh, false, true, true, nullptr );
    T r = div( sih, coh, false );
    r = rfrac( r );
    if ( debug ) std::cout << "tanh end: x_orig=" << _to_flt(x) << " tanh=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::asinh( const T& x ) const
{ 
    _log_1( asinh, x );
    T h = hypot( x, _one, false );
    T r = log( add( x, h, false ), false );
    r = rfrac( r );
    if ( debug ) std::cout << "asinh end: x_orig=" << _to_flt(x) << " asinh=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::acosh( const T& x ) const
{ 
    _log_1( acosh, x );
    T hh = hypoth( x, _one, false );
    T r = log( add( x, hh, false ), false );
    r = rfrac( r );
    if ( debug ) std::cout << "acosh end: x_orig=" << _to_flt(x) << " acosh=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::atanh( const T& x ) const
{ 
    T r = atanh2( x, _one, true, true );
    if ( debug ) std::cout << "atanh end: x_orig=" << _to_flt(x) << " atanh=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::atanh2( const T& y, const T& x ) const             
{ 
    T r = atanh2( y, x, true, false );
    if ( debug ) std::cout << "atanh2 end: y=" << _to_flt(y) << " x=" << _to_flt(x) << " atanh2=" << _to_flt(r) << "\n";
    return r;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::atanh2( const T& _y, const T& _x, bool is_final, bool x_is_one ) const             
{ 
    if ( debug ) std::cout << "atanh2 begin: y=" << _to_flt( _y ) << " x=" << _to_flt( _x ) << "\n";

    if ( is_final ) _log_2( atanh2, _y, _x );
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
    EXP_CLASS y_exp_class;
    EXP_CLASS x_exp_class;
    int32_t y_exp;
    int32_t x_exp;
    bool y_sign;
    bool x_sign;

    deconstruct( y, y_exp_class, y_exp, y_sign );

    if ( x_is_one ) {
        x = _one_fxd;
        x_exp_class = EXP_CLASS::NORMAL;
        x_exp = 0;
        x_sign = false;
    } else {
        deconstruct( x, x_exp_class, x_exp, x_sign );
    }

    // check for other special cases
    //
    bool sign = y_sign ^ x_sign;
    int32_t exp = x_exp + y_exp;
    T r;
    if ( y_exp_class == EXP_CLASS::NOT_A_NUMBER || y_exp_class == EXP_CLASS::INFINITE || exp > 0 ) {
        r = y | _quiet_NaN_fxd;
        reconstruct( r, EXP_CLASS::NOT_A_NUMBER, _exp_mask, sign );

    } else if ( x_exp_class == EXP_CLASS::NOT_A_NUMBER || x_exp_class == EXP_CLASS::ZERO ) {
        r = x | _quiet_NaN_fxd;
        reconstruct( r, EXP_CLASS::NOT_A_NUMBER, _exp_mask, sign );

    } else if ( y_exp_class == EXP_CLASS::ZERO || x_exp_class == EXP_CLASS::INFINITE ) {
        // y/x == 0 => return +/- 0
        r = sign ? _neg_zero : _zero;

    } else {
        y >>= -exp;
        exp = 0;

        T xx, yy;
        hyperbolic_vectoring( x, y, _zero, xx, yy, r );
        if ( r < 0 ) {
            r = -r;
            sign = !sign;
        }

        reconstruct( r, EXP_CLASS::NORMAL, exp, sign );
        if ( is_final ) r = rfrac( r );
    }

    if ( debug ) std::cout << "atanh2: y=" << _to_flt( _y ) << " x=" << _to_flt( _x ) << " r=" << _to_flt( r ) << "\n";
    return r;
}

template< typename T, typename FLT >
typename Cordic<T,FLT>::EXP_CLASS Cordic<T,FLT>::classify( const T& _x ) const
{
    T x = _x;
    EXP_CLASS x_exp_class;
    int32_t   x_exp;
    bool      x_sign;
    deconstruct( x, x_exp_class, x_exp, x_sign );
    return x_exp_class;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::deconstruct( T& x, EXP_CLASS& x_exp_class, int32_t& x_exp, bool& sign, bool allow_debug ) const
{
    T x_orig = x;
    if ( _is_float ) {
        // FLOAT => already in the correct format
        //
        uint32_t exp_biased = (x >> _frac_guard_w) & _exp_mask;
                 sign       = (x >> (_w-1)) & 1;
                 x         &= _frac_guard_mask;     // mantissa without implicit 1.
        if ( debug && allow_debug ) std::cout << "\ndeconstruct: sign=" << sign << " exp_biased=" << exp_biased << 
                                                 " frac_guard=" << std::hex << x << std::dec << "\n";
        if ( exp_biased == 0 ) {
            if ( x == 0 ) {
                x_exp_class = EXP_CLASS::ZERO;
            } else {
                x_exp_class = EXP_CLASS::SUBNORMAL;     // FIXIT: probably best to normalize at this point
            }
        } else if ( exp_biased == _exp_mask ) {
            if ( x == 0 ) {
                x_exp_class = EXP_CLASS::INFINITE;
            } else {
                x_exp_class = EXP_CLASS::NOT_A_NUMBER;
            }
            x_exp       = 0;
        } else {
            x_exp_class = EXP_CLASS::NORMAL;
            x_exp       = int32_t(exp_biased) - _exp_bias; 
            x          |= _one_fxd;    // add in implicit 1.
        }
    } else {
        // FIXED => normalize to look like floating-point (even if we lose some bits right-shifting)
        //
        if ( x == 0 ) {
            x_exp_class = EXP_CLASS::ZERO;
            x_exp = 0;
        } else {
            x_exp_class = EXP_CLASS::NORMAL;
            x_exp = 0;
            sign = x < 0;
            if ( sign ) x = -x;

            while( x > _one_fxd )
            {
                x_exp++;
                x = (x >> 1) | (x & 1);    // must record sticky bit
            }
            while( x < _one_fxd )
            {
                x_exp--;
                x <<= 1;
            }
        }
    }
    if ( debug && allow_debug ) std::cout << "deconstruct: x_orig_f=" << _to_flt(x_orig) << " x_orig=" << std::hex << x_orig << 
                                             " x_reduced_f=" << _to_flt(x, false, true) << " x_reduced=" << x << 
                                             " x_exp_class=" + to_str(x_exp_class) << std::dec << " x_exp=" << x_exp << 
                                             " sign=" << sign << "\n";
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reconstruct( T& x, EXP_CLASS x_exp_class, int32_t x_exp, bool sign ) const
{
    T x_orig = x;
    if ( debug ) std::cout << "reconstruct: x_orig=" << std::hex << x << std::dec << " x_orig_f=" << _to_flt(x, false, true) <<
                              " x_exp_class=" + to_str(x_exp_class) << " x_exp=" << x_exp << " sign=" << sign << "\n";
    int32_t exp = 0;
    if ( _is_float ) {
        // FLOAT
        //
        T int_part = x >> _frac_guard_w;
        x &= _frac_guard_mask;                                  // lop off the implicit 1.
        switch( x_exp_class )
        {
            case EXP_CLASS::ZERO:
                exp = 0;
                x = 0;
                break;

            case EXP_CLASS::NORMAL:
            case EXP_CLASS::SUBNORMAL:     
                cassert( int_part != 0 || x != 0, "reconstruct() normal int_part or frac_part should be non-zero" );
                while( int_part == 0 ) 
                {
                    int_part = (x >> (_frac_guard_w-1)) & 1;
                    x = (x << 1) & _frac_guard_mask;
                    x_exp--;
                }
                while( int_part > 1 ) 
                {
                    x = (x >> 1) | ((int_part&1) << (_frac_guard_w-1)) | (x & 1);
                    int_part >>= 1;
                    x_exp++;
                }
                cassert( int_part == 1, "reconstruct() normal int_part should be exactly 1, int=" + std::to_string(int_part) + 
                                        + " frac=" + std::to_string(x) + " x_orig=" + std::to_string(x_orig) );

                exp = x_exp + _exp_bias;
                if ( exp < 0 ) {
                    // subnormal
                    exp = -exp;
                    T mask = (T(1) << exp) - 1;
                    x = (x >> exp) | ((x & mask) != 0);
                    exp = 0;
                } else if ( exp >= int32_t(_exp_mask) ) {
                    // infinity
                    exp = _exp_mask;
                    x = 0;
                }
                break;

            case EXP_CLASS::INFINITE:
                exp = _exp_mask;
                x = 0;
                break;
         
            case EXP_CLASS::NOT_A_NUMBER:
                exp = _exp_mask;
                if ( x == 0 ) x = 1;                                            // to distinguish from INF
                break;

            default:
                cassert( false, "reconstruct() encountered a bad x_exp_class" );
                break;
        }
        x |= T(sign) << (_exp_w+_frac_guard_w);
        x |= T(exp)  << (       _frac_guard_w);
    } else {
        // FIXED
        //
        cassert( x_exp_class != EXP_CLASS::INFINITE,     "unexpected x_exp_class=INFINITE for fixed-point number" );
        cassert( x_exp_class != EXP_CLASS::NOT_A_NUMBER, "unexpected x_exp_class=NOT_A_NUMBER for fixed-point number" );
        if ( x_exp_class == EXP_CLASS::ZERO || x == 0 ) {
            x = 0;
        } else {
            cassert( x_exp_class == EXP_CLASS::NORMAL || x_exp_class == EXP_CLASS::SUBNORMAL, "unexpected x_exp_class for fixed-point number" );
            if ( sign ) x = -x;
            if ( x_exp > 0 ) {
                //-----------------------------------------------------
                // Crap out if we overflow.
                //-----------------------------------------------------
                T sign_mask = (x < 0) ? (T(-1) << (_w - 1)) : 0;
                cassert( (x & sign_mask) == sign_mask, "scalbn() x overflowed even before shift, x_exp=" + std::to_string(x_exp) );
                bool x_overflow = (x & sign_mask) != sign_mask;
                x <<= x_exp; 
                cassert( (x & sign_mask) == sign_mask, "scalbn() shifted x overflowed, x_exp=" + std::to_string(x_exp) );
            } else if ( x_exp < 0 ) {
                //-----------------------------------------------------
                // Keep track of the sticky bit even if not rounding.
                //-----------------------------------------------------
                x_exp = -x_exp;
                T mask = (T(1) << x_exp) - 1;
                T set_sticky = (x & mask) != 0;
                x >>= x_exp;   
                x |= set_sticky;
            }
        }
    }
    if ( debug ) std::cout << "reconstruct: x_orig=" << std::hex << x << std::dec << " x_orig_f=" << _to_flt(x, false, true) <<
                              " x_exp_class=" + to_str(x_exp_class) << " x_exp=" << x_exp << " constructed_f=" << _to_flt(x, false) << "\n";
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_add_args( T& x, T& y, EXP_CLASS& x_exp_class, int32_t& x_exp, bool& x_sign, EXP_CLASS& y_exp_class, int32_t& y_exp, bool& y_sign ) const
{
    if ( debug ) std::cout << "reduce_add_args: x_orig=" << _to_flt(x) << " y_orig=" << _to_flt(y) << "\n";
    deconstruct( x, x_exp_class, x_exp, x_sign );
    deconstruct( y, y_exp_class, y_exp, y_sign );
    if ( x_exp == y_exp ||
         x_exp_class == EXP_CLASS::ZERO || x_exp_class == EXP_CLASS::INFINITE || x_exp_class == EXP_CLASS::NOT_A_NUMBER ||
         y_exp_class == EXP_CLASS::ZERO || y_exp_class == EXP_CLASS::INFINITE || y_exp_class == EXP_CLASS::NOT_A_NUMBER ) return;
   
    // rshift fraction with smaller exponent
    if ( x_exp > y_exp ) {
        uint32_t rs = x_exp - y_exp;
        bool set_sticky = (y & ((T(1) << rs)-1)) != 0;
        y >>= rs;
        if ( set_sticky ) y |= 1;
        y_exp = x_exp;
    } else {
        uint32_t rs = y_exp - x_exp;
        bool set_sticky = (x & ((T(1) << rs)-1)) != 0;
        x >>= rs;
        if ( set_sticky ) x |= 1;
        x_exp = y_exp;
    }
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_mul_div_args( bool is_mul, T& x, T& y, EXP_CLASS& x_exp_class, int32_t& x_exp, EXP_CLASS& y_exp_class, int32_t& y_exp, bool& sign ) const
{
    if ( debug ) std::cout << "reduce_mul_div_args: is_mul=" << is_mul << " x_orig=" << _to_flt(x) << " y_orig=" << _to_flt(y) << "\n";
    bool x_sign;
    bool y_sign;
    deconstruct( x, x_exp_class, x_exp, x_sign );
    deconstruct( y, y_exp_class, y_exp, y_sign );
    sign = x_sign ^ y_sign;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_sqrt_arg( T& x, EXP_CLASS& x_exp_class, int32_t& x_exp, bool& x_sign ) const
{
    //-----------------------------------------------------
    // Identities:
    //     sqrt(x) = sqrt((x+1)^2 - (x-1)^2) / 2          
    // Strategy:
    //     Factor x=p*s where p is a power-of-2 
    //     where log2(p) is even >= 2 and s is between 0.5..1.
    //     sqrt(x) = sqrt(p) * sqrt(s) = 2^(log2(p)/2) * sqrt((s+1)^2 - (s-1)^2)/2
    //     Use hyperbolic_vectoring() for sqrt((s+1)^2 - (s-1)^2).
    //     Then exp that by log2(p)/2 - 1.
    //-----------------------------------------------------
    T x_orig = x;
    deconstruct( x, x_exp_class, x_exp, x_sign );
    if ( debug ) std::cout << "reduce_sqrt_arg mid: x=" << _to_flt(x, false) << " x_exp=" << x_exp << "\n";
    if ( x_exp_class == EXP_CLASS::NORMAL ) {
        x = (x >> 1) | (x & 1);
        x_exp++;
        if ( (x_exp & 1) ) {
            // make it even
            x_exp++;
            x = (x >> 1) | (x & 1);  // must record sticky bit
        }
        x_exp = x_exp/2 - 1;
    }
    if ( debug ) std::cout << "reduce_sqrt_arg: x_orig=" << _to_flt(x_orig) << " x_reduced=s=" << _to_flt(x, false, true) << 
                              " x_exp_class=" << to_str(x_exp_class) << " x_exp=" << x_exp << "\n";
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_exp_arg( FLT b, T& x, int32_t& i, EXP_CLASS& x_exp_class, bool& x_sign ) const
{
    //-----------------------------------------------------
    // Identities:
    //     pow(b,x)  = exp2(log2(b) * x)
    //     exp2(i+f) = exp2(i) * exp2(f)     i = integer part, f = fractional remainder
    //               = exp2(f) << i    
    //               = exp(log(2)*f) << i
    //
    // Strategy:
    //     For exp2(i+f), choose i such at f is in -1..1.  Note that i can be negative.
    //     Multiply f by log(2) and return that as the returned x.
    //     And return i as the exp.
    //-----------------------------------------------------
    T x_orig = x;
    if ( debug ) std::cout << "reduce_exp_arg begin: b=" << b << " x_orig=" << _to_flt(x_orig) << "\n";
    x_exp_class = classify( x );
    if ( x_exp_class == EXP_CLASS::ZERO || x_exp_class == EXP_CLASS::INFINITE || x_exp_class == EXP_CLASS::NOT_A_NUMBER ) {
        x_sign = signbit( x );
        i = 0;
        return;
    }

    T log2_b = to_t( std::log2( b ), false );
    x = mulc( x, log2_b, false );

    // get integer and fraction parts, still encoded;
    // convert encoded ii to int32_t;
    // multiply fraction by log(2)
    T ii;
    T f = modf( x, &ii );
    if ( signbit( f ) ) {
        ii = sub( ii, _one, false );
        f  = add( f,  _one, false );
    }
    if ( debug ) std::cout << "reduce_exp_arg mid1: x=" << _to_flt(x) << " f=" << _to_flt(f) << " ii=" << _to_flt(ii) << "\n";
    i = _to_flt( ii );
    x = mulc( f, _log2, false );
    if ( debug ) std::cout << "reduce_exp_arg mid2: f*log2=" << _to_flt(x) << "\n";

    int32_t x_exp;
    deconstruct( x, x_exp_class, x_exp, x_sign );
    while ( x_exp < 0 ) 
    {
        // turn back into un-normalized fraction (shouldn't need to shift much)
        //
        x = (x >> 1) | (x & 1);
        x_exp++;
    }

    if ( debug ) std::cout << "reduce_exp_arg end: b=" << b << " x_orig=" << _to_flt(x_orig) << 
                              " ii=" << _to_flt(ii, false) << " i=" << i << " f=" << _to_flt(f, false) << 
                              " x_reduced=log(2)*f=" << _to_flt(x, false, true) << "\n";
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_log_arg( T& x, EXP_CLASS& x_exp_class, bool& x_sign, T& addend ) const
{
    //-----------------------------------------------------
    // log(x*y)         = log(x) + log(y)
    // log(b^i)         = i*log(b)
    // log(2^i + f)     = i*log(2) + log(f)                     i=integer f=remainder
    // 
    // Normalize x so that it's in 1.00 .. 2.00.
    // Then addend = i*log(2).
    //-----------------------------------------------------
    T x_orig = x;
    int32_t x_exp;
    deconstruct( x, x_exp_class, x_exp, x_sign );
    if ( x_exp_class == EXP_CLASS::ZERO || x_exp_class == EXP_CLASS::NOT_A_NUMBER || x_exp_class == EXP_CLASS::INFINITE ) return;

    if ( x_exp_class == EXP_CLASS::NORMAL ) {
        x_exp++;
        x = (x >> 1) | (x & 1);
        if ( debug ) std::cout << "reduce_log_mid: x_shifted=0x" << std::hex << x << std::dec << "\n";
    }
    addend = to_t( FLT(x_exp) * std::log(2) );
    reconstruct( x, x_exp_class, 0, false );
    if ( debug ) std::cout << "reduce_log_arg: x_orig=" << _to_flt(x_orig) << " x_reduced=" << _to_flt(x, false) <<
                                             " (0x" << std::hex << x << ")" << std::dec <<
                                             " addend=" << _to_flt(addend, false) << "\n";
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_hypot_args( T& x, T& y, EXP_CLASS& exp_class, int32_t& exp, bool& swapped ) const
{
    //-----------------------------------------------------
    // Must shift both x and y by max( x_exp, y_exp ).
    // If x or y is negative, it's fine to negate them
    // because squaring them anyway.
    //-----------------------------------------------------
    T x_orig = x;
    T y_orig = y;
    if ( debug ) std::cout << "reduce_hypot_args begin: xy_orig=[" << _to_flt(x_orig) << ", " << _to_flt(y_orig) << "]\n";
    EXP_CLASS x_exp_class;
    EXP_CLASS y_exp_class;
    int32_t x_exp;
    int32_t y_exp;
    bool x_sign;
    bool y_sign;
    swapped = false;

    deconstruct( x, x_exp_class, x_exp, x_sign );
    deconstruct( y, y_exp_class, y_exp, y_sign );

    // check for special cases
    //
    if ( x_exp_class == EXP_CLASS::NOT_A_NUMBER ) {
        exp_class = x_exp_class;
        exp       = x_exp;

    } else if ( y_exp_class == EXP_CLASS::NOT_A_NUMBER ) {
        exp_class = y_exp_class;
        exp       = y_exp;

    } else if ( x_exp_class == EXP_CLASS::INFINITE ) {
        exp_class = x_exp_class;
        exp       = x_exp;

    } else if ( y_exp_class == EXP_CLASS::INFINITE ) {
        exp_class = y_exp_class;
        exp       = y_exp;

    } else if ( x_exp_class == EXP_CLASS::ZERO && y_exp_class == EXP_CLASS::ZERO ) {
        x = 0;
        exp_class = EXP_CLASS::ZERO;
        exp = 0;
    
    } else {  
        exp_class = EXP_CLASS::NORMAL;
        if ( x_exp_class == EXP_CLASS::ZERO ) {
            x_exp = y_exp;
            x = 0;
        } else if ( y_exp_class == EXP_CLASS::ZERO ) {
            y_exp = x_exp;
            y = 0;
        }

        int32_t diff = x_exp - y_exp;
        if ( diff > 0 ) {
            T mask = (T(1) << diff) - 1;
            y = (y >> diff) | ((y & mask) != 0);
            exp = x_exp;
            if ( debug ) std::cout << "reduce_hypot_args mid: using x_exp=" << exp << " y_rshift=" << diff << "\n";
        } else if ( diff < 0 ) {
            diff = -diff;
            T mask = (T(1) << diff) - 1;
            x = (x >> diff) | ((x & mask) != 0);
            exp = y_exp;
            if ( debug ) std::cout << "reduce_hypot_args mid: using y_exp=" << exp << " x_lshift=" << diff << "\n";
        } else {
            exp = x_exp; 
            if ( debug ) std::cout << "reduce_hypot_args mid: using x_exp=" << exp << " y_rshift=0\n";
        }

        swapped = x < y;
        if ( swapped ) {
            T tmp = x;
            x = y;
            y = tmp;
        }

        if ( x >= _one_fxd ) {
            exp++;
            x = (x >> 1) | (x & 1);
            y = (y >> 1) | (y & 1);
        }
        if ( debug ) std::cout << "reduce_hypot_args end: xy_orig=[" << _to_flt(x_orig) << "," << _to_flt(y_orig) << "]" << 
                                       " xy_reduced_f=[" << _to_flt(x, false, true, false) << "," << _to_flt(y, false, true, false) << 
                                       "] xy_reduced=[" << std::hex << x << "," << y << std::dec << "] exp_class=" << to_str(exp_class) << 
                                       " x_exp=" << x_exp << " y_exp=" << y_exp << " exp=" << exp << " swapped=" << swapped << "\n"; 
    }
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_sincos_arg( bool times_pi, T& a, uint32_t& quad, EXP_CLASS& exp_class, bool& sign, bool& did_minus_pi_div_4 ) const
{
    //-----------------------------------------------------
    // Quick check for special values.
    //-----------------------------------------------------
    switch( fpclassify( a ) )
    {
        case FP_ZERO:
        case FP_SUBNORMAL:
            exp_class = EXP_CLASS::ZERO;
            sign = signbit( a );
            break;

        case FP_NAN:
            exp_class = EXP_CLASS::NOT_A_NUMBER;
            sign = signbit( a );
            break;

        case FP_INFINITE:
            exp_class = EXP_CLASS::INFINITE;
            sign = signbit( a );
            break;

        default:
        {
            //-----------------------------------------------------
            // Normal or Subnormal
            //
            // Compute a * 4/PI, then take integer part of that.
            // Subtract int_part*PI/4 from a.
            // If we ended up subtracting an odd PI/4, then set the flag.
            //-----------------------------------------------------
            exp_class = EXP_CLASS::NORMAL;
            const T a_orig = a;
            sign = signbit( a );
            if ( sign ) a = neg( a, false );
            T m;
            T aa = 0;
            T i;
            if ( !times_pi ) {
                m = mulc( a, _four_div_pi, false );
                (void)modf( m, &i );
                aa = mulc( i, _pi_div_4, false );
                a = sub( a, aa, false );
                if ( debug ) std::cout << "reduce_sincos_arg mid: a_orig=" << _to_flt(a_orig) <<
                                          " aa_f=" << _to_flt(aa) << " aa=0x" << std::hex << aa << std::dec << 
                                          " a_reduced_f=" << _to_flt(a) << " a_reduced=0x" << std::hex << a << std::dec << "\n";
            } else {
                m = scalbn( a, -2, false );  // divide by 4
                m = modf( m, &i );
                a = mulc( m, _four_div_pi, false );
            }
            EXP_CLASS a_exp_class;
            int32_t   a_exp;
            bool      a_sign;
            deconstruct( a, a_exp_class, a_exp, a_sign );
            while( a_exp < 0 ) 
            {
                a = (a >> 1) | (a & 1);
                a_exp++;
            }

            T ii = _to_flt( i ); 
            quad = (ii >> 1) & 3;
            did_minus_pi_div_4 = ii & 1;

            if ( debug ) std::cout << "reduce_sincos_arg: times_pi=" << times_pi << " a_orig=" << _to_flt(a_orig) << 
                                      " m=" << _to_flt(m) << " aa=" << _to_flt(aa) << " i=" << _to_flt(i) << " ii=" << ii <<
                                      " a_reduced=" << _to_flt(a, false, true) << " a_exp_class=" << to_str(a_exp_class) << 
                                      " a_exp=" << a_exp << " a_sign=" << a_sign << " quadrant=" << quad << " did_minus_pi_div_4=" << did_minus_pi_div_4 << "\n"; 

            break;
        }
    }
}

template class Cordic<int64_t, double>;

#endif
