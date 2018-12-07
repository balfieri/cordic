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
    //
    // do_reduce indicates whether Cordic routines should first do argument reduction, which is usually the case.
    // If you need a mix of reduce and no-reduce, then please allocate two different Cordic objects with different do_reduce
    // settings but the same other parameters.
    //
    // n == number of iterations used during CORDIC proper (0 == default == frac_w)
    //-----------------------------------------------------
    Cordic( uint32_t int_w, uint32_t frac_w, bool do_reduce=true, uint32_t n=0 );
    ~Cordic();

    //-----------------------------------------------------
    // Explicit Conversions
    //-----------------------------------------------------
    T       to_t( FLT x ) const;                // FLT to T encoded value
    FLT     to_flt( const T& x ) const;         // T encoded value to FLT
    std::string to_string( const T& x ) const;  // T to std::string
    std::string to_bstring( const T& x ) const; // T to std::string in binary format, like "1 001 101101011010"

    T       make_fixed( bool sign, T i, T f );  // encode a fixed-point    value using sign, integer  part i, and fractional part f
    T       make_float( bool sign, T e, T f );  // encode a floating-point value using sign, exponent part 3, and fractional part f

    //-----------------------------------------------------
    // Constants 
    //-----------------------------------------------------
    uint32_t int_w( void ) const;                       // int_w  from above
    uint32_t frac_w( void ) const;                      // frac_w from above
    uint32_t n( void ) const;                           // n      from above

    T max( void ) const;                                // maximum positive encoded value 
    T min( void ) const;                                // minimum positive encoded value
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
    T    one_over( const T& x ) const;                                    // 1/x
    T    sqrt( const T& x ) const;                                        // normh( x+1, x-1 ) / 2
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
    T    hypot( const T& x, const T& y ) const;                           // norm(x, y);  (same thing)
    T    normh( const T& x, const T& y ) const;                           // sqrt(x^2 - y^2)

    T    sinh( const T& x, const T * r=nullptr ) const;                   // r*sinh(x), also r*(e^x - e^-x)/2  (default r is 1)
    T    cosh( const T& x, const T * r=nullptr ) const;                   // r*cosh(x), also r*(e^x + e^-x)/2  (default r is 1)
    void sinh_cosh( const T& x, T& sih, T& coh, const T * r=nullptr ) const;// sih=r*sinh(x), coh=r*cosh(x)    (default r is 1)
    T    tanh( const T& x ) const;                                        // sinh(x) / cosh(x)            (2)
    T    asinh( const T& x ) const;                                       // log(x + sqrt(x^2 + 1))       (2)
    T    acosh( const T& x ) const;                                       // log(x + sqrt(x^2 - 1))       (2)
    T    atanh( const T& x ) const;                                       // atanh(x)
    T    atanh2( const T& y, const T& x ) const;                          // atanh(y/x)

    //-----------------------------------------------------
    // Interesting Identities (some are used in the implementation, most are not)
    //
    // sqrt(x)          = sqrt((x+0.25)^2 - (x-0.25)^2)         allows use of hyperbolic_vectoring
    // sqrt(x)          = sqrt((x+1)^2 - (x-1)^2) / 2           ditto, but easier to make sure |atanh((x-1)/(x+1))| is small enough
    // sqrt(x*y)        = sqrt((x+y)^2 - (x-y)^2) / 2           ditto, y doesn't matter
    // sqrt(x*y)        = sqrt(|x|) * sqrt(|y|)                 assuming x*y >= 0
    // sqrt(x^2 - y^2)  = sqrt((m+d)^2 - (m-d)^2) = 2*sqrt(m*d) where m=(x+y)/2, d=(x-y)/2 
    //                                                          factor m*d=p*s so that p is a power-of-2 and s is within -1 .. 1
    //                                                          2*sqrt(m*d) = 2*sqrt(p) * sqrt(s) = sqrt(p) * sqrt((s+1)^2 - (s-1)^2)
    //                                                          if log2(p) is even and >= 2, then sqrt(p) = 2^(log2(p)/2) = some integer
    //
    // pow(b,x)         = exp(log(b) * x)
    // exp(x)           = sinh(x) + cosh(x)
    // exp(-x)          = 1/exp(x)
    // exp(x+y)         = exp(x) * exp(y)
    // exp(ix)          = cos(x) + i*sin(x)                     Euler's Formula, i = sqrt(-1)
    // exp(i*PI) + 1    = 0                                     Euler's Identity
    // log(x)           = 2*atanh2(x-1, x+1)              
    // log(x)           = atanh2((x^2 - 1, x^2 + 1) 
    // log(x*y)         = log(x) + log(y)
    // log(x/4)         = atanh2(x-0.25, x+0.25) 
    // log(x/y)         = log(x) - log(y)
    // log(x/y)         = 2*atanh2(x-y, x+y)  
    //
    // sin(-x)          = -sin(x)
    // sin(x+y)         = sin(x)*cos(y) + cos(x)*sin(y)
    // sin(x+PI/4)      = sqrt(2)/2 * (sin(x) + cos(x))
    // sin(x)           = Im(e^(ix)) = (e^(ix) + e^(-ix)) / 2
    // sin(ix)          = i*sinh(x)
    // cos(x+y)         = cos(x)*sin(y) - sin(x)*cos(y)
    // cos(x+PI/4)      = sqrt(2)/2 * (cos(x) - sin(x))
    // cos(-x)          = cos(x)
    // cos(x)           = Re(e^(ix)) = (e^(ix) - e^(-ix)) / 2i
    // cos(ix)          = cosh(x)
    // tan(-x)          = -tan(x)
    // tan(x+y)         = (tan(x) + tan(y)) / (1 - tan(x)*tan(y))
    //
    // asin(-x)         = -asin(x)
    // asin(x)          = atan2(x, sqrt(1 - x^2))
    // acos(-x)         = acos(-x)
    // acos(x)          = atan2(sqrt(1 - x^2), x)
    // acos(x+y)        = PI/2 - asin(x+y)
    // atan(x)          = atan2(x, 1)                           this is true only when second argument is > 0 (see below)
    // atan(-x)         = -atan(x)
    // atan(1/x)        = PI/2 - atan(x)                        if x > 0
    // atan(x)          = asin(x / sqrt(1 + x^2))
    // atan(x)          = 2*atan(x / (1 + sqrt(1 + x^2)))
    // atan2(y,x)       = 2*atan(y / (sqrt(x^2 + y^2) + x))     if x >  0    (x < 0 can have inflated rounding errors, so...)
    // atan2(y,x)       = 2*atan((sqrt(x^2 + y^2) - x) / y)     if x <= 0 && y != 0
    // atan2(y,x)       = PI                                    if x <  0 && y == 0
    // atan2(y,x)       = undefined                             if x == 0 && y == 0
    // PI               = 4*atan(1)                             but low 2 bits will end up as 0 for a fixed-point number
    // PI               = acos(-1)                              but that just uses atan
    // PI               = pi()                                  function in this class to return a precomputed high-precision version
    //
    // sinh(-x)         = -sinh(x)
    // sinh(x)          = (e^x - e^-x)/2
    // sinh(x+y)        = sinh(x)*cosh(y) + cosh(x)*sinh(y)
    // cosh(-x)         = cosh(x)
    // cosh(x)          = (e^x + e^-x)/2
    // cosh(x+y)        = cosh(x)*cosh(y) + sinh(x)*sinh(y)
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
    // atanh(x)         = log((1+x)/(1-x))/2 = log(1+x)/2 - log(1-x)/2    note: x must be between -1 and 1
    // atanh(x)         = asinh(x / sqrt(1 - x^2)) 
    // atanh(x)         = acosh(1 / sqrt(1 - x^2))  (+/-)
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
    void sin_cos( const T& x, T& si, T& co, bool do_reduce, bool need_si=true, bool need_co=true, const T * r=nullptr ) const;
    void sinh_cosh( const T& x, T& sih, T& coh, bool do_reduce, bool need_sih=true, bool need_coh=true, const T * r=nullptr ) const;

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
    void reduce_sin_cos_arg( T& a, uint32_t& quadrant, bool& sign, bool& did_minus_pi_div_4 ) const;
    void reduce_sinh_cosh_arg( T& x, T& sinh_i, T& cosh_i, bool& sign ) const;                  

private:
    class Impl;
    Impl * impl;
};

#ifdef DEBUG_LEVEL
static constexpr uint32_t debug = DEBUG_LEVEL;
#else
static constexpr uint32_t debug = 0;
#endif

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
    T                           max;
    T                           min;
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
    T *                         reduce_sinh_cosh_sinh_i;                // for each possible integer value, sinh(i)
    T *                         reduce_sinh_cosh_cosh_i;                // for each possible integer value, cosh(i)
    bool *                      reduce_sinh_cosh_sinh_i_oflow;          // boolean indicating if this index creates too big of a number
    bool *                      reduce_sinh_cosh_cosh_i_oflow;          // boolean indicating if this index creates too big of a number
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
    cassert( (1+int_w+frac_w) <= (sizeof( T ) * 8), "1+int_w+frac_w does not fit in T container" );
    cassert( int_w  != 0, "int_w must be > 0 currently" );
    cassert( frac_w != 0, "frac_w must be > 0 currently" );

    impl = new Impl;

    impl->int_w   = int_w;
    impl->frac_w  = frac_w;
    impl->do_reduce = do_reduce;
    impl->n       = n;

    impl->maxint        = (T(1) << int_w) - 1;
    impl->max           = (T(1) << (int_w+frac_w)) - 1;
    impl->min           = T(1);
    impl->zero          = 0;
    impl->one           = T(1) << frac_w;
    impl->half          = T(1) << (frac_w-1);
    impl->quarter       = T(1) << (frac_w-2);
    impl->sqrt2         = to_t( std::sqrt( 2.0 ) );
    impl->sqrt2_div_2   = to_t( std::sqrt( 2.0 ) / FLT(2.0) );
    impl->pi            = to_t( std::acos( FLT(-1.0) ) );
    impl->pi_div_2      = to_t( std::acos( FLT(-1.0) ) / FLT(2.0) );
    impl->pi_div_4      = to_t( std::acos( FLT(-1.0) ) / FLT(4.0) );
    impl->two_div_pi    = to_t( FLT(2.0) / std::acos( FLT(-1.0) ) );
    impl->four_div_pi   = to_t( FLT(4.0) / std::acos( FLT(-1.0) ) );
    impl->e             = to_t( std::exp( FLT(  1 ) ) );
    impl->log2          = to_t( std::log( FLT(  2 ) ) );
    impl->log10         = to_t( std::log( FLT( 10 ) ) );

    impl->circular_atan    = new T[n+1];
    impl->hyperbolic_atanh = new T[n+1];
    impl->linear_pow2      = new T[n+1];

    // compute atan/atanh table in high-resolution floating point
    //
    for( uint32_t i = 0; i <= n; i++ )
    {
        impl->linear_pow2[i]      = T(1) << (frac_w-i);
        FLT pow2                  = to_flt( impl->linear_pow2[i] );
        FLT a                     = std::atan( pow2 );
        FLT ah                    = std::atanh( pow2 );
        impl->circular_atan[i]    = to_t( a );
        impl->hyperbolic_atanh[i] = to_t( ah );

        if ( debug ) printf( "i=%2d a=%30.27g ah=%30.27g y=%30.27g\n", i, double(a), double(ah), double(pow2) );
    }

    // calculate max |z0| angle allowed
    T xx, yy, zz;
    impl->circular_angle_max   = one();   // to avoid triggering assert
    impl->hyperbolic_angle_max = zero();  // to disable assert
    circular_vectoring(   one(),     one(), zero(), xx, yy, impl->circular_angle_max );
    hyperbolic_vectoring( to_t(0.5), one(), zero(), xx, yy, impl->hyperbolic_angle_max );
    if ( debug ) std::cout << "circular_angle_max="             << std::setw(30) << to_flt(impl->circular_angle_max) << "\n";
    if ( debug ) std::cout << "hyperbolic_angle_max="           << std::setw(30) << to_flt(impl->hyperbolic_angle_max) << "\n";
    
    // calculate gain by plugging in x=1,y=0,z=0 into CORDICs
    circular_rotation(    one(), zero(), zero(), impl->circular_rotation_gain,    yy, zz );
    circular_vectoring(   one(), zero(), zero(), impl->circular_vectoring_gain,   yy, zz );
    hyperbolic_rotation(  one(), zero(), zero(), impl->hyperbolic_rotation_gain,  yy, zz );
    hyperbolic_vectoring( one(), zero(), zero(), impl->hyperbolic_vectoring_gain, yy, zz );

    // calculate 1/gain which are the multiplication factors
    impl->circular_rotation_one_over_gain    = to_t( FLT(1) / to_flt(impl->circular_rotation_gain) );
    impl->circular_vectoring_one_over_gain   = to_t( FLT(1) / to_flt(impl->circular_vectoring_gain) );
    impl->hyperbolic_rotation_one_over_gain  = to_t( FLT(1) / to_flt(impl->hyperbolic_rotation_gain) );
    impl->hyperbolic_vectoring_one_over_gain = to_t( FLT(1) / to_flt(impl->hyperbolic_vectoring_gain) );
    if ( debug ) std::cout << "circular_rotation_gain="             << std::setw(30) << to_flt(impl->circular_rotation_gain) << "\n";
    if ( debug ) std::cout << "circular_vectoring_gain="            << std::setw(30) << to_flt(impl->circular_vectoring_gain) << "\n";
    if ( debug ) std::cout << "hyperbolic_rotation_gain="           << std::setw(30) << to_flt(impl->hyperbolic_rotation_gain) << "\n";
    if ( debug ) std::cout << "hyperbolic_vectoring_gain="          << std::setw(30) << to_flt(impl->hyperbolic_vectoring_gain) << "\n";
    if ( debug ) std::cout << "circular_rotation_one_over_gain="    << std::setw(30) << to_flt(impl->circular_rotation_one_over_gain) << "\n";
    if ( debug ) std::cout << "circular_vectoring_one_over_gain="   << std::setw(30) << to_flt(impl->circular_vectoring_one_over_gain) << "\n";
    if ( debug ) std::cout << "hyperbolic_rotation_one_over_gain="  << std::setw(30) << to_flt(impl->hyperbolic_rotation_one_over_gain) << "\n";
    if ( debug ) std::cout << "hyperbolic_vectoring_one_over_gain=" << std::setw(30) << to_flt(impl->hyperbolic_vectoring_one_over_gain) << "\n";

    // construct LUT used by reduce_sinh_cosh_arg();
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
    impl->reduce_sinh_cosh_sinh_i       = sinh_i;
    impl->reduce_sinh_cosh_cosh_i       = cosh_i;
    impl->reduce_sinh_cosh_sinh_i_oflow = sinh_i_oflow;
    impl->reduce_sinh_cosh_cosh_i_oflow = cosh_i_oflow;
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
        sinh_i[i]       = sinh_i_oflow[i] ? MAX : to_t( std::sinh( i_f ) );
        cosh_i[i]       = cosh_i_oflow[i] ? MAX : to_t( std::cosh( i_f ) );
        if ( debug ) std::cout << "reduce_sinh_cosh_arg LUT: i_f=" << i_f << " sinh_i=" << to_flt(sinh_i[i]) << " cosh_i=" << to_flt(cosh_i[i]) << 
                                  " sinh_i_oflow=" << sinh_i_oflow[i] << " cosh_i_oflow=" << cosh_i_oflow[i] << "\n";
    }

    // construct LUT used by reduce_exp_arg()
    // values for negative integers come first.
    FLT * factor = new FLT[2*N];
    impl->reduce_exp_factor = factor;
    T MIN_INT = -maxint() - 1;
    for( T i = MIN_INT; i <= maxint(); i++ )
    {
        T index = i - MIN_INT;
        factor[index] = std::exp(FLT(i));
        if ( debug ) std::cout << "reduce_exp_arg LUT: factor[" << i << "]=" << factor[index] << " index=" << index << "\n";
    }

    // construct LUT used by reduce_log_arg()
    addend = new T[frac_w+int_w];
    impl->reduce_log_addend = addend;
    for( int32_t i = -frac_w; i <= int32_t(int_w); i++ )
    {
        double addend_f = std::log( std::pow( 2.0, double( i ) ) );
        addend[frac_w+i] = to_t( addend_f );
        if ( debug ) std::cout << "addend[]=0x" << std::hex << addend[frac_w+i] << "\n" << std::dec;
        if ( debug ) std::cout << "reduce_log_arg LUT: addend[" << i << "]=" << to_flt(addend[frac_w+i]) << " addend_f=" << addend_f << "\n";
    }
}

template< typename T, typename FLT >
Cordic<T,FLT>::~Cordic( void )
{
    delete impl->circular_atan;
    delete impl->hyperbolic_atanh;
    delete impl->linear_pow2;
    delete impl->reduce_sinh_cosh_sinh_i;
    delete impl->reduce_sinh_cosh_cosh_i;
    delete impl->reduce_sinh_cosh_sinh_i_oflow;
    delete impl->reduce_sinh_cosh_cosh_i_oflow;
    delete impl->reduce_exp_factor;
    delete impl->reduce_log_addend;
    delete impl;
    impl = nullptr;
}

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
    return impl->maxint;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::max( void ) const
{
    return impl->max;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::min( void ) const
{
    return impl->min;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::zero( void ) const
{
    return impl->zero;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::one( void ) const
{
    return impl->one;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::half( void ) const
{
    return impl->half;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::quarter( void ) const
{
    return impl->quarter;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sqrt2( void ) const
{
    return impl->sqrt2;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sqrt2_div_2( void ) const
{
    return impl->sqrt2_div_2;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pi( void ) const
{
    return impl->pi;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pi_div_2( void ) const
{
    return impl->pi_div_2;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pi_div_4( void ) const
{
    return impl->pi_div_4;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::two_div_pi( void ) const
{
    return impl->two_div_pi;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::four_div_pi( void ) const
{
    return impl->four_div_pi;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::e( void ) const
{
    return impl->e;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::circular_rotation_gain( void ) const
{
    return impl->circular_rotation_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::circular_vectoring_gain( void ) const
{
    return impl->circular_vectoring_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hyperbolic_rotation_gain( void ) const
{
    return impl->hyperbolic_rotation_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hyperbolic_vectoring_gain( void ) const
{
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
    return impl->circular_vectoring_one_over_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hyperbolic_rotation_one_over_gain( void ) const
{
    return impl->hyperbolic_rotation_one_over_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hyperbolic_vectoring_one_over_gain( void ) const
{
    return impl->hyperbolic_vectoring_one_over_gain;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::circular_angle_max( void ) const
{
    return impl->circular_angle_max;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::hyperbolic_angle_max( void ) const
{
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
    cassert( T(x) < (T(1) << int_w()), "to_t: integer part of |x| does not fit in int_w bits" ); 
    T x_t = x * FLT( one() );
    if ( is_neg ) x_t = -x_t;
    return x_t;
}

template< typename T, typename FLT >
inline FLT Cordic<T,FLT>::to_flt( const T& _x ) const
{
    T x = _x;
    bool is_neg = x < 0;
    if ( is_neg ) x = -x;
    FLT x_f = FLT( x ) / FLT( one() );
    if ( is_neg ) x_f = -x_f;
    return x_f;
}

template< typename T, typename FLT >
inline std::string Cordic<T,FLT>::to_string( const T& x ) const
{
    return std::to_string( to_flt( x ) );  
}

template< typename T, typename FLT >
inline std::string Cordic<T,FLT>::to_bstring( const T& _x ) const
{
    T x = _x;
    uint32_t width = 1 + int_w() + frac_w();
    std::string bs = "";
    for( uint32_t i = 0; i < width; i++ )
    {
        if ( i == (int_w()+frac_w()) || i == frac_w() ) bs = " " + bs;
        const char * b = (x & 1) ? "1" : "0";
        bs = b + bs;
        x >>= 1;
    }
    return bs;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::make_fixed( bool sign, T i, T f )
{
    cassert( i >= 0 && i <= maxint()              , "make_fixed integer part must be in range 0 .. maxint()" );
    cassert( f >= 0 && f <= ((T(1) << frac_w())-1), "make_fixed fractional part must be in range 0 .. (1 << frac_w)-1" );

    return (T(sign) << (int_w() + frac_w())) |
           (T(i)    << frac_w())             |
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
    const T ONE = one();
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
        if ( debug ) printf( "circular_rotation: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int(z >= zero()) );
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
    const T ONE = one();
    const T THREE = 3*ONE;
    const T PI  = pi();
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
        if ( debug ) printf( "circular_vectoring: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int((x < zero()) != (y < zero())) );
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
    const T ONE = one();
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
        if ( debug ) printf( "circular_vectoring_xy: i=%d xy=[%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), int(y < zero()) );
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
    const T TWO = one() << 1;
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
        if ( debug ) printf( "hyperbolic_rotation: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int(z >= zero()) );
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
    const T TWO = one() << 1;
    const T PI  = pi();
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
                             to_flt(x), to_flt(y), to_flt(z), int((x < zero()) != (y < zero())) );
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
    const T TWO = one() << 1;
    const T PI  = pi();
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
        if ( debug ) printf( "hyperbolic_vectoring_xy: i=%d xy=[%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), int(y < zero()) );
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
    const T ONE = one();
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
        if ( debug ) printf( "linear_rotation: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int(z >= zero()) );
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
    const T ONE = one();
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
                             i, to_flt(x), to_flt(y), to_flt(z), int(y < zero()) );
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
    const T TWO = one() << 1;
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
        if ( debug ) printf( "linear_vectoring_xy: i=%d xy=[%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), int(y < zero()) );
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
inline T Cordic<T,FLT>::abs( const T& x ) const
{
    T    x_abs  = x;
    bool x_sign = x_abs < T(0);
    if ( x_sign ) x_abs = -x;
    T    sign_mask = x_abs >> (int_w() + frac_w());
    cassert( (sign_mask == T(0) || sign_mask == T(-1)), "abs caused overflow" ); 
    return x_abs;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::neg( const T& x ) const
{
    bool x_sign = x < 0;
    T    x_neg  = -x;
    T    sign_mask = x_neg >> (int_w() + frac_w());
    cassert( (x == 0 || sign_mask == (x_sign ? T(0) : T(-1))), "neg caused overflow" ); 
    return x_neg;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::floor( const T& x ) const
{
    T frac_mask = (T(1) << frac_w()) - 1;
    if ( (x & frac_mask) == 0 ) {
        return x;
    } else if ( x > 0 ) {
        return x & ~frac_mask;
    } else {
        return (x & ~frac_mask) - one();
    }
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::ceil( const T& x ) const
{
    T frac_mask = (T(1) << frac_w()) - 1;
    if ( (x & frac_mask) == 0 ) {
        return x;
    } else if ( x > 0 ) {
        return (x & ~frac_mask) + one();
    } else {
        return x & ~frac_mask;
    }
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::add( const T& x, const T& y ) const
{
    bool x_sign = x < T(0);
    bool y_sign = y < T(0);
    T    sum    = x + y;
    T    sign_mask = sum >> (int_w() + frac_w());
    cassert( sign_mask == T(0) || sign_mask == T(-1), "add caused overflow" );
    return sum;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sub( const T& x, const T& y ) const
{
    bool x_sign = x < T(0);
    bool y_sign = y < T(0);
    T    sum    = x - y;
    T    sign_mask = sum >> (int_w() + frac_w());
    cassert( sign_mask == 0 || sign_mask == T(-1), "sub caused overflow" );
    return sum;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::mad( const T& _x, const T& _y, const T addend, bool do_reduce ) const
{
    T x = _x;
    T y = _y;
    if ( debug ) std::cout << "mad begin: x_orig=" << to_flt(x) << " y_orig=" << to_flt(y) << " addend=" << to_flt(addend) << " do_reduce=" << do_reduce << "\n";
    cassert( do_reduce || addend >= 0, "mad addend must be non-negative" );
    int32_t x_lshift;
    int32_t y_lshift;
    bool    sign;
    if ( do_reduce ) reduce_mul_args( x, y, x_lshift, y_lshift, sign );

    T xx, yy, zz;
    linear_rotation( x, do_reduce ? zero() : addend, y, xx, yy, zz );
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
inline T Cordic<T,FLT>::mad( const T& _x, const T& _y, const T& addend ) const
{
    return mad( _x, _y, addend, impl->do_reduce );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::fma( const T& x, const T& y, const T& addend ) const
{
    return mad( x, y, addend );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::mul( const T& x, const T& y ) const
{
    return mad( x, y, zero() );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::mul( const T& x, const T& y, bool do_reduce ) const
{
    return mad( x, y, zero(), do_reduce );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::lshift( const T& x, int ls ) const
{
    cassert( x >= 0, "lshift x should be non-negative" );
    if ( ls > 0 ) {
        //-----------------------------------------------------
        // For now, crap out if we overflow.
        // At some point, we'll have options to saturate or set a flag in the container.
        //-----------------------------------------------------
        int32_t ls_max = int_w();
        uint32_t i = x >> frac_w();
        cassert( i <= maxint(), "lshift x integer part should be <= maxint()"  );
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
    linear_vectoring( x, y, do_reduce ? zero() : addend, xx, yy, zz );
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
    return dad( y, x, zero() );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::div( const T& y, const T& x, bool do_reduce ) const
{
    return dad( y, x, zero(), do_reduce );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::one_over( const T& x ) const
{
    return div( one(), x );
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
    T x = _x;
    if ( debug ) std::cout << "sqrt begin: x_orig=" << to_flt(_x) << " do_reduce=" << impl->do_reduce << "\n";
    int32_t ls;
    if ( impl->do_reduce ) reduce_sqrt_arg( x, ls );

    T xx, yy;
    hyperbolic_vectoring_xy( x+one(), x-one(), xx, yy );  // gain*sqrt((s+1)^2 - (s-1)^2)
    xx = mul( xx, hyperbolic_vectoring_one_over_gain(), false );   // sucks that we have to do this
    if ( impl->do_reduce ) xx = lshift( xx, ls );                  // log2(p)/2 - 1

    if ( debug ) std::cout << "sqrt end: x_orig=" << to_flt(_x) << " x_reduced=s=" << to_flt(x) << " do_reduce=" << impl->do_reduce << " xx=" << to_flt(xx) << "\n";
    return xx;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::one_over_sqrt( const T& x ) const
{ 
    if ( debug ) std::cout << "one_over_sqrt begin: x_orig=" << to_flt(x) << " do_reduce=" << impl->do_reduce << "\n";
    cassert( x != 0, "one_over_sqrt x must not be 0" );

    // There might be a better way, but exp(-log(x)/2) is probably not it
    return div( one(), sqrt( x ) );
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
inline T Cordic<T,FLT>::pow( const T& b, const T& x ) const
{ 
    if ( b == zero() ) return zero();
    cassert( b >= 0, "pow base b must be non-negative" );
    return exp( mul( x, log( b, true ) ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::powc( const FLT& b, const T& x ) const
{ 
    cassert( b > 0, "powc b must be positive" );
    const FLT log_b_f = std::log( b );
    cassert( log_b_f >= 0.0, "powc log(b) must be non-negative" );
    const T   log_b   = to_t( log_b_f );
    return exp( mul( x, log_b ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pow2( const T& x ) const
{ 
    return powc( 2.0, x );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::pow10( const T& x ) const
{ 
    return powc( 10.0, x );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::log( const T& _x, bool do_reduce ) const
{ 
    T x = _x;
    cassert( x > 0, "log: x must be positive" );
    T addend;
    if ( do_reduce ) reduce_log_arg( x, addend );
    T dv = div( x-one(), x+one(), false );
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
    cassert( b > 0, "logb b must be positive" );
    return div( log(x), log(b) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::logc( const T& x, const FLT& b ) const
{ 
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
    sin_cos( x, si, co, impl->do_reduce, true, false, r );
    return si;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::cos( const T& x, const T * r ) const
{ 
    T si;
    T co;
    sin_cos( x, si, co, impl->do_reduce, false, true, r );
    return co;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::sin_cos( const T& x, T& si, T& co, const T * r ) const             
{
    sin_cos( x, si, co, impl->do_reduce, true, true, r );
}

template< typename T, typename FLT >
void Cordic<T,FLT>::sin_cos( const T& _x, T& si, T& co, bool do_reduce, bool need_si, bool need_co, const T * _r ) const             
{ 
    T x = _x;
    uint32_t quadrant;
    bool x_sign;
    bool did_minus_pi_div_4;
    if ( do_reduce ) {
        //-----------------------------------------------------
        // reduce_sin_cos_arg() will get x in the range 0 .. PI/2 and tell us the quadrant.
        // It will then check if x is still > PI/4 and, if so, subtract PI/4 and
        // set did_minus_pi_div_4.  If did_minus_pi_div_4 is true, then we need
        // to do some adjustments below after the cordic routine completes.
        //-----------------------------------------------------
        reduce_sin_cos_arg( x, quadrant, x_sign, did_minus_pi_div_4 );
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
    circular_rotation( r, zero(), x, co, si, zz );
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
    T si, co;
    sin_cos( x, si, co );
    return div( si, co );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::asin( const T& x ) const
{ 
    cassert( x >= -one() && x <= one(), "asin x must be between -1 and 1" );
    return atan2( x, normh( one(), x ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::acos( const T& x ) const
{ 
    cassert( x >= -one() && x <= one(), "acos x must be between -1 and 1" );
    return atan2( normh( one(), x ), x, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::atan( const T& x ) const
{ 
    cassert( x >= -one() && x <= one(), "atan x must be between -1 and 1" );
    return atan2( x, one(), impl->do_reduce, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::atan2( const T& y, const T& x ) const
{ 
    return atan2( y, x, impl->do_reduce );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::atan2( const T& _y, const T& _x, bool do_reduce, bool x_is_one, T * r ) const
{ 
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
                                      " zz=PI" << " r=" << ((r != nullptr) ? to_flt(*r) : to_flt(zero())) << "\n";
            return pi();
        }

        const T norm_plus_x = norm( x, y, true ) + x;
        if ( debug ) std::cout << "atan2 cordic begin: y=y=" << to_flt(y) << " x=norm_plus_x=" << to_flt(norm_plus_x) << " swapped=" << swapped << "\n";
        circular_vectoring( norm_plus_x, y, zero(), xx, yy, zz );
        zz <<= 1;
        if ( swapped ) zz = pi_div_2() - zz;
        if ( y_sign ) zz = -zz;
        if ( debug ) std::cout << "atan2 cordic end: zz=" << to_flt(zz) << "\n";
    } else {
        circular_vectoring( x, y, zero(), xx, yy, zz );
    }
    if ( debug ) std::cout << "atan2 end: y=" << to_flt(_y) << " x=" << to_flt(_x) << " do_reduce=" << do_reduce << " x_is_one=" << x_is_one << 
                              " zz=" << to_flt(zz) << " r=" << ((r != nullptr) ? to_flt(*r) : to_flt(zero())) << "\n";
    return zz;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::polar_to_rect( const T& r, const T& a, T& x, T& y ) const
{
    if ( debug ) std::cout << "polar_to_rect begin: r=" << to_flt(r) << " a=" << to_flt(a) << " do_reduce=" << impl->do_reduce << "\n";
    sin_cos( a, y, x, true, true, true, &r );
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::rect_to_polar( const T& x, const T& y, T& r, T& a ) const
{
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
    if ( debug ) std::cout << "normh begin: x=" << to_flt(x) << " y=" << to_flt(y) << " do_reduce=" << impl->do_reduce << "\n";
    cassert( x >= y, "normh x must be >= y" );
    return sqrt( mul( x+y, x-y ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::sinh( const T& x, const T * r ) const
{ 
    T sih;
    T coh;
    sinh_cosh( x, sih, coh, impl->do_reduce, true, false, r );
    return sih;
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::cosh( const T& x, const T * r ) const
{ 
    T sih;
    T coh;
    sinh_cosh( x, sih, coh, impl->do_reduce, false, true, r );
    return coh;
}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::sinh_cosh( const T& x, T& sih, T& coh, const T * r ) const
{ 
    sinh_cosh( x, sih, coh, impl->do_reduce, true, true, r );
}

template< typename T, typename FLT >
void Cordic<T,FLT>::sinh_cosh( const T& _x, T& sih, T& coh, bool do_reduce, bool need_sih, bool need_coh, const T * _r ) const
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
    T x = _x;
    T sinh_i; 
    T cosh_i;
    bool sign;
    if ( do_reduce ) reduce_sinh_cosh_arg( x, sinh_i, cosh_i, sign );  

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
    hyperbolic_rotation( r, zero(), x, cosh_f, sinh_f, zz );
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
    T sih, coh;
    sinh_cosh( x, sih, coh );
    return div( sih, coh );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::asinh( const T& x ) const
{ 
    return log( x + norm( x, one() ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::acosh( const T& x ) const
{ 
    cassert( x >= one(), "acosh x must be >= 1" );
    return log( x + normh( x, one() ) );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::atanh( const T& x ) const
{ 
    return atanh2( x, one(), impl->do_reduce, true );
}

template< typename T, typename FLT >
inline T Cordic<T,FLT>::atanh2( const T& y, const T& x ) const             
{ 
    return atanh2( y, x, impl->do_reduce );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::atanh2( const T& _y, const T& _x, bool do_reduce, bool x_is_one ) const             
{ 
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
    hyperbolic_vectoring( x, y, zero(), xx, yy, zz );
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
    T other = for_sqrt ? half() : one();
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
    const T TWO = one() << 1;
    const T MININT = -maxint() - 1;
    T x_orig = x;
    if ( debug ) std::cout << "reduce_exp_arg: b=" << b << " x_orig=" << to_flt(x_orig) << "\n";
    const FLT * factors_f = impl->reduce_exp_factor;
    T   i         = x >> frac_w();  // can be + or -
    T   index     = -MININT + i;
    FLT factor_f  = std::log(b) * factors_f[index];   // could build per-b factors_f[] LUT with multiply already done
    factor        = to_t( factor_f );
    x            -= i << frac_w();
    if ( debug ) std::cout << "reduce_exp_arg: b=" << b << " x_orig=" << to_flt(x_orig) << 
                              " i=" << i << " index=" << index << " log(b)*exp(i)=" << to_flt(factor) << 
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
    addend = addends[frac_w()+x_lshift];
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
inline void Cordic<T,FLT>::reduce_sin_cos_arg( T& a, uint32_t& quad, bool& sign, bool& did_minus_pi_div_4 ) const
{
    //-----------------------------------------------------
    // Compute a * 4/PI, then take integer part of that.
    // Subtract int_part*PI/4 from a.
    // If we ended up subtracting an odd PI/4, then set the flag.
    //-----------------------------------------------------
    const T a_orig = a;
    sign = a < 0;
    if ( sign ) a = -a;
    const T m = mul( a, four_div_pi() );
    const T i = m >> frac_w();
    const T s = i * pi_div_4();
    a        -= s;
    quad      = (i >> 1) & 3;
    did_minus_pi_div_4 = i & 1;
    if ( debug ) std::cout << "reduce_sin_cos_arg: a_orig=" << to_flt(a_orig) << " m=" << to_flt(m) << " i=" << i << 
                              " subtract=" << to_flt(s) << " a_reduced=" << to_flt(a) << 
                              " quadrant=" << quad << " did_minus_pi_div_4=" << did_minus_pi_div_4 << "\n"; 

}

template< typename T, typename FLT >
inline void Cordic<T,FLT>::reduce_sinh_cosh_arg( T& x, T& sinh_i, T& cosh_i, bool& sign ) const
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

    const T MASK = (T(1) << (int_w()+2)) - T(1);  // include 0.25 bit of fraction
    T i = (x >> (frac_w()-2)) & MASK;
    x   = x & ((T(1) << (frac_w()-2))-T(1));
    const T *    sinh_i_vals   = impl->reduce_sinh_cosh_sinh_i;
    const T *    cosh_i_vals   = impl->reduce_sinh_cosh_cosh_i;
    const bool * sinh_i_oflows = impl->reduce_sinh_cosh_sinh_i_oflow;
    const bool * cosh_i_oflows = impl->reduce_sinh_cosh_cosh_i_oflow;
    cassert( !sinh_i_oflows[i], "reduce_sinh_cosh_arg x will cause an overflow for sinh" );
    cassert( !cosh_i_oflows[i], "reduce_sinh_cosh_arg x will cause an overflow for cosh" );
    sinh_i = sinh_i_vals[i];
    cosh_i = cosh_i_vals[i];
    if ( debug ) std::cout << "reduce_sinh_cosh_arg: x_orig=" << to_flt(x_orig) << " sinh_i[" << i << "]=" << to_flt(sinh_i) << 
                              " coshh_i[" << i << "]=" << to_flt(cosh_i) << " x_reduced=" << to_flt(x) << "\n";
}

template class Cordic<int64_t, double>;

#endif
