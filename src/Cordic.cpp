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
#include "Cordic.h"

static constexpr bool debug = false;

using highres = double;

//-----------------------------------------------------
// INTERNAL IMPL STRUCTURE 
//-----------------------------------------------------
template< typename T, int FW >
struct Cordic<T,FW>::Impl
{
    uint32_t                    nc;
    uint32_t                    nh;
    uint32_t                    nl;

    std::unique_ptr<T[]>        circular_atan;                  // circular atan values
    T                           circular_one_over_gain;         // circular 1/gain

    std::unique_ptr<T[]>        hyperbolic_atanh;               // hyperbolic atanh values
    T                           hyperbolic_one_over_gain;       // hyperbolic 1/gain

    std::unique_ptr<T[]>        linear_pow2;                    // linear 2^(-i) values

    T                           log2;                           // log(2)
    T                           log10;                          // log(10)
    T                           log10_div_e;                    // log(10) / e

    T to_fp( highres x )
    {
        return T( x * highres(T(1) << T(FW)) );
    }

    highres to_flt( T x )
    {
        return highres( x ) / highres(T(1) << T(FW));
    }

};

//-----------------------------------------------------
// Constructor
//-----------------------------------------------------
template< typename T, int FW >
Cordic<T,FW>::Cordic( uint32_t nc, uint32_t nh, uint32_t nl )
{
    impl = std::make_unique<Impl>();

    impl->nc = nc;
    impl->nh = nh;
    impl->nl = nl;

    impl->circular_atan    = std::unique_ptr<T[]>( new T[nc+1] );
    impl->hyperbolic_atanh = std::unique_ptr<T[]>( new T[nh+1] );
    impl->linear_pow2      = std::unique_ptr<T[]>( new T[nl+1] );

    // compute atan/atanh table and gains in high-resolution floating point
    highres pow2  = 1.0;
    highres gain_inv  = 1.0;
    highres gainh_inv = 1.0;
    uint32_t n_max = (nh > nc)    ? nh : nc;
             n_max = (nl > n_max) ? nl : n_max;
    uint32_t next_dup_i = 4;     // for hyperbolic 
    for( uint32_t i = 0; i <= n_max; i++ )
    {
        highres a  = std::atan( pow2 );
        highres ah = std::atanh( pow2 );
        if ( i <= nc ) impl->circular_atan[i]    = T( a    * highres( T(1) << T(FW) ) );
        if ( i <= nh ) impl->hyperbolic_atanh[i] = T( ah   * highres( T(1) << T(FW) ) );
        if ( i <= nl ) impl->linear_pow2[i]      = T( pow2 * highres( T(1) << T(FW) ) );
        if ( i <= nc ) gain_inv *= std::cos( a );

        if ( i != 0 && i <= nh ) {
            gainh_inv *= std::cosh( ah );
            if ( i == next_dup_i ) {
                // for hyperbolic, we must duplicate iterations 4, 13, 40, 121, ..., 3*i+1
                gainh_inv *= std::cosh( ah );
                next_dup_i = 3*i + 1;
            }
        }

        if ( debug ) printf( "i=%2d a=%30.27g ah=%30.27g y=%30.27g gain_inv=%30.27g gainh_inv=%30.27g\n", i, double(a), double(ah), double(pow2), double(gain_inv), double(gainh_inv) );

        pow2 /= 2.0;
    }

    // now convert those last two to fixed-point
    impl->circular_one_over_gain   = T( gain_inv  * highres( T(1) << T(FW) ) );
    impl->hyperbolic_one_over_gain = T( gainh_inv * highres( T(1) << T(FW) ) );

    // constants
    impl->log2        = T( std::log( highres( 2  ) )       * highres( T(1) << T(FW) ) );
    impl->log10       = T( std::log( highres( 10 ) )       * highres( T(1) << T(FW) ) );
    impl->log10_div_e = T( std::log( highres( 10 ) / M_E ) * highres( T(1) << T(FW) ) );
}

template< typename T, int FW >
Cordic<T,FW>::~Cordic( void )
{
    impl = nullptr;
}

//-----------------------------------------------------
// Queries
//-----------------------------------------------------
template< typename T, int FW >
uint32_t Cordic<T,FW>::n_circular( void ) const
{
    return impl->nc;
}

template< typename T, int FW >
uint32_t Cordic<T,FW>::n_hyperbolic( void ) const
{
    return impl->nh;
}

template< typename T, int FW >
uint32_t Cordic<T,FW>::n_linear( void ) const
{
    return impl->nl;
}

template< typename T, int FW >
T Cordic<T,FW>::one_over_gain( void ) const
{
    return impl->circular_one_over_gain;
}

template< typename T, int FW >
T Cordic<T,FW>::one_over_gainh( void ) const
{
    return impl->hyperbolic_one_over_gain;
}

//-----------------------------------------------------
// The CORDIC Functions
//-----------------------------------------------------
template< typename T, int FW >
void Cordic<T,FW>::circular_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
{
    //-----------------------------------------------------
    // d = (z >= 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctan(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->nc;
    for( uint32_t i = 0; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "circular_rotation: i=%d xyz=[%f,%f,%f] test=%d\n", i, impl->to_flt(x), impl->to_flt(y), impl->to_flt(z), int(z >= ZERO) );
        if ( z >= ZERO ) {
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

template< typename T, int FW >
void Cordic<T,FW>::circular_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
{
    //-----------------------------------------------------
    // d = (xy < 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctan(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->nc;
    for( uint32_t i = 0; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "circular_vectoring: i=%d xyz=[%f,%f,%f] test=%d\n", i, impl->to_flt(x), impl->to_flt(y), impl->to_flt(z), int((x < ZERO) != (y < ZERO)) );
        if ( (x < ZERO) != (y < ZERO) ) {
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

template< typename T, int FW >
void Cordic<T,FW>::hyperbolic_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
{
    //-----------------------------------------------------
    // d = (z >= 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctanh(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->nh;
    uint32_t next_dup_i = 4;     
    for( uint32_t i = 1; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "hyperbolic_rotation: i=%d xyz=[%f,%f,%f] test=%d\n", i, impl->to_flt(x), impl->to_flt(y), impl->to_flt(z), int(z >= ZERO) );
        if ( z >= ZERO ) {
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
            i--;
            next_dup_i = 3*i + 1;
        }
    }
}

template< typename T, int FW >
void Cordic<T,FW>::hyperbolic_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
{
    //-----------------------------------------------------
    // d = (xy < 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctanh(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->nh;
    uint32_t next_dup_i = 4;     
    for( uint32_t i = 1; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( (x < ZERO) != (y < ZERO) ) {
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
            i--;
            next_dup_i = 3*i + 1;
        }
    }
}

template< typename T, int FW >
void Cordic<T,FW>::linear_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
{
    //-----------------------------------------------------
    // d = (z >= 0) ? 1 : -1
    // xi = x 
    // yi = y + d*(x >> i)
    // zi = z - d*2^(-i)
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->nl;
    for( uint32_t i = 0; i <= n; i++ )
    {
        if ( debug ) printf( "linear_rotation: i=%d xyz=[%f,%f,%f] test=%d\n", i, impl->to_flt(x), impl->to_flt(y), impl->to_flt(z), int(z >= ZERO) );
        T yi;
        T zi;
        if ( z >= ZERO ) {
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

template< typename T, int FW >
void Cordic<T,FW>::linear_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
{
    //-----------------------------------------------------
    // d = (y <= 0) ? 1 : -1
    // xi = x
    // yi = y + d*(x >> i)
    // zi = z - d*2^(-i)
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->nl;
    for( uint32_t i = 0; i <= n; i++ )
    {
        if ( debug ) printf( "linear_vectoring: i=%d xyz=[%f,%f,%f] test=%d\n", i, impl->to_flt(x), impl->to_flt(y), impl->to_flt(z), int((x < ZERO) != (y < ZERO)) );
        T yi;
        T zi;
        if ( (x < ZERO) != (y < ZERO) ) {
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

template< typename T, int FW >
T Cordic<T,FW>::mul( const T& x, const T& y, const T addend ) const
{
    T xx, yy, zz;
    linear_rotation( x, addend, y, xx, yy, zz );
    return yy;
}

template< typename T, int FW >
T Cordic<T,FW>::div( const T& y, const T& x, const T addend ) const
{
    T xx, yy, zz;
    linear_vectoring( x, y, addend, xx, yy, zz );
    return zz;
}

template< typename T, int FW >
T Cordic<T,FW>::sqrt( const T& x ) const
{ 
    // sqrt( (x+0.25)^2 - (x-0.25)^2 )
    return normh( x + QUARTER, x - QUARTER );
}

template< typename T, int FW >
T Cordic<T,FW>::exp( const T& x ) const
{ 
    T xx, yy, zz;
    hyperbolic_rotation( one_over_gainh(), one_over_gainh(), x, xx, yy, zz );
    return xx;
}

template< typename T, int FW >
T Cordic<T,FW>::pow( const T& b, const T& x ) const
{ 
    return exp( mul( x, log( b ) ) );
}

template< typename T, int FW >
T Cordic<T,FW>::pow2( const T& x ) const
{ 
    return exp( mul( x, impl->log2 ) );
}

template< typename T, int FW >
T Cordic<T,FW>::pow10( const T& x ) const
{ 
    // log10 is too large for mul; we use log(10/e) then add 1 after mul()
    return exp( mul( x, impl->log10_div_e ) + ONE );
}

template< typename T, int FW >
T Cordic<T,FW>::log( const T& x ) const
{ 
    return atanh2( x-ONE, x+ONE ) << 1;
}

template< typename T, int FW >
T Cordic<T,FW>::logb( const T& x, const T& b ) const
{ 
    return div( log(x), log(b) );
}

template< typename T, int FW >
T Cordic<T,FW>::log2( const T& x ) const
{ 
    return div( log(x), impl->log2 );
}

template< typename T, int FW >
T Cordic<T,FW>::log10( const T& x ) const
{ 
    return div( log(x), impl->log10 );
}

template< typename T, int FW >
T Cordic<T,FW>::sin( const T& x ) const
{ 
    T xx, yy, zz;
    circular_rotation( one_over_gain(), ZERO, x, xx, yy, zz );
    return yy;
}

template< typename T, int FW >
T Cordic<T,FW>::cos( const T& x ) const
{ 
    T xx, yy, zz;
    circular_rotation( one_over_gain(), ZERO, x, xx, yy, zz );
    return xx;
}

template< typename T, int FW >
void Cordic<T,FW>::sin_cos( const T& x, T& si, T& co ) const             
{ 
    T zz;
    circular_rotation( one_over_gain(), ZERO, x, co, si, zz );
}

template< typename T, int FW >
T Cordic<T,FW>::tan( const T& x ) const
{ 
    T si, co;
    sin_cos( x, si, co );
    return div( si, co );
}

template< typename T, int FW >
T Cordic<T,FW>::asin( const T& x ) const
{ 
    return atan2( x, normh( ONE, x ) );
}

template< typename T, int FW >
T Cordic<T,FW>::acos( const T& x ) const
{ 
    return atan2( normh( ONE, x ), x );
}

template< typename T, int FW >
T Cordic<T,FW>::atan( const T& x ) const
{ 
    T xx, yy, zz;
    circular_vectoring( ONE, x, ZERO, xx, yy, zz );
    return zz;
}

template< typename T, int FW >
T Cordic<T,FW>::atan2( const T& y, const T& x ) const
{ 
    T xx, yy, zz;
    circular_vectoring( x, y, ZERO, xx, yy, zz );
    return zz;
}

template< typename T, int FW >
void Cordic<T,FW>::polar_to_rect( const T& r, const T& a, T& x, T& y ) const
{
    T xx, yy, zz;
    circular_rotation( r, ZERO, a, xx, yy, zz );
    x = mul( xx, one_over_gain() );
    y = mul( yy, one_over_gain() );
}

template< typename T, int FW >
void Cordic<T,FW>::rect_to_polar( const T& x, const T& y, T& r, T& a ) const
{
    T rr;
    T yy;
    circular_vectoring( x, y, ZERO, rr, yy, a );
    r = mul( rr, one_over_gain() );
}

template< typename T, int FW >
T Cordic<T,FW>::norm( const T& x, const T& y ) const
{
    T xx, yy, zz;
    circular_vectoring( x, y, ZERO, xx, yy, zz );
    return mul( xx, one_over_gain() );
}

template< typename T, int FW >
T Cordic<T,FW>::normh( const T& x, const T& y ) const
{
    T xx, yy, zz;
    hyperbolic_vectoring( x, y, ZERO, xx, yy, zz );
    return mul( xx, one_over_gainh() );
}

template< typename T, int FW >
T Cordic<T,FW>::sinh( const T& x ) const
{ 
    T xx, yy, zz;
    hyperbolic_rotation( one_over_gainh(), ZERO, x, xx, yy, zz );
    return yy;
}

template< typename T, int FW >
T Cordic<T,FW>::cosh( const T& x ) const
{ 
    T xx, yy, zz;
    hyperbolic_rotation( one_over_gainh(), ZERO, x, xx, yy, zz );
    return xx;
}

template< typename T, int FW >
void Cordic<T,FW>::sinh_cosh( const T& x, T& sih, T& coh ) const
{ 
    T zz;
    hyperbolic_rotation( one_over_gainh(), ZERO, x, coh, sih, zz );
}

template< typename T, int FW >
T Cordic<T,FW>::tanh( const T& x ) const
{ 
    T sih, coh;
    sinh_cosh( x, sih, coh );
    return div( sih, coh );
}

template< typename T, int FW >
T Cordic<T,FW>::asinh( const T& x ) const
{ 
    return log( x + norm( ONE, x ) );
}

template< typename T, int FW >
T Cordic<T,FW>::acosh( const T& x ) const
{ 
    return log( x + normh( x, ONE ) );
}

template< typename T, int FW >
T Cordic<T,FW>::atanh( const T& x ) const
{ 
    T xx, yy, zz;
    hyperbolic_vectoring( ONE, x, ZERO, xx, yy, zz );
    return zz;
}

template< typename T, int FW >
T Cordic<T,FW>::atanh2( const T& y, const T& x ) const             
{ 
    T xx, yy, zz;
    hyperbolic_vectoring( x, y, ZERO, xx, yy, zz );
    return zz;
}

template class Cordic<int64_t, 56>;
