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

//-----------------------------------------------------
// INTERNAL IMPL STRUCTURE 
//-----------------------------------------------------
template< typename T, int INT_W, int FRAC_W >
struct Cordic<T,INT_W,FRAC_W>::Impl
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
        return T( x * highres(T(1) << T(FRAC_W)) );
    }

    highres to_flt( T x )
    {
        return highres( x ) / highres(T(1) << T(FRAC_W));
    }

};

//-----------------------------------------------------
// Constructor
//-----------------------------------------------------
template< typename T, int INT_W, int FRAC_W >
Cordic<T,INT_W,FRAC_W>::Cordic( uint32_t nc, uint32_t nh, uint32_t nl )
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
        if ( i <= nc ) impl->circular_atan[i]    = T( a    * highres( T(1) << T(FRAC_W) ) );
        if ( i <= nh ) impl->hyperbolic_atanh[i] = T( ah   * highres( T(1) << T(FRAC_W) ) );
        if ( i <= nl ) impl->linear_pow2[i]      = T( pow2 * highres( T(1) << T(FRAC_W) ) );
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
    impl->circular_one_over_gain   = T( gain_inv  * highres( T(1) << T(FRAC_W) ) );
    impl->hyperbolic_one_over_gain = T( gainh_inv * highres( T(1) << T(FRAC_W) ) );

    // constants
    impl->log2        = T( std::log( highres( 2  ) )       * highres( T(1) << T(FRAC_W) ) );
    impl->log10       = T( std::log( highres( 10 ) )       * highres( T(1) << T(FRAC_W) ) );
    impl->log10_div_e = T( std::log( highres( 10 ) / M_E ) * highres( T(1) << T(FRAC_W) ) );
}

template< typename T, int INT_W, int FRAC_W >
Cordic<T,INT_W,FRAC_W>::~Cordic( void )
{
    impl = nullptr;
}

//-----------------------------------------------------
// Queries
//-----------------------------------------------------
template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::to_fp( highres x ) const
{
    return T( x * highres(T(1) << T(FRAC_W)) );
}

template< typename T, int INT_W, int FRAC_W >
highres Cordic<T,INT_W,FRAC_W>::to_flt( const T& x ) const
{
    return highres( x ) / highres(T(1) << T(FRAC_W));
}

template< typename T, int INT_W, int FRAC_W >
uint32_t Cordic<T,INT_W,FRAC_W>::n_circular( void ) const
{
    return impl->nc;
}

template< typename T, int INT_W, int FRAC_W >
uint32_t Cordic<T,INT_W,FRAC_W>::n_hyperbolic( void ) const
{
    return impl->nh;
}

template< typename T, int INT_W, int FRAC_W >
uint32_t Cordic<T,INT_W,FRAC_W>::n_linear( void ) const
{
    return impl->nl;
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::one_over_gain( void ) const
{
    return impl->circular_one_over_gain;
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::one_over_gainh( void ) const
{
    return impl->hyperbolic_one_over_gain;
}

//-----------------------------------------------------
// The CORDIC Functions
//-----------------------------------------------------
template< typename T, int INT_W, int FRAC_W >
void Cordic<T,INT_W,FRAC_W>::circular_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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

template< typename T, int INT_W, int FRAC_W >
void Cordic<T,INT_W,FRAC_W>::circular_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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

template< typename T, int INT_W, int FRAC_W >
void Cordic<T,INT_W,FRAC_W>::hyperbolic_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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

template< typename T, int INT_W, int FRAC_W >
void Cordic<T,INT_W,FRAC_W>::hyperbolic_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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

template< typename T, int INT_W, int FRAC_W >
void Cordic<T,INT_W,FRAC_W>::linear_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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

template< typename T, int INT_W, int FRAC_W >
void Cordic<T,INT_W,FRAC_W>::linear_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::mul( const T& x, const T& y, const T addend, bool do_reduce ) const
{
    T xx, yy, zz;
    linear_rotation( x, addend, y, xx, yy, zz );
    return yy;
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::div( const T& y, const T& x, const T addend, bool do_reduce ) const
{
    T xx, yy, zz;
    linear_vectoring( x, y, addend, xx, yy, zz );
    return zz;
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::sqrt( const T& x, bool do_reduce ) const
{ 
    // sqrt( (x+0.25)^2 - (x-0.25)^2 )
    return normh( x + QUARTER, x - QUARTER );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::exp( const T& x, bool do_reduce ) const
{ 
    T xx, yy, zz;
    hyperbolic_rotation( one_over_gainh(), one_over_gainh(), x, xx, yy, zz );
    return xx;
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::pow( const T& b, const T& x, bool do_reduce ) const
{ 
    return exp( mul( x, log( b ), do_reduce ), false );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::pow2( const T& x, bool do_reduce ) const
{ 
    return exp( mul( x, impl->log2, do_reduce ), false );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::pow10( const T& x, bool do_reduce ) const
{ 
    // log10 is too large for mul; we use log(10/e) then add 1 after mul()
    return exp( mul( x, impl->log10_div_e, do_reduce ) + ONE, do_reduce );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::log( const T& x, bool do_reduce ) const
{ 
    return atanh2( x-ONE, x+ONE, do_reduce ) << 1;
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::logb( const T& x, const T& b, bool do_reduce ) const
{ 
    return div( log(x, do_reduce), log(b, do_reduce), false );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::log2( const T& x, bool do_reduce ) const
{ 
    return div( log(x, do_reduce), impl->log2, do_reduce );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::log10( const T& x, bool do_reduce ) const
{ 
    return div( log(x, do_reduce), impl->log10, do_reduce );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::sin( const T& x, bool do_reduce ) const
{ 
    T xx, yy, zz;
    circular_rotation( one_over_gain(), ZERO, x, xx, yy, zz );
    return yy;
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::cos( const T& x, bool do_reduce ) const
{ 
    T xx, yy, zz;
    circular_rotation( one_over_gain(), ZERO, x, xx, yy, zz );
    return xx;
}

template< typename T, int INT_W, int FRAC_W >
void Cordic<T,INT_W,FRAC_W>::sin_cos( const T& x, T& si, T& co, bool do_reduce ) const             
{ 
    T zz;
    circular_rotation( one_over_gain(), ZERO, x, co, si, zz );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::tan( const T& x, bool do_reduce ) const
{ 
    T si, co;
    sin_cos( x, si, co, do_reduce );
    return div( si, co );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::asin( const T& x, bool do_reduce ) const
{ 
    return atan2( x, normh( ONE, x, do_reduce ), do_reduce );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::acos( const T& x, bool do_reduce ) const
{ 
    return atan2( normh( ONE, x, do_reduce ), x, do_reduce );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::atan( const T& x, bool do_reduce ) const
{ 
    T xx, yy, zz;
    circular_vectoring( ONE, x, ZERO, xx, yy, zz );
    return zz;
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::atan2( const T& y, const T& x, bool do_reduce ) const
{ 
    T xx, yy, zz;
    circular_vectoring( x, y, ZERO, xx, yy, zz );
    return zz;
}

template< typename T, int INT_W, int FRAC_W >
void Cordic<T,INT_W,FRAC_W>::polar_to_rect( const T& r, const T& a, T& x, T& y, bool do_reduce ) const
{
    T xx, yy, zz;
    circular_rotation( r, ZERO, a, xx, yy, zz );
    x = mul( xx, one_over_gain(), do_reduce );
    y = mul( yy, one_over_gain(), do_reduce );
}

template< typename T, int INT_W, int FRAC_W >
void Cordic<T,INT_W,FRAC_W>::rect_to_polar( const T& x, const T& y, T& r, T& a, bool do_reduce ) const
{
    T rr;
    T yy;
    circular_vectoring( x, y, ZERO, rr, yy, a );
    r = mul( rr, one_over_gain(), do_reduce );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::norm( const T& x, const T& y, bool do_reduce ) const
{
    T xx, yy, zz;
    circular_vectoring( x, y, ZERO, xx, yy, zz );
    return mul( xx, one_over_gain(), do_reduce );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::normh( const T& x, const T& y, bool do_reduce ) const
{
    T xx, yy, zz;
    hyperbolic_vectoring( x, y, ZERO, xx, yy, zz );
    return mul( xx, one_over_gainh(), do_reduce );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::sinh( const T& x, bool do_reduce ) const
{ 
    T xx, yy, zz;
    hyperbolic_rotation( one_over_gainh(), ZERO, x, xx, yy, zz );
    return yy;
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::cosh( const T& x, bool do_reduce ) const
{ 
    T xx, yy, zz;
    hyperbolic_rotation( one_over_gainh(), ZERO, x, xx, yy, zz );
    return xx;
}

template< typename T, int INT_W, int FRAC_W >
void Cordic<T,INT_W,FRAC_W>::sinh_cosh( const T& x, T& sih, T& coh, bool do_reduce ) const
{ 
    T zz;
    hyperbolic_rotation( one_over_gainh(), ZERO, x, coh, sih, zz );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::tanh( const T& x, bool do_reduce ) const
{ 
    T sih, coh;
    sinh_cosh( x, sih, coh, do_reduce );
    return div( sih, coh );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::asinh( const T& x, bool do_reduce ) const
{ 
    return log( x + norm( ONE, x, do_reduce ), do_reduce );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::acosh( const T& x, bool do_reduce ) const
{ 
    return log( x + normh( x, ONE, do_reduce ), do_reduce );
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::atanh( const T& x, bool do_reduce ) const
{ 
    T xx, yy, zz;
    hyperbolic_vectoring( ONE, x, ZERO, xx, yy, zz );
    return zz;
}

template< typename T, int INT_W, int FRAC_W >
T Cordic<T,INT_W,FRAC_W>::atanh2( const T& y, const T& x, bool do_reduce ) const             
{ 
    T xx, yy, zz;
    hyperbolic_vectoring( x, y, ZERO, xx, yy, zz );
    return zz;
}

template class Cordic<int64_t, 7, 56>;
