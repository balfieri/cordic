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

static constexpr uint32_t debug = DEBUG_LEVEL;

//-----------------------------------------------------
// INTERNAL IMPL STRUCTURE 
//-----------------------------------------------------
template< typename T, int INT_W, int FRAC_W, typename FLT >
struct Cordic<T,INT_W,FRAC_W,FLT>::Impl
{
    uint32_t                    nc;
    uint32_t                    nh;
    uint32_t                    nl;

    std::unique_ptr<T[]>        circular_atan;                  // circular atan values
    T                           circular_gain;                  // circular gain
    T                           circular_one_over_gain;         // circular 1/gain

    std::unique_ptr<T[]>        hyperbolic_atanh;               // hyperbolic atanh values
    T                           hyperbolic_gain;                // hyperbolic gain
    T                           hyperbolic_one_over_gain;       // hyperbolic 1/gain

    std::unique_ptr<T[]>        linear_pow2;                    // linear 2^(-i) values

    T                           log2;                           // log(2)
    T                           log10;                          // log(10)
    T                           log10_div_e;                    // log(10) / e

    std::unique_ptr<T[]>        reduce_angle_addend;            // for each possible integer value, an addend to help normalize
    std::unique_ptr<uint32_t[]> reduce_angle_quadrant;          // 00,01,02,03
    std::unique_ptr<FLT[]>      reduce_exp_factor;              // for each possible integer value, exp(i)
    std::unique_ptr<T[]>        reduce_log_addend;              // for each possible lshift value, log( 1 << lshift )
};

//-----------------------------------------------------
// Constructor
//-----------------------------------------------------
template< typename T, int INT_W, int FRAC_W, typename FLT >
Cordic<T,INT_W,FRAC_W,FLT>::Cordic( uint32_t nc, uint32_t nh, uint32_t nl )
{
    impl = std::make_unique<Impl>();

    impl->nc = nc;
    impl->nh = nh;
    impl->nl = nl;

    impl->circular_atan    = std::unique_ptr<T[]>( new T[nc+1] );
    impl->hyperbolic_atanh = std::unique_ptr<T[]>( new T[nh+1] );
    impl->linear_pow2      = std::unique_ptr<T[]>( new T[nl+1] );

    // compute atan/atanh table and gains in high-resolution floating point
    FLT pow2  = 1.0;
    FLT gain_inv  = 1.0;
    FLT gainh_inv = 1.0;
    uint32_t n_max = (nh > nc)    ? nh : nc;
             n_max = (nl > n_max) ? nl : n_max;
    uint32_t next_dup_i = 4;     // for hyperbolic 
    for( uint32_t i = 0; i <= n_max; i++ )
    {
        FLT a  = std::atan( pow2 );
        if ( i <= nc ) impl->circular_atan[i]    = T( a    * FLT( T(1) << T(FRAC_W) ) );
        if ( i <= nl ) impl->linear_pow2[i]      = T( pow2 * FLT( T(1) << T(FRAC_W) ) );
        FLT ah = std::atanh( pow2 );
        //FLT ah = 0.5 * std::log( (1.0+pow2) / (1-pow2) );
        if ( i <= nh ) impl->hyperbolic_atanh[i] = T( ah   * FLT( T(1) << T(FRAC_W) ) );
        if ( i <= nc ) gain_inv *= std::cos( a );

        if ( i != 0 && i <= nh ) {
            gainh_inv *= std::cosh( ah );
            if ( i == next_dup_i ) {
                // for hyperbolic, we must duplicate iterations 4, 13, 40, 121, ..., 3*i+1
                gainh_inv *= std::cosh( ah );
                next_dup_i = 3*i + 1;
            }
        }
        pow2 /= 2.0;

        if ( debug ) printf( "i=%2d a=%30.27g ah=%30.27g y=%30.27g gain_inv=%30.27g gainh_inv=%30.27g\n", i, double(a), double(ah), double(pow2), double(gain_inv), double(gainh_inv) );

    }

    // now convert those last two to fixed-point
    impl->circular_gain            = T( 1.0/gain_inv  * FLT( T(1) << T(FRAC_W) ) );
    impl->circular_one_over_gain   = T(     gain_inv  * FLT( T(1) << T(FRAC_W) ) );
    impl->hyperbolic_gain          = T( 1.0/gainh_inv * FLT( T(1) << T(FRAC_W) ) );
    impl->hyperbolic_one_over_gain = T(     gainh_inv * FLT( T(1) << T(FRAC_W) ) );

    // constants
    impl->log2        = T( std::log( FLT( 2  ) )       * FLT( T(1) << T(FRAC_W) ) );
    impl->log10       = T( std::log( FLT( 10 ) )       * FLT( T(1) << T(FRAC_W) ) );
    impl->log10_div_e = T( std::log( FLT( 10 ) / M_E ) * FLT( T(1) << T(FRAC_W) ) );

    // construct LUT used by reduce_angle()
    dassert( INT_W < 14 && "too many cases to worry about" );
    uint32_t N = 1 << (1+INT_W);
    T *        addend   = new T[N];
    uint32_t * quadrant = new uint32_t[N];
    impl->reduce_angle_addend   = std::unique_ptr<T[]>( addend );
    impl->reduce_angle_quadrant = std::unique_ptr<uint32_t[]>( quadrant );
    const FLT PI       = M_PI;
    const FLT PI_DIV_2 = PI / 2.0;
    for( T i = 0; i <= MAX_INT; i++ )
    {
        FLT cnt = FLT(i) / PI_DIV_2;
        T   cnt_i = cnt;
        if ( debug ) std::cout << "cnt_i=" << cnt_i << "\n";
        FLT add_f = FLT(cnt_i) * PI_DIV_2;
        if ( i > 0 ) add_f = -add_f;
        addend[i]   = to_fp( add_f );
        quadrant[i] = cnt_i % 4;
        if ( debug ) std::cout << "reduce_angle_arg LUT: cnt_i=" << cnt_i << " addend[" << i << "]=" << to_flt(addend[i]) << " quadrant=" << quadrant[i] << "\n";
    }

    // construct LUT used by reduce_exp_arg()
    FLT * factor = new FLT[N];
    impl->reduce_exp_factor = std::unique_ptr<FLT[]>( factor );
    for( T i = 0; i <= MAX_INT; i++ )
    {
        factor[i] = std::exp(FLT(i));
        if ( debug ) std::cout << "reduce_exp_arg LUT: factor[" << i << "]=" << factor[i] << "\n";
    }

    // construct LUT used by reduce_log_arg()
    addend = new T[INT_W];
    impl->reduce_log_addend = std::unique_ptr<T[]>( addend );
    for( T i = 0; i <= INT_W; i++ )
    {
        addend[i] = to_fp( std::log( double( 1 << i ) ) );
        if ( debug ) std::cout << "reduce_log_arg LUT: addend[" << i << "]=" << to_flt(addend[i]) << "\n";
    }
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
Cordic<T,INT_W,FRAC_W,FLT>::~Cordic( void )
{
    impl = nullptr;
}

//-----------------------------------------------------
// Queries
//-----------------------------------------------------
template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::to_fp( FLT x ) const
{
    return T( x * FLT(T(1) << T(FRAC_W)) );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
FLT Cordic<T,INT_W,FRAC_W,FLT>::to_flt( const T& x ) const
{
    return FLT( x ) / FLT(T(1) << T(FRAC_W));
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
uint32_t Cordic<T,INT_W,FRAC_W,FLT>::n_circular( void ) const
{
    return impl->nc;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
uint32_t Cordic<T,INT_W,FRAC_W,FLT>::n_hyperbolic( void ) const
{
    return impl->nh;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
uint32_t Cordic<T,INT_W,FRAC_W,FLT>::n_linear( void ) const
{
    return impl->nl;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::gain( void ) const
{
    return impl->circular_gain;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::gainh( void ) const
{
    return impl->hyperbolic_gain;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::one_over_gain( void ) const
{
    return impl->circular_one_over_gain;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::one_over_gainh( void ) const
{
    return impl->hyperbolic_one_over_gain;
}

//-----------------------------------------------------
// The CORDIC Functions
//-----------------------------------------------------
template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::circular_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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
        if ( debug ) printf( "circular_rotation: i=%d xyz=[%f,%f,%f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int(z >= ZERO) );
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

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::circular_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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
        if ( debug ) printf( "circular_vectoring: i=%d xyz=[%f,%f,%f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int((x < ZERO) != (y < ZERO)) );
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

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::hyperbolic_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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
        if ( debug ) printf( "hyperbolic_rotation: i=%d xyz=[%f,%f,%f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int(z >= ZERO) );
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
            next_dup_i = 3*i + 1;
            i--;
        }
    }
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::hyperbolic_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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
            next_dup_i = 3*i + 1;
            i--;
        }
    }
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::linear_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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
        if ( debug ) printf( "linear_rotation: i=%d xyz=[%f,%f,%f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int(z >= ZERO) );
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

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::linear_vectoring( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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
        if ( debug ) printf( "linear_vectoring: i=%d xyz=[%f,%f,%f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int((x < ZERO) != (y < ZERO)) );
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

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::mul( const T& _x, const T& _y, const T addend, bool do_reduce ) const
{
    T x = _x;
    T y = _y;
    dassert( x >= 0 && "mul x must be non-negative" );
    dassert( y >= 0 && "mul y must be non-negative" );
    int32_t x_lshift;
    int32_t y_lshift;
    if ( do_reduce ) reduce_mul_args( x, y, x_lshift, y_lshift );

    T xx, yy, zz;
    linear_rotation( x, addend, y, xx, yy, zz );
    if ( do_reduce ) yy <<= x_lshift + y_lshift;
    return yy;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::div( const T& _y, const T& _x, const T addend, bool do_reduce ) const
{
    T x = _x;
    T y = _y;
    dassert( x >= 0 && "div x must be non-negative" );
    dassert( y > 0  && "div y must be positive" );
    int32_t x_lshift;
    int32_t y_lshift;
    if ( do_reduce ) reduce_div_args( x, y, x_lshift, y_lshift );

    T xx, yy, zz;
    linear_vectoring( x, y, addend, xx, yy, zz );
    if ( do_reduce ) zz <<= x_lshift - y_lshift;
    return zz;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::sqrt( const T& _x, bool do_reduce ) const
{ 
    T x = _x;
    dassert( x >= 0 && "sqrt x must be non-negative" );
    int32_t x_lshift;
    if ( do_reduce ) reduce_sqrt_arg( x, x_lshift );

    // sqrt( (x+0.25)^2 - (x-0.25)^2 ) = normh( x+0.25, x-0.25 )
    T n = normh( x + QUARTER, x - QUARTER );
    if ( do_reduce ) n <<= x_lshift;
    return n;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::one_over_sqrt( const T& _x, bool do_reduce ) const
{ 
    T x = _x;
    dassert( x > 0 && "one_over_sqrt x must be positive" );
    int32_t x_lshift;
    if ( do_reduce ) reduce_sqrt_arg( x, x_lshift );

    // do basic 1/sqrt for now
    // try pow( x, -0.5 ) later
    T n = div( ONE, sqrt( x, do_reduce ), false );
    if ( do_reduce ) n >>= x_lshift;
    return n;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::exp( const T& _x, bool do_reduce ) const
{ 
    T x = _x;
    dassert( x >= 0 && "exp x must be non-negative" );
    T factor;
    if ( do_reduce ) reduce_exp_arg( M_E, x, factor );

    T xx, yy, zz;
    hyperbolic_rotation( one_over_gainh(), one_over_gainh(), x, xx, yy, zz );
    if ( do_reduce ) xx = mul( xx, factor, true );
    return xx;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::pow( const T& b, const T& x, bool do_reduce ) const
{ 
    dassert( b >= 0 && "pow b must be non-negative" );
    dassert( x >= 0 && "pow x must be non-negative" );
    return exp( mul( x, log( b, true ), do_reduce ), do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::powc( const FLT& b, const T& x, bool do_reduce ) const
{ 
    dassert( b >= 0.0 && "powc b must be non-negative" );
    dassert( x >= 0   && "powc x must be non-negative" );
    const FLT log_b_f = std::log( b );
    dassert( log_b_f >= 0.0 && "powc log(b) must be non-negative" );
    const T   log_b   = to_fp( log_b_f );
    return exp( mul( x, log_b, do_reduce ), do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::pow2( const T& x, bool do_reduce ) const
{ 
    dassert( x >= 0 && "pow2 x must be non-negative" );
    return powc( 2.0, x, do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::pow10( const T& x, bool do_reduce ) const
{ 
    dassert( x >= 0 && "pow10 x must be non-negative" );
    return powc( 10.0, x, do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::log( const T& _x, bool do_reduce ) const
{ 
    T x = _x;
    dassert( x >= 0 && "log x must be non-negative" );
    T addend;
    if ( do_reduce ) reduce_log_arg( x, addend );
    T lg = atanh2( x-ONE, x+ONE, false ) << 1;
    if ( do_reduce ) lg += addend;
    return lg;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::logb( const T& x, const T& b, bool do_reduce ) const
{ 
    dassert( x >= 0 && "logb x must be non-negative" );
    dassert( b > 0  && "logb b must be positive" );
    return div( log(x, do_reduce), log(b, do_reduce), do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::logc( const T& x, const FLT& b, bool do_reduce ) const
{ 
    dassert( x >= 0  && "logc x must be non-negative" );
    dassert( b > 0.0 && "logc b must be positive" );
    const FLT one_over_log_b_f = FLT(1) / std::log( b );
    const T   one_over_log_b   = to_fp( one_over_log_b_f );
    return mul( log(x, do_reduce), one_over_log_b, do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::log2( const T& x, bool do_reduce ) const
{ 
    dassert( x >= 0 && "log2 x must be non-negative" );
    return logc( x, 2.0, do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::log10( const T& x, bool do_reduce ) const
{ 
    dassert( x >= 0 && "log10 x must be non-negative" );
    return logc( x, 10.0, do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::sin( const T& _x, bool do_reduce ) const
{ 
    T x = _x;
    dassert( x >= 0 && "sin x must be non-negative" );
    uint32_t quadrant;
    if ( do_reduce ) reduce_angle_arg( x, quadrant );

    T xx, yy, zz;
    circular_rotation( one_over_gain(), ZERO, x, xx, yy, zz );
    if ( do_reduce ) {
        if ( quadrant&1 )    yy = xx;      // use cos
        if ( quadrant >= 2 ) yy = -yy;
    }
    return yy;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::cos( const T& _x, bool do_reduce ) const
{ 
    T x = _x;
    dassert( x >= 0 && "cos x must be non-negative" );
    uint32_t quadrant;
    if ( do_reduce ) reduce_angle_arg( x, quadrant );

    T xx, yy, zz;
    circular_rotation( one_over_gain(), ZERO, x, xx, yy, zz );
    if ( do_reduce ) {
        if ( quadrant&1 )    xx = yy;      // use sin
        if ( quadrant == 1 || quadrant == 2 ) xx = -xx;
    }
    return xx;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::sin_cos( const T& _x, T& si, T& co, bool do_reduce ) const             
{ 
    T x = _x;
    dassert( x >= 0 && "sin_cos x must be non-negative" );
    uint32_t quadrant;
    if ( do_reduce ) reduce_angle_arg( x, quadrant );

    T zz;
    circular_rotation( one_over_gain(), ZERO, x, co, si, zz );
    if ( do_reduce ) {
        if ( quadrant&1 ) {
            T tmp = co;
            co = si;
            si = tmp;
        }
        if ( quadrant == 1 || quadrant == 2 ) co = -co;
        if ( quadrant >= 2 ) si = -si;
    }
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::tan( const T& x, bool do_reduce ) const
{ 
    dassert( x >= 0 && "tan x must be non-negative" );
    T si, co;
    sin_cos( x, si, co, do_reduce );
    return div( si, co, do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::asin( const T& x, bool do_reduce ) const
{ 
    dassert( x >= 0 && "asin x must be non-negative" );
    return atan2( x, normh( ONE, x, do_reduce ), do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::acos( const T& x, bool do_reduce ) const
{ 
    dassert( x >= 0 && "acos x must be non-negative" );
    return atan2( normh( ONE, x, do_reduce ), x, do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::atan( const T& x, bool do_reduce ) const
{ 
    dassert( x >= 0 && "atan x must be non-negative" );
    dassert( !do_reduce && "TODO" );
    T xx, yy, zz;
    circular_vectoring( ONE, x, ZERO, xx, yy, zz );
    return zz;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::atan2( const T& y, const T& x, bool do_reduce ) const
{ 
    dassert( y >= 0 && "atan2 y must be non-negative" );
    dassert( x > 0  && "atan2 x must be positive" );
    dassert( !do_reduce && "TODO" );
    T xx, yy, zz;
    circular_vectoring( x, y, ZERO, xx, yy, zz );
    return zz;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::polar_to_rect( const T& r, const T& a, T& x, T& y, bool do_reduce ) const
{
    dassert( r >= 0 && "polar_to_rect r must be non-negative" );
    dassert( a >= 0 && "polar_to_rect a must be non-negative" );
    dassert( !do_reduce && "TODO" );
    T xx, yy, zz;
    circular_rotation( r, ZERO, a, xx, yy, zz );
    x = mul( xx, one_over_gain(), do_reduce );
    y = mul( yy, one_over_gain(), do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::rect_to_polar( const T& _x, const T& _y, T& r, T& a, bool do_reduce ) const
{
    T x = _x;
    T y = _y;
    dassert( x >= 0 && "rect_to_polar x must be non-negative" );
    dassert( y >= 0 && "rect_to_polar y must be non-negative" );
    int32_t lshift;
    if ( do_reduce ) reduce_norm_args( x, y, lshift );

    T rr;
    T yy;
    circular_vectoring( x, y, ZERO, rr, yy, a );
    if ( do_reduce ) rr <<= lshift;
    r = mul( rr, one_over_gain(), do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::norm( const T& _x, const T& _y, bool do_reduce ) const
{
    T x = _x;
    T y = _y;
    dassert( x >= 0 && "norm x must be non-negative" );
    dassert( y >= 0 && "norm y must be non-negative" );
    int32_t lshift;
    if ( do_reduce ) reduce_norm_args( x, y, lshift );

    T xx, yy, zz;
    circular_vectoring( x, y, ZERO, xx, yy, zz );
    if ( do_reduce ) xx <<= lshift;
    return mul( xx, one_over_gain(), do_reduce );  
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::normh( const T& _x, const T& _y, bool do_reduce ) const
{
    T x = _x;
    T y = _y;
    dassert( x >= 0 && "normh x must be non-negative" );
    dassert( y >= 0 && "normh y must be non-negative" );
    dassert( x >= y && "normh x must >= y" );
    int32_t lshift;
    if ( do_reduce ) reduce_norm_args( x, y, lshift );

    T xx, yy, zz;
    hyperbolic_vectoring( x, y, ZERO, xx, yy, zz );
    if ( do_reduce ) xx <<= lshift;
    return mul( xx, one_over_gainh(), do_reduce ); 
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::sinh( const T& _x, bool do_reduce ) const
{ 
    T x = _x;
    dassert( x >= 0 && "sinh x must be non-negative" );
    uint32_t quadrant;
    if ( do_reduce ) reduce_angle_arg( x, quadrant );

    T xx, yy, zz;
    hyperbolic_rotation( one_over_gainh(), ZERO, x, xx, yy, zz );
    if ( do_reduce ) {
        if ( quadrant&1 )    yy = xx;      // use cos
        if ( quadrant >= 2 ) yy = -yy;
    }
    return yy;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::cosh( const T& _x, bool do_reduce ) const
{ 
    T x = _x;
    dassert( x >= 0 && "cosh x must be non-negative" );
    uint32_t quadrant;
    if ( do_reduce ) reduce_angle_arg( x, quadrant );

    T xx, yy, zz;
    hyperbolic_rotation( one_over_gainh(), ZERO, x, xx, yy, zz );
    if ( do_reduce ) {
        if ( quadrant&1 )    xx = yy;      // use sin
        if ( quadrant == 1 || quadrant == 2 ) xx = -xx;
    }
    return xx;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::sinh_cosh( const T& _x, T& sih, T& coh, bool do_reduce ) const
{ 
    T x = _x;
    dassert( x >= 0 && "sinh_cosh x must be non-negative" );
    uint32_t quadrant;
    if ( do_reduce ) reduce_angle_arg( x, quadrant );

    T zz;
    hyperbolic_rotation( one_over_gainh(), ZERO, x, coh, sih, zz );
    if ( do_reduce ) {
        if ( quadrant&1 ) {
            T tmp = coh;
            coh = sih;
            sih = tmp;
        }
        if ( quadrant == 1 || quadrant == 2 ) coh = -coh;
        if ( quadrant >= 2 ) sih = -sih;
    }
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::tanh( const T& x, bool do_reduce ) const
{ 
    T sih, coh;
    dassert( x >= 0 && "tanh x must be non-negative" );
    sinh_cosh( x, sih, coh, do_reduce );
    return div( sih, coh, do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::asinh( const T& x, bool do_reduce ) const
{ 
    dassert( x >= 0 && "asinh x must be non-negative" );
    return log( x + norm( ONE, x, do_reduce ), do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::acosh( const T& x, bool do_reduce ) const
{ 
    dassert( x >= 0 && "acosh x must be non-negative" );
    return log( x + normh( x, ONE, do_reduce ), do_reduce );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::atanh( const T& x, bool do_reduce ) const
{ 
    dassert( x >= 0 && "atanh x must be non-negative" );
    dassert( !do_reduce && "TODO" );
    T xx, yy, zz;
    hyperbolic_vectoring( ONE, x, ZERO, xx, yy, zz );
    return zz;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
T Cordic<T,INT_W,FRAC_W,FLT>::atanh2( const T& y, const T& x, bool do_reduce ) const             
{ 
    dassert( y >= 0 && "atanh y must be non-negative" );
    dassert( x >  0 && "atanh x must be positive" );
    dassert( !do_reduce && "TODO" );
    T xx, yy, zz;
    hyperbolic_vectoring( x, y, ZERO, xx, yy, zz );
    return zz;
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::reduce_arg( T& x, int32_t& x_lshift, bool shift_x, bool normalize ) const
{
    T x_orig = x;
    dassert( x >= 0 );
    T other = T(1) << FRAC_W;
    x_lshift = 0;
    while( x > other ) 
    {
        x_lshift++;
        if ( shift_x ) {
            x >>= 1;
        } else {
            other <<= 1;
        }
    }
    while( normalize && x < ONE )
    {
        x_lshift--;
        if ( shift_x ) {
            x <<= 1;
        } else {
            other >>= 1;
        }
    }
    if ( debug && shift_x ) std::cout << "reduce_arg: x_orig=" << to_flt(x_orig) << " x_reduced=" << to_flt(x) << " x_lshift=" << x_lshift << "\n"; 
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::reduce_mul_args( T& x, T& y, int32_t& x_lshift, int32_t& y_lshift ) const
{
    if ( debug ) std::cout << "reduce_mul_args: x_orig=" << x << " y_orig=" << y << "\n";
    reduce_arg( x, x_lshift );
    reduce_arg( y, y_lshift );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::reduce_div_args( T& x, T& y, int32_t& x_lshift, int32_t& y_lshift ) const
{
    if ( debug ) std::cout << "reduce_div_args: x_orig=" << x << " y_orig=" << y << "\n";
    reduce_arg( x, x_lshift );
    reduce_arg( y, y_lshift, true, true );
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::reduce_sqrt_arg( T& x, int32_t& x_lshift ) const
{
    //-----------------------------------------------------
    // Reduce without right-shifting x yet.
    // Then round the x_lshift up to even number.
    // And *then* shift.
    //-----------------------------------------------------
    T x_orig = x;
    reduce_arg( x, x_lshift, false );   
    if ( x_lshift & 1 ) x_lshift++;
    x >>= x_lshift;
    if ( debug ) std::cout << "reduce_sqrt_arg: x_orig=" << to_flt(x_orig) << " x_reduced=" << to_flt(x) << " x_lshift=" << x_lshift << "\n"; 
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::reduce_exp_arg( FLT b, T& x, T& factor ) const 
{
    //-----------------------------------------------------
    // Assume: x = i + f  (integer plus fraction)
    // exp(x) = exp(i) * exp(f)
    // pow(b,x) = log(b) * exp(x) = [log(b)*exp(i)] * exp(f)
    //
    // exp(i) comes for a pre-built LUT kept in FLT
    // so we can multiply it by log(b) before converting to type T.
    //-----------------------------------------------------
    T x_orig = x;
    const FLT * factors_f = impl->reduce_exp_factor.get();
    T   index    = (x >> FRAC_W) & MAX_INT;
    FLT factor_f = std::log(b) * factors_f[index];   // could build per-b factors_f[] LUT with multiply already done
    factor       = to_fp( factor_f );
    x           &= (T(1) << FRAC_W)-T(1); // fraction only
    if ( debug ) std::cout << "reduce_exp_arg: x_orig=" << to_flt(x_orig) << " x_reduced=" << to_flt(x) << " factor=" << factor << "\n"; 
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::reduce_log_arg( T& x, T& addend ) const 
{
    //-----------------------------------------------------
    // log(ab) = log(a) + log(b)
    // 
    // So right-shift x using reduce_arg().
    // Then addend = log(1 << shift).
    //-----------------------------------------------------
    T x_orig = x;
    int32_t x_lshift;
    reduce_arg( x, x_lshift, false );
    const T * addends = impl->reduce_log_addend.get();
    addend = addends[x_lshift];
    x -= addend;
    if ( debug ) std::cout << "reduce_log_arg: x_orig=" << to_flt(x_orig) << " x_reduced=" << to_flt(x) << " addend=" << to_flt(addend) << "\n"; 
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::reduce_norm_args( T& x, T& y, int32_t& lshift ) const
{
    //-----------------------------------------------------
    // Must shift both x and y by max( x_lshift, y_lshift ).
    //-----------------------------------------------------
    T x_orig = x;
    T y_orig = y;
    int32_t x_lshift;
    int32_t y_lshift;
    reduce_arg( x, x_lshift, false );   
    reduce_arg( y, y_lshift, false );   
    lshift = (x_lshift > y_lshift) ? x_lshift : y_lshift;
    x >>= lshift;
    y >>= lshift;
    if ( debug ) std::cout << "reduce_norm_arg: xy_orig=[" << to_flt(x_orig) << "," << to_flt(y_orig) << "]" << 
                                              " xy_reduced=[" << to_flt(x) << "," << to_flt(y) << "] lshift=" << lshift << "\n"; 
}

template< typename T, int INT_W, int FRAC_W, typename FLT >
void Cordic<T,INT_W,FRAC_W,FLT>::reduce_angle_arg( T& a, uint32_t& quad ) const
{
    //-----------------------------------------------------
    // Use LUT to find addend.
    //-----------------------------------------------------
    const T a_orig = a;
    dassert( a >= 0 );
    const T *  addend   = impl->reduce_angle_addend.get();
    uint32_t * quadrant = impl->reduce_angle_quadrant.get();

    T index = (a >> FRAC_W) & MAX_INT;
    quad = quadrant[index];
    a += addend[index];
    if ( debug ) std::cout << "reduce_angle_arg: a_orig=" << to_flt(a_orig) << " a_reduced=" << to_flt(a) << " quadrant=" << quad << "\n"; 
}

template class Cordic<int64_t, 7, 56>;

