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

#ifdef DEBUG_LEVEL
static constexpr uint32_t debug = DEBUG_LEVEL;
#else
static constexpr uint32_t debug = 0;
#endif

//-----------------------------------------------------
// INTERNAL IMPL STRUCTURE 
//-----------------------------------------------------
template< typename T, typename FLT >
struct Cordic<T,FLT>::Impl
{
    uint32_t                    int_w;
    uint32_t                    frac_w;
    bool                        do_reduce;
    uint32_t                    n;

    T                           maxint;
    T                           zero;
    T                           one;
    T                           quarter;
    T                           pi;
    T                           e;
    T                           log2;                           // log(2)
    T                           log10;                          // log(10)

    std::unique_ptr<T[]>        circular_atan;                  // circular atan values
    T                           circular_gain;                  // circular gain
    T                           circular_one_over_gain;         // circular 1/gain

    std::unique_ptr<T[]>        hyperbolic_atanh;               // hyperbolic atanh values
    T                           hyperbolic_gain;                // hyperbolic gain
    T                           hyperbolic_one_over_gain;       // hyperbolic 1/gain

    std::unique_ptr<T[]>        linear_pow2;                    // linear 2^(-i) values

    std::unique_ptr<T[]>        reduce_sin_cos_addend;          // for each possible integer value, an addend to help normalize
    std::unique_ptr<uint32_t[]> reduce_sin_cos_quadrant;        // 00,01,02,03
    std::unique_ptr<T[]>        reduce_sinh_cosh_sinh_i;        // for each possible integer value, sinh(i)
    std::unique_ptr<T[]>        reduce_sinh_cosh_cosh_i;        // for each possible integer value, cosh(i)
    std::unique_ptr<FLT[]>      reduce_exp_factor;              // for each possible integer value, exp(i)
    std::unique_ptr<T[]>        reduce_log_addend;              // for each possible lshift value, log( 1 << lshift )

    inline void                 do_lshift( T& x, int32_t lshift ) const
    {
        cassert( x >= 0        && "do_lshift x should be non-negative" );
        if ( lshift > 0 ) {
            //-----------------------------------------------------
            // For now, crap out if we overflow.
            // At some point, we'll have options to saturate or set a flag in the container.
            //-----------------------------------------------------
            int32_t lshift_max = int_w;
            uint32_t i = x >> frac_w;
            cassert( i <= maxint && "do_lshift x integer part should be <= maxint()"  );
            while( i != 0 ) 
            {
                lshift_max--;
                i >>= 1;
            }
            if ( lshift > lshift_max ) {
                std::cout << "do_lshift x << " << lshift << " will overflow x\n";
                exit( 1 );
            }
            x <<= lshift;
        } else if ( lshift < 0 ) {
            x >>= -lshift;
        }
    }
};

//-----------------------------------------------------
// Constructor
//-----------------------------------------------------
template< typename T, typename FLT >
Cordic<T,FLT>::Cordic( uint32_t int_w, uint32_t frac_w, bool do_reduce, uint32_t n )
{
    if ( n == 0 ) n = frac_w;
    cassert( (1+int_w+frac_w) <= (sizeof( T ) * 8) && "1+int_w+frac_w does not fit in T container" );
    cassert( int_w  != 0 && "int_w must be > 0 currently" );
    cassert( frac_w != 0 && "frac_w must be > 0 currently" );

    impl = std::make_unique<Impl>();

    impl->int_w   = int_w;
    impl->frac_w  = frac_w;
    impl->do_reduce = do_reduce;
    impl->n       = n;

    impl->maxint  = (T(1) << int_w) - 1;
    impl->zero    = 0;
    impl->one     = T(1) << frac_w;
    impl->quarter = T(1) << (frac_w-2);
    impl->pi      = to_flt( std::acos( FLT(-1.0) ) );
    impl->e       = to_flt( std::exp( FLT(1) ) );
    impl->log2    = to_t( std::log( FLT(  2 ) ) );
    impl->log10   = to_t( std::log( FLT( 10 ) ) );

    impl->circular_atan    = std::unique_ptr<T[]>( new T[n+1] );
    impl->hyperbolic_atanh = std::unique_ptr<T[]>( new T[n+1] );
    impl->linear_pow2      = std::unique_ptr<T[]>( new T[n+1] );

    // compute atan/atanh table and gains in high-resolution floating point
    // Let gain = cos(a0)*cos(a1)*... and multiply x and y 
    // Use cos(a) = 1/sqrt(1 + tan^2(a)) 
    // Let gainh = cosh(a0)*cosh(a1)*... and multiply x and y 
    // Use cosh(a) = 1/sqrt(1 - tanh^2(a)) 
    //
    FLT pow2      = 1.0;
    uint32_t next_dup_i = 4;     // for hyperbolic 
    for( uint32_t i = 0; i <= n; i++ )
    {
        impl->linear_pow2[i]      = to_t( pow2 );
        FLT a                     = std::atan( pow2 );
        FLT ah                    = std::atanh( pow2 );
        impl->circular_atan[i]    = to_t( a );
        impl->hyperbolic_atanh[i] = to_t( ah );

        if ( i == next_dup_i ) {
            // for hyperbolic, we must duplicate iterations 4, 13, 40, 121, ..., 3*i+1
            next_dup_i = 3*i + 1;
        }
        pow2 /= 2.0;

        if ( debug ) printf( "i=%2d a=%30.27g ah=%30.27g y=%30.27g\n", i, double(a), double(ah), double(pow2) );

    }

    // calculate gain and gainh by plugging x=1,y=0,z=0 into CORDICs
    T gain, gainh, yy, zz;
    circular_rotation( one(), zero(), zero(), gain, yy, zz );
    hyperbolic_rotation( one(), zero(), zero(), gainh, yy, zz );

    // calculate 1/gain and 1/gainh which are the multiplication factors 
    impl->circular_gain            = gain;
    impl->hyperbolic_gain          = gainh;
    impl->circular_one_over_gain   = to_t( FLT(1) / to_flt(gain)  );
    impl->hyperbolic_one_over_gain = to_t( FLT(1) / to_flt(gainh) );
    if ( debug ) std::cout << "circular_gain="            << std::setw(30) << to_flt(impl->circular_gain) << "\n";
    if ( debug ) std::cout << "circular_gainh="           << std::setw(30) << to_flt(impl->hyperbolic_gain) << "\n";
    if ( debug ) std::cout << "circular_one_over_gain="   << std::setw(30) << to_flt(impl->circular_one_over_gain) << "\n";
    if ( debug ) std::cout << "hyperbolic_one_over_gain=" << std::setw(30) << to_flt(impl->hyperbolic_one_over_gain) << "\n";

    // construct LUTs used by reduce_sin_cos_arg() and reduce_sinh_cosh_arg();
    // use integer part plus 0.5 bit of fraction
    cassert( int_w < 14 && "too many cases to worry about" );
    uint32_t N = 1 << (1+int_w);
    T *        addend   = new T[N];
    uint32_t * quadrant = new uint32_t[N];
    T *        sinh_i   = new T[N];
    T *        cosh_i   = new T[N];
    impl->reduce_sin_cos_addend   = std::unique_ptr<T[]>( addend );
    impl->reduce_sin_cos_quadrant = std::unique_ptr<uint32_t[]>( quadrant );
    impl->reduce_sinh_cosh_sinh_i = std::unique_ptr<T[]>( sinh_i );
    impl->reduce_sinh_cosh_cosh_i = std::unique_ptr<T[]>( cosh_i );
    const FLT PI       = M_PI;
    const FLT PI_DIV_2 = PI / 2.0;
    const T   MASK     = (T(1) << (int_w+T(1)))-T(1);  // include 0.5 bit of fraction
    for( T i = 0; i <= MASK; i++ )
    {
        FLT i_f = FLT(i) / 2.0;
        FLT cnt = i_f / PI_DIV_2; 
        T   cnt_i = cnt;
        if ( debug ) std::cout << "cnt_i=" << cnt_i << "\n";
        FLT add_f = FLT(cnt_i) * PI_DIV_2;
        if ( i > 0 ) add_f = -add_f;
        addend[i]   = to_t( add_f );
        quadrant[i] = cnt_i % 4;
        if ( debug ) std::cout << "reduce_sin_cos_arg LUT: i_f=" << i_f << " cnt_i=" << cnt_i << 
                                  " addend[" << i << "]=" << to_flt(addend[i]) << " quadrant=" << quadrant[i] << "\n";

        sinh_i[i] = std::sinh( i_f );
        cosh_i[i] = std::cosh( i_f );
        if ( debug ) std::cout << "reduce_sinh_cosh_arg LUT: i_f=" << i_f << " sinh_i=" << to_flt(sinh_i[i]) << " cosh_i=" << to_flt(cosh_i[i]) << "\n";
    }

    // construct LUT used by reduce_exp_arg()
    FLT * factor = new FLT[N];
    impl->reduce_exp_factor = std::unique_ptr<FLT[]>( factor );
    for( T i = 0; i <= maxint(); i++ )
    {
        factor[i] = std::exp(FLT(i));
        if ( debug ) std::cout << "reduce_exp_arg LUT: factor[" << i << "]=" << factor[i] << "\n";
    }

    // construct LUT used by reduce_log_arg()
    addend = new T[frac_w+int_w];
    impl->reduce_log_addend = std::unique_ptr<T[]>( addend );
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
    impl = nullptr;
}

//-----------------------------------------------------
// Constants
//-----------------------------------------------------
template< typename T, typename FLT >
uint32_t Cordic<T,FLT>::int_w( void ) const
{
    return impl->int_w;
}

template< typename T, typename FLT >
uint32_t Cordic<T,FLT>::frac_w( void ) const
{
    return impl->frac_w;
}

template< typename T, typename FLT >
uint32_t Cordic<T,FLT>::n( void ) const
{
    return impl->n;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::maxint( void ) const
{
    return impl->maxint;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::zero( void ) const
{
    return impl->zero;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::one( void ) const
{
    return impl->one;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::quarter( void ) const
{
    return impl->quarter;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::pi( void ) const
{
    return impl->pi;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::e( void ) const
{
    return impl->e;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::gain( void ) const
{
    return impl->circular_gain;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::gainh( void ) const
{
    return impl->hyperbolic_gain;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::one_over_gain( void ) const
{
    return impl->circular_one_over_gain;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::one_over_gainh( void ) const
{
    return impl->hyperbolic_one_over_gain;
}

//-----------------------------------------------------
// Conversion
//-----------------------------------------------------
template< typename T, typename FLT >
T Cordic<T,FLT>::to_t( FLT _x ) const
{
    FLT x = _x;
    bool is_neg = x < 0.0;
    if ( is_neg ) x = -x;
    cassert( T(x) < (T(1) << int_w()) && "to_t: integer part of |x| does not fit in int_w bits" ); 
    T x_t = x * FLT( one() );
    if ( is_neg ) x_t = -x_t;
    return x_t;
}

template< typename T, typename FLT >
FLT Cordic<T,FLT>::to_flt( const T& _x ) const
{
    T x = _x;
    bool is_neg = x < 0;
    if ( is_neg ) x = -x;
    FLT x_f = FLT( x ) / FLT( one() );
    if ( is_neg ) x_f = -x_f;
    return x_f;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::make_fp( bool sign, T i, T f )
{
    cassert( i >= 0 && i <= maxint()               && "make_fp integer part must be in range 0 .. maxint()" );
    cassert( f >= 0 && f <= ((T(1) << frac_w())-1) && "make_fp fractional part must be in range 0 .. (1 << frac_w)-1" );

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
    // d = (z >= 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctan(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->n;
    for( uint32_t i = 0; i < n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "circular_rotation: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int(z >= zero()) );
        if ( z >= zero() ) {
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
    // d = (xy < 0) ? 1 : -1
    // xi = x - d*(y >> i)
    // yi = y + d*(x >> i)
    // zi = z - d*arctan(2^(-i))
    //-----------------------------------------------------
    x = x0;
    y = y0;
    z = z0;
    uint32_t n = impl->n;
    for( uint32_t i = 0; i < n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "circular_vectoring: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int((x < zero()) != (y < zero())) );
        if ( (x < zero()) != (y < zero()) ) {
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
void Cordic<T,FLT>::hyperbolic_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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
    uint32_t n = impl->n;
    uint32_t next_dup_i = 4;     
    for( uint32_t i = 1; i <= n; i++ )
    {
        T xi;
        T yi;
        T zi;
        if ( debug ) printf( "hyperbolic_rotation: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int(z >= zero()) );
        if ( z >= zero() ) {
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
    // d = (xy < 0) ? 1 : -1
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
        if ( (x < zero()) != (y < zero()) ) {
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
void Cordic<T,FLT>::linear_rotation( const T& x0, const T& y0, const T& z0, T& x, T& y, T& z ) const
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
    uint32_t n = impl->n;
    for( uint32_t i = 0; i <= n; i++ )
    {
        if ( debug ) printf( "linear_rotation: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int(z >= zero()) );
        T yi;
        T zi;
        if ( z >= zero() ) {
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
    // d = (y <= 0) ? 1 : -1
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
        if ( debug ) printf( "linear_vectoring: i=%d xyz=[%.30f,%.30f,%.30f] test=%d\n", i, to_flt(x), to_flt(y), to_flt(z), int((x < zero()) != (y < zero())) );
        T yi;
        T zi;
        if ( (x < zero()) != (y < zero()) ) {
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
T Cordic<T,FLT>::mad( const T& _x, const T& _y, const T addend, bool do_reduce ) const
{
    T x = _x;
    T y = _y;
    if ( debug ) std::cout << "mad begin: x_orig=" << to_flt(x) << " y_orig=" << to_flt(y) << " addend=" << to_flt(addend) << " do_reduce=" << do_reduce << "\n";
    cassert( do_reduce || addend >= 0 && "mad addend must be non-negative" );
    int32_t x_lshift;
    int32_t y_lshift;
    bool    sign;
    if ( do_reduce ) reduce_mul_args( x, y, x_lshift, y_lshift, sign );

    T xx, yy, zz;
    linear_rotation( x, do_reduce ? zero() : addend, y, xx, yy, zz );
    if ( do_reduce ) {
        impl->do_lshift( yy, x_lshift + y_lshift );
        yy += addend;
        if ( sign ) yy = -yy;
    }
    if ( debug ) std::cout << "mad end: x_orig=" << to_flt(x) << " y_orig=" << to_flt(y) << " addend=" << to_flt(addend) << " do_reduce=" << do_reduce << " yy=" << to_flt(yy) << 
                                  " x_lshift=" << x_lshift << " y_lshift=" << y_lshift << " sign=" << sign << "\n";
    return yy;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::mad( const T& _x, const T& _y, const T addend ) const
{
    return mad( _x, _y, addend, impl->do_reduce );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::mul( const T& x, const T& y ) const
{
    return mad( x, y, zero() );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::mul( const T& x, const T& y, bool do_reduce ) const
{
    return mad( x, y, zero(), do_reduce );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::dad( const T& _y, const T& _x, const T addend, bool do_reduce ) const
{
    T x = _x;
    T y = _y;
    if ( debug ) std::cout << "dad begin: x_orig=" << to_flt(x) << " y_orig=" << to_flt(y) << " addend=" << to_flt(addend) << " do_reduce=" << do_reduce << "\n";
    cassert( y != 0  && "dad y must be non-zero" );
    cassert( do_reduce || addend >= 0 && "dad addend must be non-negative (need to fix this soon)" );
    int32_t x_lshift;
    int32_t y_lshift;
    bool    sign;
    if ( do_reduce ) reduce_div_args( y, x, y_lshift, x_lshift, sign );

    T xx, yy, zz;
    linear_vectoring( x, y, do_reduce ? zero() : addend, xx, yy, zz );
    if ( do_reduce ) {
        impl->do_lshift( zz, y_lshift-x_lshift );
        zz += addend;
        if ( sign ) zz = -zz;
    }
    if ( debug ) std::cout << "dad end: x_orig=" << to_flt(x) << " y_orig=" << to_flt(y) << " addend=" << to_flt(addend) << " do_reduce=" << do_reduce <<
                              " zz_final=" << to_flt(zz) << " x_lshift=" << x_lshift << " y_lshift=" << y_lshift << " sign=" << sign << "\n";
    return zz;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::dad( const T& _y, const T& _x, const T addend ) const
{
    return dad( _y, _x, addend, impl->do_reduce );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::div( const T& y, const T& x ) const
{
    return dad( y, x, zero() );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::div( const T& y, const T& x, bool do_reduce ) const
{
    return dad( y, x, zero(), do_reduce );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::one_over( const T& x ) const
{
    return div( one(), x );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::sqrt( const T& x ) const
{ 
    return normh( x+quarter(), x-quarter() );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::one_over_sqrt( const T& x ) const
{ 
    // later, have normh leave result normalized
    return div( one(), normh( x+quarter(), x-quarter() ) );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::exp( const T& _x ) const
{ 
    T x = _x;
    T factor;
    bool sign;
    if ( impl->do_reduce ) reduce_exp_arg( M_E, x, factor, sign );

    T xx, yy, zz;
    hyperbolic_rotation( one_over_gainh(), one_over_gainh(), x, xx, yy, zz );
    if ( impl->do_reduce ) {
        if ( !sign ) {
            xx = mul( xx, factor, true );
        } else {
            xx = div( xx, factor, true );       // could do mul() of 1/factor but not as accurate
        }
    }
    return xx;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::pow( const T& b, const T& x ) const
{ 
    cassert( b > 0 && "pow b must be positive" );
    return exp( mul( x, log( b, true ) ) );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::powc( const FLT& b, const T& x ) const
{ 
    cassert( b > 0 && "powc b must be positive" );
    const FLT log_b_f = std::log( b );
    cassert( log_b_f >= 0.0 && "powc log(b) must be non-negative" );
    const T   log_b   = to_t( log_b_f );
    return exp( mul( x, log_b ) );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::pow2( const T& x ) const
{ 
    return powc( 2.0, x );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::pow10( const T& x ) const
{ 
    return powc( 10.0, x );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::log( const T& _x, bool do_reduce ) const
{ 
    T x = _x;
    cassert( x > 0 && "log: x must be positive" );
    T addend;
    if ( do_reduce ) reduce_log_arg( x, addend );
    T lg = atanh2( x-one(), x+one(), false ) << 1;
    if ( do_reduce ) lg += addend;
    if ( debug ) std::cout << "log: x_orig=" << to_flt(_x) << " reduced_x=" << to_flt(x) << " log=" << to_flt(lg) << "\n";
    return lg;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::log( const T& _x ) const
{ 
    return log( _x, impl->do_reduce );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::logb( const T& x, const T& b ) const
{ 
    cassert( b > 0 && "logb b must be positive" );
    return div( log(x), log(b) );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::logc( const T& x, const FLT& b ) const
{ 
    cassert( b > 0.0 && "logc b must be positive" );
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
T Cordic<T,FLT>::log2( const T& x ) const
{ 
    return logc( x, 2.0 );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::log10( const T& x ) const
{ 
    return logc( x, 10.0 );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::sin( const T& x, const T * r ) const
{ 
    T si;
    T co;
    sin_cos( x, si, co, impl->do_reduce, true, false, r );
    return si;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::cos( const T& x, const T * r ) const
{ 
    T si;
    T co;
    sin_cos( x, si, co, impl->do_reduce, false, true, r );
    return co;
}

template< typename T, typename FLT >
void Cordic<T,FLT>::sin_cos( const T& x, T& si, T& co, const T * r ) const             
{
    sin_cos( x, si, co, impl->do_reduce, true, true, r );
}

template< typename T, typename FLT >
void Cordic<T,FLT>::sin_cos( const T& _x, T& si, T& co, bool do_reduce, bool need_si, bool need_co, const T * _r ) const             
{ 
    T x = _x;
    uint32_t quadrant;
    bool x_sign;
    if ( do_reduce ) reduce_sin_cos_arg( x, quadrant, x_sign );

    T r = one_over_gain();
    int32_t r_lshift;
    bool r_sign = false;
    if ( _r != nullptr ) {
        r = *_r;
        if ( do_reduce ) reduce_arg( r, r_lshift, r_sign );
        r = mul( r, one_over_gain(), false );  // should not need to reduce
    }

    T zz;
    circular_rotation( r, zero(), x, co, si, zz );
    if ( do_reduce ) {
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
T Cordic<T,FLT>::tan( const T& x ) const
{ 
    T si, co;
    sin_cos( x, si, co );
    return div( si, co );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::asin( const T& x ) const
{ 
    return atan2( x, normh( one(), x ) );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::acos( const T& x ) const
{ 
    return atan2( normh( one(), x ), x );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::atan( const T& x ) const
{ 
    return atan2( x, one(), impl->do_reduce, true );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::atan2( const T& y, const T& x ) const
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
    // Strategy:
    //     Use reduce_atan2_args() to reduce y and x and get y_sign and x_sign.
    //     Return PI if we're done.
    //     Do atan2( y, norm(x, y) + x ) << 1   if x > 0
    //     Do atan2( norm(x, y) + x, y ) << 1   if <= 0
    //-----------------------------------------------------
    if ( debug ) std::cout << "atan2 begin: y=" << to_flt(y) << " x=" << to_flt(x) << " do_reduce=" << do_reduce << " x_is_one=" << x_is_one << "\n";
    cassert( (x != 0 || y != 0) && "atan2: x or y needs to be non-zero for result to be defined" );
    T xx, yy, zz;
    if ( r != nullptr ) *r = norm( _x, _y, do_reduce );  // optimize this later with below norm() of reduced x,y
    if ( do_reduce ) {
        bool y_sign;
        bool x_sign;
        bool is_pi;
        reduce_atan2_args( y, x, y_sign, x_sign, is_pi );
        if ( is_pi ) {
            if ( debug ) std::cout << "atan2 end: y=" << to_flt(_y) << " x=" << to_flt(_x) << " do_reduce=" << do_reduce << " x_is_one=" << x_is_one << 
                                      " zz=PI" << " r=" << ((r != nullptr) ? to_flt(*r) : to_flt(zero())) << "\n";
            return pi();
        }
        const T norm_plus_x = norm( x, y, false ) + x;
        if ( x > 0 ) {
            // atan2( y, norm_plus_x );
            if ( debug ) std::cout << "atan2 cordic begin: y=y=" << to_flt(y) << " x=norm_plus_x=" << to_flt(norm_plus_x) << "\n";
            circular_vectoring( norm_plus_x, y, zero(), xx, yy, zz );
        } else {
            // atan2( norm_plus_x, y );
            if ( debug ) std::cout << "atan2 cordic begin: y=norm_plus_x=" << to_flt(norm_plus_x) << " x=y=" << to_flt(y) << "\n";
            circular_vectoring( y, norm_plus_x, zero(), xx, yy, zz );
        }
        zz <<= 1;
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
void Cordic<T,FLT>::polar_to_rect( const T& r, const T& a, T& x, T& y ) const
{
    sin_cos( a, y, x, &r );
}

template< typename T, typename FLT >
void Cordic<T,FLT>::rect_to_polar( const T& x, const T& y, T& r, T& a ) const
{
    a = atan2( y, x, true, false, &r );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::norm( const T& x, const T& y ) const
{
    return norm( x, y, impl->do_reduce );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::norm( const T& _x, const T& _y, bool do_reduce ) const
{
    T x = _x;
    T y = _y;
    if ( debug ) std::cout << "norm begin: x=" << to_flt(x) << " y=" << to_flt(y) << " do_reduce=" << do_reduce << "\n";
    int32_t lshift;
    if ( do_reduce ) reduce_norm_args( x, y, lshift );

    T xx, yy, zz;
    circular_vectoring( x, y, zero(), xx, yy, zz );
    xx = mul( xx, one_over_gain() );
    if ( do_reduce ) impl->do_lshift( xx, lshift );
    if ( debug ) std::cout << "norm end: x=" << to_flt(x) << " y=" << to_flt(y) << " do_reduce=" << do_reduce << " xx=" << to_flt(xx) << "\n";
    return xx;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::normh( const T& _x, const T& _y ) const
{
    T x = _x;
    T y = _y;
    if ( debug ) std::cout << "normh begin: x=" << to_flt(x) << " y=" << to_flt(y) << " do_reduce=" << impl->do_reduce << "\n";
    int32_t lshift;
    if ( impl->do_reduce ) reduce_norm_args( x, y, lshift );
    cassert( x >= y && "normh abs(x) must be greater than abs(y)" );

    T xx, yy, zz;
    hyperbolic_vectoring( x, y, zero(), xx, yy, zz );
    xx = mul( xx, one_over_gainh(), false );   // should not need to do reduction
    if ( impl->do_reduce ) impl->do_lshift( xx, lshift );
    if ( debug ) std::cout << "normh end: x=" << to_flt(x) << " y=" << to_flt(y) << " do_reduce=" << impl->do_reduce << " xx=" << to_flt(xx) << "\n";
    return xx;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::sinh( const T& x, const T * r ) const
{ 
    T sih;
    T coh;
    sinh_cosh( x, sih, coh, impl->do_reduce, true, false, r );
    return sih;
}

template< typename T, typename FLT >
T Cordic<T,FLT>::cosh( const T& x, const T * r ) const
{ 
    T sih;
    T coh;
    sinh_cosh( x, sih, coh, impl->do_reduce, false, true, r );
    return coh;
}

template< typename T, typename FLT >
void Cordic<T,FLT>::sinh_cosh( const T& x, T& sih, T& coh, const T * r ) const
{ 
    sinh_cosh( x, sih, coh, impl->do_reduce, true, true, r );
}

template< typename T, typename FLT >
void Cordic<T,FLT>::sinh_cosh( const T& _x, T& sih, T& coh, bool do_reduce, bool need_sih, bool need_coh, const T * _r ) const
{ 
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
    //
    T x = _x;
    T sinh_i; 
    T cosh_i;
    bool sign;
    if ( do_reduce ) reduce_sinh_cosh_arg( x, sinh_i, cosh_i, sign );  

    T r = one_over_gainh();
    int32_t r_lshift;
    bool r_sign = false;
    if ( _r != nullptr ) {
        r = *_r;
        if ( do_reduce ) reduce_arg( r, r_lshift, r_sign );
        r = mul( r, one_over_gainh(), false );  // should not need to reduce (I think)
    }

    T sinh_f;
    T cosh_f;
    T zz;
    hyperbolic_rotation( r, zero(), x, coh, sih, zz );
    if ( do_reduce ) {
        if ( need_sih ) {
            sih = mul( sinh_f, cosh_i, true ) + mul( cosh_f, sinh_i, true );
            if ( sign ) sih = -sih;
        }
        if ( need_coh ) {
            coh = mul( cosh_f, cosh_i, true ) - mul( sinh_f, sinh_i, true );
        }
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
T Cordic<T,FLT>::acosh( const T& x ) const
{ 
    return log( x + normh( x, one() ) );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::atanh( const T& x ) const
{ 
    return atanh2( x, one(), impl->do_reduce, true );
}

template< typename T, typename FLT >
T Cordic<T,FLT>::atanh2( const T& y, const T& x ) const             
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
    if ( do_reduce ) reduce_div_args( y, x, y_lshift, x_lshift, sign );
    int32_t lshift = x_lshift + y_lshift;
    cassert( lshift == 0 && "atanh2: abs(y/x) must be between 0 and 1" );

    T xx, yy, zz;
    hyperbolic_vectoring( x, y, zero(), xx, yy, zz );
    if ( sign ) zz = -zz;
    return zz;
}

template< typename T, typename FLT >
void Cordic<T,FLT>::reduce_arg( T& x, int32_t& x_lshift, bool& sign, bool shift_x, bool normalize ) const
{
    T x_orig = x;
    sign = x < 0;
    if ( sign ) x = -x;
    T other = T(1) << frac_w();
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
    while( normalize && x < one() )
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

template< typename T, typename FLT >
void Cordic<T,FLT>::reduce_mul_args( T& x, T& y, int32_t& x_lshift, int32_t& y_lshift, bool& sign ) const
{
    if ( debug ) std::cout << "reduce_mul_args: x_orig=" << to_flt(x) << " y_orig=" << to_flt(y) << "\n";
    bool x_sign;
    bool y_sign;
    reduce_arg( x, x_lshift, x_sign );
    reduce_arg( y, y_lshift, y_sign );
    sign = x_sign ^ y_sign;
}

template< typename T, typename FLT >
void Cordic<T,FLT>::reduce_div_args( T& y, T& x, int32_t& y_lshift, int32_t& x_lshift, bool& sign ) const
{
    if ( debug ) std::cout << "reduce_div_args: x_orig=" << to_flt(x) << " y_orig=" << to_flt(y) << "\n";
    bool x_sign;
    bool y_sign;
    reduce_arg( y, y_lshift, y_sign );
    reduce_arg( x, x_lshift, x_sign, true, true );
    sign = x_sign ^ y_sign;
}

template< typename T, typename FLT >
void Cordic<T,FLT>::reduce_exp_arg( FLT b, T& x, T& factor, bool& sign ) const 
{
    //-----------------------------------------------------
    // Identity:
    //     Assume: x = i + f  (integer plus fraction)
    //     exp(x) = exp(i) * exp(f)
    //     pow(b,x) = log(b) * exp(x) = [log(b)*exp(i)] * exp(f)
    //
    // If x is non-negative:
    //     exp(i) comes from a pre-built LUT kept in FLT
    //     so we can multiply it by log(b) before converting to type T and
    //     then multiplying by exp(f) in the caller.
    //
    // If x is negative:
    //     x = -x
    //     [do above but the callere will divide by factor rather than multiplying]
    //-----------------------------------------------------
    T x_orig = x;
    const FLT * factors_f = impl->reduce_exp_factor.get();
    sign         = x < 0;
    if ( sign ) x = -x;
    T   index    = (x >> frac_w()) & maxint();
    FLT factor_f = std::log(b) * factors_f[index];   // could build per-b factors_f[] LUT with multiply already done
    factor       = to_t( factor_f );
    x           &= (T(1) << frac_w())-T(1); // fraction only
    if ( debug ) std::cout << "reduce_exp_arg: x_orig=" << to_flt(x_orig) << " x_reduced=" << to_flt(x) << " factor=" << factor << "\n"; 
}

template< typename T, typename FLT >
void Cordic<T,FLT>::reduce_log_arg( T& x, T& addend ) const 
{
    //-----------------------------------------------------
    // log(ab) = log(a) + log(b)
    // 
    // Normalize x so that it's in 1.00 .. 2.00.
    // Then addend = log(1 << lshift).
    //-----------------------------------------------------
    cassert( x >= 0 && "log argument may not be negative for fixed-point numbers" );
    T x_orig = x;
    int32_t x_lshift;
    bool x_sign;
    reduce_arg( x, x_lshift, x_sign, true, true );
    const T * addends = impl->reduce_log_addend.get();
    addend = addends[frac_w()+x_lshift];
    if ( debug ) std::cout << "reduce_log_arg: x_orig=" << to_flt(x_orig) << " x_reduced=" << to_flt(x) << " addend=" << to_flt(addend) << "\n"; 
}

template< typename T, typename FLT >
void Cordic<T,FLT>::reduce_atan2_args( T& y, T& x, bool& y_sign, bool& x_sign, bool& is_pi ) const
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
    cassert( (x != 0 || y != 0) && "atan2: x or y needs to be non-zero for result to be defined" );

    x_sign = x < 0;
    y_sign = y < 0;
    if ( x_sign && y == 0 ) {
        // PI
        is_pi = true;
        if ( debug ) std::cout << "reduce_atan2_args: xy_orig=[" << to_flt(x_orig) << "," << to_flt(y_orig) << "]" << " PI returned\n";
    } else {
        is_pi = false;
        int32_t lshift;
        reduce_norm_args( x, y, lshift );
        if ( debug ) std::cout << "reduce_atan2_args: xy_orig=[" << to_flt(x_orig) << "," << to_flt(y_orig) << "]" << 
                          " xy_reduced=[" << to_flt(x) << "," << to_flt(y) << "] y_sign=" << y_sign << " x_sign=" << x_sign << "\n";
    }

}

template< typename T, typename FLT >
void Cordic<T,FLT>::reduce_norm_args( T& x, T& y, int32_t& lshift ) const
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
    if ( debug ) std::cout << "reduce_norm_args: xy_orig=[" << to_flt(x_orig) << "," << to_flt(y_orig) << "]" << 
                                               " xy_reduced=[" << to_flt(x) << "," << to_flt(y) << "] lshift=" << lshift << "\n"; 
}

template< typename T, typename FLT >
void Cordic<T,FLT>::reduce_sin_cos_arg( T& a, uint32_t& quad, bool& sign ) const
{
    //-----------------------------------------------------
    // Use LUT to find addend.
    //-----------------------------------------------------
    const T a_orig = a;
    sign = a < 0;
    if ( sign ) a = -a;
    const T *  addend   = impl->reduce_sin_cos_addend.get();
    uint32_t * quadrant = impl->reduce_sin_cos_quadrant.get();
    const T MASK = (T(1) << (int_w()+T(1))) - T(1);  // include 0.5 bit of fraction
    T i = (a >> (frac_w()-T(1))) & MASK;
    quad = quadrant[i];
    a += addend[i];
    if ( debug ) std::cout << "reduce_sin_cos_arg: a_orig=" << to_flt(a_orig) << " addend[" << i << "]=" << to_flt(addend[i]) << 
                              " a_reduced=" << to_flt(a) << " quadrant=" << quad << "\n"; 
}

template< typename T, typename FLT >
void Cordic<T,FLT>::reduce_sinh_cosh_arg( T& _x, T& sinh_i, T& cosh_i, bool& sign ) const
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
    sign = x < 0;
    if ( sign ) x = -x;

    const T MASK = (T(1) << (int_w()+T(1))) - T(1);  // include 0.5 bit of fraction
    T i = (x >> (frac_w()-T(1))) & MASK;
    x   = x & ((T(1) << (frac_w()-T(1)))-T(1));
    const T *  sinh_i_vals = impl->reduce_sinh_cosh_sinh_i.get();
    const T *  cosh_i_vals = impl->reduce_sinh_cosh_cosh_i.get();
    sinh_i = sinh_i_vals[i];
    cosh_i = cosh_i_vals[i];
}

template class Cordic<int64_t, double>;
