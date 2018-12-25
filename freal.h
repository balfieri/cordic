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
// freal.h - flexible real number 
//
// This class supports typical C++ floating-point semantics found with float and double.
// It implements freals as a user-specified fixed-point.
// All functions are implemented using Cordic.h.
//
// In the near future, it will allow the internal format to change
// dynamically based on input and output ranges, and the 
// format could be fixed-point or floating-point.
//
// Typical usage:
//
//     #include freal.h
//     using real = freal<>;
//     ...
//     real f;                                  // undefined, no type yet
//     f = real( 8, 21, 25.573822 );            // assigned using 1.8.21 fixed-point value 25.573822
//     real g = f;                              // g copies type and value of f
//
#ifndef _freal_h
#define _freal_h

#include "Cordic.h"

// T      = some signed integer type that can hold fixed-point values (default is int64_t)
// FLT    = some floating-point type that can hold constants of the desired precision (default is double)
//
template< typename T=int64_t, typename FLT=double >              
class freal
{
public:
    //-----------------------------------------------------
    // Constructors
    //-----------------------------------------------------
    freal( void );                                      // initializes value to undefined
    freal( const freal& other );                        // copy type and value of other
    freal( const freal& other, FLT f );                 // copy type of other, but value of f
    freal( const Cordic<T,FLT> * cordic, FLT f );       // use type from cordic, but value of f

    static freal make_fixed( uint32_t int_w, uint32_t frac_w, FLT init_f=FLT(0) );  // make a signed fixed-point    number 
    static freal make_float( uint32_t exp_w, uint32_t frac_w, FLT init_f=FLT(0) );  // make a signed floating-point number
    ~freal();

    //-----------------------------------------------------
    // Explicit Conversions
    //-----------------------------------------------------
    FLT         to_flt( void ) const;                   // freal to FLT
    std::string to_string( void ) const;                // freal to std::string

    //-----------------------------------------------------
    // Implicit Conversions
    //
    // IMPORTANT: implicit_to_set() must be called before 
    // using any implicit conversions TO freal.
    //
    // IMPORTANT: implicit_from_set() must be called before
    // using any implicit conversions FROM freal.
    //-----------------------------------------------------
    static void implicit_to_set( uint32_t int_w, uint32_t frac_w, bool is_fixed_point=true );
    static void implicit_to_set( const Cordic<T,FLT> * cordic );
    static void implicit_from_set( bool allow );

    freal( FLT f );
    freal( uint64_t i );
    freal( int64_t i );
    freal( uint32_t i );
    freal( int32_t i );

    operator FLT( void );
    operator uint64_t( void );
    operator int64_t( void );
    operator uint32_t( void );
    operator int32_t( void );

    //-----------------------------------------------------
    // Constants
    //-----------------------------------------------------               
    T     maxint( void );                                       // largest positive integer (just integer part, does not include fraction)
    freal max( void );                                          // maximum positive value 
    freal min( void );                                          // minimum positive value
    freal denorm_min( void );                                   // minimum positive denorm value
    freal lowest( void );                                       // most negative value
    freal epsilon( void );                                      // difference between 1 and first number above 1
    freal round_error( void );                                  // maximum rounding error
    freal zero( void );                                         // 0.0
    freal one( void );                                          // 1.0
    freal two( void );                                          // 2.0
    freal half( void );                                         // 0.5
    freal quarter( void );                                      // 0.25
    freal sqrt2( void );                                        // sqrt(2)
    freal sqrt2_div_2( void );                                  // sqrt(2)/2
    freal pi( void );                                           // PI
    freal tau( void );                                          // 2*PI
    freal pi_div_2( void );                                     // PI/2
    freal pi_div_4( void );                                     // PI/4
    freal two_div_pi( void );                                   // 2/PI
    freal four_div_pi( void );                                  // 4/PI
    freal e( void );                                            // natural exponent
    freal nan( const char * arg );                              // not-a-number (NaN)
    freal quiet_nan( void );                                    // quiet not-a-number (NaN)
    freal signaling_nan( void );                                // signaling not-a-number (NaN)
    freal inf( void );                                          // +infinity
    freal ninf( void );                                         // -infinity

    //-----------------------------------------------------
    // Standard Operators
    //-----------------------------------------------------               
    freal  operator -  ()                 const;                // -x

    freal  operator +  ( const freal& b ) const;
    freal  operator -  ( const freal& b ) const;
    freal  operator *  ( const freal& b ) const;
    freal  operator /  ( const freal& b ) const;
    freal  operator << (       int    b ) const;
    freal  operator >> (       int    b ) const;

    freal& operator =  ( const freal& b );
    freal& operator =  ( const FLT&   b );
    freal& operator += ( const freal& b );
    freal& operator -= ( const freal& b );
    freal& operator *= ( const freal& b );
    freal& operator /= ( const freal& b );
    freal& operator <<=(       int    b );
    freal& operator >>=(       int    b );

    bool   operator >  ( const freal& b ) const;
    bool   operator >= ( const freal& b ) const;
    bool   operator <  ( const freal& b ) const;
    bool   operator <= ( const freal& b ) const;
    bool   operator != ( const freal& b ) const;
    bool   operator == ( const freal& b ) const;

    //-----------------------------------------------------
    // Well-Known Math Operators and Functions 
    //
    // See Cordic.h for functionality of each of these.
    //
    // a = *this
    //-----------------------------------------------------               
    freal& assign( const freal& b );

    bool   signbit( void ) const;    
    freal  frexp( int * e ) const;
    freal  modf( freal * i ) const;
    int    ilogb( void ) const;
    freal  logb( void ) const;
    int    fpclassify( void ) const; 
    bool   isfinite( void ) const;  
    bool   isinf( void ) const;    
    bool   isnan( void ) const;   
    bool   isnormal( void ) const;

    int    fesetround( int round ) const;
    int    fegetround( void ) const;
    freal  nextafter( const freal& to ) const;
    freal  nexttoward( long double to ) const;
    freal  floor( void ) const;
    freal  ceil( void ) const;
    freal  trunc( void ) const;
    freal  round( void ) const;
    long   lround( void ) const;
    long long llround( void ) const;
    T      iround( void ) const;
    freal  rint( void ) const;
    long   lrint( void ) const;
    long long llrint( void ) const;
    T      irint( void ) const;
    freal  nearbyint( void ) const;

    freal  abs( void ) const;
    freal  neg( void ) const; 
    freal  copysign( const freal& b ) const;
    freal  add( const freal& b ) const; 
    freal  sub( const freal& b ) const; 
    freal  fma( const freal& b, const freal& c ) const;             
    freal  mul( const freal& b ) const;                             
    freal  sqr( void ) const;                             
    freal  scalbn( int b ) const;
    freal  scalbnn( int b ) const;                              // scalbn( -b )
    freal  ldexp( int b ) const;
    freal  fda( const freal& b, const freal& c ) const;      
    freal  div( const freal& b ) const;      
    freal  remainder( const freal& b ) const;
    freal  fmod( const freal& b ) const;
    freal  remquo( const freal& x, int * quo ) const;
    freal  rcp( void ) const;                                    

    bool   isgreater( const freal& b ) const;                        
    bool   isgreaterequal( const freal& b ) const;                 
    bool   isless( const freal& b ) const;                          
    bool   islessequal( const freal& b ) const;                   
    bool   islessgreater( const freal& b ) const;                
    bool   isunordered( const freal& b ) const;              
    bool   isunequal( const freal& b ) const;                 
    bool   isequal( const freal& b ) const;                   
    freal  fdim( const freal& b ) const;
    freal  fmax( const freal& b ) const;
    freal  fmin( const freal& b ) const;

    freal  sqrt( void ) const;                                        
    freal  rsqrt( void ) const;                               
    freal  cbrt( void ) const;                                        
    freal  rcbrt( void ) const;                               

    freal  exp( void ) const;                                         
    freal  expm1( void ) const;              // exp(x) - 1   (accurately)
    freal  expc( const FLT c ) const;        // c^a
    freal  exp2( void ) const;               // 2^x
    freal  exp10( void ) const;              // 10^x
    freal  pow( const freal& e ) const;      // a^e
    freal  log( void ) const;                // log base-e
    freal  log( const freal& b ) const;      // log base-b
    freal  log1p( void ) const;              // log base-e (a+1)
    freal  logc( const FLT c ) const;        // log base-c (c is a constant)
    freal  log2( void ) const;               // log base-2
    freal  log10( void ) const;              // log base-10

    freal  sin( void ) const;
    freal  sinpi( void ) const;
    freal  sin( const freal& r ) const;                                 // multiply sin by r
    freal  sinpi( const freal& r ) const;                               // multiply sin by r
    freal  cos( void ) const;
    freal  cospi( void ) const;
    freal  cos( const freal& r ) const;                                 // multiply cos by r
    freal  cospi( const freal& r ) const;                               // multiply cos by r
    void   sincos( freal& si, freal& co ) const;
    void   sinpicospi( freal& si, freal& co ) const;
    void   sincos( freal& si, freal& co, const freal& r ) const;        // multiply sin and cos by r
    void   sinpicospi( freal& si, freal& co, const freal& r ) const;    // multiply sin and cos by r
    freal  tan( void ) const;                                         
    freal  tanpi( void ) const;                                         

    freal  asin( void ) const;                                        
    freal  acos( void ) const;                                        
    freal  atan( void ) const;                                        
    freal  atan2( const freal& b ) const;    // y=a, x=b

    void   polar_to_rect( const freal& angle, freal& x, freal& y     ) const;  // a=radius
    void   rect_to_polar( const freal& b,     freal& r, freal& angle ) const;  // x=a, y=b 
    freal  hypot(  const freal& b ) const;    
    freal  hypoth( const freal& b ) const;   

    freal  sinh( void ) const;
    freal  sinh( const freal& r ) const;                                // multiply sinh by r
    freal  cosh( void ) const;
    freal  cosh( const freal& r ) const;                                // multiply cosh by r
    void   sinhcosh( freal& sih, freal& coh ) const;
    void   sinhcosh( freal& sih, freal& coh, const freal& r ) const;    // multiply sinh and cosh by r
    freal  tanh( void ) const;                                        
    freal  asinh( void ) const;                                       
    freal  acosh( void ) const;                                       
    freal  atanh( void ) const;                                       
    freal  atanh2( const freal& b ) const;   // atanh2( a, b )

    //-----------------------------------------------------
    // Introspection
    //-----------------------------------------------------
    const Cordic<T,FLT> * c( void ) const;                            // validates current cordic  and returns it
    const Cordic<T,FLT> * c( const freal& b ) const;                  // validates two     cordics and returns one to use for operation
    const Cordic<T,FLT> * c( const freal& b, const freal& _c ) const; // validates three   cordics and returns one to use for operation

    const T * raw_ptr( void ) const;                                  // useful for some gross things like manual logging by callers

private:
    static const Cordic<T,FLT> * implicit_to;
    static bool                  implicit_from;

    const Cordic<T,FLT> *        cordic;         // defines the type and most operations
    T                            v;              // this value encoded in type T

    static freal pop_value( const Cordic<T,FLT> * cordic, const T& encoded );   // pop value associated with last operation
    static bool  pop_bool(  const Cordic<T,FLT> * cordic, bool );               // pop bool  associated with last operation 
};

// Well-Known std:xxx() Functions 
//
namespace std
{

template< typename T=int64_t, typename FLT=double >              
static inline std::istream& operator >> ( std::istream &in, freal<T,FLT>& a )
{ 
    FLT a_f;
    in >> a_f; 
    const Cordic<T,FLT> * cordic = a.c();
    if ( cordic != nullptr ) {
        a = freal( cordic, a_f );     // use a's current format
    } else {
        a = a_f;                      // rely on implicit conversion, if it's currently allowed 
    }
    return in;
}

template< typename T=int64_t, typename FLT=double >              
static inline std::ostream& operator << ( std::ostream &out, const freal<T,FLT>& a )
{ 
    out << a.to_string(); 
    return out;     
}

#define _freal freal<T,FLT>

#define decl_std1( name )                                       \
    template< typename T=int64_t, typename FLT=double >         \
    static inline _freal name( const _freal& a )                \
    { return a.name(); }                                        \

#define decl_std1_ret( name, ret_type )                         \
    template< typename T=int64_t, typename FLT=double >         \
    static inline ret_type name( const _freal& a )              \
    { return a.name(); }                                        \

#define decl_std2( name )                                       \
    template< typename T=int64_t, typename FLT=double >         \
    static inline _freal name( const _freal& a, const _freal& b ) \
    { return a.name( b ); }                                     \

#define decl_std2_ret( name, ret_type )                         \
    template< typename T=int64_t, typename FLT=double >         \
    static inline ret_type name( const _freal& a, const _freal& b ) \
    { return a.name( b ); }                                     \

#define decl_std2x( name, b_type )                              \
    template< typename T=int64_t, typename FLT=double >         \
    static inline _freal name( const _freal& a, b_type b )      \
    { return a.name( b ); }                                     \

#define decl_std3( name )                                       \
    template< typename T=int64_t, typename FLT=double >         \
    static inline _freal name( const _freal& a, const _freal& b, const _freal& c ) \
    { return a.name( b, c ); }                                  \

#define decl_std3x( name, c_type )                              \
    template< typename T=int64_t, typename FLT=double >         \
    static inline _freal name( const _freal& a, const _freal& b, c_type c ) \
    { return a.name( b, c ); }                                  \

decl_std1_ret( signbit,         bool            )
decl_std2x(    frexp,           int *           )
decl_std2x(    modf,            _freal *        )
decl_std1_ret( ilogb,           int             )
decl_std1(     logb                             )
decl_std1_ret( fpclassify,      int             )
decl_std1_ret( isfinite,        bool            )
decl_std1_ret( isinf,           bool            )
decl_std1_ret( isnan,           bool            )
decl_std1_ret( isnormal,        bool            )
decl_std1_ret( to_string,       std::string     )

template< typename T=int64_t, typename FLT=double >              
static inline int fesetround( int round ) 
{ return freal<T,FLT>::implicit_to->fesetround( round ); }

template< typename T=int64_t, typename FLT=double >              
static inline int fegetround( void ) 
{ return freal<T,FLT>::implicit_to->fegetround(); }

decl_std2(     nextafter                        );
decl_std2x(    nexttoward,      long double     );
decl_std1(     floor                            );
decl_std1(     ceil                             );
decl_std1(     trunc                            );
decl_std1(     round                            );
decl_std1_ret( lround,          long            );
decl_std1_ret( llround,         long long       );
decl_std1_ret( iround,          T               );
decl_std1(     rint                             );
decl_std1_ret( lrint,           long            );
decl_std1_ret( llrint,          long long       );
decl_std1_ret( irint,           T               );
decl_std1(     nearbyint                        );
decl_std1(     abs                              );
decl_std1(     neg                              );
decl_std2(     copysign                         );
decl_std2(     add                              );
decl_std2(     sub                              );
decl_std3(     fma                              );
decl_std2(     mul                              );
decl_std1(     sqr                              );
decl_std2x(    scalbn,          int             );
decl_std2x(    ldexp,           int             );
decl_std3(     fda                              );
decl_std2(     div                              );
decl_std2(     remainder                        );
decl_std2(     fmod                             );
decl_std3x(    remquo,          int             );
decl_std1(     rcp                              );
decl_std2_ret( isgreater,       bool            );
decl_std2_ret( isgreaterequal,  bool            );
decl_std2_ret( isless,          bool            );
decl_std2_ret( islessequal,     bool            );
decl_std2_ret( islessgreater,   bool            );
decl_std2_ret( isunordered,     bool            );
decl_std2_ret( isunequal,       bool            );
decl_std2_ret( isequal,         bool            );
decl_std2(     fdim                             );
decl_std2(     fmin                             );
decl_std2(     fmax                             );
decl_std1(     sqrt                             );
decl_std1(     rsqrt                            );
decl_std1(     cqrt                             );
decl_std1(     rcqrt                            );
decl_std1(     exp                              );
decl_std1(     expm1                            );
decl_std2x(    expc,            FLT             );
decl_std1(     exp2                             );
decl_std1(     exp10                            );
decl_std2(     pow                              );
decl_std1(     log                              );
decl_std2(     log                              );
decl_std1(     log1p                            );
decl_std2x(    logc,            FLT             );
decl_std1(     log2                             );
decl_std1(     log10                            );
decl_std1(     sin                              );
decl_std1(     sinpi                            );
decl_std2(     sin                              );
decl_std2(     sinpi                            );
decl_std1(     cos                              );
decl_std1(     cospi                            );
decl_std2(     cos                              );
decl_std2(     cospi                            );

template< typename T=int64_t, typename FLT=double >              
static inline void   sincos( const freal<T,FLT>& a, freal<T,FLT>& si, freal<T,FLT>& co )               
{ a.sincos( si, co );                   }

template< typename T=int64_t, typename FLT=double >              
static inline void   sinpicospi( const freal<T,FLT>& a, freal<T,FLT>& si, freal<T,FLT>& co )               
{ a.sinpicospi( si, co );               }

template< typename T=int64_t, typename FLT=double >              
static inline void   sincos( const freal<T,FLT>& a, freal<T,FLT>& si, freal<T,FLT>& co, const freal<T,FLT>& r )               
{ a.sincos( si, co, r );                }

template< typename T=int64_t, typename FLT=double >              
static inline void   sinpicospi( const freal<T,FLT>& a, freal<T,FLT>& si, freal<T,FLT>& co, const freal<T,FLT>& r )               
{ a.sinpicospi( si, co, r );            }

decl_std1(     tan                              );
decl_std1(     tanpi                            );
decl_std1(     asin                             );
decl_std1(     acos                             );
decl_std1(     atan                             );
decl_std2(     atan2                            );

template< typename T=int64_t, typename FLT=double >              
static inline void   polar_to_rect( const freal<T,FLT>& a, const freal<T,FLT>& angle, freal<T,FLT>& x, freal<T,FLT>& y )  
{ a.polar_to_rect( angle, x, y );       }

template< typename T=int64_t, typename FLT=double >              
static inline void   rect_to_polar( const freal<T,FLT>& a, const freal<T,FLT>& b,     freal<T,FLT>& r, freal<T,FLT>& angle )  
{ a.rect_to_polar( b, r, angle );       }

decl_std2(     hypot                            );
decl_std2(     hypoth                           );
decl_std1(     sinh                             );
decl_std2(     sinh                             );
decl_std1(     cosh                             );
decl_std2(     cosh                             );

template< typename T=int64_t, typename FLT=double >              
static inline void   sinhcosh( const freal<T,FLT>& a, freal<T,FLT>& sih, freal<T,FLT>& coh )           
{ a.sinhcosh( sih, coh );               }

template< typename T=int64_t, typename FLT=double >              
static inline void   sinhcosh( const freal<T,FLT>& a, freal<T,FLT>& sih, freal<T,FLT>& coh, const freal<T,FLT>& r )           
{ a.sinhcosh( sih, coh, r );            }

decl_std1(     tanh                             );
decl_std1(     asinh                            );
decl_std1(     acosh                            );
decl_std1(     atanh                            );
decl_std2(     atanh2                           );

}

template< typename T, typename FLT >              
class std::numeric_limits<freal<T,FLT>> 
{
public:
    // any static field that depends on the current implicit_to is (re)initialized
    // when implicit_to_set() is called to change the implicit Cordic, therefore
    // these static fields are not marked const
    static bool                 is_specialized;
    static freal<T,FLT>         min() throw()           { return freal<T,FLT>::c()->min(); }
    static freal<T,FLT>         max() throw()           { return freal<T,FLT>::c()->max(); }
    static int                  digits;
    static int                  digits10;
    static const bool           is_signed = true;
    static const bool           is_integer = false;
    static const bool           is_exact = false;
    static const int            radix = 2;
    static freal<T,FLT>         epsilon() throw()       { return freal<T,FLT>::c()->epsilon(); }
    static freal<T,FLT>         round_error() throw()   { return freal<T,FLT>::c()->round_error(); }

    static int                  min_exponent;
    static int                  min_exponent10;
    static int                  max_exponent;
    static int                  max_exponent10;

    static bool                 has_infinity;
    static bool                 has_quiet_NaN;
    static bool                 has_signaling_NaN;
    static const float_denorm_style has_denorm = denorm_present;
    static const bool           has_denorm_loss = false;
    static freal<T,FLT>         infinity() throw()      { return freal<T,FLT>::c()->inf(); }
    static freal<T,FLT>         quiet_NaN() throw()     { return freal<T,FLT>::c()->quiet_nan(); }
    static freal<T,FLT>         signaling_NaN() throw() { return freal<T,FLT>::c()->signaling_nan(); }
    static freal<T,FLT>         denorm_min() throw()    { return freal<T,FLT>::c()->denorm_min(); }

    static bool                 is_iec559;
    static bool                 is_bounded;
    static bool                 is_modulo;
    static bool                 traps;
    static bool                 tinyness_before;
    static float_round_style    round_style;
};

template< typename T, typename FLT >              
bool std::numeric_limits<freal<T,FLT>>::is_specialized = false; 

//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//
// IMPLEMENTATION
//
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------

//-----------------------------------------------------
// Static Globals
//-----------------------------------------------------
template< typename T, typename FLT >              
const Cordic<T,FLT> * freal<T,FLT>::implicit_to = nullptr;  // disallow

template< typename T, typename FLT >              
bool                  freal<T,FLT>::implicit_from = false;  // disallow

//-----------------------------------------------------
// Constructors
//-----------------------------------------------------
template< typename T, typename FLT >              
inline freal<T,FLT>::freal( void )
{
    cordic = nullptr;
    v      = T(666);
}

template< typename T, typename FLT >              
inline freal<T,FLT>::~freal()
{
    if ( cordic != nullptr ) {
        cordic->destructed( v );
        cordic = nullptr;
    }
    v = T(668);
}

template< typename T, typename FLT >              
inline freal<T,FLT>::freal( const Cordic<T,FLT> * _cordic, FLT f )
{
    cassert( _cordic != nullptr, "freal(cordic, f) cordic argument must be non-null" );
    cordic = _cordic;
    cordic->constructed( v );
    cordic->pop_value( v, cordic->to_t( f, true ) );
}

template< typename T, typename FLT >              
inline freal<T,FLT>::freal( const freal& other )
{
    cordic = other.c();
    v      = other.v;
    cordic->constructed( v );
    cordic->assign( v, other.v ); 
}

template< typename T, typename FLT >              
inline freal<T,FLT>::freal( const freal& other, FLT f )
{
    freal( other.cordic, f );
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::make_fixed( uint32_t int_w, uint32_t frac_w, FLT init_f )
{
    return freal( new Cordic<T,FLT>( int_w, frac_w ), init_f );
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::make_float( uint32_t exp_w, uint32_t frac_w, FLT init_f )
{
    cassert( false, "can't encode floating-point values right now" );
    (void)exp_w;   // unused
    (void)frac_w;
    (void)init_f;
    return freal();
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::pop_value( const Cordic<T,FLT> * cordic, const T& encoded )
{
    cassert( cordic != nullptr, "pop_value(cordic, encoded) called with null cordic" );
    freal r;
    r.cordic = cordic;
    cordic->constructed( r.v );
    cordic->pop_value( r.v, encoded );
    return r;
}

template< typename T, typename FLT >              
inline bool freal<T,FLT>::pop_bool( const Cordic<T,FLT> * cordic, bool b )
{
    cassert( cordic != nullptr, "pop_bool(cordic, b) called with null cordic" );
    cordic->pop_bool( b );
    return b;
}

//-----------------------------------------------------
// Explicit Conversions
//-----------------------------------------------------
template< typename T, typename FLT >              
FLT    freal<T,FLT>::to_flt( void ) const                                                               
{ return c()->to_flt( v );              }

template< typename T, typename FLT >              
std::string freal<T,FLT>::to_string( void ) const                                                               
{ return c()->to_string( v );           }

//-----------------------------------------------------
// Implicit Conversions
//-----------------------------------------------------
template< typename T, typename FLT >              
inline void freal<T,FLT>::implicit_to_set( const Cordic<T,FLT> * cordic )
{ 
    implicit_to = cordic;

    if ( cordic != nullptr ) {
        std::numeric_limits<freal<T,FLT>>::is_specialized       = true;
        std::numeric_limits<freal<T,FLT>>::digits               = cordic->int_w() + cordic->frac_w();
        std::numeric_limits<freal<T,FLT>>::digits10             = std::numeric_limits<freal<T,FLT>>::digits / std::log2( 10 );
        std::numeric_limits<freal<T,FLT>>::min_exponent         = 0;
        std::numeric_limits<freal<T,FLT>>::min_exponent10       = 0;
        std::numeric_limits<freal<T,FLT>>::max_exponent         = 0;
        std::numeric_limits<freal<T,FLT>>::max_exponent10       = 0;
        std::numeric_limits<freal<T,FLT>>::has_inifinity        = false;
        std::numeric_limits<freal<T,FLT>>::has_quiet_NaN        = false;
        std::numeric_limits<freal<T,FLT>>::has_signaling_NaN    = false;
        std::numeric_limits<freal<T,FLT>>::is_iec559            = false;
        std::numeric_limits<freal<T,FLT>>::is_bounded           = true;
        std::numeric_limits<freal<T,FLT>>::is_modulo            = true;
        std::numeric_limits<freal<T,FLT>>::traps                = false;
        std::numeric_limits<freal<T,FLT>>::tinyness_before      = true;
        std::numeric_limits<freal<T,FLT>>::round_style          = std::numeric_limits<freal<T,FLT>>::round_to_nearest;
    } else {
        std::numeric_limits<freal<T,FLT>>::is_specialized       = false;
    }
}

template< typename T, typename FLT >              
inline void freal<T,FLT>::implicit_to_set( uint32_t int_w, uint32_t frac_w, bool is_fixed_point )
{ 
    cassert( is_fixed_point, "implicit_to_set() floating-point is not implemented yet" );
    implicit_to = new Cordic<T,FLT>( int_w, frac_w );
}

template< typename T, typename FLT >              
inline void freal<T,FLT>::implicit_from_set( bool allow )
{ 
    implicit_from = allow;
}

template< typename T, typename FLT >              
inline freal<T,FLT>::freal( FLT f )
{
    cassert( implicit_to != nullptr, "implicit_to_set() must be called before relying on any implicit from FLT to freal<>" );
    *this = freal( implicit_to, f );
}

template< typename T, typename FLT >              
inline freal<T,FLT>::freal( uint64_t i )
{
    cassert( implicit_to != nullptr, "implicit_to_set() must be called before relying on any implicit from uint64_t to freal<>" );
    *this = freal( implicit_to, FLT(i) );
}

template< typename T, typename FLT >              
inline freal<T,FLT>::freal( int64_t i )
{
    cassert( implicit_to != nullptr, "implicit_to_set() must be called before relying on any implicit from int64_t to freal<>" );
    *this = freal( implicit_to, FLT(i) );
}

template< typename T, typename FLT >              
inline freal<T,FLT>::freal( uint32_t i )
{
    cassert( implicit_to != nullptr, "implicit_to_set() must be called before relying on any implicit from uint32_t to freal<>" );
    *this = freal( implicit_to, FLT(i) );
}

template< typename T, typename FLT >              
inline freal<T,FLT>::freal( int32_t i )
{
    cassert( implicit_to != nullptr, "implicit_to_set() must be called before relying on any implicit from int32_t to freal<>" );
    *this = freal( implicit_to, FLT(i) );
}

template< typename T, typename FLT >              
inline freal<T,FLT>::operator FLT( void )
{ 
    cassert( implicit_from, "implicit_from_set( true ) must be called before relying on any implicit from freal<> to FLT" );
    return to_flt();
}

template< typename T, typename FLT >              
inline freal<T,FLT>::operator uint64_t( void )
{ 
    cassert( implicit_from, "implicit_from_set( true ) must be called before relying on any implicit from freal<> to uint64_t" );
    return uint64_t( to_flt() );
}

template< typename T, typename FLT >              
inline freal<T,FLT>::operator int64_t( void )
{ 
    cassert( implicit_from, "implicit_from_set( true ) must be called before relying on any implicit from freal<> to int64_t" );
    return int64_t( to_flt() );
}

template< typename T, typename FLT >              
inline freal<T,FLT>::operator uint32_t( void )
{ 
    cassert( implicit_from, "implicit_from_set( true ) must be called before relying on any implicit from freal<> to uint32_t" );
    return uint32_t( to_flt() );
}

template< typename T, typename FLT >              
inline freal<T,FLT>::operator int32_t( void )
{ 
    cassert( implicit_from, "implicit_from_set( true ) must be called before relying on any implicit from freal<> to int32_t" );
    return int32_t( to_flt() );
}

//-----------------------------------------------------
// Check that cordic(s) defined and return it.
//-----------------------------------------------------
template< typename T, typename FLT >              
inline const Cordic<T,FLT> * freal<T,FLT>::c( void ) const
{
    cassert( cordic != nullptr, "undefined type" );
    return cordic;
}

template< typename T, typename FLT >              
inline const Cordic<T,FLT> * freal<T,FLT>::c( const freal<T,FLT>& b ) const
{
    cassert( cordic   != nullptr, "a has undefined type" );
    cassert( b.cordic != nullptr, "b has undefined type" );
    cassert( cordic == b.cordic, "a and b must have same type currently" );
    return cordic;
}

template< typename T, typename FLT >              
inline const Cordic<T,FLT> * freal<T,FLT>::c( const freal<T,FLT>& b, const freal<T,FLT>& _c ) const
{
    cassert( cordic    != nullptr, "a has undefined type" );
    cassert( b.cordic  != nullptr, "b has undefined type" );
    cassert( _c.cordic != nullptr, "c has undefined type" );
    cassert( cordic == b.cordic && cordic == _c.cordic, "a and b and c must have same type currently" );
    return cordic;
}

template< typename T, typename FLT >              
inline const T * freal<T,FLT>::raw_ptr( void ) const 
{ 
    return &v; 
}

//-----------------------------------------------------
// Constants
//-----------------------------------------------------               
#define decl_const( name )                      \
    template< typename T, typename FLT >        \
    inline _freal freal<T,FLT>::name( void )    \
    { return( c(), cordic->name() ); }          \

#define decl_const1x( name, a_type )            \
    template< typename T, typename FLT >        \
    inline _freal freal<T,FLT>::name( a_type a ) \
    { return( c(), cordic->name( a ) ); }       \

decl_const( max )
decl_const( min )
decl_const( denorm_min )
decl_const( lowest )
decl_const( epsilon )
decl_const( round_error )
decl_const( zero )
decl_const( one )
decl_const( two )
decl_const( half )
decl_const( quarter )
decl_const( sqrt2 )
decl_const( sqrt2_div_2 )
decl_const( pi )
decl_const( tau )
decl_const( pi_div_2 )
decl_const( pi_div_4 )
decl_const( two_div_pi )
decl_const( four_div_pi )
decl_const( e )
decl_const1x( nan, const char * )
decl_const( quiet_nan )
decl_const( signaling_nan )
decl_const( inf )
decl_const( ninf )

//-----------------------------------------------------
// Standard Operators 
//-----------------------------------------------------               
#define decl_op1( op, name )                                    \
    template< typename T, typename FLT >                        \
    inline _freal  _freal::operator op () const                 \
    { return name();                    }                       \

#define decl_op1_code( op, code )                               \
    template< typename T, typename FLT >                        \
    inline _freal  _freal::operator op () const                 \
    { return code;                      }                       \

#define decl_op2( op, name )                                    \
    template< typename T, typename FLT >                        \
    inline _freal  _freal::operator op ( const _freal& b ) const \
    { return name( b );                 }                       \

#define decl_op2_ret( op, name, ret_type )                      \
    template< typename T, typename FLT >                        \
    inline ret_type _freal::operator op ( const _freal& b ) const \
    { return name( b );                 }                       \

#define decl_op2x( op, name, b_type )                           \
    template< typename T, typename FLT >                        \
    inline _freal  _freal::operator op ( b_type b ) const       \
    { return name( b );                 }

#define decl_op2a( op, name )                                   \
    template< typename T, typename FLT >                        \
    inline _freal& _freal::operator op ( const _freal& b )      \
    { return assign( name( b ) );       }                       \

#define decl_op2ax( op, name, b_type )                          \
    template< typename T, typename FLT >                        \
    inline _freal& _freal::operator op ( b_type b )             \
    { return assign( name( b ) );       }                       \

decl_op1(     -,        neg             )
decl_op2(     +,        add             )
decl_op2(     -,        sub             )
decl_op2(     *,        mul             )
decl_op2(     /,        div             )
decl_op2x(    <<,       scalbn,  int    )
decl_op2x(    >>,       scalbnn, int    )
decl_op2a(    =,        _freal          )
decl_op2ax(   =,        _freal,  const FLT& )
decl_op2a(    +=,       add             )
decl_op2a(    -=,       sub             )
decl_op2a(    *=,       mul             )
decl_op2a(    /=,       div             )
decl_op2ax(   <<=,      scalbn,  int    )
decl_op2ax(   >>=,      scalbnn, int    )
decl_op2_ret( >,        isgreater,      bool )
decl_op2_ret( >=,       isgreaterequal, bool )
decl_op2_ret( <,        isless,         bool )
decl_op2_ret( <=,       islessequal,    bool )
decl_op2_ret( !=,       isunequal,      bool )
decl_op2_ret( ==,       isequal,        bool )

//-----------------------------------------------------
// Well-Known Math Operators and Functions 
//
// a = *this
//-----------------------------------------------------               

template< typename T, typename FLT >              
inline freal<T,FLT>& freal<T,FLT>::assign( const freal<T,FLT>& b ) 
{ 
    cassert( b.cordic != nullptr, "assigning from undefined value b" );
    cordic = b.cordic;
    cordic->assign( v, b.v ); 
    return *this;
}

#define decl_pop1( name )                               \
    template< typename T, typename FLT >                \
    inline _freal _freal::name( void ) const            \
    { return( c(), pop_value( cordic, cordic->name( v ) ) ); } \

#define decl_pop2( name )                               \
    template< typename T, typename FLT >                \
    inline _freal _freal::name( const _freal& b ) const \
    { return( c( b ), pop_value( cordic, cordic->name( v, b.v ) ) ); } \

#define decl_popb2( name )                              \
    template< typename T, typename FLT >                \
    inline bool _freal::name( const _freal& b ) const   \
    { return( c( b ), pop_bool( cordic, cordic->name( v, b.v ) ) ); } \

#define decl_pop2p( name )                              \
    template< typename T, typename FLT >                \
    inline _freal _freal::name( _freal * b ) const      \
    { return( c(), pop_value( cordic, cordic->name( v, &b->v ) ) ); } \

#define decl_pop2x( name, b_type )                      \
    template< typename T, typename FLT >                \
    inline _freal _freal::name( b_type b ) const        \
    { return( c(), pop_value( cordic, cordic->name( v, b ) ) ); } \

#define decl_pop3( name )                               \
    template< typename T, typename FLT >                \
    inline _freal _freal::name( const _freal& b, const _freal& c ) const \
    { return( c( b, c ), pop_value( cordic, cordic->name( v, b.v, c.v ) ) ); } \

#define decl_pop3x( name, c_type )                      \
    template< typename T, typename FLT >                \
    inline _freal _freal::name( const _freal& b, c_type c ) const \
    { return( c( b ), pop_value( cordic, cordic->name( v, b.v, c ) ) ); } \

#define decl_nopop0( name, ret_type )                   \
    template< typename T, typename FLT >                \
    inline ret_type _freal::name( void ) const          \
    { return( c(), cordic->name() ); }                  \

#define decl_nopop1( name, ret_type )                   \
    template< typename T, typename FLT >                \
    inline ret_type _freal::name( void ) const          \
    { return( c(), cordic->name( v ) ); }               \

#define decl_nopop1x( name, ret_type, b_type )          \
    template< typename T, typename FLT >                \
    inline ret_type _freal::name( b_type b ) const      \
    { return( c(), cordic->name( b ) ); }               \

decl_nopop1(    signbit,        bool                    )
decl_pop2x(     frexp,          int *                   )
decl_pop2p(     modf                                    )
decl_nopop1(    ilogb,          int                     )
decl_pop1(      logb                                    )
decl_nopop1(    fpclassify,     int                     )
decl_nopop1(    isfinite,       bool                    )
decl_nopop1(    isinf,          bool                    )
decl_nopop1(    isnan,          bool                    )
decl_nopop1(    isnormal,       bool                    )
decl_nopop1x(   fesetround,     int,    int             )
decl_nopop0(    fegetround,     int                     )
decl_pop2(      nextafter                               )
decl_pop2x(     nexttoward,     long double             )
decl_pop1(      floor                                   )
decl_pop1(      ceil                                    )
decl_pop1(      trunc                                   )
decl_pop1(      round                                   )
decl_nopop1(    lround,         long                    )
decl_nopop1(    llround,        long long               )
decl_nopop1(    iround,         T                       )
decl_pop1(      rint                                    )
decl_nopop1(    lrint,          long                    )
decl_nopop1(    llrint,         long long               )
decl_nopop1(    irint,          T                       )
decl_pop1(      nearbyint                               )
decl_pop1(      abs                                     )
decl_pop1(      neg                                     )
decl_pop2(      copysign                                )
decl_pop2(      add                                     )
decl_pop2(      sub                                     )
decl_pop3(      fma                                     )
decl_pop2(      mul                                     )
decl_pop1(      sqr                                     )
decl_pop3(      fda                                     )
decl_pop2(      div                                     )
decl_pop2(      remainder                               )
decl_pop2(      fmod                                    )
decl_pop3x(     remquo,         int *                   )
decl_pop1(      rcp                                     )
decl_popb2(     isgreater                               )
decl_popb2(     isgreaterequal                          )
decl_popb2(     isless                                  )
decl_popb2(     islessequal                             )
decl_popb2(     islessgreater                           )
decl_popb2(     isunordered                             )
decl_popb2(     isunequal                               )
decl_popb2(     isequal                                 )
decl_pop2(      fdim                                    )
decl_pop2(      fmax                                    )
decl_pop2(      fmin                                    )
decl_pop1(      sqrt                                    )
decl_pop1(      rsqrt                                   )
decl_pop1(      cbrt                                    )
decl_pop1(      rcbrt                                   )
decl_pop1(      exp                                     )
decl_pop1(      expm1                                   )
decl_pop2x(     expc,           FLT                     )
decl_pop1(      exp2                                    )
decl_pop1(      exp10                                   )
decl_pop2(      pow                                     )
decl_pop1(      log                                     )
decl_pop2(      log                                     )
decl_pop1(      log1p                                   )
decl_pop2x(     logc,           FLT                     )
decl_pop1(      log2                                    )
decl_pop1(      log10                                   )
decl_pop1(      sin                                     )
decl_pop1(      sinpi                                   )
decl_pop2(      sin                                     )
decl_pop2(      sinpi                                   )
decl_pop1(      cos                                     )
decl_pop1(      cospi                                   )
decl_pop2(      cos                                     )
decl_pop2(      cospi                                   )

template< typename T, typename FLT >              
inline void freal<T,FLT>::sincos( freal<T,FLT>& si, freal<T,FLT>& co ) const                           
{ 
    T si_t, co_t;
    c()->sincos( v, si_t, co_t );
    si = pop_value( cordic, si_t );
    co = pop_value( cordic, co_t );
}

template< typename T, typename FLT >              
inline void freal<T,FLT>::sinpicospi( freal<T,FLT>& si, freal<T,FLT>& co ) const                           
{ 
    T si_t, co_t;
    c()->sinpicospi( v, si_t, co_t );
    si = pop_value( cordic, si_t );
    co = pop_value( cordic, co_t );
}

template< typename T, typename FLT >              
inline void freal<T,FLT>::sincos( freal<T,FLT>& si, freal<T,FLT>& co, const freal<T,FLT>& r ) const                           
{ 
    T si_t, co_t;
    c()->sincos( v, si_t, co_t, &r.v );
    si = pop_value( cordic, si_t );
    co = pop_value( cordic, co_t );
}

template< typename T, typename FLT >              
inline void freal<T,FLT>::sinpicospi( freal<T,FLT>& si, freal<T,FLT>& co, const freal<T,FLT>& r ) const                           
{ 
    T si_t, co_t;
    c()->sinpicospi( v, si_t, co_t, &r.v );
    si = pop_value( cordic, si_t );
    co = pop_value( cordic, co_t );
}

decl_pop1(      tan                                     )
decl_pop1(      tanpi                                   )
decl_pop1(      asin                                    )
decl_pop1(      acos                                    )
decl_pop1(      atan                                    )
decl_pop2(      atan2                                   )

template< typename T, typename FLT >              
inline void freal<T,FLT>::polar_to_rect( const freal<T,FLT>& angle, freal<T,FLT>& x, freal<T,FLT>& y ) const    
{ 
    T x_t, y_t;
    c( angle )->polar_to_rect( v, angle.v, x_t, y_t );
    x = pop_value( cordic, x_t );
    y = pop_value( cordic, y_t );
}

template< typename T, typename FLT >              
inline void freal<T,FLT>::rect_to_polar( const freal<T,FLT>& b,     freal<T,FLT>& r, freal<T,FLT>& angle ) const    
{ 
    T r_t, a_t;
    c( b )->rect_to_polar( v, b.v, r_t, a_t );
    r     = pop_value( cordic, r_t );
    angle = pop_value( cordic, a_t );
}

decl_pop2(      hypot                                   )
decl_pop2(      hypoth                                  )
decl_pop1(      sinh                                    )
decl_pop2(      sinh                                    )
decl_pop1(      cosh                                    )
decl_pop2(      cosh                                    )

template< typename T, typename FLT >              
inline void freal<T,FLT>::sinhcosh( freal<T,FLT>& sih, freal<T,FLT>& coh ) const                       
{ 
    T sih_t, coh_t;
    c()->sinhcosh( v, sih_t, coh_t );
    sih = pop_value( cordic, sih_t );
    coh = pop_value( cordic, coh_t );
}

template< typename T, typename FLT >              
inline void freal<T,FLT>::sinhcosh( freal<T,FLT>& sih, freal<T,FLT>& coh, const freal<T,FLT>& r ) const                       
{ 
    T sih_t, coh_t;
    c()->sinhcosh( v, sih_t, coh_t, &r.v );
    sih = pop_value( cordic, sih_t );
    coh = pop_value( cordic, coh_t );
}

decl_pop1(      tanh                                    )
decl_pop1(      asinh                                   )
decl_pop1(      acosh                                   )
decl_pop1(      atanh                                   )
decl_pop2(      atanh2                                  )

#endif // _freal_h
