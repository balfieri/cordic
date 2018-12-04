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
// freal.h - flexible real number class
//
// This class supports typical C++ floating-point semantics.
// It implements freals as a user-specified fixed-point.
// Most functions are implemented using Cordic.h.
//
// In the near future, it will allow the internal format to change
// dynamically based on input and output ranges, and the 
// format could be fixed-point or floating-point.
//
// Typical usage:
//
//     #include freal.h
//     using real = freal<int64_t, double>;
//     [use "real" in the rest of your program]
//
#ifndef _freal_h
#define _freal_h

#include "Misc.h"

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
    freal( void );                              // initializes value to undefined
    freal( const freal& other );                // copy type and value of other
    freal( const freal& other, FLT f );         // copy type of other, but value of f

    static freal make_fixed( uint32_t int_w, uint32_t frac_w, FLT init_f=FLT(0) );  // make a signed fixed-point    number with initial value
    static freal make_float( uint32_t exp_w, uint32_t frac_w, FLT init_f=FLT(0) );  // make a signed floating-point number with initial value
    ~freal();

    //-----------------------------------------------------
    // Conversions
    //-----------------------------------------------------
    //FLT    FLT( void ) const;                   // freal to FLT

    //-----------------------------------------------------
    // Standard Operators
    //-----------------------------------------------------               
    freal  operator =  ( const freal& b ) const;
    freal  operator =  ( const FLT&   b ) const;
    freal  operator +  ( const freal& b ) const;
    freal  operator -  ( const freal& b ) const;
    freal  operator *  ( const freal& b ) const;
    freal  operator /  ( const freal& b ) const;
    freal  operator << (       int    b ) const;
    freal  operator >> (       int    b ) const;

    freal& operator += ( const freal& b ) const;
    freal& operator -= ( const freal& b ) const;
    freal& operator *= ( const freal& b ) const;
    freal& operator /= ( const freal& b ) const;
    freal& operator <<=(       int    b ) const;
    freal& operator >>=(       int    b ) const;

    bool   operator == ( const freal& b ) const;
    bool   operator != ( const freal& b ) const;
    bool   operator <  ( const freal& b ) const;
    bool   operator <= ( const freal& b ) const;
    bool   operator >  ( const freal& b ) const;
    bool   operator >= ( const freal& b ) const;

    //-----------------------------------------------------
    // Well-Known Math Operators and Functions 
    //
    // a = *this
    //-----------------------------------------------------               
    freal  add( const freal& b ) const; 
    freal  sub( const freal& b ) const; 
    freal  mad( const freal& b, const freal& c ) const;             
    freal  mul( const freal& b ) const;                             
    freal  sqr( void ) const;
    freal  dad( const freal& b, const freal& c ) const;      
    freal  div( const freal& b ) const;      // a/b
    freal  one_over( void ) const;                                    
    freal  sqrt( void ) const;                                        
    freal  one_over_sqrt( void ) const;                               

    freal  exp( void ) const;                                         
    freal  pow( const freal& e ) const;      // a^e
    freal  powc( const FLT c ) const;        // c^a
    freal  pow2( void ) const;               // 2^a
    freal  pow10( void ) const;              // 10^a
    freal  log( void ) const;                                         
    freal  logb( const freal& b ) const;     // log-base-b(a)
    freal  logc( const FLT c ) const;        // log-base-c(a)                  
    freal  log2( void ) const;                                        
    freal  log10( void ) const;                                       

    freal  sin( void ) const;
    freal  cos( void ) const;
    void   sin_cos( freal& si, freal& co ) const;
    freal  tan( void ) const;                                         
    freal  asin( void ) const;                                        
    freal  acos( void ) const;                                        
    freal  atan( void ) const;                                        
    freal  atan2( const freal& b ) const;    // y=a, x=b

    void   polar_to_rect( const freal& angle, freal& x, freal& y     ) const;  // a=radius
    void   rect_to_polar( const freal& b,     freal& r, freal& angle ) const;  // x=a, y=b 
    freal  norm(  const freal& b ) const;    // x=a, y=b
    freal  normh( const freal& b ) const;    // x=a, y=b

    freal  sinh( void ) const;
    freal  cosh( void ) const;
    void   sinh_cosh( freal& sih, freal& coh, const freal * r=nullptr ) const;
    freal  tanh( void ) const;                                        
    freal  asinh( void ) const;                                       
    freal  acosh( void ) const;                                       
    freal  atanh( void ) const;                                       
    freal  atanh2( const freal& b ) const;   // atanh2( a, b )

private:
};

// std:xxx() functions should pick these up 
//
template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  add( const freal<T,FLT>& a, const freal<T,FLT>& b ); 

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  sub( const freal<T,FLT>& a, const freal<T,FLT>& b ); 

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  mad( const freal<T,FLT>& a, const freal<T,FLT>& b, const freal<T,FLT>& c );             

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  mul( const freal<T,FLT>& a, const freal<T,FLT>& b );                             

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  sqr( const freal<T,FLT>& a );

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  dad( const freal<T,FLT>& a, const freal<T,FLT>& b, const freal<T,FLT>& c );      

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  div( const freal<T,FLT>& a, const freal<T,FLT>& b );      // a/b

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  one_over( const freal<T,FLT>& a );                                    

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  sqrt( const freal<T,FLT>& a );                                        

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  one_over_sqrt( const freal<T,FLT>& a );

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  exp( const freal<T,FLT>& a );

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  pow( const freal<T,FLT>& a, const freal<T,FLT>& e );      // a^e

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  powc( const freal<T,FLT>& a, const FLT c );  // c^a

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  pow2( const freal<T,FLT>& a );               // 2^a

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  pow10( const freal<T,FLT>& a );              // 10^a

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  log( const freal<T,FLT>& a );                                         

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  logb( const freal<T,FLT>& a, const freal<T,FLT>& b );     // log-base-b(a)

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  logc( const freal<T,FLT>& a, const FLT c );        // log-base-c(a)                  

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  log2( const freal<T,FLT>& a );                                        

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  log10( const freal<T,FLT>& a );                                       

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  sin( const freal<T,FLT>& a );

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  cos( const freal<T,FLT>& a );

template< typename T=int64_t, typename FLT=double >              
static void   sin_cos( const freal<T,FLT>& a, freal<T,FLT>& si, freal<T,FLT>& co );

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  tan( const freal<T,FLT>& a );                                         

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  asin( const freal<T,FLT>& a );                                        

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  acos( const freal<T,FLT>& a );                                        

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  atan( const freal<T,FLT>& a );                                        

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  atan2( const freal<T,FLT>& a, const freal<T,FLT>& b );    // y=a, x=b

template< typename T=int64_t, typename FLT=double >              
static void   polar_to_rect( const freal<T,FLT>& a, const freal<T,FLT>& angle, freal<T,FLT>& x, freal<T,FLT>& y     );  // a=radius

template< typename T=int64_t, typename FLT=double >              
static void   rect_to_polar( const freal<T,FLT>& a, const freal<T,FLT>& b,     freal<T,FLT>& r, freal<T,FLT>& angle );  // x=a, y=b 

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  norm( const freal<T,FLT>& a,  const freal<T,FLT>& b );    // x=a, y=b

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  normh( const freal<T,FLT>& a, const freal<T,FLT>& b );    // x=a, y=b

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  sinh( const freal<T,FLT>& a );

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  cosh( const freal<T,FLT>& a );

template< typename T=int64_t, typename FLT=double >              
static void   sinh_cosh( const freal<T,FLT>& a, freal<T,FLT>& sih, freal<T,FLT>& coh, const freal<T,FLT> * r=nullptr );

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  tanh( const freal<T,FLT>& a );                                        

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  asinh( const freal<T,FLT>& a );                                       

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  acosh( const freal<T,FLT>& a );                                       

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  atanh( const freal<T,FLT>& a );                                       

template< typename T=int64_t, typename FLT=double >              
static freal<T,FLT>  atanh2( const freal<T,FLT>& a, const freal<T,FLT>& b );                          

#endif
