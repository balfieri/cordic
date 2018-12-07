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
    static freal make_raw(   const Cordic<T,FLT> * cordic, const T& encoded );      // use encoded as the raw bits (dangerous)
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
    static void implicit_to_set( const Cordic<T,FLT> * cordic );
    static void implicit_from_set( bool allow );

    freal( float f );
    freal( double f );
    freal( uint64_t i );
    freal( int64_t i );
    freal( uint32_t i );
    freal( int32_t i );

    operator float( void );
    operator double( void );
    operator uint64_t( void );
    operator int64_t( void );
    operator uint32_t( void );
    operator int32_t( void );

    //-----------------------------------------------------
    // Constants
    //-----------------------------------------------------               
    freal max( void );                                       // maximum positive value 
    freal min( void );                                       // minimum positive value
    T     maxint( void );                                       // largest positive integer (just integer part, does not include fraction)
    freal zero( void );                                         // 0.0
    freal one( void );                                          // 1.0
    freal half( void );                                         // 0.5
    freal quarter( void );                                      // 0.25
    freal sqrt2( void );                                        // sqrt(2)
    freal sqrt2_div_2( void );                                  // sqrt(2)/2
    freal pi( void );                                           // PI
    freal pi_div_2( void );                                     // PI/2
    freal pi_div_4( void );                                     // PI/4
    freal two_div_pi( void );                                   // 2/PI
    freal four_div_pi( void );                                  // 4/PI
    freal e( void );                                            // natural exponent

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

    bool   isgreater( const freal& b ) const;                        
    bool   isgreaterequal( const freal& b ) const;                 
    bool   isless( const freal& b ) const;                          
    bool   islessequal( const freal& b ) const;                   
    bool   islessgreater( const freal& b ) const;                
    bool   isunordered( const freal& b ) const;              
    bool   isunequal( const freal& b ) const;                 
    bool   isequal( const freal& b ) const;                   

    freal  abs( void ) const;
    freal  neg( void ) const; 
    freal  floor( void ) const;
    freal  ceil( void ) const;

    freal  add( const freal& b ) const; 
    freal  sub( const freal& b ) const; 
    freal  mad( const freal& b, const freal& c ) const;             
    freal  fma( const freal& b, const freal& c ) const;         // same as mad()
    freal  mul( const freal& b ) const;                             
    freal  lshift( int b ) const;
    freal  rshift( int b ) const;
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
    freal  hypot( const freal& b ) const;    // same as norm()
    freal  normh( const freal& b ) const;    // x=a, y=b

    freal  sinh( void ) const;
    freal  cosh( void ) const;
    void   sinh_cosh( freal& sih, freal& coh ) const;
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

private:
    static const Cordic<T,FLT> * implicit_to;
    static bool                  implicit_from;

    const Cordic<T,FLT> *       cordic;         // defines the type and most operations
    T                           v;              // this value encoded in type T
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

template< typename T=int64_t, typename FLT=double >              
static inline std::string   to_string( const freal<T,FLT>& a )
{ return a.to_string();                 }

template< typename T=int64_t, typename FLT=double >              
static inline bool          isgreater( const freal<T,FLT>& a, const freal<T,FLT>& b )
{ return a.isgreater( b );              }

template< typename T=int64_t, typename FLT=double >              
static inline bool          isgreaterequal( const freal<T,FLT>& a, const freal<T,FLT>& b )
{ return a.isgreaterequal( b );         }

template< typename T=int64_t, typename FLT=double >              
static inline bool          isless( const freal<T,FLT>& a, const freal<T,FLT>& b )
{ return a.isless( b );                 }

template< typename T=int64_t, typename FLT=double >              
static inline bool          islessequal( const freal<T,FLT>& a, const freal<T,FLT>& b )
{ return a.islessequal( b );            }

template< typename T=int64_t, typename FLT=double >              
static inline bool          islessgreater( const freal<T,FLT>& a, const freal<T,FLT>& b )
{ return a.islessgreater( b );          }

template< typename T=int64_t, typename FLT=double >              
static inline bool          isunordered( const freal<T,FLT>& a, const freal<T,FLT>& b )
{ return a.isunordered( b );            }

template< typename T=int64_t, typename FLT=double >              
static inline bool          isunequal( const freal<T,FLT>& a, const freal<T,FLT>& b )
{ return a.isunequal( b );              }

template< typename T=int64_t, typename FLT=double >              
static inline bool          isequal( const freal<T,FLT>& a, const freal<T,FLT>& b )
{ return a.isequal( b );                }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  abs( const freal<T,FLT>& a )
{ return a.abs();                       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  neg( const freal<T,FLT>& a )
{ return a.neg();                       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  floor( const freal<T,FLT>& a )
{ return a.floor();                     }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  ceil( const freal<T,FLT>& a )
{ return a.ceil();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  add( const freal<T,FLT>& a, const freal<T,FLT>& b )                         
{ return a.add( b );                    }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  sub( const freal<T,FLT>& a, const freal<T,FLT>& b )                         
{ return a.sub( b );                    }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  mad( const freal<T,FLT>& a, const freal<T,FLT>& b, const freal<T,FLT>& c )  
{ return a.mad( b, c );                 }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  fma( const freal<T,FLT>& a, const freal<T,FLT>& b, const freal<T,FLT>& c )  
{ return a.fma( b, c );                 }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  mul( const freal<T,FLT>& a, const freal<T,FLT>& b )                         
{ return a.mul( b );                    }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  sqr( const freal<T,FLT>& a )                                                
{ return a.sqr();                       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  dad( const freal<T,FLT>& a, const freal<T,FLT>& b, const freal<T,FLT>& c )  
{ return a.dad( b, c );                 }      

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  div( const freal<T,FLT>& a, const freal<T,FLT>& b )                         
{ return a.div( b );                    }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  one_over( const freal<T,FLT>& a )                                           
{ return a.one_over();                  }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  sqrt( const freal<T,FLT>& a )                                               
{ return a.sqrt();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  one_over_sqrt( const freal<T,FLT>& a )                                      
{ return a.one_over_sqrt();             }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  exp( const freal<T,FLT>& a )                                                
{ return a.exp();                       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  pow( const freal<T,FLT>& a, const freal<T,FLT>& e )                         
{ return a.pow( e );                    }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  powc( const freal<T,FLT>& a, const FLT c )                                  
{ return a.powc( c );                   }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  pow2( const freal<T,FLT>& a )                                               
{ return a.pow2();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  pow10( const freal<T,FLT>& a )                                              
{ return a.pow10();                     }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  log( const freal<T,FLT>& a )                                                
{ return a.log();                       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  logb( const freal<T,FLT>& a, const freal<T,FLT>& b )                        
{ return a.logb( b );                   }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  logc( const freal<T,FLT>& a, const FLT c )                                  
{ return a.logc( c );                   }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  log2( const freal<T,FLT>& a )                                               
{ return a.log2();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  log10( const freal<T,FLT>& a )                                              
{ return a.log10();                     }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  sin( const freal<T,FLT>& a )                                                
{ return a.sin();                       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  cos( const freal<T,FLT>& a )                                                
{ return a.cos();                       }

template< typename T=int64_t, typename FLT=double >              
static inline void   sin_cos( const freal<T,FLT>& a, freal<T,FLT>& si, freal<T,FLT>& co )               
{ a.sin_cos( si, co );                  }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  tan( const freal<T,FLT>& a )                                                
{ return a.tan();                       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  asin( const freal<T,FLT>& a )                                               
{ return a.asin();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  acos( const freal<T,FLT>& a )                                               
{ return a.acos();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  atan( const freal<T,FLT>& a )                                               
{ return a.atan();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  atan2( const freal<T,FLT>& a, const freal<T,FLT>& b )                       
{ return a.atan2( b );                  }

template< typename T=int64_t, typename FLT=double >              
static inline void   polar_to_rect( const freal<T,FLT>& a, const freal<T,FLT>& angle, freal<T,FLT>& x, freal<T,FLT>& y )  
{ a.polar_to_rect( angle, x, y );       }

template< typename T=int64_t, typename FLT=double >              
static inline void   rect_to_polar( const freal<T,FLT>& a, const freal<T,FLT>& b,     freal<T,FLT>& r, freal<T,FLT>& angle )  
{ a.rect_to_polar( b, r, angle );       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  norm( const freal<T,FLT>& a,  const freal<T,FLT>& b )                       
{ return a.norm( b );                   }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  hypot( const freal<T,FLT>& a, const freal<T,FLT>& b )                       
{ return a.norm( b );                   }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  normh( const freal<T,FLT>& a, const freal<T,FLT>& b )                       
{ return a.normh( b );                  }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  sinh( const freal<T,FLT>& a )                                               
{ return a.sinh();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  cosh( const freal<T,FLT>& a )                                               
{ return a.cosh();                      }

template< typename T=int64_t, typename FLT=double >              
static inline void   sinh_cosh( const freal<T,FLT>& a, freal<T,FLT>& sih, freal<T,FLT>& coh )           
{ a.sinh_cosh( sih, coh );              }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  tanh( const freal<T,FLT>& a )                                               
{ return a.tanh();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  asinh( const freal<T,FLT>& a )                                              
{ return a.asinh();                     }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  acosh( const freal<T,FLT>& a )                                              
{ return a.acosh();                     }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  atanh( const freal<T,FLT>& a )                                              
{ return a.atanh();                     }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  atanh2( const freal<T,FLT>& a, const freal<T,FLT>& b )                      
{ return a.atanh2( b );                 }

}

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
inline bool                  freal<T,FLT>::implicit_from = false;      // disallow

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
inline freal<T,FLT>::freal( const freal& other )
{
    cordic = other.c();
    v      = other.v;
}

template< typename T, typename FLT >              
inline freal<T,FLT>::freal( const freal& other, FLT f )
{
    cordic = other.c();
    v      = cordic->to_t( f );
}

template< typename T, typename FLT >              
inline freal<T,FLT>::freal( const Cordic<T,FLT> * _cordic, FLT f )
{
    cassert( _cordic != nullptr, "freal(cordic, f) cordic argument must be non-null" );
    cordic = _cordic;
    v      = cordic->to_t( f );
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
inline freal<T,FLT> freal<T,FLT>::make_raw( const Cordic<T,FLT> * cordic, const T& encoded )
{
    cassert( cordic != nullptr, "make_raw(cordic, encoded) called with null cordic" );
    freal r;
    r.cordic = cordic;
    r.v      = encoded;
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT>::~freal()
{
    cordic = nullptr;
    v      = T(668);
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
}

template< typename T, typename FLT >              
inline void freal<T,FLT>::implicit_from_set( bool allow )
{ 
    implicit_from = allow;
}

template< typename T, typename FLT >              
inline freal<T,FLT>::freal( double f )
{
    cassert( implicit_to != nullptr, "implicit_to_set() must be called before relying on any implicit from double to freal<>" );
    *this = freal( implicit_to, FLT(f) );
}

template< typename T, typename FLT >              
inline freal<T,FLT>::freal( float f )
{
    cassert( implicit_to != nullptr, "implicit_to_set() must be called before relying on any implicit from float to freal<>" );
    *this = freal( implicit_to, FLT(f) );
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
inline freal<T,FLT>::operator double( void )
{ 
    cassert( implicit_from, "implicit_from_set( true ) must be called before relying on any implicit from freal<> to double" );
    return to_flt();
}

template< typename T, typename FLT >              
inline freal<T,FLT>::operator float( void )
{ 
    cassert( implicit_from, "implicit_from_set( true ) must be called before relying on any implicit from freal<> to float" );
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

//-----------------------------------------------------
// Constants
//-----------------------------------------------------               
template< typename T, typename FLT >              
inline T            freal<T,FLT>::maxint( void ) 
{ return( c(), cordic->maxint() );      } 

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::max( void ) 
{ return( c(), make_raw( cordic, cordic->max() ) ); }

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::min( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->min();
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::zero( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->zero();
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::one( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->one();
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::half( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->half();
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::quarter( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->quarter();
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::sqrt2( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->sqrt2();
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::sqrt2_div_2( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->sqrt2_div_2();
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::pi( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->pi();
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::pi_div_2( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->pi_div_2();
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::pi_div_4( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->pi_div_4();
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::two_div_pi( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->two_div_pi();
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::four_div_pi( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->four_div_pi();
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::e( void ) 
{ 
    freal r;
    r.cordic = c();
    r.v      = cordic->e();
    return r;
}

//-----------------------------------------------------
// Standard Operators 
//-----------------------------------------------------               
template< typename T, typename FLT >              
inline freal<T,FLT>  freal<T,FLT>::operator -  () const 
{ return neg();                         }

template< typename T, typename FLT >              
inline freal<T,FLT>  freal<T,FLT>::operator +  ( const freal<T,FLT>& b ) const                          
{ return add( b );                      }

template< typename T, typename FLT >              
inline freal<T,FLT>  freal<T,FLT>::operator -  ( const freal<T,FLT>& b ) const                          
{ return sub( b );                      }

template< typename T, typename FLT >              
inline freal<T,FLT>  freal<T,FLT>::operator *  ( const freal<T,FLT>& b ) const                          
{ return mul( b );                      }

template< typename T, typename FLT >              
inline freal<T,FLT>  freal<T,FLT>::operator /  ( const freal<T,FLT>& b ) const                          
{ return div( b );                      }

template< typename T, typename FLT >              
inline freal<T,FLT>  freal<T,FLT>::operator << (       int    b ) const                                 
{ return lshift( b );                   }

template< typename T, typename FLT >              
inline freal<T,FLT>  freal<T,FLT>::operator >> (       int    b ) const                                 
{ return rshift( b );                   }

template< typename T, typename FLT >              
inline freal<T,FLT>& freal<T,FLT>::operator =  ( const freal<T,FLT>& b ) 
{ return assign( b );                   }

template< typename T, typename FLT >              
inline freal<T,FLT>& freal<T,FLT>::operator =  ( const FLT&   b )
{ return assign( freal<T,FLT>( b ) );   }

template< typename T, typename FLT >              
inline freal<T,FLT>& freal<T,FLT>::operator += ( const freal<T,FLT>& b )
{ return assign( add( b ) );            }

template< typename T, typename FLT >              
inline freal<T,FLT>& freal<T,FLT>::operator -= ( const freal<T,FLT>& b )
{ return assign( sub( b ) );            }

template< typename T, typename FLT >              
inline freal<T,FLT>& freal<T,FLT>::operator *= ( const freal<T,FLT>& b )
{ return assign( mul( b ) );            }

template< typename T, typename FLT >              
inline freal<T,FLT>& freal<T,FLT>::operator /= ( const freal<T,FLT>& b )
{ return assign( div( b ) );            }

template< typename T, typename FLT >              
inline freal<T,FLT>& freal<T,FLT>::operator <<=(       int    b )
{ return assign( lshift( b ) );         }

template< typename T, typename FLT >              
inline freal<T,FLT>& freal<T,FLT>::operator >>=(       int    b )
{ return assign( rshift( b ) );         }

template< typename T, typename FLT >             
inline bool   freal<T,FLT>::operator >  ( const freal<T,FLT>& b ) const                                        
{ return isgreater( b );                }

template< typename T, typename FLT >              
inline bool   freal<T,FLT>::operator >= ( const freal<T,FLT>& b ) const                                        
{ return isgreaterequal( b );           }

template< typename T, typename FLT >              
inline bool   freal<T,FLT>::operator <  ( const freal<T,FLT>& b ) const                                        
{ return isless( b );                   }

template< typename T, typename FLT >              
inline bool   freal<T,FLT>::operator <= ( const freal<T,FLT>& b ) const                                        
{ return islessequal( b );              }

template< typename T, typename FLT >              
inline bool   freal<T,FLT>::operator != ( const freal<T,FLT>& b ) const                                        
{ return isunequal( b );                }

template< typename T, typename FLT >              
inline bool   freal<T,FLT>::operator == ( const freal<T,FLT>& b ) const                                        
{ return isequal( b );                  }

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

template< typename T, typename FLT >              
inline bool freal<T,FLT>::isgreater( const freal<T,FLT>& b ) const                                    
{ return c( b )->isgreater( v, b.v );   }

template< typename T, typename FLT >              
inline bool freal<T,FLT>::isgreaterequal( const freal<T,FLT>& b ) const                                    
{ return c( b )->isgreaterequal( v, b.v );   }

template< typename T, typename FLT >              
inline bool freal<T,FLT>::isless( const freal<T,FLT>& b ) const                                    
{ return c( b )->isless( v, b.v );   }

template< typename T, typename FLT >              
inline bool freal<T,FLT>::islessequal( const freal<T,FLT>& b ) const                                    
{ return c( b )->islessequal( v, b.v );   }

template< typename T, typename FLT >              
inline bool freal<T,FLT>::islessgreater( const freal<T,FLT>& b ) const                                    
{ return c( b )->islessgreater( v, b.v );   }

template< typename T, typename FLT >              
inline bool freal<T,FLT>::isunordered( const freal<T,FLT>& b ) const                                    
{ return c( b )->isunordered( v, b.v );   }

template< typename T, typename FLT >              
inline bool freal<T,FLT>::isunequal( const freal<T,FLT>& b ) const                                    
{ return c( b )->isunequal( v, b.v );   }

template< typename T, typename FLT >              
inline bool freal<T,FLT>::isequal( const freal<T,FLT>& b ) const                                    
{ return c( b )->isequal( v, b.v );   }

template< typename T, typename FLT >              
inline freal<T,FLT>  freal<T,FLT>::abs( void ) const
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->abs( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT>  freal<T,FLT>::neg( void ) const
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->neg( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT>  freal<T,FLT>::floor( void ) const
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->floor( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT>  freal<T,FLT>::ceil( void ) const
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->ceil( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::add( const freal<T,FLT>& b ) const                                    
{ 
    freal r;
    r.cordic = c( b );
    r.v      = r.cordic->add( v, b.v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::sub( const freal<T,FLT>& b ) const                                    
{ 
    freal r;
    r.cordic = c( b );
    r.v      = r.cordic->sub( v, b.v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::mad( const freal<T,FLT>& b, const freal<T,FLT>& c ) const             
{ 
    freal r;
    r.cordic = c( b, c );
    r.v      = r.cordic->mad( v, b.v, c.v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::mul( const freal<T,FLT>& b ) const                                    
{ 
    freal r;
    r.cordic = c( b );
    r.v      = r.cordic->mul( v, b.v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::sqr( void ) const                                                     
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->sqr( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::dad( const freal<T,FLT>& b, const freal<T,FLT>& c ) const             
{ 
    freal r;
    r.cordic = c( b, c );
    r.v      = r.cordic->dad( v, b.v, c.v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::div( const freal<T,FLT>& b ) const                                    
{ 
    freal r;
    r.cordic = c( b );
    r.v      = r.cordic->div( v, b.v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::one_over( void ) const                                                
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->one_over( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::sqrt( void ) const                                                    
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->sqrt( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::one_over_sqrt( void ) const                                           
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->one_over_sqrt( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::exp( void ) const                                                     
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->exp( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::pow( const freal<T,FLT>& e ) const                                    
{ 
    freal r;
    r.cordic = c( e );
    r.v      = r.cordic->pow( v, e.v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::powc( const FLT c ) const                                             
{ 
    freal r;
    r.cordic = c( c );
    r.v      = r.cordic->powc( v, c.v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::pow2( void ) const                                                    
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->pow2( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::pow10( void ) const                                                   
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->pow10( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::log( void ) const                                                     
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->log( v );
    return r;
}

template< typename T, typename FLT >             
inline freal<T,FLT> freal<T,FLT>::logb( const freal<T,FLT>& b ) const                                   
{ 
    freal r;
    r.cordic = c( b );
    r.v      = r.cordic->logb( b.v, v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::logc( const FLT c ) const                                             
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->logc( c, v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::log2( void ) const                                                    
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->log2( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::log10( void ) const                                                   
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->log10( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::sin( void ) const                                                     
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->sin( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::cos( void ) const                                                     
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->cos( v );
    return r;
}

template< typename T, typename FLT >              
inline void freal<T,FLT>::sin_cos( freal<T,FLT>& si, freal<T,FLT>& co ) const                           
{ 
    si.cordic = cordic; co.cordic = cordic;
    c()->sin_cos( v, si.v, co.v );                
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::tan( void ) const                                                     
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->tan( v );
    return r;
}

template< typename T, typename FLT >                                                     
inline freal<T,FLT> freal<T,FLT>::asin( void ) const                                                    
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->asin( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::acos( void ) const                                                    
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->acos( v );
    return r;
}

template< typename T, typename FLT >                                     
inline freal<T,FLT> freal<T,FLT>::atan( void ) const                                                    
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->atan( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::atan2( const freal<T,FLT>& b ) const                                  
{ 
    freal r;
    r.cordic = c( b );
    r.v      = r.cordic->atan2( v, b.v );
    return r;
}

template< typename T, typename FLT >              
inline void freal<T,FLT>::polar_to_rect( const freal<T,FLT>& angle, freal<T,FLT>& x, freal<T,FLT>& y ) const    
{ 
    x.cordic = cordic; y.cordic = cordic;
    c( angle )->polar_to_rect( v, x.v, y.v );     
}

template< typename T, typename FLT >              
inline void freal<T,FLT>::rect_to_polar( const freal<T,FLT>& b,     freal<T,FLT>& r, freal<T,FLT>& angle ) const    
{ 
    r.cordic = cordic; angle.cordic = cordic;     
    c( b )->rect_to_polar( v, b.v, r.v, angle.v );
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::norm(  const freal<T,FLT>& b ) const                                  
{ 
    freal r;
    r.cordic = c( b );
    r.v      = r.cordic->norm( v, b.v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::hypot(  const freal<T,FLT>& b ) const                                  
{ 
    freal r;
    r.cordic = c( b );
    r.v      = r.cordic->hypot( v, b.v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::normh( const freal<T,FLT>& b ) const                                  
{ 
    freal r;
    r.cordic = c( b );
    r.v      = r.cordic->normh( v, b.v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::sinh( void ) const                                                    
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->sinh( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::cosh( void ) const                                                    
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->cosh( v );
    return r;
}

template< typename T, typename FLT >              
inline void freal<T,FLT>::sinh_cosh( freal<T,FLT>& sih, freal<T,FLT>& coh ) const                       
{ 
    sih.cordic = cordic; coh.cordic = cordic;
    c()->sinh_cosh( v, sih.v, coh.v );            
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::tanh( void ) const                                                    
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->tanh( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::asinh( void ) const                                                   
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->asinh( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::acosh( void ) const                                                   
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->acosh( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::atanh( void ) const                                                   
{ 
    freal r;
    r.cordic = c();
    r.v      = r.cordic->atanh( v );
    return r;
}

template< typename T, typename FLT >              
inline freal<T,FLT> freal<T,FLT>::atanh2( const freal<T,FLT>& b ) const                                 
{ 
    freal r;
    r.cordic = c( b );
    r.v      = r.cordic->atanh2( v, b.v );
    return r;
}

#endif // _freal_h
