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
    freal( void );                                      // initializes value to undefined
    freal( const freal& other );                        // copy type and value of other
    freal( const freal& other, FLT f );                 // copy type of other, but value of f
    freal( const Cordic<T,FLT> * cordic, FLT f );       // use type from cordic, but value of f

    static freal make_fixed( uint32_t int_w, uint32_t frac_w, FLT init_f=FLT(0) );  // make a signed fixed-point    number with initial value
    static freal make_float( uint32_t exp_w, uint32_t frac_w, FLT init_f=FLT(0) );  // make a signed floating-point number with initial value
    ~freal();

    //-----------------------------------------------------
    // Conversions
    //-----------------------------------------------------
    FLT    to_flt( void ) const;                        // freal to FLT

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
    void   sinh_cosh( freal& sih, freal& coh ) const;
    freal  tanh( void ) const;                                        
    freal  asinh( void ) const;                                       
    freal  acosh( void ) const;                                       
    freal  atanh( void ) const;                                       
    freal  atanh2( const freal& b ) const;   // atanh2( a, b )

private:
    const Cordic<T,FLT> *       cordic;         // defines the type and most operations
    T                           a;              // this value encoded in type T
};

// std:xxx() calls should pick these up automatically
//
template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  add( const freal<T,FLT>& a, const freal<T,FLT>& b )                         { return a.add( b );                    }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  sub( const freal<T,FLT>& a, const freal<T,FLT>& b )                         { return a.sub( b );                    }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  mad( const freal<T,FLT>& a, const freal<T,FLT>& b, const freal<T,FLT>& c )  { return a.mad( b, c );                 }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  mul( const freal<T,FLT>& a, const freal<T,FLT>& b )                         { return a.mul( b );                    }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  sqr( const freal<T,FLT>& a )                                                { return a.sqr();                       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  dad( const freal<T,FLT>& a, const freal<T,FLT>& b, const freal<T,FLT>& c )  { return a.dad( b, c );                 }      

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  div( const freal<T,FLT>& a, const freal<T,FLT>& b )                         { return a.div( b );                    }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  one_over( const freal<T,FLT>& a )                                           { return a.one_over();                  }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  sqrt( const freal<T,FLT>& a )                                               { return a.sqrt();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  one_over_sqrt( const freal<T,FLT>& a )                                      { return a.one_over_sqrt();             }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  exp( const freal<T,FLT>& a )                                                { return a.exp();                       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  pow( const freal<T,FLT>& a, const freal<T,FLT>& e )                         { return a.pow( e );                    }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  powc( const freal<T,FLT>& a, const FLT c )                                  { return a.powc( c );                   }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  pow2( const freal<T,FLT>& a )                                               { return a.pow2();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  pow10( const freal<T,FLT>& a )                                              { return a.pow10();                     }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  log( const freal<T,FLT>& a )                                                { return a.log();                       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  logb( const freal<T,FLT>& a, const freal<T,FLT>& b )                        { return a.logb( b );                   }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  logc( const freal<T,FLT>& a, const FLT c )                                  { return a.logc( c );                   }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  log2( const freal<T,FLT>& a )                                               { return a.log2();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  log10( const freal<T,FLT>& a )                                              { return a.log10();                     }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  sin( const freal<T,FLT>& a )                                                { return a.sin();                       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  cos( const freal<T,FLT>& a )                                                { return a.cos();                       }

template< typename T=int64_t, typename FLT=double >              
static inline void   sin_cos( const freal<T,FLT>& a, freal<T,FLT>& si, freal<T,FLT>& co )               { a.sin_cos( si, co );                  }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  tan( const freal<T,FLT>& a )                                                { return a.tan();                       }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  asin( const freal<T,FLT>& a )                                               { return a.asin();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  acos( const freal<T,FLT>& a )                                               { return a.acos();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  atan( const freal<T,FLT>& a )                                               { return a.atan();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  atan2( const freal<T,FLT>& a, const freal<T,FLT>& b )                       { return a.atan2( b );                  }

template< typename T=int64_t, typename FLT=double >              
static inline void   polar_to_rect( const freal<T,FLT>& a, const freal<T,FLT>& angle, freal<T,FLT>& x, freal<T,FLT>& y )  
                                                                                                        { a.polar_to_rect( angle, x, y );       }
template< typename T=int64_t, typename FLT=double >              
static inline void   rect_to_polar( const freal<T,FLT>& a, const freal<T,FLT>& b,     freal<T,FLT>& r, freal<T,FLT>& angle )  
                                                                                                        { rect_to_polar( b, r, angle );         }
template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  norm( const freal<T,FLT>& a,  const freal<T,FLT>& b )                       { return a.norm( b );                   }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  normh( const freal<T,FLT>& a, const freal<T,FLT>& b )                       { return a.normh( b );                  }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  sinh( const freal<T,FLT>& a )                                               { return a.sinh();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  cosh( const freal<T,FLT>& a )                                               { return a.cosh();                      }

template< typename T=int64_t, typename FLT=double >              
static inline void   sinh_cosh( const freal<T,FLT>& a, freal<T,FLT>& sih, freal<T,FLT>& coh )           { a.sinh_cosh( sih, coh );              }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  tanh( const freal<T,FLT>& a )                                               { return a.tanh();                      }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  asinh( const freal<T,FLT>& a )                                              { return a.asinh();                     }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  acosh( const freal<T,FLT>& a )                                              { return a.acosh();                     }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  atanh( const freal<T,FLT>& a )                                              { return a.atanh();                     }

template< typename T=int64_t, typename FLT=double >              
static inline freal<T,FLT>  atanh2( const freal<T,FLT>& a, const freal<T,FLT>& b )                      { return a.atanh2( b );                 }

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

#endif // _freal_h
