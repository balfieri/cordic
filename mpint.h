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
// mpint.h - very limited multi-precision signed integer class for C++ 
//
// This class is a very simplisitic version of arbitrary precision integers, or "big integers".
// It provides only a bare minimum set of operations needed by Cordic.h.
//
// Typical usage:
//
//     #include mpint.h
//     mpint::implicit_int_w_set( 128 );        // change default int_w to 128 bits (from 64)
//     mpint i;                                 // will get allocated 128 bits and initialized to 0
//     mpint j( 12 );                           // will get allocated 128 bits and initialized to 12
//     mpint k( 10, 56 );                       // will get allocated 56  bits and initialized to 10
//
#ifndef _mpint_h
#define _mpint_h

#include <cmath>
#include <iostream>
#include <string>

class mpint
{
public:
    static void implicit_int_w_set( uint32_t int_w );
    
    mpint( void );
    mpint( int64_t i, uint32_t int_w=0 );
    ~mpint();
    
    // minimum set of operators needed by Cordic.h:
    mpint& operator =  ( const mpint& other );
    mpint  operator -  () const;
    mpint  operator +  ( const mpint& other ) const;
    mpint  operator -  ( const mpint& other ) const;
    mpint  operator << ( int shift ) const;
    mpint  operator >> ( int shift ) const;

    static mpint to_mpint( std::string );
    std::string  to_string( void ) const;                

private:
    static uint32_t     implicit_int_w;
    uint32_t            int_w;
    uint32_t            word_cnt;
    union
    {
        uint64_t   w0;          // if fits in 64 bits
        uint64_t * w;           // if doesn't fit in 64 bits
    } u;
};

// Well-Known std:xxx() Functions 
//
namespace std
{

template< typename T=int64_t, typename FLT=double >              
static inline std::string to_string( const mpint& a ) 
{ 
    return a.to_string();
}

template< typename T=int64_t, typename FLT=double >              
static inline std::istream& operator >> ( std::istream &in, mpint& a )
{ 
    std::string s = "";
    // TODO: collect characters
    a = mpint::to_mpint( s );
    return in;
}

template< typename T=int64_t, typename FLT=double >              
static inline std::ostream& operator << ( std::ostream &out, const mpint& a )
{ 
    out << a.to_string(); 
    return out;     
}

}

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

#define iassert(expr, msg) if ( !(expr) ) \
                { std::cout << "ERROR: assertion failure: " << (msg) << " at " << __FILE__ << ":" << __LINE__ << "\n"; exit( 1 ); }

uint32_t mpint::implicit_int_w = 64;

inline void mpint::implicit_int_w_set( uint32_t int_w )
{
    implicit_int_w = int_w;
}

inline mpint::mpint( void )
{
    // mark it undefined
    int_w    = 0;
    word_cnt = 0;
}

inline mpint::mpint( int64_t init, uint32_t _int_w )
{
    int_w = (_int_w == 0) ? implicit_int_w : _int_w;
    iassert( int_w > 0, "int_w must be > 0" );
    word_cnt = (int_w + 63) / 64;
    if ( word_cnt == 1 ) {
        u.w0 = init;
    } else {
        u.w = new uint64_t[word_cnt];

        uint64_t sign = init < 0;
        uint64_t sign_mask = uint64_t( -sign );

        for( uint32_t i = 0; i < word_cnt-1; i++ )
        {
            u.w[i] = sign_mask;
        }
        u.w[word_cnt-1] = init;
    }
}

inline mpint::~mpint()
{
    if ( word_cnt > 1 ) {
        delete u.w;
        u.w = nullptr;
    }
}

inline mpint& mpint::operator = ( const mpint& other )
{
    iassert( other.int_w > 0, "rhs int_w must be > 0" );
    if ( int_w == 0 ) {
        int_w    = other.int_w;
        word_cnt = other.word_cnt;
        if ( word_cnt > 1 ) u.w = new uint64_t[word_cnt];
    }

    if ( word_cnt == 1 ) {
        // need to truncate and sign-extend
        u.w0 = other.u.w0;
    } else {
        // need to truncate and sign-extend
        for( int32_t i = (word_cnt-1); i >= 0; i-- ) u.w[i] = other.u.w[i];
    }

    return *this;
}

inline mpint mpint::operator - () const
{
    // negate = 2's complement
    iassert( int_w > 0, "trying to negate in undefined mpint" );
    mpint r = *this;
    if ( word_cnt == 1 ) {
        r.u.w0 = ~u.w0 + 1;
    } else {
        int64_t cin = 1;
        for( uint32_t i = 0; i < word_cnt; i++ )
        {
            r.u.w[i] = ~u.w[i] + cin;
            cin = r.u.w[i] < u.w[i];
        }
    }
    return r;
}

inline mpint mpint::operator + ( const mpint& other ) const
{
    mpint r;

    // pick larger of the two for result
    r.int_w    = (int_w > other.int_w) ? int_w : other.int_w;
    r.word_cnt = (int_w + 63) / 64;
    if ( r.word_cnt > 1 ) r.u.w = new uint64_t[r.word_cnt];
    
    if ( r.word_cnt == 1 ) {
        r.u.w0 = u.w0 + other.u.w0;
    } else {
        uint64_t cin = 0;
        for( uint32_t i = 0; i < word_cnt; i++ ) 
        {
            r.u.w[i] = u.w[i] + other.u.w[i] + cin;
            cin = r.u.w[i] < u.w[i];
        }
        // TODO: fix any overflow
    }
    return r;
}

inline mpint mpint::operator - ( const mpint& other ) const
{
    return *this + -other;
}

inline mpint mpint::operator << ( int shift ) const
{
    if ( shift == 0 ) return *this;

    if ( shift < 0 ) return *this >> -shift;
    
    mpint r;
    if ( word_cnt == 1 ) {
        r.u.w0 = u.w0 << shift;
        return r;
    }

    for( uint32_t i = 0; i < word_cnt; i++ )
    {
        r.u.w[i] = 0;
    }

    for( uint32_t tb = 0; tb < word_cnt*8; tb++ )
    {
        uint32_t fb  = tb + shift;
        uint32_t fw  = fb / 64;
        uint32_t fwb = fb % 64;
        bool     b   = (fb >= word_cnt*8) ? 0 : ((u.w[fw] >> (63-fwb)) & 1);

        uint32_t tw  = tb / 64;
        uint32_t twb = tb % 64;
        r.u.w[tw] |= b << (63-twb);
    }

    return r;
}

inline mpint mpint::operator >> ( int shift ) const
{
    if ( shift == 0 ) return *this;

    if ( shift < 0 ) return *this << -shift;
    
    uint64_t sign;

    mpint r;
    if ( word_cnt == 1 ) {
        // easy
        sign = u.w0 < 1;
        r.u.w0 = uint64_t( -sign ) << (64-shift);
        r.u.w0 |= u.w0 >> shift;
        return r;
    }

    sign = (u.w[0] >> 63) & 1;

    for( uint32_t i = 0; i < word_cnt; i++ )
    {
        r.u.w[i] = 0;
    }

    for( int32_t tb = word_cnt*8-1; tb >= 0; tb-- )
    {
        uint32_t fb  = tb - shift;
        uint32_t fw  = fb / 64;
        uint32_t fwb = fb % 64;
        bool     b   = (fb < 0) ? sign : ((u.w[fw] >> fwb) & 1);

        uint32_t tw  = tb / 64;
        uint32_t twb = tb % 64;
        r.u.w[tw] |= b << twb;
    }

    return r;
}

inline mpint mpint::to_mpint( std::string s )
{
    mpint r( 0 );

    if ( r.word_cnt == 1 ) {
        r.u.w0 = std::atoi( s.c_str() );
    } else {
        // TODO: there must be code somewhere to do this
    }

    return r;
}

inline std::string mpint::to_string( void ) const
{
    if ( word_cnt == 1 ) {
        return std::to_string( int64_t(u.w0) );
    } else {
        // TODO: convert to decimal
        std::string s = "";
        uint32_t k = 0;
        for( uint32_t i = 0; i < word_cnt; i++ )
        {
            for( uint32_t j = 0; j < 64 && k < int_w; j++, k++ )
            {
                char b = ((u.w[i] >> j) & 1) ? '1' : '0';
                s = b + s;
            }
        }
        s = "0b" + s;
        return s;
    }
}

#endif
