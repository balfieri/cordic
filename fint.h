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
// fint.h - flexible signed integer class
//
// This class is a very simplisitic version arbitrary precision integers, or "big integers".
// It provides a bare minimum set of operations needed by Cordic.h.
//
// Typical usage:
//
//     #include fint.h
//     using myint = fint<128>;                 // myint is a flexible integer with 128 bits
//
#ifndef _fint_h
#define _fint_h

#include <ctypes>

template< uint32_t int_w >
class fint
{
    static constexpr word_cnt = (2*int_w - 1) / 64;
    
    fint( void );
    fint( int64_t other );
    ~fint();

    // bare minimum set of operators
    fint& operator =  ( const fint& other );
    fint  operator +  ( const fint& other ) const;
    fint  operator -  ( const fint& other ) const;
    fint  operator << ( int shift ) const;
    fint  operator >> ( int shift ) const;

private:
    uint64_t    w[word_cnt];
};

template< uint32_t int_w >
inline fint::fint( void )
{
}

template< uint32_t int_w >
inline fint::fint( int64_t other )
{
    uint64_t sign = other < 0;
    uint64_t sign_mask = uint64_t( -sign );

    for( uint32_t i = 0; i < word_cnt-1; i++ )
    {
        w[i] = sign_mask;
    }
    w[word_cnt-1] = other;
}

template< uint32_t int_w >
inline fint::~fint()
{
}

template< uint32_t int_w >
inline fint& fint::operator = ( const fint& other )
{
    for( int32_t i = (word_cnt-1); i >= 0; i-- )
    {
        w[i] = other.w[i];
    }
    return *this;
}

template< uint32_t int_w >
inline fint fint::operator + ( const fint& other )
{
    fint r;
    uint64_t cin = 0;
    for( int32_t i = (word_cnt-1); i >= 0; i-- )
    {
        r.w[i] = w[i] + other.w[i] + cin;
        cin = r.w[i] < w[i];
    }
    return r;
}

template< uint32_t int_w >
inline fint fint::operator - ( int shift )
{
    fint r;
    uint64_t cin = 0;
    for( int32_t i = (word_cnt-1); i >= 0; i-- )
    {
        r.w[i] = w[i] - other.w[i] - cin;
        cin = w[i] < other.w[i];
    }
    return r;
}

template< uint32_t int_w >
inline fint fint::operator << ( int shift )
{
    if ( shift == 0 ) return *this;

    if ( shift < 0 ) return *this >> -shift;
    
    fint r;
    if ( word_cnt == 1 ) {
        // easy
        r.w[0] = w[0] << shift;
        return r;
    }

    for( uint32_t i = 0; i < word_cnt; i++ )
    {
        r.w[i] = 0;
    }

    for( uint32_t tb = 0; tb < word_cnt*8; tb++ )
    {
        uint32_t fb  = tb + shift;
        uint32_t fw  = fb / 64;
        uint32_t fwb = fb % 64;
        bool     b   = (fb >= word_cnt*8) ? 0 : ((w[fw] >> (63-fwb)) & 1);

        uint32_t tw  = tb / 64;
        uint32_t twb = tb % 64;
        r.w[tw] |= b << (63-twb);
    }

    return r;
}

template< uint32_t int_w >
inline fint fint::operator >> ( int shift )
{
    if ( shift == 0 ) return *this;

    if ( shift < 0 ) return *this << -shift;
    
    uint64_t sign = (w[0] >> 63) & 1;

    fint r;
    if ( word_cnt == 1 ) {
        // easy
        r.w[0] = uint64_t( -sign ) << (64-shift);
        r.w[0] |= w[0] >> shift;
        return r;
    }

    for( uint32_t i = 0; i < word_cnt; i++ )
    {
        r.w[i] = 0;
    }

    for( int32_t tb = word_cnt*8-1; tb >= 0; tb-- )
    {
        uint32_t fb  = tb - shift;
        uint32_t fw  = fb / 64;
        uint32_t fwb = fb % 64;
        bool     b   = (fb < 0) ? sign : ((w[fw] >> fwb) & 1);

        uint32_t tw  = tb / 64;
        uint32_t twb = tb % 64;
        r.w[tw] |= b << twb;
    }

    return r;
}
