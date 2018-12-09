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
// mpint.h - multi-precision signed integer class for C++ 
//
// This class is a very simplisitic version of arbitrary precision integers, or "big integers".
// It provides only a bare minimum set of operations needed by Cordic.h.
//
// Typical usage:
//
//     #include mpint.h
//     mpint::implicit_int_w = 128;             // change default int_w to 128 bits
//     mpint i;                                 // will get allocated 128 bits and initialized to 0
//     mpint j( 12 );                           // will get allocated 128 bits and initialized to 12
//     mpint k( 10, 56 );                       // will get allocated 56  bits and initialized to 10
//
#ifndef _mpint_h
#define _mpint_h

#include <cmath>

class mpint
{
public:
    static void implicit_int_w_set( uint32_t int_w );
    
    mpint( int64_t init=0, uint32_t int_w=0 );
    ~mpint();
    
    static mpint make_int( uint32_t int_w );

    // bare minimum set of operators
    mpint& operator =  ( const mpint& other );
    mpint  operator +  ( const mpint& other ) const;
    mpint  operator -  ( const mpint& other ) const;
    mpint  operator << ( int shift ) const;
    mpint  operator >> ( int shift ) const;

private:
    static uint32_t     implicit_int_w;
    uint32_t            int_w;
    uint32_t            word_cnt;
    union
    {
        uint64_t   w0;
        uint64_t * w;
    } u;
};

uint32_t mpint::implicit_int_w = 64;

inline mpint::mpint( int64_t init, uint32_t _int_w )
{
    int_w    = (_int_w == 0) ? implicit_int_w : _int_w;
    word_cnt = int_w / 64;
    if ( word_cnt > 1 ) {
        u.w = new uint64_t[word_cnt];

        uint64_t sign = init < 0;
        uint64_t sign_mask = uint64_t( -sign );

        for( uint32_t i = 0; i < word_cnt-1; i++ )
        {
            u.w[i] = sign_mask;
        }
        u.w[word_cnt-1] = init;
    } else {
        u.w0 = init;
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
    if ( word_cnt == 1 ) {
        u.w0 = other.u.w0;
        return *this;
    }

    for( int32_t i = (word_cnt-1); i >= 0; i-- )
    {
        u.w[i] = other.u.w[i];
    }
    return *this;
}

inline mpint mpint::operator + ( const mpint& other ) const
{
    mpint r;

    if ( word_cnt == 1 ) {
        r.u.w0 = u.w0 + other.u.w0;
        return r;
    }

    uint64_t cin = 0;
    for( int32_t i = (word_cnt-1); i >= 0; i-- )
    {
        r.u.w[i] = u.w[i] + other.u.w[i] + cin;
        cin = r.u.w[i] < u.w[i];
    }
    return r;
}

inline mpint mpint::operator - ( const mpint& other ) const
{
    mpint r;

    if ( word_cnt == 1 ) {
        r.u.w0 = u.w0 - other.u.w0;
        return r;
    }

    uint64_t cin = 0;
    for( int32_t i = (word_cnt-1); i >= 0; i-- )
    {
        r.u.w[i] = u.w[i] - other.u.w[i] - cin;
        cin = u.w[i] < other.u.w[i];
    }
    return r;
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

#endif
