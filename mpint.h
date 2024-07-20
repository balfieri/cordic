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
// This class is a very simplisitic version of arbitrary precision integers, or "big integers."
// It provides only a bare minimum set of operations needed by Cordic.h.
// I didn't use an official class because I wanted to make sure Cordic wasn't
// using any integer math operations besides adding and shifting.
//
// Typical usage:
//
//     #include mpint.h
//     mpint::implicit_int_w_set( 128 );        // change default int_w to 128 bits (from 64)
//     mpint i;                                 // will get allocated 128 bits and initialized to 0
//     mpint j( 12 );                           // will get allocated 128 bits and initialized to 12
//     mpint k( 12, 56 );                       // will get allocated 56  bits and initialized to 12
//
#ifndef _mpint_h
#define _mpint_h

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>

class mpint
{
public:
    static void implicit_int_w_set( size_t int_w );
    
    mpint( void );
    mpint( int64_t i );
    mpint( int64_t i, size_t int_w );
    mpint( const mpint& b );
    ~mpint();
    
    // minimum set of operators needed by Cordic.h:
    bool   signbit     ( void ) const;
    mpint  neg         ( void ) const;
    mpint  operator -  () const;
    mpint& operator =  ( const mpint& b );
    mpint  operator +  ( const mpint& b ) const;
    mpint  operator -  ( const mpint& b ) const;
    mpint  operator << ( int shift ) const;
    mpint  operator >> ( int shift ) const;

    bool   operator >  ( const mpint& b ) const;
    bool   operator >= ( const mpint& b ) const;
    bool   operator <  ( const mpint& b ) const;
    bool   operator <= ( const mpint& b ) const;
    bool   operator != ( const mpint& b ) const;
    bool   operator == ( const mpint& b ) const;

    mpint& operator +=  ( const mpint& b );
    mpint& operator -=  ( const mpint& b );
    mpint& operator <<= ( int shift );
    mpint& operator >>= ( int shift );

    static mpint to_mpint( std::string, bool allow_no_conversion=false, int base=10, size_t * pos=nullptr );  
    std::string  to_string( int base=10, int width=0 ) const;                

private:
    static size_t     implicit_int_w;
    size_t            int_w;
    size_t            word_cnt;
    union
    {
        uint64_t   w0;          // if fits in 64 bits
        uint64_t * w;           // if doesn't fit in 64 bits
    } u;

    bool bit( size_t i ) const;         // returns bit i
    void fixsign( void );               // re-extend the sign after possible overflow
    int  compare( const mpint& b ) const; // -1 is <, 0 is ==, 1 is >
};

// Well-Known std:xxx() Functions 
//
namespace std
{

static inline bool signbit( const mpint& a ) 
{ 
    return a.signbit();
}

template< typename T=int64_t, typename FLT=double >              
static inline std::string to_string( const mpint& a, int base=10, int width=0 )
{ 
    return a.to_string( base, width );
}

static inline mpint stoi( std::string str, size_t * pos=nullptr, int base=10 )
{ 
    return mpint::to_mpint( str, true, base, pos );
}

static inline std::istream& operator >> ( std::istream &in, mpint& a )
{ 
    int base = 10;              // need to query base 
    std::string s = "";
    in >> std::ws;          // eat up whitespace
    for( bool is_first = true; ; is_first = false )
    {
        int c = in.peek();
        if ( (!is_first && c == '-') || (c < '0' && c > '9') ) break;
        in >> c; // consume it
        s += c;
    }
    a = mpint::to_mpint( s, true, base ); // quietly produce 0 for bad input
    return in;
}

static inline std::ostream& operator << ( std::ostream &out, const mpint& a )
{
    int base = 10;              // need to query base and width
    int width = 0;
    out << a.to_string( base, width );       
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

size_t mpint::implicit_int_w = 64;

inline void mpint::implicit_int_w_set( size_t int_w )
{
    implicit_int_w = int_w;
}

inline mpint::mpint( void )
{
    // mark it undefined
    int_w    = 0;
    word_cnt = 0;
    u.w      = nullptr;
}

inline mpint::mpint( int64_t init, size_t _int_w )
{
    int_w = _int_w;
    iassert( int_w > 0, "int_w must be > 0" );
    word_cnt = (int_w + 63) / 64;
    if ( word_cnt == 1 ) {
        u.w0 = init;
    } else {
        u.w = new uint64_t[word_cnt];

        uint64_t sign = init < 0;
        uint64_t sign_mask = uint64_t( -sign );

        u.w[0] = init;
        for( size_t i = 1; i < word_cnt; i++ )
        {
            u.w[i] = sign_mask;
        }
    }
}

inline mpint::mpint( int64_t init ) : mpint( init, implicit_int_w )
{
}

inline mpint::mpint( const mpint& b ) : mpint() 
{
    *this = b;
}

inline mpint::~mpint()
{
    if ( word_cnt > 1 ) {
        delete u.w;
        u.w = nullptr;
    }
}

inline bool mpint::bit( size_t i ) const
{
    iassert( int_w > 0, "mpint is undefined" );
    iassert( i < int_w, "mpint bit i is out of range" );

    if ( word_cnt == 1 ) {
        return (u.w0 >> i) & 1;
    } else {
        size_t w = i / 64;
        size_t b = i % 64;
        return (u.w[w] >> b) & 1;
    }
}

inline bool mpint::signbit( void ) const
{
    return bit( int_w-1 );
}

inline mpint mpint::to_mpint( std::string s, bool allow_no_conversion, int base, size_t * pos )
{
    iassert( base == 10, "to_mpint() currently supports only base 10" );

    //--------------------------------------------------------------
    // Do this conversion without doing any multiplies.
    // When we have to multiply by 10 we simply do shifts and adds.
    // Then add in the new digit. [This can be extended to b bases up
    // to the allowed base of 36.]
    //--------------------------------------------------------------
    mpint r( 0, implicit_int_w+1 );  // to hold the largest negative integer (as positive integer)
    bool is_neg = false;
    bool got_digit = false;
    size_t len = s.length();
    size_t i;
    for( i = 0; i < len; i++ )
    {
        char c = s[i];
        if ( !is_neg && !got_digit && (c == ' ' || c == '\t' || c == '\n') ) continue; // skip whitespace

        if ( c == '-' ) {
            if ( is_neg || got_digit ) {
                iassert( allow_no_conversion, "to_mpint '-' is not allowed after first sign or digit" );
                break;
            }
            is_neg = true;
        } else if ( c >= '0' && c <= '9' ) {
            mpint r1 = r << 3;
            mpint r2 = r << 1;
            mpint r3( c - '0' );
            r = r1 + r2 + r3;
            iassert( !r.signbit(), "to_mpint string does not fit: " + s );
            got_digit = true;
        } else {
            break;
        }
    }
    iassert( got_digit || allow_no_conversion, "to_mpint did not find any digits in '" + s + "'" ); 
    if ( pos != nullptr ) *pos = is_neg ? (i - 1) : i;
    if ( is_neg ) r = r.neg();
    mpint rr( 0 );
    rr = r;    // will cause it to truncate
    iassert( rr.signbit() == r.signbit(), "to_mpint string does not fit: " + s );
    return rr;
}

inline std::string mpint::to_string( int base, int width ) const
{
    iassert( int_w > 0, "to_string: this mpint is undefined" );
    iassert( base >= 0 && base <= 36, "base must be between 0 and 36" );
    if ( base == 0 ) base = 10;
   
    std::string s;
    if ( base == 2 ) {
        //--------------------------------------------------------------
        // Fast path for base-2.
        //--------------------------------------------------------------
        s = "";
        for( size_t i = 0; i < int_w; i++ )
        {
            char b = bit( i ) ? '1' : '0';
            s = b + s;
        }
    } else {
        //--------------------------------------------------------------
        // General Path - like elementary school addition
        //
        // Maintain a power-of-2 as a character string in
        // the proper base.  The string starts off as "1" and 
        // gets wider every time we multiply by 2 
        // by adding the current power-of-2 to itself (no need to multiply).
        //
        // If mpint bit i is set, then we add the power of two into a similar
        // character string that starts out as 0.  We calculate
        // the next power-of-2 at the same time using a similar add.
        //--------------------------------------------------------------
        bool is_neg = signbit();
        mpint a( 0, int_w+1 );  // extra bit to deal with most negative integer
        a = *this;
        if ( is_neg ) a = -a;
        s = "0";
        std::string pow2 = "1";
        for( size_t i = 0; i < (a.int_w-1); i++ )
        {
            std::string  cpow2     = pow2;  // current pow-of-2
            const char * cpow2_c   = cpow2.c_str();
            size_t       cpow2_len = cpow2.length();
            std::string  cs        = s;     // current s
            const char * cs_c      = cs.c_str();
            size_t       cs_len    = s.length();

            s = "";                         // gonna re-create these
            pow2 = "";
            uint32_t cins = 0;              // carry ins
            uint32_t cin2 = 0;
            bool     a_bit = a.bit( i );    // don't add pow2 to s if !a_bit
            for( size_t j = 0; j <= cpow2_len; j++ )
            {
                size_t   k2  = cpow2_len-1 - j;
                size_t   ks  = cs_len-1 - j;
                uint32_t dc2 = (j >= cpow2_len) ? 0 : ((cpow2_c[k2] <= '9') ? (cpow2_c[k2] - '0') : (cpow2_c[k2] - 'a'));
                uint32_t dcs = (j >= cs_len)    ? 0 : ((cs_c[ks]    <= '9') ? (cs_c[ks]    - '0') : (cs_c[ks]    - 'a'));
                uint32_t ds  = (a_bit ? dc2 : 0) + dcs + cins;  // will end up simply with no change in s if !a_bit
                uint32_t d2  = dc2 + dc2 + cin2;                // must always update pow2
                cins = ds / base;  
                ds   = ds % base;
                cin2 = d2 / base;
                d2   = d2 % base;
                char chs = (ds <= 9) ? ('0' + ds) : ('a' + ds);
                char ch2 = (d2 <= 9) ? ('0' + d2) : ('a' + d2);
                s    = chs + s;
                pow2 = ch2 + pow2;
            }
        }
        while( s[0] == '0' ) 
        {
            s = s.substr( 1 );
        }
        if ( is_neg ) s = "-" + s;
    }
    while( size_t(width) > s.length() )
    {
        s = " " + s;  // pad
    }
    return s;
}

inline mpint& mpint::operator = ( const mpint& b )
{
    iassert( b.int_w > 0, "rhs int_w must be > 0" );
    if ( int_w == 0 ) {
        // inherit b's int_w
        int_w    = b.int_w;
        word_cnt = b.word_cnt;
        if ( word_cnt > 1 ) u.w = new uint64_t[word_cnt];
    }

    if ( word_cnt == 1 ) {
        u.w0 = (b.word_cnt == 1) ? b.u.w0 : b.u.w[0];
    } else {
        bool b_sign = b.signbit();
        for( size_t i = 0; i < word_cnt; i++ )
        {
            u.w[i] = (i < b.word_cnt) ? ((b.word_cnt == 1) ? b.u.w0 : b.u.w[i]) : (b_sign ? uint64_t(-1) : 0);
        }
    }

    if ( int_w != b.int_w ) fixsign();

    return *this;
}

inline mpint mpint::neg() const
{
    // negate = 2's complement
    iassert( int_w > 0, "trying to negate in undefined mpint" );
    mpint r( 0, int_w );
    if ( word_cnt == 1 ) {
        r.u.w0 = ~u.w0 + 1;
    } else {
        int64_t cin = 1;
        for( size_t i = 0; i < word_cnt; i++ )
        {
            r.u.w[i] = ~u.w[i] + cin;
            cin = u.w[i] == 0 && cin;
        }
    }
    return r;
}

inline mpint mpint::operator -() const
{
    return neg();
}

inline void mpint::fixsign( void ) 
{
    //-------------------------------------------------------
    // If int_w is not an integral multiple of word_cnt, then
    // we need to re-extend the new sign bit in the top word.
    //-------------------------------------------------------
    size_t sign_pos = (int_w-1) % 64;  // in the top word
    if ( sign_pos != 63 ) {
        bool       sign      = signbit();
        uint64_t   sign_mask = uint64_t(-1) << sign_pos;
        uint64_t * word_ptr = (word_cnt == 1) ? &u.w0 : &u.w[word_cnt-1];
        if ( sign ) {
            *word_ptr |= sign_mask;             // propagate 1
        } else {
            *word_ptr &= ~sign_mask;            // propagate 0
        }
    }
}

inline mpint mpint::operator + ( const mpint& b ) const
{
    // pick larger of the two for result
    mpint r( 0, (int_w > b.int_w) ? int_w : b.int_w );
    
    if ( r.word_cnt == 1 ) {
        r.u.w0 = u.w0 + b.u.w0;
    } else {
        uint64_t cin = 0;
        for( size_t i = 0; i < r.word_cnt; i++ ) 
        {
            uint64_t wt = (i < word_cnt)   ? ((word_cnt > 1)   ? u.w[i]   : u.w0)   : 0;
            uint64_t wo = (i < b.word_cnt) ? ((b.word_cnt > 1) ? b.u.w[i] : b.u.w0) : 0;
            r.u.w[i] = wt + wo + cin;
            cin = r.u.w[i] < wt || r.u.w[i] < wo;
        }
    }
    
    r.fixsign();  // after possible overflow corruption of sign bits
    return r;
}

inline mpint mpint::operator - ( const mpint& b ) const
{
    return *this + b.neg();
}

inline mpint mpint::operator << ( int shift ) const
{
    if ( shift == 0 ) return *this;

    if ( shift < 0 ) return *this >> -shift;
    
    mpint r( 0, int_w );
    if ( word_cnt == 1 ) {
        r.u.w0 = u.w0 << shift;
        return r;
    }

    for( size_t tb = 0; tb < int_w; tb++ )
    {
        int64_t  fb  = int64_t(tb) - shift;
        int64_t  fw  = fb / 64;
        int64_t  fwb = fb % 64;
        uint64_t b   = (fb < 0) ? 0 : ((u.w[fw] >> fwb) & 1);

        size_t   tw  = tb / 64;
        size_t   twb = tb % 64;
        r.u.w[tw] |= b << twb;
    }

    return r;
}

inline mpint mpint::operator >> ( int shift ) const
{
    if ( shift == 0 ) return *this;

    if ( shift < 0 ) return *this << -shift;
    
    bool sign = signbit();

    mpint r( 0, int_w );
    if ( word_cnt == 1 ) {
        // easy
        r.u.w0 = uint64_t( -sign ) << (64-shift);
        r.u.w0 |= u.w0 >> shift;
        return r;
    }

    for( size_t tb = 0; tb < int_w; tb++ )
    {
        size_t   fb  = tb + shift;
        size_t   fw  = fb / 64;
        size_t   fwb = fb % 64;
        uint64_t b   = (fb >= int_w) ? sign : ((u.w[fw] >> fwb) & 1);

        size_t   tw  = tb / 64;
        size_t   twb = tb % 64;
        r.u.w[tw] |= b << twb;
    }

    return r;
}

inline mpint& mpint::operator +=  ( const mpint& b ) 
{ 
    *this = *this + b;      
    return *this;
}

inline mpint& mpint::operator -=  ( const mpint& b ) 
{ 
    *this = *this - b;      
    return *this;
}

inline mpint& mpint::operator <<= ( int shift )      
{ 
    *this = *this << shift; 
    return *this;
}

inline mpint& mpint::operator >>= ( int shift )      
{ 
    *this = *this >> shift; 
    return *this;
}

int mpint::compare( const mpint& b ) const
{
    int a_sign = signbit()   ? -1 : 1;
    int b_sign = b.signbit() ? -1 : 1;
    if ( a_sign != b_sign ) {
        // easy case
        return a_sign;
    }

    // need to start looking at words
    // be careful about different numbers of words
    size_t cnt = (word_cnt > b.word_cnt) ? word_cnt : b.word_cnt;
    for( size_t i = 0; i < cnt; i++ )
    {
        size_t   k = cnt-i;
        uint64_t wa = (k < word_cnt)   ? ((word_cnt   == 1) ?   u.w0 :   u.w[k]) : (a_sign ? uint64_t(-1) : 0);
        uint64_t wb = (k < b.word_cnt) ? ((b.word_cnt == 1) ? b.u.w0 : b.u.w[k]) : (b_sign ? uint64_t(-1) : 0);
        if ( wa != wb ) {
            return (wa < wb) ? b_sign : a_sign;
        }
    }
    return 0;  // equal
}

bool mpint::operator >  ( const mpint& b ) const
{
    return compare( b ) > 0;
}

bool mpint::operator >= ( const mpint& b ) const
{
    return compare( b ) >= 0;
}

bool mpint::operator <  ( const mpint& b ) const
{
    return compare( b ) < 0;
}

bool mpint::operator <= ( const mpint& b ) const
{
    return compare( b ) <= 0;
}

bool mpint::operator != ( const mpint& b ) const
{
    return compare( b ) != 0;
}

bool mpint::operator == ( const mpint& b ) const
{
    return compare( b ) == 0;
}

#endif
