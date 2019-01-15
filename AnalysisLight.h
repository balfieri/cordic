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
// AnalysisLight.h - class intended to just count op totals and be very fast
//
#ifndef _AnalysisLight_h
#define _AnalysisLight_h

#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "Cordic.h"
#include "Logger.h"
#include "Analysis.h"

template< typename T=int64_t, typename FLT=double >
class AnalysisLight : public Logger<T,FLT>
{
public:
    AnalysisLight( std::string base_name = "log" );     
    ~AnalysisLight();

    // Call this if you want to be thread-safe
    virtual void tid_set( uint32_t t );      

    // Logger Overrides
    //
    virtual void cordic_constructed( const void * cordic, uint32_t int_exp_w, uint32_t frac_w, 
                                     bool is_float, uint32_t guard_w, uint32_t n );
    virtual void cordic_destructed(  const void * cordic );

    virtual void enter( const char * name );
    virtual void leave( const char * name );

    virtual void constructed( const T * v, const void * cordic );
    virtual void destructed(  const T * v, const void * cordic );

    virtual void op(  uint16_t op, uint32_t opnd_cnt, const T * opnd[] );
    virtual void op1( uint16_t op, const T *  opnd1 );
    virtual void op1( uint16_t op, const T&   opnd1 );
    virtual void op1( uint16_t op, const bool opnd1 );
    virtual void op1( uint16_t op, const FLT& opnd1 );
    virtual void op2( uint16_t op, const T *  opnd1, const T * opnd2 );
    virtual void op2( uint16_t op, const T *  opnd1, const T&  opnd2 );
    virtual void op2( uint16_t op, const T *  opnd1, const FLT&opnd2 );
    virtual void op3( uint16_t op, const T *  opnd1, const T * opnd2, const T * opnd3 );
    virtual void op4( uint16_t op, const T *  opnd1, const T * opnd2, const T * opnd3, const T * opnd4 );

    using OP                         = typename Cordic<T,FLT>::OP;
    static constexpr uint64_t OP_cnt = Cordic<T,FLT>::OP_cnt;

    virtual void inc_op_cnt( OP op, uint32_t by=1 );

    virtual void parse( void );
    virtual void clear_stats( void );
    virtual void print_stats( std::string basename="", double scale_factor=1.0, const std::vector<std::string>& ignore_funcs=std::vector<std::string>() ) const;

private:
    std::string         base_name;

    static constexpr uint32_t INT_W_MAX = 32;           
    static constexpr uint32_t THREAD_CNT_MAX = 64;
    uint64_t op_cnt[THREAD_CNT_MAX][OP_cnt];                            // keep totals for each thread
};

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

template< typename T, typename FLT >
AnalysisLight<T,FLT>::AnalysisLight( std::string _base_name ) : Logger<T,FLT>( Cordic<T,FLT>::op_to_str )
{
    base_name = _base_name;

    clear_stats();
}

template< typename T, typename FLT >
AnalysisLight<T,FLT>::~AnalysisLight()
{
}

template< typename T, typename FLT >
void AnalysisLight<T,FLT>::tid_set( uint32_t t )
{
    cassert( t < THREAD_CNT_MAX, "tid_set: tid must be < THREAD_CNT_MAX" );
    tid = t;
}

//-----------------------------------------------------
// Logger Method Overrides
//-----------------------------------------------------
template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::cordic_constructed( const void * cordic_ptr, uint32_t int_exp_w, uint32_t frac_w, 
                                          bool is_float, uint32_t guard_w, uint32_t n )
{
    (void)cordic_ptr;
    (void)int_exp_w;
    (void)frac_w;
    (void)is_float;
    (void)guard_w;
    (void)n;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::cordic_destructed( const void * cordic_ptr )
{
    (void)cordic_ptr;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::enter( const char * _name )
{
    (void)_name;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::leave( const char * _name )
{
    (void)_name;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::constructed( const T * v, const void * cordic_ptr )
{
    (void)v;
    (void)cordic_ptr;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::destructed(  const T * v, const void * cordic_ptr )
{
    (void)v;
    (void)cordic_ptr;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op( uint16_t _op, uint32_t opnd_cnt, const T * opnd[] )
{
    (void)opnd_cnt;
    (void)opnd;
    op_cnt[tid][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op1( uint16_t _op, const T * opnd1 )
{
    (void)opnd1;
    op_cnt[tid][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op1( uint16_t _op, const T& opnd1 )
{
    (void)opnd1;
    op_cnt[tid][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op1( uint16_t _op, bool opnd1 )
{
    (void)opnd1;
    op_cnt[tid][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op1( uint16_t _op, const FLT& opnd1 )
{
    (void)opnd1;
    op_cnt[tid][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op2( uint16_t _op, const T * opnd1, const T * opnd2 )
{
    (void)opnd1;
    (void)opnd2;
    op_cnt[tid][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op2( uint16_t _op, const T * opnd1, const T& opnd2 )
{
    (void)opnd1;
    (void)opnd2;
    op_cnt[tid][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op2( uint16_t _op, const T * opnd1, const FLT& opnd2 ) 
{
    (void)opnd1;
    (void)opnd2;
    op_cnt[tid][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op3( uint16_t _op, const T * opnd1, const T * opnd2, const T * opnd3 )
{
    (void)opnd1;
    (void)opnd2;
    (void)opnd3;
    op_cnt[tid][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op4( uint16_t _op, const T * opnd1, const T * opnd2, const T * opnd3, const T * opnd4 )
{
    (void)opnd1;
    (void)opnd2;
    (void)opnd3;
    (void)opnd4;
    op_cnt[tid][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::inc_op_cnt( OP _op, uint32_t by )
{
    op_cnt[tid][uint16_t(_op)] += by;
}

//-----------------------------------------------------
// We don't parse.
//-----------------------------------------------------
template< typename T, typename FLT >
void AnalysisLight<T,FLT>::parse( void )
{
    std::cout << "ERROR: AnalysisLight doesn't parse\n";
    exit( 1 );
}

template< typename T, typename FLT >
void AnalysisLight<T,FLT>::clear_stats( void )
{
    for( uint32_t t = 0; t < THREAD_CNT_MAX; t++ )
    {
        for( uint32_t o = 0; o < Cordic<T,FLT>::OP_cnt; o++ )
        {
            op_cnt[t][o] = 0;
        }
    }
}

template< typename T, typename FLT >
void AnalysisLight<T,FLT>::print_stats( std::string basename, double scale_factor, const std::vector<std::string>& ignore_funcs ) const
{
    (void)ignore_funcs;

    std::string out_name = basename + ".out";
    FILE * out = fopen( out_name.c_str(), "w" );
    std::ofstream csv( basename + ".csv", std::ofstream::out );

    fprintf( out, "\n\nOP Totals:\n" );
    csv << "\n\n\"Totals:\"" << "\n";
    for( uint32_t i = 0; i < OP_cnt; i++ )
    {
        uint64_t cnt = 0;
        for( uint32_t t = 0; t < THREAD_CNT_MAX; t++ )
        {
            cnt += op_cnt[t][i];
        }
        if ( cnt == 0 ) continue;

        uint64_t scaled_cnt = double(cnt) * scale_factor + 0.5;
        fprintf( out, "    %-40s:  %10lld   %10lld\n", Cordic<T,FLT>::op_to_str( i ).c_str(), cnt, scaled_cnt );
        csv << "\"" << Cordic<T,FLT>::op_to_str( i ) << "\", " << cnt << ", " << scaled_cnt << "\n";
    }

    fclose( out );
    csv.close();
    std::cout << "\nWrote stats to " + basename + ".{out,csv}\n";
}

#endif
