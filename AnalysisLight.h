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
#ifdef issubnormal
#undef issubnormal
#endif
#include <cstdint>
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

    virtual void enter( uint16_t func_id );
    virtual void leave( uint16_t func_id );

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
    virtual void print_stats( std::string basename, double scale_factor,
                              const std::vector<std::string>& func_names, const std::vector<uint16_t>& ignore_funcs=std::vector<uint16_t>() ) const;

private:
    std::string         base_name;

    static constexpr uint32_t INT_W_MAX = 32;           
    static constexpr uint32_t FUNC_CNT_MAX = 64;
    static constexpr uint32_t THREAD_CNT_MAX = 64;
    static constexpr uint32_t STACK_CNT_MAX = 1024;
    uint64_t                  op_cnt[THREAD_CNT_MAX][FUNC_CNT_MAX][OP_cnt];             // keep totals for each thread-function
    uint16_t                  stack[THREAD_CNT_MAX][STACK_CNT_MAX];                     // func call stack
    uint32_t                  stack_cnt[THREAD_CNT_MAX];                                // func call stack depth

    void                stack_push( uint16_t func_id );
    uint16_t            stack_top( void );
    void                stack_pop( void );
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

    for( uint32_t t = 0; t < THREAD_CNT_MAX; t++ )
    {
        stack_cnt[t] = 0;
    }

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
inline void AnalysisLight<T,FLT>::stack_push( uint16_t func_id )
{
    cassert( stack_cnt[tid] < STACK_CNT_MAX, "depth of call stack exceeded" );
    stack[tid][stack_cnt[tid]++] = func_id;
}

template< typename T, typename FLT >
inline uint16_t AnalysisLight<T,FLT>::stack_top( void )
{
    cassert( stack_cnt[tid] > 0, "can't get top of an empty call stack" );
    return stack[tid][stack_cnt[tid]-1];
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::stack_pop( void )
{
    cassert( stack_cnt[tid] > 0, "can't pop an empty call stack" );
    stack_cnt[tid]--;
}

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
inline void AnalysisLight<T,FLT>::enter( uint16_t func_id )
{
    cassert( func_id < FUNC_CNT_MAX, "func_id is too large" );
    stack_push( func_id );
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::leave( uint16_t func_id )
{
    cassert( stack_top() == func_id , "trying to leave a routine that's not at the top of the stack: entered " + 
                                      std::to_string(stack_top()) + " leaving " + std::to_string(func_id) );
    stack_pop();
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
    op_cnt[tid][stack_top()][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op1( uint16_t _op, const T * opnd1 )
{
    (void)opnd1;
    op_cnt[tid][stack_top()][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op1( uint16_t _op, const T& opnd1 )
{
    (void)opnd1;
    op_cnt[tid][stack_top()][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op1( uint16_t _op, bool opnd1 )
{
    (void)opnd1;
    op_cnt[tid][stack_top()][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op1( uint16_t _op, const FLT& opnd1 )
{
    (void)opnd1;
    op_cnt[tid][stack_top()][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op2( uint16_t _op, const T * opnd1, const T * opnd2 )
{
    (void)opnd1;
    (void)opnd2;
    op_cnt[tid][stack_top()][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op2( uint16_t _op, const T * opnd1, const T& opnd2 )
{
    (void)opnd1;
    (void)opnd2;
    op_cnt[tid][stack_top()][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op2( uint16_t _op, const T * opnd1, const FLT& opnd2 ) 
{
    (void)opnd1;
    (void)opnd2;
    op_cnt[tid][stack_top()][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op3( uint16_t _op, const T * opnd1, const T * opnd2, const T * opnd3 )
{
    (void)opnd1;
    (void)opnd2;
    (void)opnd3;
    op_cnt[tid][stack_top()][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::op4( uint16_t _op, const T * opnd1, const T * opnd2, const T * opnd3, const T * opnd4 )
{
    (void)opnd1;
    (void)opnd2;
    (void)opnd3;
    (void)opnd4;
    op_cnt[tid][stack_top()][_op]++;
}

template< typename T, typename FLT >
inline void AnalysisLight<T,FLT>::inc_op_cnt( OP _op, uint32_t by )
{
    op_cnt[tid][stack_top()][uint16_t(_op)] += by;
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
        for( uint32_t f = 0; f < FUNC_CNT_MAX; f++ )
        {
            for( uint32_t o = 0; o < Cordic<T,FLT>::OP_cnt; o++ )
            {
                op_cnt[t][f][o] = 0;
            }
        }
    }
}

template< typename T, typename FLT >
void AnalysisLight<T,FLT>::print_stats( std::string basename, double scale_factor, 
                                        const std::vector<std::string>& func_names, const std::vector<uint16_t>& ignore_funcs ) const
{
    std::string out_name = basename + ".out";
    FILE * out = fopen( out_name.c_str(), "w" );
    std::ofstream csv( basename + ".csv", std::ofstream::out );

    uint64_t total_op_cnt[OP_cnt];
    for( uint32_t i = 0; i < OP_cnt; i++ )
    {
        total_op_cnt[i] = 0;
    }
    for( uint32_t f = 0; f < FUNC_CNT_MAX; f++ )
    {
        bool ignored = false;
        for( uint32_t i = 0; !ignored && i < ignore_funcs.size(); i++ )
        {
            ignored = f == ignore_funcs[i];
        }
        if ( ignored ) continue;

        bool have_any = false;
        for( uint32_t i = 0; !have_any && i < OP_cnt; i++ )
        {
            for( uint32_t t = 0; t < THREAD_CNT_MAX; t++ )
            {
                have_any |= op_cnt[t][f][i];
            }
        }
        if ( !have_any ) continue;
        cassert( f < func_names.size(), "func_names doesn't have enough names" );

        fprintf( out, "\n\n%s OP Totals:\n", func_names[f].c_str() );
        csv << "\n\n\"" << func_names[f] << " OP Totals:\"\n";
        for( uint32_t i = 0; i < OP_cnt; i++ )
        {
            uint64_t cnt = 0;
            for( uint32_t t = 0; t < THREAD_CNT_MAX; t++ )
            {
                cnt += op_cnt[t][f][i];
            }
            if ( cnt == 0 ) continue;
            total_op_cnt[i] += cnt;

            uint64_t scaled_cnt = double(cnt) * scale_factor + 0.5;
            fprintf( out, "    %-40s:  %10" FMT_LLU "   %10" FMT_LLU "\n", Cordic<T,FLT>::op_to_str( i ).c_str(), cnt, scaled_cnt );
            csv << "\"" << Cordic<T,FLT>::op_to_str( i ) << "\", " << cnt << ", " << scaled_cnt << "\n";
        }
    }

    fprintf( out, "\n\nGrand OP Totals:\n" );
    csv << "\n\n\"Grand OP Totals:\"\n";
    for( uint32_t i = 0; i < OP_cnt; i++ )
    {
        uint64_t cnt = total_op_cnt[i];
        if ( cnt == 0 ) continue;
        uint64_t scaled_cnt = double(cnt) * scale_factor + 0.5;
        fprintf( out, "    %-40s:  %10" FMT_LLU "   %10" FMT_LLU "\n", Cordic<T,FLT>::op_to_str( i ).c_str(), cnt, scaled_cnt );
        csv << "\"" << Cordic<T,FLT>::op_to_str( i ) << "\", " << cnt << ", " << scaled_cnt << "\n";
    }

    fclose( out );
    csv.close();
    std::cout << "\nWrote stats to " + basename + ".{out,csv}\n";
}

#endif
