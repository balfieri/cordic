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
// Analysis.h - class for analyzing Logger.h text output,
//              but it's also a derived class from Logger and
//              best used directly in a Cordic program to
//              perform the analysis on-the-fly.
//
#ifndef _Analysis_h
#define _Analysis_h

#include <string>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <mutex>

#include "Cordic.h"
#include "Logger.h"

template< typename T=int64_t, typename FLT=double >
class Analysis : public Logger<T,FLT>
{
public:
    Analysis( std::string base_name = "log" );     
    ~Analysis();

    // Call this if you want this to be thread-safe
    //
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

    std::istream *      in;
    bool                in_text;

    static constexpr uint32_t INT_W_MAX = 32;           // maximum int_w we expect

    struct FuncInfo
    {
        uint64_t call_cnt;                                      // number of enters for this function
        uint64_t op_cnt[OP_cnt];                                // total op counts from all calls
        uint64_t opnd_cnt[OP_cnt];                              // total number of operands per OP
        uint64_t opnd_is_const_cnt[OP_cnt];                     // total number of operands per OP that are constants
        uint64_t opnd_all_are_const_cnt[OP_cnt];                // same but only if all operands to OP are constants
        uint64_t opnd_int_w_used_cnt[OP_cnt][INT_W_MAX+1];      // for each op, total number of operands that use each int_w
        uint64_t opnd_max_int_w_used_cnt[OP_cnt][INT_W_MAX+1];  // same but using maximum int_w
    };

    struct FrameInfo
    {
        uint16_t    func_id;
    };

    struct CordicInfo
    {
        size_t   cordic_i;
        bool     is_alive;
        bool     is_float;
        uint32_t int_exp_w;
        uint32_t frac_w;
        uint32_t guard_w;
        uint32_t n;
    };

    struct ValInfo
    {
        bool     is_alive;
        bool     is_assigned;
        size_t   cordic_i;
        bool     is_float;
        uint32_t int_exp_w;
        uint32_t frac_w;
        uint32_t guard_w;
        size_t   opnd_i[3];
        T        encoded;
        uint32_t encoded_int_w_used;
        bool     is_constant;
        FLT      constant;
        FLT      min;
        FLT      max;
    };

    enum class KIND
    {
        cordic_constructed,
        cordic_destructed,
        enter,
        leave,
        constructed,
        destructed,
        op1, 
        op2, 
        op3, 
        op4, 
        op1i, 
        op1b, 
        op1f, 
        op2i, 
        op2f, 
    };

    std::mutex                                  lock;                   // to make this thread-safe

    static constexpr uint32_t KIND_cnt = 12;

    std::map<std::string, KIND>                 kinds;
    std::map<std::string, OP>                   ops;
    std::vector<FuncInfo>                       funcs;
    std::map<uint64_t, CordicInfo>              cordics;
    std::map<uint64_t, ValInfo>                 vals;

    static constexpr uint32_t                   THREAD_CNT_MAX = 64;
    static constexpr uint32_t                   STACK_CNT_MAX = 1024;
    FrameInfo                                   stack[THREAD_CNT_MAX][STACK_CNT_MAX];   // func call stack
    uint32_t                                    stack_cnt[THREAD_CNT_MAX];              // func call stack depth

    static constexpr uint32_t                   VAL_STACK_CNT_MAX = 2;
    ValInfo                                     val_stack[THREAD_CNT_MAX][VAL_STACK_CNT_MAX];
    uint32_t                                    val_stack_cnt[THREAD_CNT_MAX];

    void                stack_push( const FrameInfo& info );
    FrameInfo&          stack_top( void );
    void                stack_pop( void );

    void                val_stack_push( const ValInfo& val );
    ValInfo             val_stack_pop( void );

    void                calc_int_w_used( ValInfo& val );
    void                inc_op_cnt_nolock( OP op, uint32_t by=1 );
    void                inc_opnd_cnt( OP op, const ValInfo& val, uint32_t by=1 );
    void                inc_all_opnd_cnt( OP op, bool all_are_const, uint32_t max_int_w_used, uint32_t by=1 );

    static void         _skip_junk( char *& c );
    static std::string  parse_name( char *& c );
    static KIND         parse_kind( char *& c );
    static const void * parse_addr( char *& c );
    static const T *    parse_val_addr( char *& c );
    static T            parse_int( char *& c );
    static FLT          parse_flt( char *& c );
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

static inline void _die( std::string msg )
{
    std::cout << "ERROR: " << msg << "\n";
    exit( 1 );
}

template< typename T, typename FLT >
Analysis<T,FLT>::Analysis( std::string _base_name ) : Logger<T,FLT>( Cordic<T,FLT>::op_to_str )
{
    base_name = _base_name;

    in_text = true;  // for now
    if ( in_text ) {
        in = &std::cin;
    }

    // set up ops map
    for( uint32_t o = 0; o < Cordic<T,FLT>::OP_cnt; o++ )
    {
        std::string name = Cordic<T,FLT>::op_to_str( o );
        ops[name] = OP(o);
    }

    // set up kinds map
    kinds["cordic_constructed"] = KIND::cordic_constructed;
    kinds["cordic_destructed"]  = KIND::cordic_destructed;
    kinds["enter"]              = KIND::enter;
    kinds["leave"]              = KIND::leave;
    kinds["constructed"]        = KIND::constructed;
    kinds["destructed"]         = KIND::destructed;
    kinds["op1"]                = KIND::op1;
    kinds["op1i"]               = KIND::op1i;
    kinds["op1b"]               = KIND::op1b;
    kinds["op1f"]               = KIND::op1f;
    kinds["op2"]                = KIND::op2;
    kinds["op2i"]               = KIND::op2i;
    kinds["op2f"]               = KIND::op2f;
    kinds["op3"]                = KIND::op3;
    kinds["op4"]                = KIND::op4;

    for( uint32_t t = 0; t < THREAD_CNT_MAX; t++ )
    {
        stack_cnt[t]     = 0;
        val_stack_cnt[t] = 0;
    }
}

template< typename T, typename FLT >
Analysis<T,FLT>::~Analysis()
{
}

static thread_local uint32_t tid = 0;

template< typename T, typename FLT >
void Analysis<T,FLT>::tid_set( uint32_t t )
{
    cassert( t < THREAD_CNT_MAX, "tid_set: tid must be < THREAD_CNT_MAX" );
    tid = t;
}

//-----------------------------------------------------
// Logger Method Overrides
//-----------------------------------------------------
template< typename T, typename FLT >
void Analysis<T,FLT>::cordic_constructed( const void * cordic_ptr, uint32_t int_exp_w, uint32_t frac_w, 
                                          bool is_float, uint32_t guard_w, uint32_t n )
{
    std::lock_guard<std::mutex> guard(lock);
    CordicInfo info;
    uint64_t cordic = uint64_t(cordic_ptr); 
    info.is_alive   = true;
    info.is_float   = is_float;
    info.int_exp_w  = int_exp_w;
    info.frac_w     = frac_w;
    info.guard_w    = guard_w;
    info.n          = n;

    auto it = cordics.find( cordic );
    cassert( it == cordics.end() || !it->second.is_alive, "Cordic reconstructed before previous was destructed" );
    cordics[cordic] = info;
}

template< typename T, typename FLT >
void Analysis<T,FLT>::cordic_destructed( const void * cordic_ptr )
{
    std::lock_guard<std::mutex> guard(lock);
    uint64_t cordic = uint64_t(cordic_ptr); 
    auto it = cordics.find( cordic );
    cassert( it != cordics.end() && it->second.is_alive, "Cordic destructed before being constructed" );
    it->second.is_alive = false;
}

template< typename T, typename FLT >
void Analysis<T,FLT>::enter( uint16_t func_id )
{
    std::lock_guard<std::mutex> guard(lock);
    size_t size = funcs.size();
    if ( size <= func_id ) {
        funcs.resize( func_id+1 );
        for( size_t i = size; i <= func_id; i++ )
        {
            memzero( &funcs[i], sizeof(funcs[0]) );
        }
    }
    funcs[func_id].call_cnt++; 
    FrameInfo frame;
    frame.func_id = func_id;
    stack_push( frame );
}

template< typename T, typename FLT >
void Analysis<T,FLT>::leave( uint16_t func_id )
{
    std::lock_guard<std::mutex> guard(lock);
    cassert( funcs.size() > func_id, "func_id " + std::to_string(func_id) + " does not exist" );
    FrameInfo& frame = stack_top();
    cassert( frame.func_id == func_id , "trying to leave a routine that's not at the top of the stack: entered " + 
                                        frame.func_id + " leaving " + func_id );
    stack_pop();
}

template< typename T, typename FLT >
void Analysis<T,FLT>::constructed( const T * v, const void * cordic_ptr )
{
    std::lock_guard<std::mutex> guard(lock);
    uint64_t val    = reinterpret_cast<uint64_t>( v );
    uint64_t cordic = reinterpret_cast<uint64_t>( cordic_ptr );
    ValInfo info;
    info.is_alive    = true;
    info.is_assigned = false;
    info.is_constant = false;
    if ( cordic != 0 ) {
        auto cit = cordics.find( cordic );
        cassert( cit != cordics.end() && cit->second.is_alive, "val constructed using unknown cordic" );
        info.cordic_i  = cit->second.cordic_i;
        info.is_float  = cit->second.is_float;
        info.int_exp_w = cit->second.int_exp_w;
        info.frac_w    = cit->second.frac_w;
        info.guard_w   = cit->second.guard_w;
    } else {
        info.cordic_i  = size_t(-1);
        info.is_float  = true;
        info.int_exp_w = 0;
        info.frac_w    = 0;
        info.guard_w   = 0;
    }
    auto vit = vals.find( val );
    if ( vit == vals.end() ) {
        vals[val] = info;
    } else {
        //cassert( !vit->second.is_alive, "val constructed before previous was desctructed" );
        vit->second = info;
    }
}

template< typename T, typename FLT >
void Analysis<T,FLT>::destructed(  const T * v, const void * cordic_ptr )
{
    std::lock_guard<std::mutex> guard(lock);
    uint64_t val    = reinterpret_cast<uint64_t>( v );
    uint64_t cordic = reinterpret_cast<uint64_t>( cordic_ptr );
    auto it = vals.find( val );
    cassert( it != vals.end() && it->second.is_alive, "val destructed before being constructed" );
    it->second.is_alive = false;
}

template< typename T, typename FLT >
void Analysis<T,FLT>::op( uint16_t _op, uint32_t opnd_cnt, const T * opnd[] )
{
    std::lock_guard<std::mutex> guard(lock);
    OP op = OP(_op);
    inc_op_cnt_nolock( op );
    uint32_t max_int_w_used = 0;
    bool     all_are_const = true;
    for( uint32_t i = 0; i < opnd_cnt; i++ )
    {
        if ( !(i == 0 && op == OP::assign) &&
             !(i == 1 && op == OP::sincos) &&
             !(i == 2 && op == OP::sincos) &&
             !(i == 1 && op == OP::sinhcosh) &&
             !(i == 2 && op == OP::sinhcosh) ) {
            auto it = vals.find( reinterpret_cast<uint64_t>( opnd[i] ) );
            cassert( it != vals.end() && it->second.is_alive, "opnd[" + std::to_string(i) + "] does not exist" );
            cassert( it->second.is_assigned, "opnd[" + std::to_string(i) + "] used when not previously assigned" );
            inc_opnd_cnt( op, it->second );
            if ( it->second.encoded_int_w_used > max_int_w_used ) max_int_w_used = it->second.encoded_int_w_used;
            all_are_const &= it->second.is_constant;
            if ( i == 1 && op == OP::assign ) vals[reinterpret_cast<uint64_t>(opnd[0])] = it->second;
            if ( debug && it->second.is_constant ) {
                std::cout << "    opnd[" + std::to_string(i) + "] is constant " << it->second.constant << "\n";
            }
        }
    }

    inc_all_opnd_cnt( op, all_are_const, max_int_w_used );

    // push result if not assign
    uint32_t cnt = (op == OP::sincos || op == OP::sinhcosh) ? 2 : 
                   (op == OP::assign)                       ? 0 : 1;
    ValInfo val;
    val.is_alive    = true;
    val.is_assigned = true;
    val.is_constant = false;
    for ( uint32_t i = 0; i < cnt; i++ ) val_stack_push( val );
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::op1( uint16_t _op, const T * opnd1 )
{
    op( _op, 1, &opnd1 );
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::op1( uint16_t _op, const T& opnd1 )
{
    (void)_op;
    (void)opnd1;
    _die( "op1i should not be used right now" );
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::op1( uint16_t _op, bool opnd1 )
{
    std::lock_guard<std::mutex> guard(lock);
    (void)opnd1;
    OP op = OP(_op);
    cassert( op == OP::pop_bool, "op1b allowed only for pop_bool right now" );
    inc_op_cnt_nolock( op );
    (void)val_stack_pop();
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::op1( uint16_t _op, const FLT& opnd1 )
{
    std::lock_guard<std::mutex> guard(lock);
    OP op = OP(_op);
    cassert( op == OP::push_constant, "op1f allowed only for make_constant" );
    inc_op_cnt_nolock( op );
    ValInfo val;
    val.is_alive    = true;
    val.is_assigned = true;
    val.is_constant = true;
    val.constant    = opnd1;
    val_stack_push( val );
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::op2( uint16_t _op, const T * opnd1, const T * opnd2 )
{
    const T * opnds[] = { opnd1, opnd2 };
    op( _op, 2, opnds );
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::calc_int_w_used( ValInfo& val )
{
    // calculate number of bits needed to hold integer part of abs(x)
    cassert( (val.int_exp_w + val.frac_w + val.guard_w) != 0, "calc_int_w_used: int_exp_w/frac_w/guard_w are not defined" );
    T x = val.encoded;
    if ( x < T(0) ) x = -x;
    x >>= val.frac_w + val.guard_w;

    // ceil(log2(x))
    uint32_t lg2 = 0;
    while( x != T(0) )
    {
        lg2++;
        x >>= 1;
    }
    val.encoded_int_w_used = lg2;
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::op2( uint16_t _op, const T * opnd1, const T& opnd2 )
{
    std::lock_guard<std::mutex> guard(lock);
    OP op = OP(_op);
    cassert( op == OP::scalbn || op == OP::pop_value, "op2i allowed only for scalbn/pop_value" );
    inc_op_cnt_nolock( op );
    auto it = vals.find( reinterpret_cast<uint64_t>( opnd1 ) );
    cassert( it != vals.end() && it->second.is_alive, "opnd[0] does not exist" );
    switch( op )
    {
        case OP::pop_value:
        {
            // pop result
            ValInfo pval = val_stack_pop();
            it->second.is_assigned = true;
            it->second.is_constant = pval.is_constant;
            it->second.encoded     = opnd2;
            calc_int_w_used( it->second );
            break;
        }

        default:
        {
            // push result
            ValInfo val;
            val.is_alive    = true;
            val.is_assigned = true;
            val.is_constant = false;
            val.encoded     = opnd2;
            val_stack_push( val );
            break;
        }
    }
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::op2( uint16_t op, const T * opnd1, const FLT& opnd2 ) 
{
    std::lock_guard<std::mutex> guard(lock);
    inc_op_cnt_nolock( OP(op) );
    auto it = vals.find( reinterpret_cast<uint64_t>( opnd1 ) );
    cassert( it != vals.end() && it->second.is_alive, "opnd1 does not exist" );
    cassert( it->second.is_assigned,                  "opnd1 is used before being assigned" );
    ValInfo val;
    val.is_alive    = true;
    val.is_assigned = true;
    val.is_constant = false;
    (void)opnd2;
//  val.constant    = opnd2;   // save conversion to FLT
    val_stack_push( val );
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::op3( uint16_t _op, const T * opnd1, const T * opnd2, const T * opnd3 )
{
    const T * opnds[] = { opnd1, opnd2, opnd3 };
    op( _op, 3, opnds );
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::op4( uint16_t _op, const T * opnd1, const T * opnd2, const T * opnd3, const T * opnd4 )
{
    const T * opnds[] = { opnd1, opnd2, opnd3, opnd4 };
    op( _op, 4, opnds );
}

//-----------------------------------------------------
// Parsing Stuff
//-----------------------------------------------------
template< typename T, typename FLT >
inline void Analysis<T,FLT>::_skip_junk( char *& c )
{
    // skip spaces, '(' and ','
    for( ;; )
    {
        char ch = *c;
        if ( ch == '\0' ) return;
        if ( ch == ' ' || ch == ',' || ch == '(' ) {
            c++;
        } else {
            return;
        }
    }
}

template< typename T, typename FLT >
inline std::string Analysis<T,FLT>::parse_name( char *& c )
{
    // read string of letters, "::" and numbers
    _skip_junk( c );
    std::string name = "";
    for( ;; )
    {
        char ch = *c;
        if ( ch == ':' || ch == '_' || ch == '-' || ch == '.' || 
             (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') || (ch >= '0' && ch <= '9') ) {
            name += ch;
            c++;
        } else {
            break;
        }
    }

    cassert( name.size() > 0, "could not parse a name" );
    return name;
}

template< typename T, typename FLT >
inline typename Analysis<T,FLT>::KIND Analysis<T,FLT>::parse_kind( char *& c )
{
    std::string name = parse_name( c );

    if ( name == "cordic_constructed" ) return KIND::cordic_constructed;
    if ( name == "cordic_destructed" )  return KIND::cordic_destructed;
    if ( name == "enter" )              return KIND::enter;
    if ( name == "leave" )              return KIND::leave;
    if ( name == "constructed" )        return KIND::constructed;
    if ( name == "destructed" )         return KIND::destructed;
    if ( name == "op1" )                return KIND::op1;
    if ( name == "op1i" )               return KIND::op1i;
    if ( name == "op1f" )               return KIND::op1f;
    if ( name == "op2" )                return KIND::op2;
    if ( name == "op2i" )               return KIND::op2i;
    if ( name == "op2f" )               return KIND::op2f;
    if ( name == "op3" )                return KIND::op3;
    if ( name == "op4" )                return KIND::op4;
    return KIND(-1);
}

template< typename T, typename FLT >
inline const void * Analysis<T,FLT>::parse_addr( char *& c )
{
    std::string addr_s = parse_name( c );
    char * end;
    return reinterpret_cast<const void *>( std::strtoull( addr_s.c_str(), &end, 16 ) );
}

template< typename T, typename FLT >
inline const T * Analysis<T,FLT>::parse_val_addr( char *& c )
{
    return reinterpret_cast<const T *>( parse_addr( c ) );
}

template< typename T, typename FLT >
inline T Analysis<T,FLT>::parse_int( char *& c )
{
    std::string int_s = parse_name( c );
    return std::atoi( int_s.c_str() );
}

template< typename T, typename FLT >
inline FLT Analysis<T,FLT>::parse_flt( char *& c )
{
    std::string flt_s = parse_name( c );
    return std::atof( flt_s.c_str() );
}

template< typename T, typename FLT > inline void Analysis<T,FLT>::stack_push( const FrameInfo& info )
{
    cassert( stack_cnt[tid] < STACK_CNT_MAX, "depth of call stack exceeded" );
    stack[tid][stack_cnt[tid]++] = info;
}

template< typename T, typename FLT >
inline typename Analysis<T,FLT>::FrameInfo& Analysis<T,FLT>::stack_top( void )
{
    cassert( stack_cnt[tid] > 0, "can't get top of an empty call stack" );
    return stack[tid][stack_cnt[tid]-1];
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::stack_pop( void )
{
    cassert( stack_cnt[tid] > 0, "can't pop an empty call stack" );
    stack_cnt[tid]--;
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::inc_op_cnt_nolock( OP op, uint32_t by )
{
    FrameInfo& frame = stack_top();
    FuncInfo& func = funcs[frame.func_id];
    func.op_cnt[uint32_t(op)] += by;
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::inc_op_cnt( OP op, uint32_t by )
{
    std::lock_guard<std::mutex> guard(lock);
    inc_op_cnt_nolock( op, by );
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::inc_opnd_cnt( OP op, const ValInfo& val, uint32_t by )
{
    FrameInfo& frame = stack_top();
    FuncInfo& func = funcs[frame.func_id];
    uint16_t op_i = uint16_t(op);
    func.opnd_cnt[op_i] += by;
    if ( val.is_constant ) func.opnd_is_const_cnt[op_i] += by;
    uint32_t int_w_used = val.encoded_int_w_used;
    if ( int_w_used > INT_W_MAX ) int_w_used = INT_W_MAX;
    func.opnd_int_w_used_cnt[op_i][int_w_used] += by;
}

template< typename T, typename FLT > 
inline void Analysis<T,FLT>::inc_all_opnd_cnt( OP op, bool all_are_const, uint32_t max_int_w_used, uint32_t by )
{
    FrameInfo& frame = stack_top();
    FuncInfo& func = funcs[frame.func_id];
    uint16_t op_i = uint16_t(op);
    if ( all_are_const ) func.opnd_all_are_const_cnt[op_i] += by;
    if ( max_int_w_used > INT_W_MAX ) max_int_w_used = INT_W_MAX;
    func.opnd_max_int_w_used_cnt[op_i][max_int_w_used] += by;
}

template< typename T, typename FLT > 
inline void Analysis<T,FLT>::val_stack_push( const ValInfo& info )
{
    cassert( val_stack_cnt[tid] < VAL_STACK_CNT_MAX, "depth of val_stack exceeded" );
    val_stack[tid][val_stack_cnt[tid]++] = info;
}

template< typename T, typename FLT >
inline typename Analysis<T,FLT>::ValInfo Analysis<T,FLT>::val_stack_pop( void )
{
    cassert( val_stack_cnt[tid] > 0, "can't pop an empty val_stack" );
    return val_stack[tid][--val_stack_cnt[tid]];
}

template< typename T, typename FLT >
void Analysis<T,FLT>::parse( void )
{
    // assume parsing text at this point
    std::string line;
    char cs[1024];
    for( uint32_t t = 0; t < THREAD_CNT_MAX; t++ )
    {
        stack_cnt[t] = 0;
        val_stack_cnt[t] = 0;
    }
    while( std::getline( *in, line ) )
    {
        if ( debug ) std::cout << line << "\n";
        strcpy( cs, line.c_str() );
        char * c = cs;
        const KIND kind = parse_kind( c );
        switch( kind )
        {
            case KIND::cordic_constructed:
            {
                const void * cordic    = parse_addr( c );
                bool         is_float  = parse_int( c ) != 0;
                T            int_exp_w = parse_int( c );
                T            frac_w    = parse_int( c );
                T            guard_w   = parse_int( c );
                T            n         = parse_int( c );
                cordic_constructed( cordic, int_exp_w, frac_w, is_float, guard_w, n );
                break;
            }

            case KIND::cordic_destructed:
            {
                const void * cordic  = parse_addr( c );
                cordic_destructed( cordic );
                break;
            }

            case KIND::enter:
            {
                std::string name = parse_name( c );
                enter( name.c_str() );
                break;
            }

            case KIND::leave:
            {
                std::string name = parse_name( c );
                leave( name.c_str() );
                break;
            }

            case KIND::constructed:
            {
                const T    * val     = parse_val_addr( c );
                const void * cordic  = parse_addr( c );
                constructed( val, cordic );
                break;
            }

            case KIND::destructed:
            {
                const T    * val     = parse_val_addr( c );
                const void * cordic  = parse_addr( c );
                destructed( val, cordic );
                break;
            }

            case KIND::op1:
            case KIND::op2:
            case KIND::op3:
            case KIND::op4:
            {
                std::string name = parse_name( c );
                uint16_t o = uint16_t( ops[name] );
                uint32_t opnd_cnt = uint32_t(kind) - uint32_t(KIND::op1) + 1;
                const T * opnd[4];
                for( uint32_t i = 0; i < opnd_cnt; i++ )
                {
                    opnd[i] = parse_val_addr( c );
                }
                op( uint16_t(o), opnd_cnt, opnd );
                break;
            }

            case KIND::op1i:
            {
                std::string name  = parse_name( c );
                T           opnd0 = parse_int( c );
                OP op = ops[name];
                op1( uint16_t(op), opnd0 );
                break;
            }

            case KIND::op1b:
            {
                std::string name  = parse_name( c );
                bool        opnd0 = parse_int( c );
                OP op = ops[name];
                op1( uint16_t(op), opnd0 );
                break;
            }

            case KIND::op1f:
            {
                std::string name  = parse_name( c );
                FLT         opnd0 = parse_flt( c );
                OP op = ops[name];
                op1( uint16_t(op), opnd0 );
                break;
            }

            case KIND::op2i:
            {
                std::string name = parse_name( c );
                OP op = ops[name];
                const T * opnd0 = parse_val_addr( c );
                T         opnd1 = parse_int( c );
                op2( uint16_t(op), opnd0, opnd1 );
                break;
            }

            case KIND::op2f:
            {
                std::string name = parse_name( c );
                OP op = ops[name];
                const T * opnd0 = parse_val_addr( c );   
                FLT       opnd1 = parse_flt( c );   
                op2( uint16_t(op), opnd0, opnd1 );
                break;
            }

            default:
            {
                continue;
            }
        }
    }
}

template< typename T, typename FLT >
void Analysis<T,FLT>::clear_stats( void )
{
    //--------------------------------------------------------
    // Print only the non-zero counts from non-ignored functions.
    //--------------------------------------------------------
    for( size_t i = 0; i < funcs.size(); i++ )
    {
        FuncInfo& func = &funcs[i];
        func.call_cnt = 0;
        for( uint32_t j = 0; j < OP_cnt; j++ )
        {
            func.op_cnt[j] = 0;
        }
    }
}

template< typename T, typename FLT >
void Analysis<T,FLT>::print_stats( std::string basename, double scale_factor, 
                                   const std::vector<std::string>& func_names, const std::vector<uint16_t>& ignore_funcs ) const
{
    //--------------------------------------------------------
    // Print only the non-zero counts from non-ignored functions.
    //--------------------------------------------------------
    if ( basename == "" ) basename = base_name;
    std::map<uint16_t, bool> func_ignored;
    for( auto it = ignore_funcs.begin(); it != ignore_funcs.end(); it++ )
    {
        func_ignored[*it] = true;
    }
    std::string out_name = basename + ".out";
    FILE * out = fopen( out_name.c_str(), "w" );
    std::ofstream csv( basename + ".csv", std::ofstream::out );
    uint64_t total_op_cnt[OP_cnt];
    uint64_t total_opnd_cnt[OP_cnt];
    uint64_t total_opnd_is_const_cnt[OP_cnt];                     
    uint64_t total_opnd_all_are_const_cnt[OP_cnt];               
    uint64_t total_opnd_int_w_used_cnt[OP_cnt][INT_W_MAX+1];    
    uint64_t total_opnd_max_int_w_used_cnt[OP_cnt][INT_W_MAX+1];
    for( uint32_t i = 0; i < OP_cnt; i++ ) 
    { 
        total_op_cnt[i] = 0;
        total_opnd_cnt[i] = 0;
        total_opnd_is_const_cnt[i] = 0;
        total_opnd_all_are_const_cnt[i] = 0;
        for( uint32_t j = 0; j < (INT_W_MAX+1); j++ )
        {
            total_opnd_int_w_used_cnt[i][j] = 0;
            total_opnd_max_int_w_used_cnt[i][j] = 0;
        }
    }
    cassert( funcs.size() <= func_names.size(), "func_names doesn't have enough names" );
    for( uint32_t for_opnds = 0; for_opnds < 2; for_opnds++ )
    {
        fprintf( out, for_opnds ? "\nOPERAND COUNTS:\n" : "\nOP COUNTS:\n" );
        for( size_t i = 0; i < funcs.size(); i++ )
        {
            if ( func_ignored.find( i ) != func_ignored.end() ) continue;
            const FuncInfo& func = funcs[i];
            if ( !for_opnds ) {
                fprintf( out, "\n%-20s: %8lld calls\n", func_names[i].c_str(), func.call_cnt );
                csv << "\n\"" << func_names[i] << "\", " << func.call_cnt << "\n";
            } else {
                fprintf( out, "\n%4d:\n", i );
            }
            for( uint32_t j = 0; j < OP_cnt; j++ )
            {
                OP op = OP(j);
                if ( op == OP::push_constant || op == OP::assign || op == OP::pop_value || op == OP::pop_bool ) continue; // consume no hardware

                uint64_t cnt = func.op_cnt[i];
                if ( cnt == 0 ) continue;

                if ( !for_opnds ) {
                    total_op_cnt[i] += cnt;
                    double avg = double(cnt) / double(func.call_cnt);
                    uint64_t scaled_cnt = double(cnt) * scale_factor + 0.5;
                    fprintf( out, "    %-40s: %8.1f/call   %10lld total   %10lld scaled_total\n", Cordic<T,FLT>::op_to_str( i ).c_str(), avg, cnt, scaled_cnt );
                    csv << "\"" << Cordic<T,FLT>::op_to_str( i ) << "\", " << avg << ", " << cnt << ", " << scaled_cnt << "\n";
                } else {
                    fprintf( out, "    %s:\n", Cordic<T,FLT>::op_to_str( i ).c_str() );
                    fprintf( out, "        %-50s: %lld\n", "Total op count", cnt );
                    fprintf( out, "        %-50s: %lld\n", "Total operand count", func.opnd_cnt[i] );
                    fprintf( out, "        %-50s: %lld\n", "Total operands that were constants", func.opnd_is_const_cnt[i] );
                    fprintf( out, "        %-50s: %lld\n", "Total times all operands were constants", func.opnd_all_are_const_cnt[i] );
                    total_opnd_cnt[i] += func.opnd_cnt[i];
                    total_opnd_is_const_cnt[i] += func.opnd_is_const_cnt[i];
                    total_opnd_all_are_const_cnt[i] += func.opnd_all_are_const_cnt[i];
                    for( uint32_t w = 0; w <= INT_W_MAX; w++ )
                    {
                        uint64_t wcnt = func.opnd_int_w_used_cnt[i][w]; 
                        if ( wcnt == 0 ) continue;
                        std::string s = "Total operands that fit into " + std::to_string(w) + " integer bits";
                        fprintf( out, "        %-50s: %lld\n", s.c_str(), wcnt );
                        total_opnd_int_w_used_cnt[i][w] += wcnt;
                    }
                }
            }
        }

        //--------------------------------------------------------
        // And the totals.
        //--------------------------------------------------------
        if ( !for_opnds ) {
            fprintf( out, "\n\nOP Totals:\n" );
            csv << "\n\n\"Totals:\"" << "\n";
            for( uint32_t i = 0; i < OP_cnt; i++ )
            {
                if ( total_op_cnt[i] == 0 ) continue;

                uint64_t cnt = total_op_cnt[i];
                uint64_t scaled_cnt = double(cnt) * scale_factor + 0.5;
                fprintf( out, "    %-40s:  %10lld   %10lld\n", Cordic<T,FLT>::op_to_str( i ).c_str(), cnt, scaled_cnt );
                csv << "\"" << Cordic<T,FLT>::op_to_str( i ) << "\", " << cnt << ", " << scaled_cnt << "\n";
            }
        } else {
            fprintf( out, "\n\nOPND Totals:\n" );
            for( uint32_t i = 0; i < OP_cnt; i++ )
            {
                if ( total_op_cnt[i] == 0 ) continue;

                uint64_t cnt = total_op_cnt[i];
                uint64_t scaled_cnt = double(cnt) * scale_factor + 0.5;
                fprintf( out, "    %s:\n", Cordic<T,FLT>::op_to_str( i ).c_str() );
                fprintf( out, "        %-50s: %lld\n", "Total op count", total_op_cnt[i] );
                fprintf( out, "        %-50s: %lld\n", "Total operand count", total_opnd_cnt[i] );
                fprintf( out, "        %-50s: %lld\n", "Total operands that were constants", total_opnd_is_const_cnt[i] );
                fprintf( out, "        %-50s: %lld\n", "Total times all operands were constants", total_opnd_all_are_const_cnt[i] );
                for( uint32_t w = 0; w <= INT_W_MAX; w++ )
                {
                    uint64_t wcnt = total_opnd_int_w_used_cnt[i][w]; 
                    if ( wcnt == 0 ) continue;
                    std::string s = "Total operands that fit into " + std::to_string(w) + " integer bits";
                    fprintf( out, "        %-50s: %lld\n", s.c_str(), wcnt );
                }
            }
        }
    }

    fclose( out );
    csv.close();
    std::cout << "\nWrote stats to " + basename + ".{out,csv}\n";
}

#endif
