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
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "Cordic.h"
#include "Logger.h"

template< typename T=int64_t, typename FLT=double >
class Analysis : public Logger<T,FLT>
{
public:
    Analysis( std::string base_name = "log" );     
    ~Analysis();

    // Logger Overrides
    //
    virtual void cordic_constructed( const void * cordic, uint32_t int_w, uint32_t frac_w, uint32_t guard_w, uint32_t n );
    virtual void cordic_destructed(  const void * cordic );

    virtual void enter( const char * name );
    virtual void leave( const char * name );

    virtual void constructed( const T * v, const void * cordic );
    virtual void destructed(  const T * v, const void * cordic );

    virtual void op(  uint16_t op, uint32_t opnd_cnt, const T * opnd[] );
    virtual void op1( uint16_t op, const T * opnd1 );
    virtual void op1( uint16_t op, const T&  opnd1 );
    virtual void op1( uint16_t op, const FLT&opnd1 );
    virtual void op2( uint16_t op, const T * opnd1, const T * opnd2 );
    virtual void op2( uint16_t op, const T * opnd1, const T&  opnd2 );
    virtual void op2( uint16_t op, const T * opnd1, const FLT&opnd2 );
    virtual void op3( uint16_t op, const T * opnd1, const T * opnd2, const T * opnd3 );
    virtual void op4( uint16_t op, const T * opnd1, const T * opnd2, const T * opnd3, const T * opnd4 );

    void         parse( void );
    void         print_stats( double scale_factor=1.0, const std::vector<std::string>& ignore_funcs=std::vector<std::string>() ) const;

private:
    std::string         base_name;

    std::istream *      in;
    bool                in_text;

    using OP                         = typename Cordic<T,FLT>::OP;
    static constexpr uint64_t OP_cnt = Cordic<T,FLT>::OP_cnt;

    struct FuncInfo
    {
        uint64_t call_cnt;                              // number of enters for this function
        uint64_t op_cnt[OP_cnt];                        // total op counts from all calls
    };

    struct FrameInfo
    {
        std::string func_name;                          
    };

    struct CordicInfo
    {
        size_t   cordic_i;
        bool     is_alive;
        uint32_t int_w;
        uint32_t frac_w;
        uint32_t guard_w;
        uint32_t n;
    };

    struct ValInfo
    {
        bool    is_alive;
        bool    is_assigned;
        size_t  cordic_i;
        size_t  opnd_i[3];
        T       encoded;
        bool    is_constant;
        FLT     constant;
        FLT     min;
        FLT     max;
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
        op1f, 
        op2i, 
        op2f, 
    };

    static constexpr uint32_t KIND_cnt = 12;

    std::map<std::string, KIND>                 kinds;
    std::map<std::string, OP>                   ops;
    std::vector<std::string>                    func_names;             // in order of appearance
    std::map<std::string, FuncInfo>             funcs;
    std::map<uint64_t, CordicInfo>              cordics;
    std::map<uint64_t, ValInfo>                 vals;

    static constexpr uint32_t                   STACK_CNT_MAX = 1024;
    FrameInfo                                   stack[STACK_CNT_MAX];   // func call stack
    uint32_t                                    stack_cnt;              // func call stack depth

    static constexpr uint32_t                   VAL_STACK_CNT_MAX = 2;
    ValInfo                                     val_stack[VAL_STACK_CNT_MAX];
    uint32_t                                    val_stack_cnt;

    void                stack_push( const FrameInfo& info );
    FrameInfo&          stack_top( void );
    void                stack_pop( void );
    void                inc_op_cnt( OP op );

    void                val_stack_push( const ValInfo& val );
    ValInfo             val_stack_pop( void );

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

//-----------------------------------------------------
// Logger Method Overrides
//-----------------------------------------------------
template< typename T, typename FLT >
void Analysis<T,FLT>::cordic_constructed( const void * cordic_ptr, uint32_t int_w, uint32_t frac_w, uint32_t guard_w, uint32_t n )
{
    CordicInfo info;
    uint64_t cordic = uint64_t(cordic_ptr); 
    info.is_alive   = true;
    info.int_w      = int_w;
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
    uint64_t cordic = uint64_t(cordic_ptr); 
    auto it = cordics.find( cordic );
    cassert( it != cordics.end() && it->second.is_alive, "Cordic destructed before being constructed" );
    it->second.is_alive = false;
}

template< typename T, typename FLT >
void Analysis<T,FLT>::enter( const char * _name )
{
    std::string name = _name;
    auto it = funcs.find( name );
    if ( it == funcs.end() ) {
        func_names.push_back( name );       // for printouts
        FuncInfo info;
        funcs[name] = info;
        it = funcs.find( name );
        it->second.call_cnt = 0;
        for( uint32_t i = 0; i < OP_cnt; i++ )
        {
            it->second.op_cnt[i] = 0;
        }
    } 
    it->second.call_cnt++;
    FrameInfo frame;
    frame.func_name = name;
    stack_push( frame );
}

template< typename T, typename FLT >
void Analysis<T,FLT>::leave( const char * _name )
{
    std::string name = _name;
    auto it = funcs.find( name );
    cassert( it != funcs.end(), "leave should have found function " + name );
    FrameInfo& frame = stack_top();
    cassert( frame.func_name == name, "trying to leave a routine that's not at the top of the stack: entered " + 
                                      frame.func_name + " leaving " + name );
    stack_pop();
}

template< typename T, typename FLT >
void Analysis<T,FLT>::constructed( const T * v, const void * cordic_ptr )
{
    uint64_t val    = reinterpret_cast<uint64_t>( v );
    uint64_t cordic = reinterpret_cast<uint64_t>( cordic_ptr );
    ValInfo info;
    info.is_alive    = true;
    info.is_assigned = false;
    info.is_constant = false;
    if ( cordic != 0 ) {
        auto cit = cordics.find( cordic );
        cassert( cit != cordics.end() && cit->second.is_alive, "val constructed using unknown cordic" );
        info.cordic_i = cit->second.cordic_i;
    } else {
        info.cordic_i = size_t(-1);
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
    uint64_t val    = reinterpret_cast<uint64_t>( v );
    uint64_t cordic = reinterpret_cast<uint64_t>( cordic_ptr );
    auto it = vals.find( val );
    cassert( it != vals.end() && it->second.is_alive, "val destructed before being constructed" );
    it->second.is_alive = false;
}

template< typename T, typename FLT >
void Analysis<T,FLT>::op( uint16_t _op, uint32_t opnd_cnt, const T * opnd[] )
{
    OP op = OP(_op);
    inc_op_cnt( op );
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
            if ( i == 1 && op == OP::assign ) vals[reinterpret_cast<uint64_t>(opnd[0])] = it->second;
            if ( debug && it->second.is_constant ) {
                std::cout << "    opnd[" + std::to_string(i) + "] is constant " << it->second.constant << "\n";
            }
        }
    }
    
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
inline void Analysis<T,FLT>::op1( uint16_t _op, const FLT& opnd1 )
{
    OP op = OP(_op);
    cassert( op == OP::push_constant, "op1f allowed only for make_constant" );
    inc_op_cnt( op );
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
inline void Analysis<T,FLT>::op2( uint16_t _op, const T * opnd1, const T& opnd2 )
{
    OP op = OP(_op);
    cassert( op == OP::lshift || op == OP::rshift || op == OP::pop_value, "op2i allowed only for lshift/rshift/pop_value" );
    inc_op_cnt( op );
    auto it = vals.find( reinterpret_cast<uint64_t>( opnd1 ) );
    cassert( it != vals.end() && it->second.is_alive, "opnd[0] does not exist" );
    switch( op )
    {
        case OP::pop_value:
        {
            // pop result
            it->second  = val_stack_pop();
            it->second.encoded = opnd2;
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
    inc_op_cnt( OP(op) );
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
    cassert( stack_cnt < STACK_CNT_MAX, "depth of call stack exceeded" );
    stack[stack_cnt++] = info;
}

template< typename T, typename FLT >
inline typename Analysis<T,FLT>::FrameInfo& Analysis<T,FLT>::stack_top( void )
{
    cassert( stack_cnt > 0, "can't get top of an empty call stack" );
    return stack[stack_cnt-1];
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::stack_pop( void )
{
    cassert( stack_cnt > 0, "can't pop an empty call stack" );
    stack_cnt--;
}

template< typename T, typename FLT >
inline void Analysis<T,FLT>::inc_op_cnt( OP op )
{
    FrameInfo& frame = stack_top();
    //std::cout << "incr " << frame.func_name << " " << Cordic<T,FLT>::op_to_str( uint16_t(op) ) << "\n";
    FuncInfo& func = funcs[frame.func_name];
    func.op_cnt[uint32_t(op)]++;
}

template< typename T, typename FLT > 
inline void Analysis<T,FLT>::val_stack_push( const ValInfo& info )
{
    cassert( val_stack_cnt < VAL_STACK_CNT_MAX, "depth of val_stack exceeded" );
    val_stack[val_stack_cnt++] = info;
    //if ( 0 && debug ) std::cout << "    new stack depth=" << val_stack_cnt << "\n";
}

template< typename T, typename FLT >
inline typename Analysis<T,FLT>::ValInfo Analysis<T,FLT>::val_stack_pop( void )
{
    cassert( val_stack_cnt > 0, "can't pop an empty val_stack" );
    //if ( debug ) std::cout << "    new stack depth=" << (val_stack_cnt-1) << "\n";
    return val_stack[--val_stack_cnt];
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
    kinds["op2"]                = KIND::op2;
    kinds["op2i"]               = KIND::op2i;
    kinds["op2f"]               = KIND::op2f;
    kinds["op3"]                = KIND::op3;
    kinds["op4"]                = KIND::op4;

}

template< typename T, typename FLT >
Analysis<T,FLT>::~Analysis()
{
}

template< typename T, typename FLT >
void Analysis<T,FLT>::parse( void )
{
    // assume parsing text at this point
    std::string line;
    char cs[1024];
    stack_cnt = 0;
    val_stack_cnt = 0;
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
                const void * cordic  = parse_addr( c );
                T            int_w   = parse_int( c );
                T            frac_w  = parse_int( c );
                T            guard_w = parse_int( c );
                T            n       = parse_int( c );
                cordic_constructed( cordic, int_w, frac_w, guard_w, n );
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
void Analysis<T,FLT>::print_stats( double scale_factor, const std::vector<std::string>& ignore_funcs ) const
{
    //--------------------------------------------------------
    // Print only the non-zero counts from non-ignored functions.
    //--------------------------------------------------------
    std::map<std::string, bool> func_ignored;
    for( auto it = ignore_funcs.begin(); it != ignore_funcs.end(); it++ )
    {
        func_ignored[*it] = true;
    }
    std::string out_name = base_name + ".out";
    FILE * out = fopen( out_name.c_str(), "w" );
    std::ofstream csv( base_name + ".csv", std::ofstream::out );
    uint64_t total_cnt[OP_cnt];
    for( uint32_t i = 0; i < OP_cnt; i++ ) total_cnt[i] = 0;

    for( auto nit = func_names.begin(); nit != func_names.end(); nit++ )
    {
        if ( func_ignored.find( *nit ) != func_ignored.end() ) continue;
        auto it = funcs.find( *nit );
        const FuncInfo& func = it->second;
        fprintf( out, "\n%-44s: %8lld calls\n", it->first.c_str(), it->second.call_cnt );
        csv << "\n\"" << it->first + "\", " << it->second.call_cnt << "\n";
        for( uint32_t i = 0; i < OP_cnt; i++ )
        {
            OP op = OP(i);
            if ( op == OP::push_constant || op == OP::assign || op == OP::pop_value ) continue; // consume no hardware

            uint64_t cnt = func.op_cnt[i];
            if ( cnt != 0 ) {
                total_cnt[i] += cnt;
                double avg = double(cnt) / double(it->second.call_cnt);
                uint64_t scaled_cnt = double(cnt) * scale_factor + 0.5;
                fprintf( out, "    %-40s: %8.1f/call   %10lld total   %10lld scaled_total\n", Cordic<T,FLT>::op_to_str( i ).c_str(), avg, cnt, scaled_cnt );
                csv << "\"" << Cordic<T,FLT>::op_to_str( i ) << "\", " << avg << ", " << cnt << ", " << scaled_cnt << "\n";
            }
        }
    }

    //--------------------------------------------------------
    // And the totals.
    //--------------------------------------------------------
    fprintf( out, "\n\nTotals:\n" );
    for( uint32_t i = 0; i < OP_cnt; i++ )
    {
        if ( total_cnt[i] != 0 ) {
            uint64_t cnt = total_cnt[i];
            uint64_t scaled_cnt = double(cnt) * scale_factor + 0.5;
            fprintf( out, "    %-40s:  %10lld   %10lld\n", Cordic<T,FLT>::op_to_str( i ).c_str(), cnt, scaled_cnt );
            csv << "\"" << Cordic<T,FLT>::op_to_str( i ) << "\", " << cnt << ", " << scaled_cnt << "\n";
        }
    }
    
    fclose( out );
    csv.close();
    std::cout << "\nWrote stats to " + base_name + ".{out,csv}\n";
}

#endif
