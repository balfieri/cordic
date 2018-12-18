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
// analyze.cpp - simple main program that uses Analysis.h
//
//      zcat xxx.log.gz | analyze
//
#include "Analysis.h"

using T   = int64_t;
using FLT = double;

int main( int argc, const char * argv[] )
{
    if ( argc < 3 ) {
        std::cout << "usage: analyze <base_name> <scale_factor> <funcs to ignore>\n";
        exit( 1 );
    }
    std::string base_name = argv[1];
    double      scale_factor = std::atof( argv[2] );
    std::vector<std::string> ignore_funcs;
    for( int i = 3; i < argc; i++ ) 
    {
        std::string ignore_name = argv[i];
        ignore_funcs.push_back( ignore_name );
    }
    auto a = new Analysis<T,FLT>( base_name );
    a->print_stats( scale_factor, ignore_funcs );
}
