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
// test_cordic.cpp - test Cordic math
//
#include "Cordic.h"

using FLT = double;
using FP  = int64_t;
constexpr int64_t FW = 56;
constexpr FLT     TOL = 1.0 / FLT( 1LL << (FW/4) ); 

FP to_fp( FLT x )
{
    return FP( x * FLT(FP(1) << FP(FW)) );
}

FLT to_flt( FP x )
{
    return FLT( x ) / FLT(FP(1) << FP(FW));
}

int main( int argc, const char * argv[] )
{
    Cordic<FP, FW> cordic;
    FLT fltx;
    FLT flty;
    FP  fpx;
    FP  fpy;
    FP  fpz;
    FP  fpz0;
    FP  fpz1;
    FLT fltz;
    FLT fltz0;
    FLT fltz1;
    FLT flte;
    FLT flte0;
    FLT flte1;
    FLT flterr;
    FLT flterr0;
    FLT flterr1;

    std::cout << "tol: " << TOL << "\n";

    // x*y
    //
    fltx = 0.681807431807431031;
    flty = 0.810431798013170871;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    fpz  = cordic.mul( fpx, fpy );
    fltz = to_flt( fpz );
    flte = fltx*flty;
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nx*y\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(x), " << flty << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // y/x
    //
    fltx = 0.681807431807431031;
    flty = 0.810431798013170871;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    fpz  = cordic.div( fpy, fpx );
    fltz = to_flt( fpz );
    flte = flty/fltx;
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\ny/x\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(x), " << flty << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // sqrt( x )
    fltx = 0.756728943106177373;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    fpz  = cordic.sqrt( fpx );
    fltz = to_flt( fpz );
    flte = std::sqrt( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nsqrt(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // e^x 
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.exp( fpx );
    fltz = to_flt( fpz );
    flte = std::exp( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nexp(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // x^y
    fltx = 0.456728943106177373;
    flty = 0.790843170754708943;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    fpz  = cordic.pow( fpx, fpy );
    fltz = to_flt( fpz );
    flte = std::pow( fltx, flty );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\npow(x, y)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(x), " << flty << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // 2^x 
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.pow2( fpx );
    fltz = to_flt( fpz );
    flte = std::pow( 2.0, fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\npow2(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // 10^x 
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.pow10( fpx );
    fltz = to_flt( fpz );
    flte = std::pow( 10.0, fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\npow10(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    //dassert( flterr <= TOL ); // TODO

    // log(x)
    fltx = 1.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.log( fpx );
    fltz = to_flt( fpz );
    flte = std::log( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nlog(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // logb(x, b)
    fltx = 0.456728943106177373;
    flty = 1.790843170754708943;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    fpz  = cordic.logb( fpx, fpy );
    fltz = to_flt( fpz );
    flte = std::log( fltx ) / std::log( flty );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nlogb(x, b)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(x), " << flty << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // log2(x)
    fltx = 1.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.log2( fpx );
    fltz = to_flt( fpz );
    flte = std::log( fltx ) / std::log( 2.0 );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nlog2(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // log10(x)
    fltx = 1.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.log10( fpx );
    fltz = to_flt( fpz );
    flte = std::log( fltx ) / std::log( 10.0 );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nlog10(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // sin() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.sin( fpx );
    fltz = to_flt( fpz );
    flte = std::sin( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nsin(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // cos() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.cos( fpx );
    fltz = to_flt( fpz );
    flte = std::cos( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\ncos(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // sin( x ), cos( x )
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    cordic.sin_cos( fpx, fpz0, fpz1 );
    fltz0 = to_flt( fpz0 );
    fltz1 = to_flt( fpz1 );
    flte0 = std::sin( fltx );
    flte1 = std::cos( fltx );
    flterr0 = std::abs( flte0-fltz0 );
    flterr1 = std::abs( flte1-fltz1 );

    std::cout.precision(24);
    std::cout << "\nsin_cos(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx  << "\n";
    std::cout << "Expected: " << std::setw(30) << flte0 << ", " << flte1 << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz0 << ", " << fltz1 << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr0 << ", " << flterr1 << "\n";
    dassert( flterr0 <= TOL );
    dassert( flterr1 <= TOL );

    // tan() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.tan( fpx );
    fltz = to_flt( fpz );
    flte = std::tan( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\ntan(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // asin() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.asin( fpx );
    fltz = to_flt( fpz );
    flte = std::asin( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nasin(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // acos() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.acos( fpx );
    fltz = to_flt( fpz );
    flte = std::acos( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nacos(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // atan() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.atan( fpx );
    fltz = to_flt( fpz );
    flte = std::atan( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\natan(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );
    
    // atan2() 
    //
    fltx = 0.456728943106177373;
    flty = 0.709831990704326039;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    fpz  = cordic.atan2( fpy, fpx );
    fltz = to_flt( fpz );
    flte = std::atan2( flty, fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\natan2(y, x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(x), " << flty << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // norm( x, y ) = sqrt( x^2 + y^2 )
    fltx = 0.956728943106177373;
    flty = 0.708473170947310947;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    fpz  = cordic.norm( fpx, fpy );
    fltz = to_flt( fpz );
    flte = std::sqrt(fltx*fltx + flty*flty);
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nnorm(x, y)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(x), " << flty << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // normh( x, y ) = sqrt( x^2 - y^2 )
    fltx = 0.708473170947310947;
    flty = 0.556728943106177373;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    fpz  = cordic.normh( fpx, fpy );
    fltz = to_flt( fpz );
    flte = std::sqrt(fltx*fltx - flty*flty);
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nnormh(x, y)\n";
    std::cout << "Input:    " << std::setw(30) << fltx  << "(x), " << flty  << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // rect_to_polar
    fltx = 0.456728943106177373;
    flty = 0.708473170947310947;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    cordic.rect_to_polar( fpx, fpy, fpz0, fpz1 );
    fltz0 = to_flt( fpz0 );
    fltz1 = to_flt( fpz1 );
    flte0 = sqrt(fltx*fltx + flty*flty);
    flte1 = atan2( flty, fltx );
    flterr0 = std::abs( flte0-fltz0 );
    flterr1 = std::abs( flte1-fltz1 );

    std::cout.precision(24);
    std::cout << "\nrect_to_polar(x, y)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(x), " << flty << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte0 << ", " << flte1 << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz0 << ", " << fltz1 << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr0 << ", " << flterr1 << "\n";
    dassert( flterr0 <= TOL );
    dassert( flterr1 <= TOL );

    // polar_to_rect
    fltx = 0.456728943106177373;
    flty = 0.708473170947310947;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    cordic.polar_to_rect( fpx, fpy, fpz0, fpz1 );
    fltz0 = to_flt( fpz0 );
    fltz1 = to_flt( fpz1 );
    flte0 = fltx*cos(flty);
    flte1 = fltx*sin(flty);
    flterr0 = std::abs( flte0-fltz0 );
    flterr1 = std::abs( flte1-fltz1 );

    std::cout.precision(24);
    std::cout << "\npolar_to_rect(r, a)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(r), " << flty << "(a)\n";
    std::cout << "Expected: " << std::setw(30) << flte0 << ", " << flte1 << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz0 << ", " << fltz1 << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr0 << ", " << flterr1 << "\n";
    dassert( flterr0 <= TOL );
    dassert( flterr1 <= TOL );

    // sinh() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.sinh( fpx );
    fltz = to_flt( fpz );
    flte = std::sinh( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nsinh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // cosh() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.cosh( fpx );
    fltz = to_flt( fpz );
    flte = std::cosh( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\ncosh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // sinh( x ), cosh( x )
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    cordic.sinh_cosh( fpx, fpz0, fpz1 );
    fltz0 = to_flt( fpz0 );
    fltz1 = to_flt( fpz1 );
    flte0 = std::sinh( fltx );
    flte1 = std::cosh( fltx );
    flterr0 = std::abs( flte0-fltz0 );
    flterr1 = std::abs( flte1-fltz1 );

    std::cout.precision(24);
    std::cout << "\nsinh_cosh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx  << "\n";
    std::cout << "Expected: " << std::setw(30) << flte0 << ", " << flte1 << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz0 << ", " << fltz1 << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr0 << ", " << flterr1 << "\n";
    dassert( flterr0 <= TOL );
    dassert( flterr1 <= TOL );

    // tanh() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.tanh( fpx );
    fltz = to_flt( fpz );
    flte = std::tanh( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\ntanh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // asinh() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.asinh( fpx );
    fltz = to_flt( fpz );
    flte = std::asinh( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nasinh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // acosh() 
    //
    fltx = 1.591370341781322;
    fpx  = to_fp( fltx );
    fpz  = cordic.acosh( fpx );
    fltz = to_flt( fpz );
    flte = std::acosh( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\nacosh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // atanh() 
    //
    fltx = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpz  = cordic.atanh( fpx );
    fltz = to_flt( fpz );
    flte = std::atanh( fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\natanh(x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    // atanh2() 
    //
    fltx = 0.709831990704326039;
    flty = 0.456728943106177373;
    fpx  = to_fp( fltx );
    fpy  = to_fp( flty );
    fpz  = cordic.atanh2( fpy, fpx );
    fltz = to_flt( fpz );
    flte = std::atanh( flty/fltx );
    flterr = std::abs( flte-fltz );

    std::cout.precision(24);
    std::cout << "\natanh2(y, x)\n";
    std::cout << "Input:    " << std::setw(30) << fltx << "(x), " << flty << "(y)\n";
    std::cout << "Expected: " << std::setw(30) << flte << "\n";
    std::cout << "Actual:   " << std::setw(30) << fltz << "\n";
    std::cout << "Error:    " << std::setw(30) << flterr << "\n";
    dassert( flterr <= TOL );

    return 0;
}
