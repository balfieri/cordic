<p>
This repository contains some C++ code that shows how to implement CORDIC math. <b>See Cordic.h.</b>
</p>

<p>
<b>This is all open-source.  Refer to the LICENSE.md for licensing details.
This code is intended for tutorial purposes only. You should not assume that the library is bug-free or accurate enough
for use as a production reference model.</b>
</p>

<p>
Addition, subtraction, and multiplication are relatively easy tasks to implement in a computer chip.  Divide, sqrt(), and other
transcendental functions are challenging even when high precision is not required.  CORDIC math makes it easy and cheap to implement
these complicated math functions by using a simple sequence of shifts and adds.  The only complication is that inputs need to 
be reduced to a small range of values, typically 0.0 .. 2.0.
</p>

<p>
CORDIC was invented in 1956 by Jack Volder using mathematics from 1624 and 1771.  CORDIC, which stands for for 
<b>CO</b>ordinate <b>R</b>otation <b>DI</b>gital <b>C</b>omputer
was first used in the navigation system of the B-58 bomber, which used an early digital computer
that did even have a multiply instruction.  Prior to that, navigation systems were done using
analog circuitry.  So CORDIC is over 60 years old and all been done before.  The hope of
this library is to create a small, yet complete, tutorial package.
</p>

<p>
The library assumes that values are stored as fixed-point with user-defined integer width (int_w) and fraction width (frac_w).  
The fixed-point container type T must be a signed integer at least as wide as 1+int_w+frac_w. A fixed-point number stores
the sign in the most-significant bit, followed by the int_w binary integer bits, followed by the frac_w binary fraction bits 
in the least-significant
bits.  If T is larger than the required number of bits, the extra upper bits are assumed to contain replications of the sign bit
(1=negative, 0=non-negative).  In other words, fixed-point values are stored in 2's-complement integer containers, so -A == ~A + 1.
</p>

<p>
Here are some examples of fixed-point numbers.  1.3.8 means 1 sign bit, 3 integer bits, and 8 fraction bits.
</p>
<pre>
Format         value            binary (spaces added for readability)
---------------------------------------------------------------------
1.3.8           0.0             0 000 00000000
1.3.8           1.0             0 001 00000000
1.3.8           2.0             0 010 00000000
1.3.8           0.5             0 000 10000000
1.3.8           2<sup>-8</sup>             0 000 00000001      
1.3.8           -1.0            1 111 00000000          ~(0 001 00000000) + 1  == (1 110 11111111) + 1
</pre>

<p>
This code automatically performs appropriate argument range reductions and post-CORDIC adjustments, 
but it can be turned off if you know that the inputs to the math functions are in the proper range already.
</p>

<p>
The CORDIC routines perform frac_w iterations in order to arrive at frac_w precision.  You may, however,
reduce the number of iterations by passing a different value for n to the constructor.  The numeric encoding
will be the same, but the result will be less precise.  Again, two different Cordic() instances with idential parameters
except for do_reduce and n can be used on an operation-by-operation basis.
</p>

<p>
At some point in the near future, the library will be enhanced to allow T to hold IEEE floating-point numbers with
arbitrary exponent width (exp_w). So the library will allow T (still an integer container) to hold 
fixed-point OR floating-point encodings including, but not
limited to, standard IEEE floating-point formats such as float (fp32), double (fp64), quadruple (fp128), half (fp16), and quarter (fp8).
Only one of int_w or exp_w may be non-zero.  In IEEE floating-point numbers, int_w is 0 because the binary fraction is 
assumed to have an implied '1' before the 
binary point, which is known as a normalized number (e.g., 1.110110<sub>2</sub> * 2<sup>24</sup>).  One exception is when
the value is less than the smallest normalized number 1.0 * 2<sup>MIN_EXP</sup>, which makes it a denormalized number or "denorm." Other special
numbers include +Infinity (e.g., from 1/0), -Infinity (e.g., from -1/0), or NaN (not a number, e.g., from sqrt(-1)), which 
are identified using special encodings of the exponent.
</p>

<p>
Fixed-point numbers naturally support denorms (i.e., they are all denorms), but have no way to indicate a value 
outside their allowed range.  The library needs
an option to mark a number as +Infinity, -Infinity, or NaN.  An additional needed option is
to gracefully flush large numbers to +/- "max value" and NaNs to zero.
</p>

<p>
Also not handled are the different IEEE rounding modes: round-to-nearest (current behavior 
and what one normally expects), 
round-toward-zero, round-toward-plus-infinity, round-toward-minus-infinity.  There's another one worth doing
called round-away-from-zero (aka round-toward-plus-or-minus-infinity). The rounding mode would be set
during the constructor and used for all operations.  These modes will be implemented for both fixed-point and floating-point
encodings.
</p>

<p>
There are a couple other options supported in most floating-point libraries that we'll need to add for
floating-point encodings only.
Flush-To-Zero (FTZ) means that any time an operation would produce a denom (again, a number smaller than the smallest normalized
number), it is changed to zero.
Denorm-As-Zero (DAZ) means that any denorm input to an operation is first changed to zero.
</p>

<p>
To install this on your computer, you'll need git and a C++ compiler, then:
</p>
<pre>
git clone https://github.com/balfieri/cordic
cd cordic
</pre>

<p>
To build and run the basic "smoke" test, <b>test_basic.cpp</b>, on Linux, Cygwin, or macOS:
</p>
<pre>
doit.test
doit.test 1                             - run with debug spew 
doit.test 0 test_basic                  - same as doit.test with no args
doit.test 1 test_basic -new_bugs        - run new bugs too (those not fixed yet)
doit.test 0 test_basic -int_w 8         - change int_w from default to 8 bits
</pre>

<p>
test_basic.cpp does its own checking using macros in test_helpers.h.  In the near future, 
Cordic.cpp should do optional checking of computations so that test_helpers.h can be deleted or greatly simplified.
</p>

<p>
The library needs an "freal" flexible real number class that follows all the rules of any C++ floating-point number but uses Cordic 
as its underlying implementation.  freal should keep track of ranges and automatically change the underlying fixed-point and floating-point 
representations, ideally mixing fixed-point and floating-point.

<p>
There should never be a need for this library to do something special for complex numbers. The C++ complex<> 
class and associated math functions fall out naturally by using complex&lt;freal&gt; or whatever.
</p>

<p>
Bob Alfieri<br>
Chapel Hill, NC
</p>
