<p>
Addition, subtraction, and multiplication are relatively easy tasks to implement in a computer chip.  Divide, sqrt(), and other
transcendental functions are challenging even when high precision is not required.  CORDIC math makes it easy and cheap to implement
these complicated math functions by using a simple sequence of shifts and adds.  The only complication is that inputs need to 
be reduced to a small range of values, typically 1.0 .. 2.0.  
</p>

<p>
This repository contains some C++ code that shows how to implement CORDIC math. <b>See Cordic.h.</b>
The library currently assumes that values are stored as fixed-point with user-defined integer width (int_w) and fraction width (frac_w).  
The fixed-point container type T must be a signed integer at least as wide as 1+int_w+frac_w.  [A fixed-point number stores
the sign in the most-significant bit, followed by the int_w integer bits, followed by the frac_w fraction bits in the least-significant
bits.  If T is larger than the required number of bits, the extra upper bits are assumed to contain replications of the sign 
(1=negative, 0=non-negative).]
</p>

<p>
By default, this code automatically performs appropriate argument range reductions and post-CORDIC adjustments, 
but it can be turned off if you know that the inputs to the math functions are in the proper range already.
</p>

<p>
At some point in the near future, the library will be enhanced to allow T to hold IEEE floating-point numbers with
arbitrary exponent width (exp_w). So the library will allow T (still an integer container) to hold 
fixed-point OR floating-point encodings including, but not
limited to, standard IEEE floating-point formats such as float (fp32), double (fp64), half (fp16), and quarter (fp8).
Only one of int_w or exp_w may be non-zero.  In IEEE floating-point numbers, the integer part is assumed to be an implicit '1' unless
the entire value is 0 or a special number such as +Infinity, -Infinity, or NaN (not a number). So int_w is 0 for floating-point numbers.
Conversely, fixed-point numbers naturally have no exponent.
</p>

<p>
CORDIC was invented in the 1950's and first used in the navigation system of the B-38 bomber.  Prior to that,
navigation systems used analog computers.
There is nothing novel or proprietary here.  This is intended
for tutorial purposes only. You should not assume that the library is bug-free or accurate enough
for use as a production reference model.
</p>

<p>
<b>This is all open-source.  Refer to the LICENSE.md for licensing details.</b>
</p>

<p>
To build and run the test, <b>test_cordic.cpp</b>, on Linux, Cygwin, or macOS:
</p>
<pre>
doit.test
</pre>

<p>
Currently, test_cordic.cpp does its own checking within a small tolerance.  In the near future, I'd like to have
the Cordic.cpp library itself do optional checking of computations.
</p>

<p>
We need to add a "cordreal" real number class that follows all the rules of any C++ floating-point number but uses Cordic 
as its underlying implementation. Note that there will be no need to do anything special for complex numbers.  C++ complex<> 
class and associated math functions fall out naturally by using complex<cordreal>.
</p>

<p>
Bob Alfieri<br>
Chapel Hill, NC
</p>
