<p>
Addition, subtraction, and multiplication are relatively easy tasks to implement in a computer chip.  Divide, sqrt(), and other
transcendental functions are challenging even when high-precision is not required.  CORDIC math makes it easy and cheap to implement
these complicated math functions using a simple sequence of shifts and adds. 
</p>

<p>
This repository contains some C++ code that shows how to implement CORDIC math.  
The library currently assumes that values are stored as fixed-point with user-defined integer width (int_w) and fraction width (frac_w).  
The fixed-point container type T must be a signed integer at least as wide as 1+int_w+frac_w.
</p>

<p>
At some point in the near future, the library will be enhanced to allow T to represent IEEE floating-point numbers with
arbitrary exponent width (exp_W). So the library will allow T to be fixed-point OR floating-point.
However only one of int_w or exp_w may be non-zero.  In floating-point numbers, the integer part is assumed to be 1 unless
the entire value is 0 or a special number.
</p>

<p>
CORDIC was invented in the 1950's.  There is nothing novel or proprietary here.  This is intended
for tutorial purposes only. You should not assume that the library is bug-free or accurate enough
for use as a production reference model.
</p>

<p>
This is all open-source.  Refer to the LICENSE.md for licensing details.
</p>

<p>
To build and run the test (test_cordic.cpp) on Linux, Cygwin, or macOS:
</p>
<pre>
doit.test
</pre>

<p>
Besides thorough testing, I would like to add support for complex numbers and associated math operations,
taking advantage of known identities where possible to optimize the functions.
</p>

<p>
Bob Alfieri<br>
Chapel Hill, NC
</p>
