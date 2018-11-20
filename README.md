<p>
This repository contains some C++ code that shows how to implement CORDIC math.  
The library currently assumes that values are stored as fixed-point with user-defined integer width (INT_W) and fraction width (FRAC_W).  
The fixed-point container type T must be a signed integer at least as wide as 1+INT_W+FRAC_W.
At some point in the future, the library will be enhanced to allow T to represent IEEE floating-point numbers with
arbitrary exponent width (EXP_W) and mantissa width (MANT_W).  So the library will allow
T to be fixed-point OR floating-point.
</p>

<p>
CORDIC was invented in the 1950's.  There is nothing novel or proprietary here.  This is intended
for tutorial purposes only.  That said, as the library gets tested thoroughly, it might be good
enough for production code.  It's not even close to ready for that currently.
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
Bob Alfieri<br>
Chapel Hill, NC
</p>
