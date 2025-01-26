# Machine Factor

Relatively fast factorisation library for n < 2^128. Approximately 3 times faster than GNU Factor.

Permits constant-time evaluation, although this is impractical for 128-bit integers. 

Provides a C-style api. 

Note that this library uses Pollard-rho with deterministic parameters, so factorisation may loop
indefinitely although this is extremely improbable.

The factors are output in the form of an array of elements with the format prime,power,prime2,power2, 
and the length of the array.

This library primarily exists as a system library, so fancy formatting is ignored.
