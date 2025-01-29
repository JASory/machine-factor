# Machine Factor

Relatively fast factorisation library for n < 2^128. Approximately 3 times faster than GNU Factor.

Permits constant-time evaluation, although this is impractical for 128-bit integers. 

Provides a C-style api. 

Note that this library uses Pollard-rho with deterministic parameters, so factorisation may loop
indefinitely although this is extremely improbable. A semiprime N will fail at a rate of 0.5^x where x is the period of a xorshift seeded with N, 
random computer errors and hardware failure are far more likely.

The factors are output in the form of an array of prime factors and an array of powers with the length. 

This library primarily exists as a system library, so fancy formatting is ignored.
