# Approximation
Implemented 1D function approximation in C++ 

Overall, files have functions for interpolation(Lagrange), approximation(with polinome, chebishow and exp()) 
and natural cubic spline interpolation.

APPROXIMATION

Let's say we have a function, but its shape is too inconvenient(hard to find derivative or to count). 
Than you can try to compute a simple function Q{n}(x), that has similar values to f(x). In simple terms:

f(x) ≈ Q{n}(x) = sum(c{i}*phi{i}(x)), i from 0 to n

Approximation allowes to find coeficients c{i} for different set of phi{i}(x)

INTERPOLATION

Allowes to compute a value of x on a table-given function. Interpolation is more time-consuming than approximation, and can give only a meaning of a function in one dot at a time - not a simplified version of a function. However, it gives much more accurate results.

More detailed about concepts: https://en.wikipedia.org/wiki/Approximation, https://en.wikipedia.org/wiki/Interpolation
