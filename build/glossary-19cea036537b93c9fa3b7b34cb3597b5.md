# Glossary

```{glossary}
automatic differentiation
: Accurate numerical computation of a derivative through application of standard differentiation rules, as implemented by `ForwardDiff` and other packages.

comprehension
: Construction of an array by an implied loop.

destructure
: Assign multiple names to a tuple value in order to break it into its constituent parts.

export
: Module declaration that makes a name available without {term}`namespace` qualification.

generator
: Object that can generate iterable values on demand, like the interior of a {term}`comprehension`.

interpolate
: Insert a defined value, as in a string; invoked by the use of `$`.

macro
: Transformation of code applied at interpretation/compile time; macro names all start with `@`.

mutating
: Changing the value or state of one or more of input arguments; by convention, mutating functions have `!` at the end of their name.

namespace
: Compartmentalization of defined names and symbols into modules, in order to avoid conflicts.

short-circuiting
: Terminating the evaluation of a logical AND `&&` or OR `||` as soon as its final result can be deduced.

splatting
: The syntax `x...`, used to indicate that all the elements of the iterable `x` are to be included at that spot, as though they had been typed out in place.

symbol
: A type of named constant value starting with `:`.

ternary
: The syntax `<predicate> ? <if true> : <if false>` for conditional evaluation.

tuple
: An immutable sequence of values, typically created by enclosing the values in parentheses and separating them by commas.

vectorize
: Apply a function that expects a scalar argument to all the elements of an array.
```