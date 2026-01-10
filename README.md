# PeterMath

![License](https://img.shields.io/github/license/peter-hunt/petermath)
![GitHub repo size](https://img.shields.io/github/repo-size/peter-hunt/petermath)
![Last Commit](https://img.shields.io/github/last-commit/peter-hunt/petermath)
![GitHub Repo stars](https://img.shields.io/github/stars/peter-hunt/petermath)


A math library for Python functions like Computer Algebra System (CAS) for math expressions, set theory and number domains, equation/inequality solvers, and future contents to come.

## Usage

The code in the project is not organized to be executed directly. There are some that require dependency libraries listed in the `requirements.txt`. Example usages might be added in the future. See documentations and test code for example usages for now.

If used in a folder, put `.` in front of the relative imports, e.g. `utils.py`, to have it work properly. The `__init__.py` is already implemented to contain all content and work with `import *`.

## Table of Content

- `__init__.py`: Compilation of all the content in the library. Note that `polynomial` is not added to the `__init__.py` and must be imported separately like `import petermath.polynomial` to avoid confusion with the `Expr` structures;
- `expr.py`: Python computer algebra system (CAS) for mathematical expressions;
- `sets.py`: (Work in Progress) Set theory logic for number domains, set operations, solution sets and so on;
- `statement.py`: (Work in Progress) Types for equations, inequalities, diffeq, etc.;
- `polynomial.py`: (Work in Progress) Data structure optimized for single-variable natural-number-power preferably-rational-coefficients polynomials for solving and factoring;
- `solve.py`: (Work in Progress) Solvers for equations, inequalities, diffeq, etc.;
- `xxx_test.py`: Relevant test code for each module and for example usage;
- `xxx_utils.py`: Utility code based on the original implement but not required for the original;
- `utils.py`: Common code for all other usages.

## Computer Algebra System

This is an experimental/educational project to recreate a computer algebra system to perform algebra and mathematical tasks automatically, like on graphing calculators. The CAS part is motivated by a potential project in generating and checking math practice problems of up to around high school level with a CLI Engine and profile managing code.

Here are the current structures and functionalities supported or planned:

### Current Structures

- `Expr` superclass for CAS classes;
- Number literals. Complex numbers are supported but not much suggested;
- Single-letter `Var`. See designer notes below for more information;
- Mathematical `Constant`s:
  - Regular constants like `pi` and `e`, both available kept as themselves and evaluation for opened options and usage;
  - Special constants like `inf`, `-inf`, and `undefined` for intervals, invalid operations, and so on. (Work in Progress)
- `Add`, `Mul`, and `Pow` nodes, where subtraction and division are represented with those only but printed properly. See designer notes below for more information;
- `Function` structures for one and two variable functions, where currently implemented are abs, trigs, arctrigs, ln/log.
- `Limit`, `Derivative`, and `Integral` structures, although mostly not yet functional;
- `Set` superclass;
- Set operation structures: `Union`, `Intersection`, `Difference`, and `Complement`;
- `UniversalSet` and `EmptySet`;
- `FiniteSet` for finite element evaluations;
- `NumberDomain` superclass for natural, integer, rational, real, and complex sets;
- `ConstInterval` for constant interval to enable denoting intervals for solutions.

### Current Functionalities

- Constructing corresponding classes when python operation are used on Expr subclass, done with operator methods;
- Supporting Fraction interactions but not Decimal. See designer notes below for more information;
- Automatically casting python values for `inf`, `-inf`, and `nan` to functional corresponding `Constant` values;
- Template for recursive functionality implementation, implementation/structure-planning on following functions:
- - `doit(expr)`: Evaluating evaluatable limit/derivative/integral structures and plugging in values for known constant radiants in trigs (this part not implemented yet);
- - `subs(expr, expr_map)`: Applying substitution without modifying the corresponding values or the result expression;
- - `evalf(expr, value_map)`: Evaluating the number value of the expression with the given values of the related variables, or partially so if not all variables are given. See designer notes below for more information;
- - `expand(expr)/expand_xxx(expr)`: Basic expansion of `Mul`, `Pow`, `Ln`/`Log`, trigs with both separated and all-in-one functions;
- - `diff(expr, var, order=0, evaluate=True)`: Differentiation either in evaluating or constructing a `Derivative` instance with given variable and order. Equivalent to `Derivative(expr, var, order).doit()`. This doesn't work for `Abs` yet since piecewise function class is missing.
- - `integrate(expr, var, a=None, b=None)`: Integration by evaluating an `Integral` instance with given variable and optional bounds. Equivalent to `Integral(expr, var, a, b).doit()`. This only works for monomial/polynomial with basic rules;
- Expanded `exprhash` as `sethash` with each `Set`s added.

### Planned Structures

- An implementation of sets to support and organize statements and solutions in a much more systematic manner;
- `Statement` class with equality and inequalities:
  - Equality, not to be confused with Equation.See designer notes below for more information ;
  - Inequality structures. See designer notes below for more information.
- `SetBuilder` with separate code since it requires the `Statement` class, which requires `Set` class for the number domains.
- Contraints like domain sets and intervals to feed into statements;
- Functions: factorial, ceil, floor, sgn, modulo, exp (construction only).
- Piecewise functions with binary search optimization for given values;
- `Sum` and `Prod` for discrete operations.

### Planned Functionalities

- Variable class, expanded or separate, for sets;
- Simplifying nested set operators;
- Fixing the hash tuple sorting in a more systematic way;
- Checking for subset and superset relationships;
- Evaluating set operations on constant intervals;
- Backward supporting inputting `Set` in some `Expr` functions that would make sense;
- Supporting set theory concepts;
- Solving equations and inequalities;
- Limit and integral evaluation with function to construct or evaluate;
- The series of factor/combine functions to factor `Add`, `Mul` and combine `Ln`/`Log`, trigs and so on;
- Fraction operations:
- - `cancel()`: Canceling factors in divisions. Optionally ignoring the x/x -> 0 where x might be 0.
- - `together()`: Combining factors with common denominator.
- - `apart()`: Partial fraction decomposition.
- - `collect()`: Grouping the expression, or mostly polynomials, by the powers of one variable.
- Implementing `doit()` functionalities for trig constant values.
- Implementing `simplify()` after implementing the factor/combine and all the fraction ones.
- Use binomial coefficients or similar combinatorics to optimize power of sum distribution.
- Make multiplication print out more nicely for letters and constants. Potentially flag options.

### Designer Notes

- For the variable names, an option is to allow especially subscripts to support variable names and literals like A_1 for more variable spaces, but that will make the command line printing more complicated and not support compact multiplication like writing `abc` in place of `a * b * c` as easily. And especially as an experimental/educational project, this will only support single letter variables from latin and greek letters that are not mathematically used yet for now;
- For the sake of avoiding confusion and reserving the use of python comparison operators, and since python also doesn't have an internal `===` operation other than `is` (which doesn't really work with artificially immutable instances), builtin operations like `==` and `>` will be reserved to compare exprhash values for ordering (which will not call on methods like `simplify()`) instead of creating (in)equality instances;
- As a statement, the equal sign structure will be called `Equality` instead of `Equation` as it serves more as a equality comparison for the set theory logic instead of a solvable construction, which will be implemented later with different solving algorithms and the basics of set theory implemented to be able to have solution sets;
- More on inequalities, this will be implemented later with domains, like the integer, rational, real, and complex sets, and constraints, like `x` being positive or negative, that would simplify the process of solving equations and allow solving or simplifying inequalities.
- Although I haven't found an instance that leads to this, the currently planned logic for `expand()` with applying all the basic expansion functions might lead to different results for different expanding order, which could be avoided by changing the order or ensuring that the expansion logic won't lead to separate local results; the same applies for other methods that might return unexpectedly different results. Just for future reference;
- With the current `evalf()` implementation, the value map allows missing variables or even with empty value. This is by design since this would potentially allow faster evaluation specified on the variable(s) left out, and it will also support the passing of some partial value map that isn't related to this specific expression;
- For the factoring functionality, I would love to try to design and implement a general usage one that works for all polynomials and so on. However, with computers not being able to recognize general patterns in a humane way, factoring looks more like solving for the roots and writing it as the solution form. Therefore, the planned implementation is to use quadratic up to quartic formulas, simplified after plugging in, from the original polynomial or arithmetic series of powers with up to four terms, to arrive to the factoring of polynomials. Note that this would require a lot of progress in `simplify` and more to work for even simpler polynomials to return sane expressions, but this approach will ensure a rather more compatible factoring formula, even at the cost of not recognizing patterns like `x^n-1` always factorable by `x-1` (which is still implementable as exceptions but the logic wouldn't be learned);
- Python builtin library types for only Fraction and not Decimal is supported, because the interactions within them are will still lead to TypeErrors as python doesn't support them by default. This is overcome by casting decimal instances to fractions automatically, the same process where the python values for infinities and NaNs are casted to expr constants.

# License: MIT

[MIT License](./LICENSE.txt)
