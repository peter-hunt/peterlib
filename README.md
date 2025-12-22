# PeterLib

A common framework library for Python functions like Computer Algebra System (CAS), CLI engine for CLI interface and RPG, RPG profile manager template, Pydantic/dataclass styled data type base class for easier data management and JSON loading/dumping, and future contents to come.

## Usage

The code in the project is not organized to be executed directly. There are some that require dependency libraries listed in the `requirements.txt`. Example usages might be added in the future.

## Table of Content

- `algebra.py`: Python computer algebra system (CAS).
- `cliengine.py`: CLI engine for CLI interface and RPG.
- `datatype.py`: Pydantic/dataclass styled data type base class for easier data management and JSON loading/dumping.
- `profile_manage.py` and `profile_template`: profile manager with working folder management with `saves` and `settings.json`. The profile class template is for functionality integration with the management code.

## Computer Algebra System

This is an experimental/educational project to recreate a computer algebra system to perform algebra and mathematical tasks automatically, like on graphing calculators. The CAS part is motivated by a potential project in generating and checking math practice problems of up to around high school level with the CLIEngine and profile managing code.

Here are the current structures and functionalities supported or planned:

### Current Structures
- Expr class to have all other CAS classes as subclasses;
- Number literals. Complex numbers are supported but not much suggested;
- Single-letter `Var`. See designer notes below for more information;
- Mathematical `Constant`s, like `pi` and `e`, both available kept as themselves and evaluation for opened options and usage;
- `Add`, `Mul`, and `Pow` nodes, where subtraction and division are represented with those only but printed properly. See designer notes below for more information;
- `Function` structures for one and two variable functions, where currently implemented are abs, trigs, arctrigs, ln/log.
- `Limit`, `Derivative`, and `Integral` structures, although mostly not yet functional.

### Current Functionalities
- Constructing corresponding classes when python operation are used on Expr subclass, done with operator methods;
- Template for recursive functionality implementation, implementation/structure-planning on following functions:
- - `doit(expr)`: Evaluating evaluatable limit/derivative/integral structures and plugging in values for known constant radiants in trigs (this part not implemented yet);
- - `subs(expr, expr_map)`: Applying substitution without modifying the corresponding values or the result expression;
- - `evalf(expr, value_map)`: Evaluating the number value of the expression with the given values of the related variables, or partially so if not all variables are given. See designer notes below for more information;
- - `expand(expr)/expand_xxx(expr)`: Basic expansion of `Mul`, `Pow`, `Ln`/`Log`, trigs with both separated and all-in-one functions;
- - `diff(expr, var, order=0, evaluate=True)`: Differentiation either in evaluating or constructing a `Derivative` instance with given variable and order. Equivalent to `Derivative(expr, var, order).doit()`. This doesn't work for `Abs` yet since piecewise function class is missing.
- - `integrate(expr, var, a=None, b=None)`: Integration by evaluating an `Integral` instance with given variable and optional bounds. Equivalent to `Integral(expr, var, a, b).doit()`. This only works for monomial/polynomial with basic rules.

### Planned Structures
- Functions: factorial, ceil, floor, sgn, modulo, exp (construction only).
- Piecewise functions with binary search optimization for given values;
- Equations and inequality structures. See designer notes below for more information.

### Planned Functionalities
- Limit and integral evaluation with function to construct or evaluate;
- The series of factor/combine functions to factor `Add`, `Mul` and combine `Ln`/`Log`, trigs and so on;
- Fraction operations:
- - `cancel()`: Canceling factors in divisions. Optionally ignoring the x/x -> 0 where x might be 0.
- - `together()`: Combining factors with common denominator.
- - `apart()`: Partial fraction decomposition.
- - `collect()`: Grouping the expression, or mostly polynomials, by the powers of one variable.
- Solving equations and inequalities;
- Implementing `doit()` functionalities for trig constant values.
- Implementing `simplify()` after implementing the factor/combine and all the fraction ones.
- Use binomial coefficients or similar combinatorics to optimize power of sum distribution.
- Make multiplication print out more nicely for letters and constants. Potentially flag options.

### Designer Notes
- For the variable names, an option is to allow especially subscripts to support variable names and literals like A_1 for more variable spaces, but that will make the command line printing more complicated and not support compact multiplication like writing `abc` in place of `a * b * c` as easily. And especially as an experimental/educational project, this will only support single letter variables from latin and greek letters that are not mathematically used yet for now;
- For implementing and supporting equations and inequalities, it would make sense to make the python comparison symbols trigger creating corresponding instance. However, for the sake of refering to checking expression equality and python not having internal `===` other than `is`, which doesn't really work with artificially immutable instances, the default `==` will be reserved for checking for direct equality after simplication, for example. In addition, to reduce computational burden and avoid potential confusion, the equality check by `==` will not automatically call on `simplify` or related methods;
- Although I haven't found an instance that leads to this, the currently planned logic for `expand()` with applying all the basic expansion functions might lead to different results for different expanding order, which could be avoided by changing the order or ensuring that the expansion logic won't lead to separate local results; the same applies for other methods that might return unexpectedly different results. Just for future reference;
- With the current `evalf()` implementation, the value map allows missing variables or even with empty value. This is by design since this would potentially allow faster evaluation specified on the variable(s) left out, and it will also support the passing of some partial value map that isn't related to this specific expression;
- For the factoring functionality, I would love to try to design and implement a general usage one that works for all polynomials and so on. However, with computers not being able to recognize general patterns in a humane way, factoring looks more like solving for the roots and writing it as the solution form. Therefore, the planned implementation is to use quadratic up to quartic formulas, simplified after plugging in, from the original polynomial or arithmetic series of powers with up to four terms, to arrive to the factoring of polynomials. Note that this would require a lot of progress in `simplify` and more to work for even simpler polynomials to return sane expressions, but this approach will ensure a rather more compatible factoring formula, even at the cost of not recognizing patterns like `x^n-1` always factorable by `x-1` (which is still implementable as exceptions but the logic wouldn't be learned).

## Command Line Engine (CLIEngine)
A work-in-progress engine for command line interactions for applications. The engine handles the command parsing, basic argument types, and fitting the command to execute the corresponding functions. More documentation will be added as the functionalities refine.
