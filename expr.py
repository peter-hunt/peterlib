from decimal import Decimal
from fractions import Fraction
from itertools import product as iterprod
from math import e as math_e, pi as math_pi, tau as math_tau, \
    sqrt, prod, sin, cos, tan, log, asin, acos, atan, \
    inf as math_inf, nan as math_nan, \
    isinf as math_isinf, isnan as math_isnan, isfinite as math_isfinite
from numbers import Number, Rational
from struct import pack
from types import FunctionType
from typing import Iterable, Union

from typeguard import typechecked

from tags import *
from utils import inherit_docstrings


__all__ = [
    "VARIABLE_LETTERS", "RESERVED_LETTERS",
    "CanonicalKey",

    "frdiv", "frpow",

    "float_hash",
    "exprhash", "exprsorted",

    "cast_number", "number_casted",
    "cast_map", "map_casted",
    "Expr", "ExprLike",
    "Var", "ExprMap", "ExprValueMap", "symbols",

    "apply",
    "doit", "subs", "evalf",

    "expand_mul", "expand_pow", "expand_distribute", "expand_trig", "expand_log",
    "expand",
    "factor_poly",
    "factor",
    "cancel", "together", "apart", "collect",
    "simplify",

    "_diff", "diff",

    "Constant",
    "e", "pi", "tau", "phi",
    "inf", "undefined",
    "is_special_const", "is_literal", "is_infinite", "is_finite", "get_value", "get_sign",
    "format_term", "Add",
    "format_factor", "Mul",
    "Pow",
    "Function", "UnaryFunction", "BinaryFunction",
    "Abs",
    "Sin", "Cos", "Tan", "Sec", "Csc", "Cot",
    "Arcsin", "Arccos", "Arctan", "Arcsec", "Arccsc", "Arccot",
    "Ln", "Log",
    "Limit", "Derivative", "Integral",

    "is_constant", "is_rat_constant",
    "is_mono", "split_mono",
    "is_poly", "is_pos_poly", "is_rat_poly", "is_perfect_poly", "split_poly",

    "integrate",
]


VARIABLE_LETTERS = (
    "abcdefghijklmnopqrstuvwxyz"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "αβγδεζηθικλμνξοπρστυφχψω"
    "ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ"
)

RESERVED_LETTERS = "eπτφΣΠ"


CanonicalKey = tuple[int | float | tuple, ...]


@typechecked
def number_casted(func: FunctionType) -> FunctionType:
    def wrapper(*args, **kwargs):
        return func(*(cast_number(arg) for arg in args),
                    **{key: cast_number(value) for key, value in kwargs.items()})
    return wrapper


@typechecked
def map_casted(func: FunctionType) -> FunctionType:
    def wrapper(self, expr_map):
        return func(self, cast_map(expr_map))
    return wrapper


@typechecked
class Expr:
    def __init__(self, *args, **kwargs):
        raise Exception(
            "Cannot create an instance of base Expr class. "
            "To use a subclass, redefine the init code."
        )

    def sqrt(self) -> Union["Expr", Number]:
        """
        Returns the square root of expression,
        equivalent to expr ** (1/2)

        :param self: The expression.
        :return: The square root of the expression.
        :rtype: Expr | Number
        """
        return Pow(self, Fraction(1, 2))

    def cbrt(self) -> Union["Expr", Number]:
        """
        Returns the cube root of expression,
        equivalent to expr ** (1/3)

        :param self: The expression.
        :return: The cube root of the expression.
        :rtype: Expr | Number
        """
        return Pow(self, Fraction(1, 3))

    def __repr__(self):
        raise NotImplementedError

    def __str__(self):
        raise NotImplementedError

    def exprhash(self) -> CanonicalKey:
        """
        Hash tuple generator for expressions.
        Use the global function for supported int/float/Fraction numbers.
        This is used for sorting and hashing expressions.

        :param self: The expression instance.
        :return: The tuple of relevant information.
        :rtype: tuple[int | tuple, ...]
        """
        raise NotImplementedError

    def __hash__(self):
        return hash(self.exprhash())

    def __eq__(self, other: Union["Expr", Number]) -> bool:
        if isinstance(other, (Expr, Number)):
            return exprhash(self) == exprhash(other)
        else:
            return NotImplemented

    @number_casted
    def __add__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if self == undefined or other == undefined:
            return undefined
        elif other == 0:
            return self
        elif is_literal(self) and is_literal(other):
            if (self == inf and is_finite(other) or other == inf and is_finite(self)
                    or self == inf and other == inf):
                return inf
            elif (self == -inf and is_finite(other) or other == -inf and is_finite(self)
                    or self == -inf and other == -inf):
                return -inf
            elif self == inf and other == -inf or other == inf and self == -inf:
                return undefined
            else:
                return Add(self, other)
        elif isinstance(self, Add) and isinstance(other, Add):
            return Add(self.const + other.const, *self.terms, *other.terms)
        elif isinstance(self, Add) and isinstance(other, Number):
            return Add(self.const + other, *self.terms)
        elif isinstance(self, Add):
            return Add(self.const, *self.terms, other)
        elif isinstance(other, Add):
            return Add(self, other.const, *other.terms)
        else:
            return Add(self, other)

    @number_casted
    def __radd__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if self == undefined or other == undefined:
            return undefined
        elif other == 0:
            return self
        elif is_literal(self) and is_literal(other):
            if self == inf:
                return inf
            elif self == -inf:
                return -inf
            else:
                return Add(other, self)
        elif isinstance(self, Add) and isinstance(other, Number):
            return Add(self.const + other, *self.terms)
        elif isinstance(self, Add):
            return Add(other, self.const, *self.terms)
        else:
            return Add(other, self)

    @number_casted
    def __sub__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if self == undefined or other == undefined:
            return undefined
        elif other == 0:
            return self
        elif is_literal(self) and is_literal(other):
            if (self == inf and is_finite(other) or other == -inf and is_finite(self)
                    or self == inf and other == -inf):
                return inf
            elif (self == -inf and is_finite(other) or other == inf and is_finite(self)
                    or self == -inf and other == inf):
                return -inf
            elif self == inf and other == inf or other == -inf and self == -inf:
                return undefined
            else:
                return Add(self, -other)
        elif isinstance(self, Add) and isinstance(other, Number):
            return Add(self.const - other, *self.terms)
        elif isinstance(self, Add):
            return Add(-other, self.const, *self.terms)
        else:
            return Add(self, -other)

    @number_casted
    def __rsub__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if self == Undefined:
            return undefined
        elif other == 0:
            return -self
        elif is_literal(self) and is_literal(other):
            if self == -inf:
                return inf
            elif self == inf:
                return -inf
            else:
                return Add(other, -self)
        else:
            return Add(other, -self)

    @number_casted
    def __mul__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if self == undefined or other == undefined:
            return undefined
        elif other == 1:
            return self
        elif is_literal(self) and is_literal(other):
            if is_infinite(self) and is_infinite(other):
                return inf if self == other else -inf
            elif is_infinite(self) and is_finite(other):
                sign = get_sign(other)
                return self if sign == 1 else undefined if sign == 0 else -self
            elif is_infinite(other) and is_finite(self):
                sign = get_sign(self)
                return other if sign == 1 else undefined if sign == 0 else -other
            elif other == 0:
                return 0
            else:
                return Mul(self, other)
        elif isinstance(self, Mul) and isinstance(other, Number):
            return Mul(self.coef * other, *self.factors)
        elif isinstance(self, Mul) and isinstance(other, Mul):
            return Mul(self.coef * other.coef, *self.factors, *other.factors)
        elif isinstance(self, Mul):
            return Mul(self.coef, *self.factors, other)
        elif isinstance(other, Mul):
            return Mul(self, other.coef, *other.factors)
        else:
            return Mul(self, other)

    @number_casted
    def __rmul__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if self == undefined:
            return undefined
        elif other == 1:
            return self
        elif is_literal(self) and is_literal(other):
            if is_infinite(self):
                sign = get_sign(other)
                return self if sign == 1 else undefined if sign == 0 else -self
            elif other == 0:
                return 0
            else:
                return Mul(other, self)
        elif isinstance(self, Mul) and isinstance(other, Number):
            return Mul(self.coef * other, *self.factors)
        elif isinstance(self, Mul):
            return Mul(other, self.coef, *self.factors)
        else:
            return Mul(other, self)

    @number_casted
    def __truediv__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if self == undefined or other == undefined:
            return undefined
        elif is_literal(self) and is_literal(other):
            if is_infinite(self) and is_infinite(other):
                return undefined
            elif is_infinite(self) and is_finite(other):
                sign = get_sign(other)
                return self if sign > 0 else undefined if sign == 0 else -self
            elif is_finite(self) and is_infinite(other):
                return 0
            elif other == 1:
                return self
            elif other == 0:
                return undefined
            else:
                return Mul(self, frdiv(1, other))
        elif isinstance(self, Mul) and isinstance(other, Number):
            return Mul(frdiv(self.coef, other), *self.factors)
        elif isinstance(self, Mul):
            return Mul(self.coef, *self.factors, frdiv(1, other))
        else:
            return Mul(self, frdiv(1, other))

    @number_casted
    def __rtruediv__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if self == undefined or other == undefined:
            return undefined
        elif is_literal(other) and is_literal(self):
            if is_infinite(self):
                return 0
            elif other == 1:
                return Pow(self, -1)
            else:
                return Mul(other, Pow(self, -1))
        elif other == 1:
            return Pow(self, -1)
        else:
            return Mul(other, Pow(self, -1))

    @number_casted
    def __pow__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if self == undefined or other == undefined:
            return undefined
        elif is_literal(self) and is_literal(other):
            if is_infinite(self) and is_infinite(other):
                if self == inf and other == inf:
                    return inf
                elif self == inf and other == -inf:
                    return 0
                else:
                    return undefined
            elif is_infinite(self) and is_finite(other):
                if self == inf and isinstance(other, int):
                    return inf if other > 0 else undefined if other == 0 else 0
                elif self == -inf and isinstance(other, int):
                    return undefined if other == 0 else inf if other % 2 == 0 else -inf
                else:
                    return Pow(self, other)
            elif is_finite(self) and is_infinite(other):
                if self.value < 0:
                    return undefined
                elif other == inf:
                    return inf if self.value > 1 else 1 if self.value == 1 else 0
                else:
                    return 0 if self.value > 1 else 1 if self.value == 1 else inf
            elif other == 0:
                return Pow(self, 0)
            elif other == 1:
                return self
            else:
                return Pow(self, other)
        elif isinstance(self, Pow):
            return Pow(self.base, self.expo * other)
        else:
            return Pow(self, other)

    @number_casted
    def __rpow__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if self == undefined or other == undefined:
            return undefined
        elif is_literal(other) and is_literal(self):
            if is_infinite(self):
                if isinstance(other, Constant):
                    other = other.value
                if other < 0:
                    return undefined
                elif other == 0:
                    return 0 if self == inf else undefined
                elif self == inf:
                    return inf if other > 1 else 1 if other == 1 else 0
                else:
                    return 0 if other > 1 else 1 if other == 1 else inf
            elif other == 1:
                return 1
            elif other == 0:
                return 0 if get_sign(self) > 0 else undefined
            else:
                return Pow(other, self)
        elif other == 0:
            return Pow(0, self)
        elif other == 1:
            return self
        else:
            return Pow(other, self)

    def __neg__(self) -> Union["Expr", Number]:
        return Mul(self, -1)

    def __pos__(self) -> Union["Expr", Number]:
        return self

    def __abs__(self) -> Union["Expr", Number]:
        return undefined if self == undefined else Abs(self)

    def apply(self, func: FunctionType, *args) -> any:
        """
        Recursively apply some functionality to the expression instance.
        The function usually refers to some internally implemented
        method different for each subclass.
        Used for evaluation, simplifcation, and so on.

        :param self: The expression.
        :param func: The function to recursively apply.
        :type func: FunctionType
        :param args: Related arguments to pass to the function if any.
        :return: The result of the recursive function call.
        :rtype: Any
        """
        return func(self, *args)

    # expansions
    def expand_mul(self) -> Union["Expr", Number]:
        return self.apply(lambda x: x if isinstance(x, Number) else x._expand_mul())

    def _expand_mul(self) -> Union["Expr", Number]:
        return self

    def expand_dist(self) -> Union["Expr", Number]:
        return self.apply(lambda x: x if isinstance(x, Number) else x._expand_dist())

    def _expand_dist(self) -> Union["Expr", Number]:
        return self

    def expand_pow(self) -> Union["Expr", Number]:
        return self.apply(lambda x: x if isinstance(x, Number) else x._expand_pow())

    def _expand_pow(self) -> Union["Expr", Number]:
        return self

    def expand_trig(self) -> Union["Expr", Number]:
        return self.apply(lambda x: x if isinstance(x, Number) else x._expand_trig())

    def _expand_trig(self) -> Union["Expr", Number]:
        return self

    def expand_log(self) -> Union["Expr", Number]:
        return self.apply(lambda x: x if isinstance(x, Number) else x._expand_log())

    def _expand_log(self) -> Union["Expr", Number]:
        return self

    def expand(self) -> Union["Expr", Number]:
        expr = self
        while True:
            original = expr
            expr = expand_trig(expr)
            expr = expand_log(expr)
            expr = expand_pow(expr)
            expr = expand_distribute(expr)
            expr = expand_mul(expr)
            if original == expr:
                break
        return expr

    # substitutions and evaluations
    def doit(self) -> Union["Expr", Number]:
        return self.apply(lambda x: x if isinstance(x, Number) else x._doit())

    def _doit(self) -> Union["Expr", Number]:
        return self

    @map_casted
    def subs(self, expr_map: dict["Var", Union["Expr", Number]] | None = None, /) -> Union["Expr", Number]:
        return self.apply(lambda x: x if isinstance(x, Number) else x._subs(expr_map))

    def _subs(self, expr_map: dict["Var", Union["Expr", Number]] | None = None, /) -> Union["Expr", Number]:
        return self

    @map_casted
    def evalf(self, value_map: dict["Var", Number] | None = None, /) -> Union["Expr", Number]:
        return self.apply(lambda x: x if isinstance(x, Number) else x._evalf(value_map))

    def _evalf(self, value_map: dict["Var", Number] | None = None, /) -> Union["Expr", Number]:
        return self

    def diff(self, var: "Var", /, order: int = 1,
             *, evaluate: bool = True) -> Union["Expr", Number]:
        """
        Public differentiation method for either evaluating the derivative (default)
        or constructing a Derivative instance with the given variable and order.

        :param self: The expression instance.
        :param var: The variable to take the derivative with.
        :type var: "Var"
        :param order: The order of the derivative.
        :type order: int
        :param evaluate: Whether to evaluate the derivative instead of constructing
                         the Derivative class.
        :type evaluate: bool
        :return: Either the evaluated derivative or the Derivative class
                 of the given variable and order.
        :rtype: Expr | Number
        """
        if evaluate:
            result = self
            for _ in range(order):
                result = _diff(result, var)
                if result == 0:
                    break
            return result
        else:
            return Derivative(self, var, order)

    def _diff(self, var: "Var", /) -> Union["Expr", Number]:
        """
        Internal differentiation method for the classes to define derivative
        rules with. Should always be order 1 and evaluated.

        :param self: The expression instance.
        :param var: The variable to take the derivative with.
        :type var: "Var"
        :return: The evaluated derivative.
        :rtype: Expr | Number
        """
        raise NotImplementedError


ExprLike = Expr | Number


@typechecked
def frdiv(a: ExprLike, b: ExprLike, /) -> ExprLike:
    """Fraction-safe version for division."""
    if isinstance(a, int) and isinstance(b, int):
        return undefined if b == 0 else Fraction(a, b)
    elif b == 0:
        return undefined
    else:
        return a / b


@typechecked
def frpow(a: ExprLike, b: ExprLike, /) -> ExprLike:
    """Fraction-safe version for power."""
    if a == 0:
        if isinstance(b, Number):
            return 0 if b > 0 else undefined
        else:
            return a ** b
    elif isinstance(a, int) and isinstance(b, int) and b < 0:
        return Fraction(1, a ** -b)
    else:
        return a ** b


# Convert float -> 8 bytes big-endian → integer
@typechecked
def float_hash(x: float, /) -> int:
    """Convert float to an integer hash."""
    return int.from_bytes(pack(">d", x), "big")


@typechecked
def exprhash(expr: ExprLike, /) -> tuple[int | tuple, ...]:
    """
    Hash tuple generator for expressions and numbers.
    Only int/float/Fraction are supported as non-Expr numbers.
    This is used for sorting and hashing expressions.

    :param expr: The expression or number to hash.
    :type expr: ExprLike
    :return: The tuple of relevant information.
    :rtype: tuple[int | tuple, ...]
    """
    if isinstance(expr, int):
        return (KEY_INT, expr)
    elif isinstance(expr, Fraction):
        return (KEY_FRC, float(expr), expr.numerator)
    elif isinstance(expr, float):
        return (KEY_FLT, float_hash(expr))
    elif isinstance(expr, Expr):
        return expr.exprhash()
    elif isinstance(expr, Number):
        raise TypeError(f"unsupported number type for"
                        f" exprhash: {type(expr).__name__}")
    else:
        raise TypeError(f"unsupported value type for"
                        f" exprhash: {type(expr).__name__}")


@typechecked
def exprsorted(iterable: Iterable[ExprLike], /) -> list[ExprLike]:
    """
    Sort iterable of expressions by descending exprhash.

    :param iterable: The iterable of expressions to sort.
    :type iterable: Iterable[ExprLike]
    :return: The list of sorted expressions in descending order.
    :rtype: list[ExprLike]
    """
    return sorted(iterable, key=exprhash, reverse=True)


@typechecked
@inherit_docstrings
class Var(Expr):
    name: str

    def __init__(self, name: str, /):
        """
        Initialize a variable instance with given letter name.

        :param self: Instance to be initialized.
        :param name: The one-letter name from latin/greek alphabet.
                     Note that letters with existing constant
                     defined cannot be used: e, pi, tau, phi.
        :type name: str
        """
        if len(name) != 1:
            raise ValueError(
                "only single letters are allowed as variable names")
        if name not in VARIABLE_LETTERS:
            raise ValueError(
                f"variable name not in allowed list of "
                f"english and greek alphabets: {name}")
        if name in RESERVED_LETTERS:
            raise ValueError(
                f"letter is reserved for constants "
                f"or expressions: {name}")
        self.name = name

    def __repr__(self):
        return f"Var({self.name})"

    def __str__(self):
        return self.name

    def exprhash(self) -> CanonicalKey:
        # return (KEY_VAR, VARIABLE_LETTERS.index(self.name))
        return (KEY_VAR, VARIABLE_LETTERS[::-1].index(self.name))

    def _subs(self, expr_map: dict["Var", ExprLike] | None = None, /) -> ExprLike:
        return (expr_map or {}).get(self, self)

    def _evalf(self, value_map: dict["Var", Number] | None = None, /) -> ExprLike:
        return (value_map or {}).get(self, self)

    def _diff(self, var: "Var", /) -> ExprLike:
        return 1 if var == self else 0


ExprMap = dict[Var, ExprLike]
ExprValueMap = dict[Var, Number]


@typechecked
def symbols(letters: str, /) -> tuple[Var, ...]:
    """
    Creates variable instance(s) supporting multiple letters.

    :param letters: The letters of the variables.
    :type letters: str
    :return: A tuple of the variable instance(s).
    :rtype: tuple[Var]
    """
    return tuple(Var(letter) for letter in letters)


@typechecked
def apply(expr: ExprLike, func: FunctionType, *args) -> ExprLike:
    """
    Recursively apply some function with ExprLike instance.
    Any children node inside the expression is applied with the function
    before creating a new class, as passed arbitrarily or defined
    otherwise in the corresponding class method to recursively
    achieve functionalities like evaluations and simplifications.

    In each level, apply() is called to the children nodes before
    applying the function to the expr itself to achieve the recursion.

    :param expr: The expression to recursively apply function on.
    :type expr: ExprLike (Expr | Number)
    :param func: The function to recursively apply to the expression.
    :type func: FunctionType
    :rtype: ExprLike (Expr | Number)
    :raises TypeError: If any argument doesn't match expected type.
    """
    return expr.apply(func, *args) if isinstance(expr, Expr) else expr


@typechecked
def doit(expr: ExprLike, /) -> ExprLike:
    """
    Recursively apply Expr.doit() method with ExprLike instance.
    This function executes symbolic operations, like derivatives and
    integrals, and evaluates closed-form expressions, like Sin(π).

    :param expr: The expression to execute symbolic operations on.
    :type expr: ExprLike (Expr | Number)
    :rtype: ExprLike (Expr | Number)
    :raises TypeError: If any argument doesn't match expected type.
    """
    return apply(expr, lambda x: x._doit() if isinstance(x, Expr) else x)


@typechecked
def subs(expr: ExprLike, expr_map: ExprMap | None = None, /) -> ExprLike:
    """
    Recursively apply Expr.subs(mapping) method with ExprLike instance.
    This function replaces symbols exactly as written, doing no
    simplification, differentiation, integration, algebra, or
    arithmetic other than the basic structural flattening.
    Should be used when exact values of trigs are expected follwed by doit().
    Should not be used when numerical values are expected to be evaluated,
    which is done by evalf() instead.

    :param expr: The expression to perform substitution on.
    :type expr: ExprLike (Expr | Number)
    :param mapping: The dictionary with variables and matching expressions/values.
    :type mapping: ExprMap (dict[Var, ExprLike])
    :rtype: ExprLike
    :raises TypeError: If any argument doesn't match expected type.
    """
    return expr.subs(expr_map) if isinstance(expr, Expr) else expr


@typechecked
def evalf(expr: ExprLike, value_map: ExprValueMap | None = None, /) -> ExprLike:
    """
    Recursively apply Expr.evalf(mapping) method with ExprLike instance.
    This function calculates numeral approximation with floating-point
    computation. Small margins of error are to be expected when used on
    trignometry functions in place of doit(). Evaluates constants into values.
    If not all variables are given, all the given ones will be plugged in.
    Should be used when numerical values are expected to be evaluated.
    Should not be used when exact values of trigs are expected.
    which is done by evalf() instead, follwed by doit().

    :param expr: The expression to evaluate floating-point value on.
    :type expr: ExprLike (Expr | Number)
    :param mapping: The dictionary with variables and matching values.
    :type mapping: ExprValueMap (dict[Var, Number])
    :rtype: ExprLike
    :raises TypeError: If any argument doesn't match expected type.
    :raises NotImplementedError: If any subclass doesn't have
                                 .evalf() method implemented.
    """
    return expr.evalf(value_map) if isinstance(expr, Expr) else expr


@typechecked
def expand_mul(expr: ExprLike, /) -> ExprLike:
    """Distribute compact algebra in multiplication."""
    return expr.expand_mul() if isinstance(expr, Expr) else expr


@typechecked
def expand_pow(expr: ExprLike, /) -> ExprLike:
    """Distribute compact algebra in power."""
    return expr.expand_pow() if isinstance(expr, Expr) else expr


@typechecked
def expand_distribute(expr: ExprLike, /) -> ExprLike:
    """Distribute polynomial in natural number power."""
    return expr.expand_dist() if isinstance(expr, Expr) else expr


@typechecked
def expand_trig(expr: ExprLike, /) -> ExprLike:
    """Distribute compact algebra in trigonometry functions."""
    return expr.expand_trig() if isinstance(expr, Expr) else expr


@typechecked
def expand_log(expr: ExprLike, /) -> ExprLike:
    """Distribute compact algebra in logs."""
    return expr.expand_log() if isinstance(expr, Expr) else expr


# TODO: make this consistent for results considering expanding order.
@typechecked
def expand(expr: ExprLike, /) -> ExprLike:
    """Expand compact algebra in multiplication, power, trig, and log."""
    if isinstance(expr, Number):
        return expr
    while True:
        original = expr
        expr = expand_trig(expr)
        expr = expand_log(expr)
        expr = expand_pow(expr)
        expr = expand_distribute(expr)
        expr = expand_mul(expr)
        if original == expr:
            break
    return expr


# factor_poly is basically using linear/quadratic/cubic/quartic formulas
# on one variable for the solutions for zero value and writing the
# polynomial in factor form.
# quintic and up are still doable if polynomial is single variable
# but at some point it's pointless to factor
@typechecked
def factor_poly(expr: ExprLike, /) -> ExprLike:
    """Factor polynomial-like expressions."""
    return expr if isinstance(expr, Expr) else expr


# TODO: expand to add, mul, pow, trig, log, func
@typechecked
def factor(expr: ExprLike, /) -> ExprLike:
    """Combine/factor flattened algebra in multiplication, power, and so on."""
    return expr if isinstance(expr, Expr) else expr


@typechecked
def cancel(expr: ExprLike, /) -> ExprLike:
    """Reduce rational with an expression."""
    return expr if isinstance(expr, Expr) else expr


@typechecked
def together(expr: ExprLike, /) -> ExprLike:
    """Combine over common denominator."""
    return expr if isinstance(expr, Expr) else expr


@typechecked
def apart(expr: ExprLike, /) -> ExprLike:
    """Perform partial fraction decomposition."""
    return expr if isinstance(expr, Expr) else expr


@typechecked
def collect(expr: ExprLike, var: Var, /) -> ExprLike:
    """Group the expression by given variable."""
    return expr if isinstance(expr, Expr) else expr


@typechecked
def simplify(expr: ExprLike, /) -> ExprLike:
    """Simplify the given expression with basic math rules."""
    return expr if isinstance(expr, Expr) else expr


@typechecked
def _diff(expr: ExprLike, var: Var, /) -> ExprLike:
    """
    Calculates the first order derivative of the expression with
    respect to the given variable with.

    :param expr: The expression to take the derivative of.
    :type expr: ExprLike
    :param var: The variable to take the derivative with.
    :type var: Var
    :return: The derivative of the expression.
    :rtype: ExprLike
    """
    return expr._diff(var) if isinstance(expr, Expr) else 0


@typechecked
def diff(expr: ExprLike, var: Var, /, order: int = 1, *, evaluate=True) -> ExprLike:
    """
    Either evaluate the derivative (default) or construct a Derivative instance
    with the given variable and order. Calls on Expr.diff() for the implementation.

    :param self: The expression instance.
    :param var: The variable to take the derivative with.
    :type var: "Var"
    :param order: The order of the derivative.
    :type order: int
    :param evaluate: Whether to evaluate the derivative instead of constructing
                        the Derivative class.
    :type evaluate: bool
    :return: Either the evaluated derivative or the Derivative class
                of the given variable and order.
    :rtype: Expr | Number
    """
    return expr.diff(var, order, evaluate=evaluate) if isinstance(expr, Expr) else 0


@typechecked
@inherit_docstrings
class Constant(Expr):
    name: str
    value: Number

    def __init__(self, name: str, value: Number, /):
        """
        Initialize a mathematical constant instance with given name
        and value.

        :param self: Instance to be initialized.
        :param name: The one-letter name from latin/greek alphabet.
        :type name: str
        :param value: The constant value for evaluation purposes.
        :type value: Number
        """
        if value == 0:
            # This is to avoid invalid results and simplifications in Expr operators
            # since Constant.value is never checked to be 0.
            raise ValueError("Constant cannot be of value 0")
        elif not math_isfinite(value):
            # Avoiding manual construction of inf/-inf/nan constants
            raise ValueError("Constant cannot be of value inf/-inf/nan,"
                             " use special constants instead.")
        self.name = name
        self.value = value

    def __repr__(self):
        return f"Constant({self.name})"

    def __str__(self):
        return self.name

    def exprhash(self) -> CanonicalKey:
        return (KEY_CNS, VARIABLE_LETTERS[::-1].index(self.name) + 3)

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        return self.value

    def _diff(self, var: Var, /) -> ExprLike:
        return 0


@typechecked
@inherit_docstrings
class Infinity(Constant):
    name: str = "Infinity"

    def __init__(self, /):
        """
        Initialize an infinity constant.
        """
        pass

    def _evalf(self, value_map: ExprValueMap | None = None, /) -> ExprLike:
        return math_inf

    def exprhash(self) -> CanonicalKey:
        return (KEY_CNS, CONSTKEY_INF)

    def __neg__(self) -> ExprLike:
        return NegativeInfinity()

    def __pos__(self) -> ExprLike:
        return self


@typechecked
@inherit_docstrings
class NegativeInfinity(Constant):
    name: str = "-Infinity"

    def __init__(self, /):
        """
        Initialize a negative infinity constant.
        """
        pass

    def _evalf(self, value_map: ExprValueMap | None = None, /) -> ExprLike:
        return -math_inf

    def exprhash(self) -> CanonicalKey:
        return (KEY_CNS, CONSTKEY_NEGINF)

    def __neg__(self) -> ExprLike:
        return Infinity()

    def __pos__(self) -> ExprLike:
        return self


@typechecked
@inherit_docstrings
class Undefined(Constant):
    name: str = "Undefined"

    def __init__(self, /):
        """
        Initialize an undefined constant.
        """
        pass

    def _evalf(self, value_map: ExprValueMap | None = None, /) -> ExprLike:
        return math_nan

    def exprhash(self) -> CanonicalKey:
        return (KEY_CNS, CONSTKEY_UNDEF)

    def __neg__(self) -> ExprLike:
        return self

    def __pos__(self) -> ExprLike:
        return self


e = Constant('e', math_e)
pi = π = Constant('π', math_pi)
tau = τ = Constant('τ', math_tau)
phi = φ = Constant('φ', (1 + sqrt(5)) / 2)
inf = Infinity()
undefined = Undefined()


@typechecked
@number_casted
def is_special_const(expr: ExprLike, /) -> bool:
    """
    Determine if the expression or number is inf, -inf, or undefined.

    :param expr: The expression to check.
    :type expr: ExprLike
    :return: Whether if the expression is a special constant value.
    :rtype: bool
    """
    return expr in (inf, -inf, undefined)


@typechecked
@number_casted
def is_literal(expr: ExprLike, /) -> bool:
    """
    Determine if the expression or number is a literal constant, i.e.,
    a python number or a constant class. Does not count constant
    structures like Add or Mul.

    :param expr: The expression to check.
    :type expr: ExprLike
    :return: Whether if the expression is a direct constant value.
    :rtype: bool
    """
    return isinstance(expr, Number | Constant)


@typechecked
@number_casted
def is_infinite(expr: ExprLike, /) -> bool:
    """
    Determine if the expression is inf or -inf.

    :param expr: The expression to determine.
    :type expr: ExprLike
    :return: Whether if the expression is inf or -inf
    :rtype: bool
    """
    return cast_number(expr) in (inf, -inf)


@typechecked
@number_casted
def is_finite(expr: ExprLike, /) -> bool:
    """
    Determine if the expression or number is a literal constant that's not
    inf, -inf, or undefined.

    :param expr: The expression to check.
    :type expr: ExprLike
    :return: Whether if the expression is a finite constant value.
    :rtype: bool
    """
    return is_literal(expr) and not is_special_const(expr)


@typechecked
@number_casted
def get_value(expr: ExprLike, /) -> Number:
    """
    Obtain the number value of the constant or number.

    :param expr: The number or constant to get the value of.
    :type expr: ExprLike
    :return: The number value.
    :rtype: Number
    """
    if not is_finite(expr):
        raise ValueError(f"value must be finite, got {expr}")
    return expr.value if isinstance(expr, Constant) else expr


@typechecked
@number_casted
def get_sign(expr: ExprLike, /) -> int | None:
    """
    Determine the sign of the constant in terms of -1/0/1,
    return None if undefined or isn't a number or constant.

    :param expr: The literal to determine sign of.
    :type expr: ExprLike
    :return: The sign in -1/0/1 format or None.
    :rtype: int | None
    """
    if is_literal(expr):
        if expr == undefined:
            return
        elif expr == inf:
            return 1
        elif expr == -inf:
            return 1
        elif isinstance(expr, Constant):
            return 1 if expr.value > 0 else 0 if expr.value == 0 else -1
        else:
            return 1 if expr > 0 else 0 if expr == 0 else -1
    else:
        return


@typechecked
def cast_number(expr: ExprLike | Decimal, /) -> ExprLike:
    """
    Normalize an expression to replace python values of inf, -inf, and nan
    to `Constant` expressions and `Decimal` to `Fraction`. Only replaces top
    level value and does not search recursively.

    :param expr: The expression to cast.
    :type expr: ExprLike
    :return: The expression replaced by constants.
    :rtype: ExprLike
    """
    if isinstance(expr, float):
        if math_isinf(expr):
            return inf if expr > 0 else -inf
        elif math_isnan(expr):
            return undefined
        else:
            return expr
    elif isinstance(expr, Decimal):
        return Fraction(expr)
    else:
        return expr


@typechecked
def cast_map(subs: ExprMap, /) -> ExprMap:
    """
    Normalize substitution map with `cast_number` to replace special constants
    with `Constant` expressions and `Decimal` to `Fraction`. Only replaces top
    level value and does not search recursively.

    :param subs: The substitution map to normalize
    :type subs: ExprMap
    :return: The normalized substitution map.
    :rtype: ExprMap
    """
    return {var: cast_number(expr) for var, expr in subs.items()}


@typechecked
@number_casted
def format_term(term: ExprLike, /) -> tuple[bool, str]:
    """Format term for Add instance printing."""
    if isinstance(term, Mul):
        if len(term.factors) == 1 and term.coef == -1:
            return (False, f"{term.factors[0]}")
        elif len(term.factors) == 1 and term.coef < 0:
            return (False, f"{Mul(-term.coef, term.factors[0])}")
    return (True, f"{term}")


@typechecked
@inherit_docstrings
class Add(Expr):
    terms: list[ExprLike]
    const: Number = 0

    @number_casted
    def __init__(self, *terms: list[ExprLike]):
        """
        Initialize an addition instance with given terms.
        Python numbers are combined and terms are sorted.

        :param self: Instance to be initialized.
        :param terms: The elements of the addition.
        :type terms: list[ExprLike]
        """
        self.terms = []
        self.const = 0
        for term in terms:
            if isinstance(term, Number):
                self.const += term
            else:
                self.terms.append(term)
        self.terms = exprsorted(self.terms)

    def __repr__(self):
        return (
            "Add(" + ", ".join(f"{term!r}" for term in self.terms) +
            (f", {self.const}" if self.const != 0 else "") + ")"
        )

    def __str__(self):
        result = ''
        for term in self.terms:
            is_pos, content = format_term(term)
            if result == '':
                result += ('' if is_pos else '-') + content
            else:
                result += (' + ' if is_pos else ' - ') + content
        if result != '':
            if self.const != 0:
                sign = '+' if self.const >= 0 else '-'
                return f"{result} {sign} {abs(self.const)}"
            else:
                return result
        else:
            return f"{self.const}"

    def exprhash(self) -> CanonicalKey:
        return (KEY_ADD,) + tuple(
            term.exprhash() for term in self.terms
        ) + (exprhash(self.const),)

    def apply(self, func: FunctionType, *args) -> any:
        return func(
            self.const + sum(apply(term, func, *args)
                             for term in self.terms), *args,
        )

    def _diff(self, var: Var, /) -> ExprLike:
        return sum(_diff(term, var) for term in self.terms)


@typechecked
@number_casted
def format_factor(term: ExprLike, /) -> tuple[bool, str]:
    """Format factor for Mul instance printing."""
    if isinstance(term, Pow):
        if is_constant(term.expo):
            expo_value = evalf(term.expo)
            if expo_value == -1:
                return (False, f"{term.base}")
            elif expo_value < 0:
                return (False, f"{term.base ** (-term.expo)}")
        elif isinstance(term.expo, Mul):
            if term.expo.coef < 0:
                return (False, f"{term.base ** (-term.expo)}")
    if isinstance(term, Add):
        return (True, f"({term})")
    else:
        return (True, f"{term}")


@typechecked
@inherit_docstrings
class Mul(Expr):
    factors: list[ExprLike]
    coef: Number = 1

    @number_casted
    def __init__(self, *factors: list[ExprLike]):
        """
        Initialize a multiplication instance with given factors.
        Python numbers are combined and factors are sorted.

        :param self: Instance to be initialized.
        :param factors: The elements of the multiplication.
        :type factors: list[ExprLike]
        """
        self.factors = []
        self.coef = 1
        for factor in factors:
            if isinstance(factor, Number):
                if factor == 0:
                    self.coef == 0
                    self.factors.clear()
                    return
                self.coef *= factor
            else:
                self.factors.append(factor)
        self.factors = exprsorted(self.factors)

    def __repr__(self):
        return (
            "Mul(" + (f"{self.coef}, " if self.coef != 1 else "") +
            ", ".join(f"{factor!r}" for factor in self.factors) + ")"
        )

    def __str__(self):
        sign = '' if self.coef >= 0 else '-'
        numer_factors = [f"{abs(self.coef)}"]
        denom_factors = []
        for factor in self.factors:
            is_mul, content = format_factor(factor)
            if is_mul:
                numer_factors.append(content)
            else:
                denom_factors.append(content)
        numer = sign + " * ".join(numer_factors)
        if numer.startswith("1 * "):
            numer = numer[4:]
        denom = " * ".join(denom_factors)
        if denom == '':
            return numer
        else:
            if '*' in numer and len(numer_factors) > 1 + (abs(self.coef) != 1):
                numer = f"({numer})"
            if '*' in denom and len(denom_factors) > 1:
                denom = f"({denom})"
            return f"{numer} / {denom}"

    def exprhash(self) -> CanonicalKey:
        return (KEY_MUL,) + tuple(
            factor.exprhash() for factor in self.factors
        ) + (exprhash(self.coef),)

    def apply(self, func: FunctionType, *args) -> any:
        return func(
            self.coef * prod(apply(factor, func, *args)
                             for factor in self.factors), *args,
        )

    def _expand_mul(self) -> ExprLike:
        return sum(prod(group) for group in iterprod(
            [self.coef],
            *([factor.const] + factor.terms if isinstance(factor, Add)
              else [factor] for factor in self.factors),
        ))

    def _diff(self, var: Var, /) -> ExprLike:
        return self.coef * sum(
            prod(self.factors[:i]) * diff(factor, var)
            * prod(self.factors[i+1:])
            for i, factor in enumerate(self.factors)
        )


@typechecked
@inherit_docstrings
class Pow(Expr):
    base: ExprLike
    expo: ExprLike

    @number_casted
    def __init__(self, base: ExprLike, expo: ExprLike, /):
        """
        Initialize a power instance with given base and exponent.

        :param self: Instance to be initialized.
        :param base: The base of the power.
        :type base: ExprLike
        :param expo: The exponent of the power.
        :type expo: ExprLike
        """
        self.base = base
        self.expo = expo

    def __repr__(self):
        return f"Pow({self.base!r}, {self.expo!r})"

    def __str__(self):
        base_str = f"{self.base}"
        if isinstance(self.base, Add | Mul):
            base_str = f"({base_str})"
        expo_str = f"{self.expo}"
        if is_constant(self.expo):
            expo_value = evalf(self.expo)
            if expo_value == -1:
                return f"1 / {base_str}"
            elif expo_value < 0:
                return f"1 / {base_str}^{-self.expo}"
        elif isinstance(self.expo, Mul):
            if self.expo.coef < 0:
                expo_str = f"{-self.expo}"
                if '*' in expo_str:
                    expo_str = f"({expo_str})"
                return f"1 / {base_str}^{expo_str}"
        if isinstance(self.expo, Add | Mul):
            expo_str = f"({expo_str})"
        return f"{base_str}^{expo_str}"

    def exprhash(self) -> CanonicalKey:
        return (KEY_POW, exprhash(self.base), exprhash(self.expo))

    def apply(self, func: FunctionType, *args) -> any:
        return func(
            frpow(apply(self.base, func, *args),
                  apply(self.expo, func, *args)), *args,
        )

    def _expand_dist(self) -> ExprLike:
        if not isinstance(self.base, Add):
            return self
        elif not isinstance(self.expo, Number):
            return self
        elif self.expo % 1 != 0:
            return self
        # to maximize floating point safety in cost of efficiency
        expo = round(self.expo)

        all_terms = [self.base.const] + self.base.terms
        return sum(prod(group) for group in iterprod(
            *([all_terms] * expo)
        ))

    def _expand_pow(self) -> ExprLike:
        if isinstance(self.base, Mul):
            base_objs = [self.base.coef] + self.base.factors
        else:
            base_objs = [self.base]
        if isinstance(self.expo, Add):
            expo_objs = [self.expo.const] + self.expo.terms
        else:
            expo_objs = [self.expo]
        return prod(frpow(base, expo) for base, expo in iterprod(base_objs, expo_objs))

    def _diff(self, var: Var, /) -> ExprLike:
        if self.base == e:
            return self.base ** self.expo * self.expo._diff(var)
        elif is_constant(self.expo, var):
            return self.expo * frpow(self.base, (self.expo - 1)) * self.base._diff(var)
        else:
            return (e ** (Ln(self.base) * self.expo)) * (Ln(self.base) * self.expo)._diff(var)


@typechecked
@inherit_docstrings
class Function(Expr):
    name: str


@typechecked
@inherit_docstrings
class UnaryFunction(Function):
    name: str
    arg: ExprLike

    @number_casted
    def __init__(self, arg: ExprLike, /):
        """
        Initialize a unary function instance with given argument.

        :param self: Instance to be initialized.
        :param arg: The argument for the function.
        :type arg: ExprLike
        """
        self.arg = arg

    def __repr__(self):
        return f"{self.name}({self.arg!r})"

    def __str__(self):
        return f"{self.name}({self.arg})"

    def exprhash(self) -> CanonicalKey:
        return (KEY_FNC, self.func_key, exprhash(self.arg))

    def apply(self, func: FunctionType, *args) -> any:
        return func(
            type(self)(apply(self.arg, func, *args)), *args,
        )

    def _diff(self, var: Var, /) -> ExprLike:
        return self._basediff(self.arg) * self.arg._diff(var)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        """
        Base derivative for an unary function. Used by `_diff()`
        method for chain rule.

        :param self: The unary function instance.
        :param expr: The expression for the derivative function/expression.
                     Usually `self.arg`.
        :type expr: ExprLike
        :return: The expression plugged into the derivative function/expression.
        :rtype: ExprLike
        """
        raise NotImplementedError


@typechecked
@inherit_docstrings
class BinaryFunction(Function):
    name: str
    arg1: ExprLike
    arg2: ExprLike

    @number_casted
    def __init__(self, arg1: ExprLike, arg2: ExprLike, /):
        """
        Initialize a binary function instance with given arguments.

        :param self: Instance to be initialized.
        :param arg1: First argument for the function.
        :type arg1: ExprLike
        :param arg2: Second argument for the function.
        :type arg2: ExprLike
        """
        self.arg1 = arg1
        self.arg2 = arg2

    def __repr__(self):
        return f"{self.name}({self.arg1!r}, {self.arg2!r})"

    def __str__(self):
        return f"{self.name}({self.arg1}, {self.arg2})"

    def exprhash(self) -> CanonicalKey:
        return (KEY_FNC, self.func_key,
                exprhash(self.arg1), exprhash(self.arg2))

    def apply(self, func: FunctionType, *args) -> any:
        return func(
            type(self)(apply(self.arg1, func, *args),
                       apply(self.arg2, func, *args)), *args,
        )


@typechecked
@inherit_docstrings
class Abs(UnaryFunction):
    name = "Abs"
    func_key = FUNC_ABS

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return abs(arg_value) if isinstance(arg_value, Number) else Abs(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        # TODO: implement piecewise and this
        # return Piecewise(expr)
        raise NotImplementedError


@typechecked
@inherit_docstrings
class Sin(UnaryFunction):
    name = "Sin"
    func_key = FUNC_SIN

    # def _doit(self) -> ExprLike:
    #     arg = doit(self.arg)
    #     # TODO: simplify

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return sin(arg_value) if isinstance(arg_value, Number) else Sin(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return Cos(expr)


@typechecked
@inherit_docstrings
class Cos(UnaryFunction):
    name = "Cos"
    func_key = FUNC_COS

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return cos(arg_value) if isinstance(arg_value, Number) else Cos(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return -Sin(expr)


@typechecked
@inherit_docstrings
class Tan(UnaryFunction):
    name = "Tan"
    func_key = FUNC_TAN

    def _expand_trig(self) -> ExprLike:
        return Sin(self.arg) / Cos(self.arg)

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return tan(arg_value) if isinstance(arg_value, Number) else Tan(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return Sec(expr) ** 2


@typechecked
@inherit_docstrings
class Sec(UnaryFunction):
    name = "Sec"
    func_key = FUNC_SEC

    def _expand_trig(self) -> ExprLike:
        return 1 / Cos(self.arg)

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return (1 / cos(arg_value)) if isinstance(arg_value, Number) else Sec(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return Sec(expr) * Tan(expr)


@typechecked
@inherit_docstrings
class Csc(UnaryFunction):
    name = "Csc"
    func_key = FUNC_CSC

    def _expand_trig(self) -> ExprLike:
        return 1 / Sin(self.arg)

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return (1 / sin(arg_value)) if isinstance(arg_value, Number) else Csc(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return -Csc(expr) * Cot(expr)


@typechecked
@inherit_docstrings
class Cot(UnaryFunction):
    name = "Cot"
    func_key = FUNC_COT

    def _expand_trig(self) -> ExprLike:
        return Cos(self.arg) / Sin(self.arg)

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return (1 / tan(arg_value)) if isinstance(arg_value, Number) else Cot(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return -Csc(expr) ** 2


@typechecked
@inherit_docstrings
class Arcsin(UnaryFunction):
    name = "Arcsin"
    func_key = FUNC_ARCSIN

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return asin(arg_value) if isinstance(arg_value, Number) else Arcsin(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return 1 / (1 - expr ** 2).sqrt()


@typechecked
@inherit_docstrings
class Arccos(UnaryFunction):
    name = "Arccos"
    func_key = FUNC_ARCCOS

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return acos(arg_value) if isinstance(arg_value, Number) else Arccos(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return -1 / (1 - expr ** 2).sqrt()


@typechecked
@inherit_docstrings
class Arctan(UnaryFunction):
    name = "Arctan"
    func_key = FUNC_ARCTAN

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return atan(arg_value) if isinstance(arg_value, Number) else Arctan(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return 1 / (1 + expr ** 2)


@typechecked
@inherit_docstrings
class Arcsec(UnaryFunction):
    name = "Arcsec"
    func_key = FUNC_ARCSEC

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return acos(1 / arg_value) if isinstance(arg_value, Number) else Arcsec(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return 1 / Abs(expr) / (expr ** 2 - 1).sqrt()


@typechecked
@inherit_docstrings
class Arccsc(UnaryFunction):
    name = "Arccsc"
    func_key = FUNC_ARCCSC

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return asin(1 / arg_value) if isinstance(arg_value, Number) else Arccsc(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return -1 / Abs(expr) / (expr ** 2 - 1).sqrt()


@typechecked
@inherit_docstrings
class Arccot(UnaryFunction):
    name = "Arccot"
    func_key = FUNC_ARCCOT

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return atan(1 / arg_value) if isinstance(arg_value, Number) else Arccot(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return -1 / (1 + expr ** 2)


@typechecked
@inherit_docstrings
class Ln(UnaryFunction):
    name = "Ln"
    func_key = FUNC_LN

    def _expand_log(self) -> ExprLike:
        if isinstance(self.arg, Mul):
            return sum(
                ([] if self.arg.coef == 1 else [Ln(self.arg.coef)]) +
                [Ln(factor).expand_log() for factor in self.arg.factors]
            )
        elif isinstance(self.arg, Pow):
            return self.arg.expo * Ln(self.arg.base)
        else:
            return self

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return log(arg_value) if isinstance(arg_value, Number) else Ln(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return 1 / expr


@typechecked
@inherit_docstrings
class Log(BinaryFunction):
    name = "Log"
    func_key = FUNC_LOG

    def _expand_log(self) -> ExprLike:
        if isinstance(self.arg1, Mul):
            return sum(
                ([] if self.arg1.coef == 1 else [Log(self.arg1.coef, self.arg2)]) +
                [Log(factor, self.arg2).expand_log()
                 for factor in self.arg1.factors]
            )
        elif isinstance(self.arg1, Pow):
            return self.arg1.expo * Log(self.arg1.base, self.arg2)
        else:
            return self

    def _evalf(self, value_map: ExprValueMap, /) -> ExprLike:
        arg1_value = evalf(self.arg1, value_map)
        arg2_value = evalf(self.arg2, value_map)
        if isinstance(arg1_value, Number) and isinstance(arg2_value, Number):
            return log(arg1_value, arg2_value)
        else:
            return Log(arg1_value, arg2_value)

    def _basediff(self, var: Var, /) -> ExprLike:
        return (Ln(self.arg1) / Ln(self.arg2))._diff()


@typechecked
@inherit_docstrings
class Limit(Expr):
    expr: ExprLike
    var: Var
    point: ExprLike
    direction: int = 0

    @number_casted
    def __init__(self, expr: ExprLike, var: Var,
                 point: ExprLike, direction: int = 0, /):
        """
        Initialize a limit instance with given expression, variable,
        limit value point, and optional direction.

        :param self: Instance to be initialized.
        :param expr: The expression for the limit.
        :type expr: ExprLike
        :param var: The variable of the limit.
        :type var: Var
        :param point: The value point of the limit.
        :type point: ExprLike
        :param direction: The optional direction for the value point.
                          -1 for left (-), 0 for both ( ), and 1 for right (+).
                          Default to both sides.
        :type direction: int
        """
        if direction not in (-1, 0, 1):
            raise TypeError(f"direction must be -1 or 0 or 1,"
                            f" not {direction}")
        self.expr = expr
        self.var = var
        self.point = point
        self.direction = direction

    def __str__(self):
        dir_str = ('-', '', '+')[self.direction + 1]
        return f"limit({self.var} -> {self.point}{dir_str}) ({self.expr})"

    def __repr__(self):
        dir_str = ('-', '', '+')[self.direction + 1]
        return f"limit({self.var!r} -> {self.point!r}{dir_str}) ({self.expr!r})"

    def exprhash(self) -> CanonicalKey:
        return (KEY_LIM, exprhash(self.expr), exprhash(self.var),
                exprhash(self.point), exprhash(self.direction))

    def apply(self, func: FunctionType, *args) -> any:
        return func(
            Limit(apply(self.expr, func, *args),
                  apply(self.var, func, *args),
                  apply(self.point, func, *args),
                  self.direction), *args,
        )

    def _doit(self):
        if is_constant(self.expr, self.var):
            return self.expr
        else:
            return self


@typechecked
@inherit_docstrings
class Derivative(Expr):
    expr: ExprLike
    var: Var
    order: int = 1

    @number_casted
    def __init__(self, expr: ExprLike, var: Var,
                 order: int = 1, /):
        """
        Initialize a derivative instance with given expression, variable,
        and optional order (default to first order).

        :param self: Instance to be initialized.
        :param expr: The expression for the derivative.
        :type expr: ExprLike
        :param var: The variable to take the derivative with.
        :type var: Var
        :param order: The order of the derivative, default to 1.
        :type order: int
        """
        self.expr = expr
        self.var = var
        if order <= 0:
            raise ValueError(f"the order of a derivative must be"
                             f" a positive integer, got {order}")
        self.order = order

    def __str__(self):
        order_str = '' if self.order == 1 else f"^{self.order}"
        return f"d{order_str}/d{self.var}{order_str} ({self.expr})"

    def __repr__(self):
        order_str = '' if self.order == 1 else f"^{self.order}"
        return f"d{order_str}/d{self.var}{order_str} ({self.expr!r})"

    def exprhash(self) -> CanonicalKey:
        return (KEY_DER, exprhash(self.expr),
                exprhash(self.var), exprhash(self.order))

    def apply(self, func: FunctionType, *args) -> any:
        return func(
            Derivative(apply(self.expr, func, *args),
                       apply(self.var, func, *args),
                       self.order), *args,
        )

    def _doit(self):
        return diff(self.expr, self.var, self.order)

    def diff(self, var: Var, order: int = 1, /, evaluate: bool = True) -> ExprLike:
        if evaluate:
            result = self.doit()
            for _ in range(order):
                result = _diff(result, var)
                if result == 0:
                    break
            return result
        elif var == self.var:
            return Derivative(self.expr, self.var, self.order + 1)
        else:
            return Derivative(self, var, self.order + 1)

    def _diff(self, var: Var, /) -> ExprLike:
        return _diff(self.doit(), var)


@typechecked
@inherit_docstrings
class Integral(Expr):
    expr: ExprLike
    var: Var
    lower: ExprLike | None = None
    upper: ExprLike | None = None

    @number_casted
    def __init__(self, expr: ExprLike, var: Var,
                 lower: ExprLike | None = None, upper: ExprLike | None = None, /):
        """
        Initialize an integral instance with given expression, variable,
        and optional bounds (default to be both None). Bounds must either
        be both ExprLike values for a definite integral or both None
        for an indefinite integral.

        :param self: Instance to be initialized.
        :param expr: The expression for the integral.
        :type expr: ExprLike
        :param var: The variable to take the integral with.
        :type var: Var
        :param lower: Optional lower bound.
        :type lower: ExprLike | None
        :param upper: Optional upper bound.
        :type upper: ExprLike | None
        """
        self.expr = expr
        self.var = var
        if (lower is None) != (upper is None):
            raise ValueError("Integrals must either have both"
                             " bounds defined or both undefined")
        self.lower = lower
        self.upper = upper

    def __str__(self):
        bound_str = '' if self.lower is None else f"[{self.lower}, {self.upper}]"
        return f"∫{bound_str} {self.expr} d{self.var}"

    def __repr__(self):
        bound_str = '' if self.lower is None else f"[{self.lower!r}, {self.upper!r}]"
        return f"∫{bound_str} {self.expr!r} d{self.var}"

    def exprhash(self) -> CanonicalKey:
        if self.lower is None:
            return (KEY_INT, exprhash(self.expr), exprhash(self.var))
        else:
            return (KEY_INT, exprhash(self.expr), exprhash(self.var),
                    exprhash(self.lower), exprhash(self.upper))

    def apply(self, func: FunctionType, *args) -> any:
        return func(
            Integral(apply(self.expr, func, *args),
                     apply(self.var, func, *args),
                     None if self.lower is None else apply(
                         self.lower, func, *args),
                     None if self.upper is None else apply(self.upper, func, *args)), *args,
        )

    def indef_doit(self) -> ExprLike | None:
        if is_constant(self.expr, self.var):
            return self.expr * self.var
        elif is_mono(self.expr, self.var):  # ! could be buggy
            coef, expo = split_mono(self.expr, self.var)
            if expo == -1:
                return coef * Ln(self.var)
            else:
                return coef / (expo + 1) * self.var ** (expo + 1)
        elif is_poly(expanded := expand(self.expr), self.var):
            return sum(integrate(mono, self.var)
                       for mono in split_poly(expanded))

    # TODO: Risch/Risch-Norman algorithm
    def _doit(self) -> ExprLike:
        indef_result = self.indef_doit()
        if indef_result is None:
            return self
        elif self.lower is None:
            return indef_result
        else:
            return indef_result.subs({self.var: self.upper}) - indef_result.subs({self.var: self.lower})

    def _diff(self, var: Var, /) -> ExprLike:
        if self.var == var and self.lower is None and self.upper is None:
            return self.expr
        else:
            result = doit(self)
            return result if isinstance(result, Integral) else result._diff(var)


@typechecked
def is_constant(expr: ExprLike, var: Var | None = None, /) -> bool:
    """
    Determines if the expression is a constant value, aka. a structure
    with only known functions and numbers and constants. An optional
    variable variable can be passed to check if the variable is included
    in the expression.

    :param expr: The expression to check.
    :type expr: ExprLike
    :param var: The optional variable to check relationship with.
    :type var: Var | None
    :return: Whether if the expression is constant
             or unrelated to the given variable
    :rtype: bool
    """
    if isinstance(expr, Number | Constant):
        return True
    elif isinstance(expr, Var):
        return var is not None and expr != var
    elif isinstance(expr, Add):
        return all(is_constant(term, var) for term in expr.terms)
    elif isinstance(expr, Mul):
        return all(is_constant(factor, var) for factor in expr.factors)
    elif isinstance(expr, Pow):
        return is_constant(expr.base, var) and is_constant(expr.expo, var)
    elif isinstance(expr, UnaryFunction):
        return is_constant(expr.arg, var)
    elif isinstance(expr, BinaryFunction):
        return is_constant(expr.arg1, var) and is_constant(expr.arg2, var)
    elif isinstance(expr, Limit):
        return (is_constant(expr.expr, var) and is_constant(expr.var, var)
                and is_constant(expr.point, var))
    elif isinstance(expr, Derivative):
        return is_constant(expr.expr, var) and is_constant(expr.var, var)
    elif isinstance(expr, Integral):
        return (is_constant(expr.expr, var) and is_constant(expr.var, var)
                and is_constant(expr.lower, var) and is_constant(expr.upper, var))
    else:
        return False


@typechecked
def is_rat_constant(expr: ExprLike, var: Var | None = None, /) -> bool:
    """
    Determines if the expression is a rational constant value, aka.
    a structure with only known functions and rational numbers. An optional
    variable variable can be passed to check if the variable is included
    in the expression.

    :param expr: The expression to check.
    :type expr: ExprLike
    :param var: The optional variable to check relationship with.
    :type var: Var | None
    :return: Whether if the expression is rational constant
             or unrelated to the given variable
    :rtype: bool
    """
    if isinstance(expr, Number):
        return isinstance(expr, Rational)
    elif isinstance(expr, Var):
        return var is not None and expr != var
    elif isinstance(expr, Add):
        return all(is_rat_constant(term, var) for term in expr.terms)
    elif isinstance(expr, Mul):
        return all(is_rat_constant(factor, var) for factor in expr.factors)
    elif isinstance(expr, Pow):
        return is_rat_constant(expr.base, var) and is_rat_constant(expr.expo, var)
    elif isinstance(expr, UnaryFunction):
        return is_rat_constant(expr.arg, var)
    elif isinstance(expr, BinaryFunction):
        return is_rat_constant(expr.arg1, var) and is_rat_constant(expr.arg2, var)
    elif isinstance(expr, Limit):
        return (is_rat_constant(expr.expr, var) and is_rat_constant(expr.var, var)
                and is_rat_constant(expr.point, var))
    elif isinstance(expr, Derivative):
        return is_rat_constant(expr.expr, var) and is_rat_constant(expr.var, var)
    elif isinstance(expr, Integral):
        return (is_rat_constant(expr.expr, var) and is_rat_constant(expr.var, var)
                and is_rat_constant(expr.lower, var) and is_rat_constant(expr.upper, var))
    else:
        return False


@typechecked
def is_mono(expr: ExprLike, var: Var | None = None, /) -> bool:
    """
    Determines if the expression is a monomial only consisting of
    constants, variables, and mul/pow nodes. Exponents are required
    to be constants. If a variable is given, all other variables
    are seen as constants.

    :param expr: The expression to check.
    :type expr: ExprLike
    :return: Whether if the expression has a monomial structure.
    :rtype: bool
    """
    if isinstance(expr, Number | Constant | Var):
        return True
    elif isinstance(expr, Mul):
        return all(is_mono(factor, var) for factor in expr.factors)
    elif isinstance(expr, Pow):
        return is_mono(expr.base, var) and is_constant(expr.expo, var)
    else:
        return False


@typechecked
def split_mono(expr: ExprLike, var: Var, /) -> tuple[ExprLike, ExprLike]:
    """
    Organize the monomial into ax^b form with x being the given variable.
    Return the tuple of (a, b).

    :param expr: Description
    :type expr: ExprLike
    :param var: Description
    :type var: Var
    :return: Description
    :rtype: tuple[ExprLike, ExprLike]
    """

    if not is_mono(expr, var):
        raise ValueError(f"a monomial of {var} expected, got {expr}")
    elif is_constant(expr, var):
        return (expr, 0)

    expanded = expand(expr)
    factors = ([expanded.coef] + expanded.factors
               if isinstance(expanded, Mul) else [expanded])
    coef = 1
    expo = 0
    for factor in factors:
        if is_constant(factor, var):
            coef *= factor
        elif isinstance(factor, Var):
            expo += 1
        elif isinstance(factor, Mul):
            raise ValueError(
                "Mul instance not expected as monomial is expanded")
        elif isinstance(factor, Pow):
            if factor.base == var:
                expo += factor.expo
            else:
                raise ValueError(f"unexpected unconstant power base with"
                                 f" expanded monomial: {factor.base}")
        else:
            raise ValueError(f"unexpected factor: {factor}")
    return (coef, expo)


@typechecked
def is_poly(expr: ExprLike, var: Var | None = None, /) -> bool:
    """
    Determines if the expression is a polynomial only consisting of
    constants, variables, add/mul/pow nodes. Exponents are required
    to be constants. Does NOT check if constants are int/frac or if
    the exponents are integers or positive. If a variable is given,
    all other variables are seen as constants

    :param expr: The expression to check.
    :type expr: ExprLike
    :return: Whether if the expression has a polynomial structure.
    :rtype: bool
    """
    if isinstance(expr, Var) or is_constant(expr, var):
        return True
    elif isinstance(expr, Add):
        return all(is_poly(term, var) for term in expr.terms)
    elif isinstance(expr, Mul):
        return all(is_poly(factor, var) for factor in expr.factors)
    elif isinstance(expr, Pow):
        return is_poly(expr.base, var) and is_constant(expr.expo, var)
    else:
        return False


@typechecked
def is_pos_poly(expr: ExprLike, var: Var | None = None, /) -> bool:
    """
    Determines if the expression is a polynomial only consisting of
    constants, variables, add/mul/pow nodes. Exponents are required
    to be natural numbers. Does NOT check if constants are int/frac
    If a variable is given, all other variables are seen as constants.

    :param expr: The expression to check.
    :type expr: ExprLike
    :return: Whether if the expression has a natural number exponent
             only polynomial structure.
    :rtype: bool
    """
    if isinstance(expr, Var) or is_constant(expr, var):
        return True
    elif isinstance(expr, Add):
        return all(is_pos_poly(term, var) for term in expr.terms)
    elif isinstance(expr, Mul):
        return all(is_pos_poly(factor, var) for factor in expr.factors)
    elif isinstance(expr, Pow):
        return (is_pos_poly(expr.base, var) and isinstance(expr.expo, Number)
                and expr.expo % 1 == 0 and expr.expo >= 0)
    else:
        return False


@typechecked
def is_rat_poly(expr: ExprLike, var: Var | None = None, /) -> bool:
    """
    Determines if the expression is a polynomial only consisting of
    rational constants, variables, add/mul/pow nodes. Does NOT check
    if the exponents are integers or positive. If a variable is given,
    all other variables are seen as constants

    :param expr: The expression to check.
    :type expr: ExprLike
    :return: Whether if the expression has a rational constant only
             polynomial structure.
    :rtype: bool
    """
    if isinstance(expr, Var) or is_rat_constant(expr, var):
        return True
    elif isinstance(expr, Number):
        return isinstance(expr, Rational)
    elif is_constant(expr, var):
        return True
    elif isinstance(expr, Add):
        return isinstance(expr.const, Rational) and all(is_rat_poly(term, var) for term in expr.terms)
    elif isinstance(expr, Mul):
        return isinstance(expr.coef, Rational) and all(is_rat_poly(factor, var) for factor in expr.factors)
    elif isinstance(expr, Pow):
        return is_rat_poly(expr.base, var) and isinstance(expr.expo, Rational)
    else:
        return False


@typechecked
def is_perfect_poly(expr: ExprLike, var: Var | None = None, /) -> bool:
    """
    Determines if the expression is a polynomial only consisting of
    rational constants, variables, add/mul/pow nodes. Exponents are
    required to be natural numbers. If a variable is given, all other
    variables are seen as constants

    :param expr: The expression to check.
    :type expr: ExprLike
    :return: Whether if the expression has a natural number exponent
             and rational constants only polynomial structure.
    :rtype: bool
    """
    if isinstance(expr, Var) or is_rat_constant(expr, var):
        return True
    elif isinstance(expr, Add):
        return isinstance(expr.const, Rational) and all(is_perfect_poly(term, var) for term in expr.terms)
    elif isinstance(expr, Mul):
        return isinstance(expr.coef, Rational) and all(is_perfect_poly(factor, var) for factor in expr.factors)
    elif isinstance(expr, Pow):
        return is_perfect_poly(expr.base, var) and isinstance(expr.expo, int) and expr.expo >= 0
    else:
        return False


@typechecked
def split_poly(expr: ExprLike, /) -> list[ExprLike]:
    """
    Return a list of addition operants or the expression itself if not
    an add node. Useful for processing polynomials into monomials when
    checked to be polynomials.

    :param expr: The expression to split.
    :type expr: ExprLike
    :return: The list of addition operants or
             a list of the expressiion itself.
    :rtype: list[ExprLike]
    """
    return [expr.const] + expr.terms if isinstance(expr, Add) else [expr]


@typechecked
def integrate(expr: ExprLike, var: Var,
              lower: ExprLike | None = None, upper: ExprLike | None = None, /):
    return Integral(expr, var, lower, upper).doit()


def main():
    x, y, z = symbols("xyz")
    print(Constant('b', ''))


if __name__ == "__main__":
    main()
