from itertools import product as iterprod
from math import e as math_e, pi as math_pi, \
    sqrt, prod, sin, cos, tan, log, asin, acos, atan
from numbers import Number
from struct import pack
from types import FunctionType
from typing import Iterable, Union

from typeguard import typechecked


VARIABLE_LETTERS = (
    "abcdefghijklmnopqrstuvwxyz"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "αβγδεζηθικλμνξοπρστυφχψω"
    "ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ"
)

RESERVED_LETTERS = "eπφΣΠ"


def inherit_docstrings(cls: type) -> type:
    """
    Decorator for the subclass to inherit the class methods
    docstrings if it doesn't override it.

    :param cls: The parent class
    :type cls: type
    :return: The class with the method docstrings filled.
    :rtype: type
    """
    for name, attr in cls.__dict__.items():
        func = getattr(cls, name)
        if not getattr(func, "__doc__", None):
            for base in cls.__mro__[1:]:
                parent = getattr(base, name, None)
                if parent and getattr(parent, "__doc__", None):
                    func.__doc__ = parent.__doc__
                    break
    return cls


CanonicalKey = tuple[int | float | tuple, ...]

KEY_INT = 0
KEY_FLT = 1
KEY_CNS = 2
KEY_VAR = 3
KEY_ADD = 4
KEY_MUL = 5
KEY_POW = 6
KEY_FNC = 7
KEY_LIM = 8
KEY_DER = 9
KEY_ITG = 10

FUNC_ABS = 0
FUNC_SIN = 0
FUNC_COS = 0
FUNC_TAN = 0
FUNC_SEC = 0
FUNC_CSC = 0
FUNC_COT = 0
FUNC_ARCSIN = 0
FUNC_ARCCOS = 0
FUNC_ARCTAN = 0
FUNC_ARCSEC = 0
FUNC_ARCCSC = 0
FUNC_ARCCOT = 0
FUNC_LN = 0
FUNC_LOG = 0


@typechecked
class Expr:
    def __init__(self, *args):
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
        return Pow(self, 1 / 2)

    def cbrt(self) -> Union["Expr", Number]:
        """
        Returns the cube root of expression,
        equivalent to expr ** (1/3)

        :param self: The expression.
        :return: The cube root of the expression.
        :rtype: Expr | Number
        """
        return Pow(self, 1 / 3)

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
            expr = expand_mul(expr)
            if original == expr:
                break
        return expr

    # substitutions and evaluations
    def doit(self) -> Union["Expr", Number]:
        return self.apply(lambda x: x if isinstance(x, Number) else x._doit())

    def _doit(self) -> Union["Expr", Number]:
        return self

    def subs(self, expr_map: dict["Var", Union["Expr", Number]] | None = None, /) -> Union["Expr", Number]:
        return self.apply(lambda x: x if isinstance(x, Number) else x._subs(expr_map))

    def _subs(self, expr_map: dict["Var", Union["Expr", Number]] | None = None, /) -> Union["Expr", Number]:
        return self

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

    def __repr__(self):
        raise NotImplementedError

    def __str__(self):
        raise NotImplementedError

    def __hash__(self):
        return hash(self.exprhash())

    def exprhash(self) -> CanonicalKey:
        raise NotImplementedError

    def __eq__(self, other: Union["Expr", Number]) -> bool:
        return isinstance(other, (Expr, Number)) and exprhash(self) == exprhash(other)

    def __add__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if other == 0:
            return self
        elif isinstance(self, Add) and isinstance(other, Add):
            return Add(self.const + other.const, *self.terms, *other.terms)
        elif isinstance(self, Add):
            return Add(self.const, *self.terms, other)
        elif isinstance(other, Add):
            return Add(self, other.const, *other.terms)
        else:
            return Add(self, other)

    def __radd__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if other == 0:
            return self
        elif isinstance(self, Add):
            return Add(other, self.const, *self.terms)
        else:
            return Add(other, self)

    def __sub__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if other == 0:
            return self
        elif isinstance(self, Add):
            return Add(Mul(-1, other), self.const, *self.terms)
        else:
            return Add(self, Mul(-1, other))

    def __rsub__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if other == 0:
            return Mul(-1, self)
        else:
            return Add(other, Mul(-1, self))

    def __mul__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if other == 1:
            return self
        elif isinstance(self, Mul) and isinstance(other, Mul):
            return Mul(self.coef * other.coef, *self.factors, *other.factors)
        elif isinstance(self, Mul):
            return Mul(self.coef, *self.factors, other)
        elif isinstance(other, Mul):
            return Mul(self, other.coef, *other.factors)
        else:
            return Mul(self, other)

    def __rmul__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if other == 1:
            return self
        elif isinstance(self, Mul):
            return Mul(other, self.coef, *self.factors)
        else:
            return Mul(other, self)

    def __truediv__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if other == 1:
            return self
        elif isinstance(self, Mul):
            return Mul(self.coef, *self.factors, Pow(other, -1))
        else:
            return Mul(self, Pow(other, -1))

    def __rtruediv__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if other == 1:
            return Pow(self, -1)
        elif isinstance(other, Mul):
            return Mul(other.coef, *other.factors, Pow(self, -1))
        else:
            return Mul(other, Pow(self, -1))

    def __pow__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        if other == 0:
            # TODO: check for 0^0
            return 1
        elif other == 1:
            return self
        elif isinstance(self, Pow):
            return Pow(self.base, self.exp * other)
        else:
            return Pow(self, other)

    def __rpow__(self, other: Union["Expr", Number]) -> Union["Expr", Number]:
        # TODO: check for 0^0
        if other == 0 or other == 1:
            return other
        else:
            return Pow(other, self)

    def __neg__(self) -> Union["Expr", Number]:
        return self * -1

    def __pos__(self) -> Union["Expr", Number]:
        return self

    def __abs__(self) -> Union["Expr", Number]:
        return Abs(self)


ExprLike = Expr | Number


def float_bits(x: float, /) -> int:
    # Convert float -> 8 bytes big-endian → integer
    return int.from_bytes(pack(">d", x), "big")


@typechecked
def exprhash(expr: ExprLike, /) -> tuple[int | tuple, ...]:
    if isinstance(expr, int):
        return (KEY_INT, expr)
    elif isinstance(expr, float):
        return (KEY_FLT, float_bits(expr))
    elif isinstance(expr, Expr):
        return expr.exprhash()


@typechecked
def exprsorted(iterable: Iterable[ExprLike], /) -> list[ExprLike]:
    return sorted(iterable, key=exprhash, reverse=True)


@inherit_docstrings
class Var(Expr):
    name: str

    def __init__(self, name: str, /):
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
ValueMap = dict[Var, Number]


@inherit_docstrings
class Constant(Expr):
    name: str
    value: Number

    def __init__(self, name: str, value: Number, /):
        self.name = name
        self.value = value

    def __repr__(self):
        return f"Constant({self.name})"

    def __str__(self):
        return self.name

    def exprhash(self) -> CanonicalKey:
        # return (KEY_CNS, VARIABLE_LETTERS.index(self.name))
        return (KEY_CNS, VARIABLE_LETTERS[::-1].index(self.name))

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        return self.value

    def _diff(self, var: Var, /) -> ExprLike:
        return 0


def format_term(term: ExprLike, /) -> tuple[bool, str]:
    if isinstance(term, Mul):
        if len(term.factors) == 1 and term.coef == -1:
            return (False, f"{term.factors[0]}")
        elif len(term.factors) == 1 and term.coef < 0:
            return (False, f"{Mul(-term.coef, term.factors[0])}")
    return (True, f"{term}")


e = Constant('e', math_e)
pi = π = Constant('π', math_pi)
phi = φ = Constant('φ', (1 + sqrt(5)) / 2)


@inherit_docstrings
class Add(Expr):
    terms: list[Expr]
    const: Number = 0

    def __init__(self, *terms: list[ExprLike]):
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


def format_factor(term: ExprLike, /) -> tuple[bool, str]:
    if isinstance(term, Pow):
        if is_constant(term.exp):
            exp_value = evalf(term.exp)
            if exp_value == -1:
                return (False, f"{term.base}")
            elif exp_value < 0:
                return (False, f"{term.base ** (-term.exp)}")
        elif isinstance(term.exp, Mul):
            if term.exp.coef < 0:
                return (False, f"{term.base ** (-term.exp)}")
    if isinstance(term, Add):
        return (True, f"({term})")
    else:
        return (True, f"{term}")


@inherit_docstrings
class Mul(Expr):
    factors: list[Expr]
    coef: Number = 1

    def __init__(self, *factors: list[ExprLike]):
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

    def _expand_mul(self) -> Union["Expr", Number]:
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


@inherit_docstrings
class Pow(Expr):
    base: ExprLike
    exp: ExprLike

    def __init__(self, base: ExprLike, exp: ExprLike, /):
        self.base = base
        self.exp = exp

    def __repr__(self):
        return f"Pow({self.base!r}, {self.exp!r})"

    def __str__(self):
        base_str = f"{self.base}"
        if isinstance(self.base, Add | Mul):
            base_str = f"({base_str})"
        exp_str = f"{self.exp}"
        if is_constant(self.exp):
            exp_value = evalf(self.exp)
            if exp_value == -1:
                return f"1 / {base_str}"
            elif exp_value < 0:
                return f"1 / {base_str}^{-self.exp}"
        elif isinstance(self.exp, Mul):
            if self.exp.coef < 0:
                exp_str = f"{-self.exp}"
                if '*' in exp_str:
                    exp_str = f"({exp_str})"
                return f"1 / {base_str}^{exp_str}"
        if isinstance(self.exp, Add | Mul):
            exp_str = f"({exp_str})"
        return f"{base_str}^{exp_str}"

    def exprhash(self) -> CanonicalKey:
        return (KEY_POW, exprhash(self.base), exprhash(self.exp))

    def apply(self, func: FunctionType, *args) -> any:
        return func(
            apply(self.base, func, *args) ** apply(self.exp, func, *args), *args,
        )

    def _expand_pow(self) -> Union["Expr", Number]:
        if isinstance(self.base, Mul):
            base_objs = [self.base.coef] + self.base.factors
        else:
            base_objs = [self.base]
        if isinstance(self.exp, Add):
            exp_objs = [self.exp.const] + self.exp.terms
        else:
            exp_objs = [self.exp]
        return prod(base ** exp for base, exp in iterprod(base_objs, exp_objs))

    def _diff(self, var: Var, /) -> ExprLike:
        if self.base == e:
            return self.base ** self.exp * self.exp._diff(var)
        elif is_constant(self.exp, var):
            return self.exp * self.base ** (self.exp - 1) * self.base._diff(var)
        else:
            return (e ** (Ln(self.base) * self.exp)) * (Ln(self.base) * self.exp)._diff(var)


@inherit_docstrings
class Function(Expr):
    name: str


@inherit_docstrings
class UnaryFunction(Function):
    name: str
    arg: ExprLike

    def __init__(self, arg: ExprLike, /):
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


@inherit_docstrings
class BinaryFunction(Function):
    name: str
    arg1: ExprLike
    arg2: ExprLike

    def __init__(self, arg1: ExprLike, arg2: ExprLike, /):
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


@inherit_docstrings
class Abs(UnaryFunction):
    name = "Abs"
    func_key = FUNC_ABS

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return abs(arg_value) if isinstance(arg_value, Number) else Abs(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        # TODO: implement piecewise and this
        # return Piecewise(expr)
        raise NotImplementedError


@inherit_docstrings
class Sin(UnaryFunction):
    name = "Sin"
    func_key = FUNC_SIN

    # def _doit(self) -> ExprLike:
    #     arg = doit(self.arg)
    #     # TODO: simplify

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return sin(arg_value) if isinstance(arg_value, Number) else Sin(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return Cos(expr)


@inherit_docstrings
class Cos(UnaryFunction):
    name = "Cos"
    func_key = FUNC_COS

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return cos(arg_value) if isinstance(arg_value, Number) else Cos(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return -Sin(expr)


@inherit_docstrings
class Tan(UnaryFunction):
    name = "Tan"
    func_key = FUNC_TAN

    def _expand_trig(self) -> Union["Expr", Number]:
        return Sin(self.arg) / Cos(self.arg)

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return tan(arg_value) if isinstance(arg_value, Number) else Tan(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return Sec(expr) ** 2


@inherit_docstrings
class Sec(UnaryFunction):
    name = "Sec"
    func_key = FUNC_SEC

    def _expand_trig(self) -> Union["Expr", Number]:
        return 1 / Cos(self.arg)

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return (1 / cos(arg_value)) if isinstance(arg_value, Number) else Sec(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return Sec(expr) * Tan(expr)


@inherit_docstrings
class Csc(UnaryFunction):
    name = "Csc"
    func_key = FUNC_CSC

    def _expand_trig(self) -> Union["Expr", Number]:
        return 1 / Sin(self.arg)

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return (1 / sin(arg_value)) if isinstance(arg_value, Number) else Csc(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return -Csc(expr) * Cot(expr)


@inherit_docstrings
class Cot(UnaryFunction):
    name = "Cot"
    func_key = FUNC_COT

    def _expand_trig(self) -> Union["Expr", Number]:
        return Cos(self.arg) / Sin(self.arg)

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return (1 / tan(arg_value)) if isinstance(arg_value, Number) else Cot(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return -Csc(expr) ** 2


@inherit_docstrings
class Arcsin(UnaryFunction):
    name = "Arcsin"
    func_key = FUNC_ARCSIN

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return asin(arg_value) if isinstance(arg_value, Number) else Arcsin(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return 1 / (1 - expr ** 2).sqrt()


@inherit_docstrings
class Arccos(UnaryFunction):
    name = "Arccos"
    func_key = FUNC_ARCCOS

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return acos(arg_value) if isinstance(arg_value, Number) else Arccos(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return -1 / (1 - expr ** 2).sqrt()


@inherit_docstrings
class Arctan(UnaryFunction):
    name = "Arctan"
    func_key = FUNC_ARCTAN

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return atan(arg_value) if isinstance(arg_value, Number) else Arctan(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return 1 / (1 + expr ** 2)


@inherit_docstrings
class Arcsec(UnaryFunction):
    name = "Arcsec"
    func_key = FUNC_ARCSEC

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return acos(1 / arg_value) if isinstance(arg_value, Number) else Arcsec(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return 1 / Abs(expr) / (expr ** 2 - 1).sqrt()


@inherit_docstrings
class Arccsc(UnaryFunction):
    name = "Arccsc"
    func_key = FUNC_ARCCSC

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return asin(1 / arg_value) if isinstance(arg_value, Number) else Arccsc(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return -1 / Abs(expr) / (expr ** 2 - 1).sqrt()


@inherit_docstrings
class Arccot(UnaryFunction):
    name = "Arccot"
    func_key = FUNC_ARCCOT

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return atan(1 / arg_value) if isinstance(arg_value, Number) else Arccot(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return -1 / (1 + expr ** 2)


@inherit_docstrings
class Ln(UnaryFunction):
    name = "Ln"
    func_key = FUNC_LN

    def _expand_log(self) -> Union["Expr", Number]:
        if isinstance(self.arg, Mul):
            return sum(
                ([] if self.arg.coef == 1 else [Ln(self.arg.coef)]) +
                [Ln(factor).expand_log() for factor in self.arg.factors]
            )
        elif isinstance(self.arg, Pow):
            return self.arg.exp * Ln(self.arg.base)
        else:
            return self

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return log(arg_value) if isinstance(arg_value, Number) else Ln(arg_value)

    def _basediff(self, expr: ExprLike, /) -> ExprLike:
        return 1 / expr


@inherit_docstrings
class Log(BinaryFunction):
    name = "Log"
    func_key = FUNC_LOG

    def _expand_log(self) -> Union["Expr", Number]:
        if isinstance(self.arg1, Mul):
            return sum(
                ([] if self.arg1.coef == 1 else [Log(self.arg1.coef, self.arg2)]) +
                [Log(factor, self.arg2).expand_log()
                 for factor in self.arg1.factors]
            )
        elif isinstance(self.arg1, Pow):
            return self.arg1.exp * Log(self.arg1.base, self.arg2)
        else:
            return self

    def _evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg1_value = evalf(self.arg1, value_map)
        arg2_value = evalf(self.arg2, value_map)
        if isinstance(arg1_value, Number) and isinstance(arg2_value, Number):
            return log(arg1_value, arg2_value)
        else:
            return Log(arg1_value, arg2_value)

    def _basediff(self, var: Var, /) -> ExprLike:
        return (Ln(self.arg1) / Ln(self.arg2))._diff()


@inherit_docstrings
class Limit(Expr):
    expr: ExprLike
    var: Var
    point: ExprLike
    direction: int = 0

    def __init__(self, expr: ExprLike, var: Var,
                 point: ExprLike, direction: int = 0, /):
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


@inherit_docstrings
class Derivative(Expr):
    expr: ExprLike
    var: Var
    order: int = 1

    def __init__(self, expr: ExprLike, var: Var,
                 order: int = 1, /):
        self.expr = expr
        self.var = var
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


@inherit_docstrings
class Integral(Expr):
    expr: ExprLike
    var: Var
    a: ExprLike | None = None
    b: ExprLike | None = None

    def __init__(self, expr: ExprLike, var: Var,
                 a: ExprLike | None = None, b: ExprLike | None = None, /):
        self.expr = expr
        self.var = var
        if (a is None) != (b is None):
            raise ValueError("Integrals must either have both"
                             " bounds defined or undefined")
        self.a = a
        self.b = b

    def __str__(self):
        bound_str = '' if self.a is None else f"[{self.a}, {self.b}]"
        return f"∫{bound_str} {self.expr} d{self.var}"

    def __repr__(self):
        bound_str = '' if self.a is None else f"[{self.a!r}, {self.b!r}]"
        return f"∫{bound_str} {self.expr!r} d{self.var}"

    def exprhash(self) -> CanonicalKey:
        if self.a is None:
            return (KEY_INT, exprhash(self.expr), exprhash(self.var))
        else:
            return (KEY_INT, exprhash(self.expr), exprhash(self.var),
                    exprhash(self.a), exprhash(self.b))

    def apply(self, func: FunctionType, *args) -> any:
        return func(
            Integral(apply(self.expr, func, *args),
                     apply(self.var, func, *args),
                     apply(self.a, func, *args),
                     apply(self.b, func, *args)), *args,
        )

    def _diff(self, var: Var, /) -> ExprLike:
        if self.var == var and self.a is None and self.b is None:
            return self.expr
        else:
            result = doit(self)
            return result if isinstance(result, Integral) else result._diff(var)


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
def is_constant(expr: ExprLike, var: Var | None = None, /) -> bool:
    """
    Checks of the expression is a constant value, aka. a structure
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
        return all(is_constant(term) for term in expr.terms)
    elif isinstance(expr, Mul):
        return all(is_constant(factor) for factor in expr.factors)
    elif isinstance(expr, Pow):
        return is_constant(expr.base) and is_constant(expr.exp)
    elif isinstance(expr, UnaryFunction):
        return is_constant(expr.arg)
    elif isinstance(expr, BinaryFunction):
        return is_constant(expr.arg1) and is_constant(expr.arg2)
    elif isinstance(expr, Limit):
        return is_constant(expr.expr) and is_constant(expr.var) and is_constant(expr.point)
    elif isinstance(expr, Derivative):
        return is_constant(expr.expr) and is_constant(expr.var)
    elif isinstance(expr, Integral):
        return is_constant(expr.expr) and is_constant(expr.var) and is_constant(expr.a) and is_constant(expr.b)
    else:
        return False


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
    if isinstance(expr, Number):
        return expr
    elif isinstance(expr, Expr):
        return expr.apply(func, *args)


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
    if isinstance(expr, Number):
        return expr
    elif isinstance(expr, Expr):
        return expr.subs(expr_map)


@typechecked
def evalf(expr: ExprLike, value_map: ValueMap | None = None, /) -> ExprLike:
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
    :type mapping: ValueMap (dict[Var, Number])
    :rtype: ExprLike
    :raises TypeError: If any argument doesn't match expected type.
    :raises NotImplementedError: If any subclass doesn't have
                                 .evalf() method implemented.
    """
    if isinstance(expr, Number):
        return expr
    elif isinstance(expr, Expr):
        return expr.evalf(value_map)


# factor_mul is basically using linear/quadratic/cubic/quartic formulas
# on one variable for the solutions for zero value and writing the
# polynomial in factor form.
# quintic and up are still doable if polynomial is single variable
# but at some point it's pointless to factor
@typechecked
def factor(expr: ExprLike, /) -> ExprLike:
    if isinstance(expr, Number):
        return expr
    return expr


@typechecked
def expand_mul(expr: ExprLike, /) -> ExprLike:
    """Distribute compact algebra in multiplication."""
    return expr.expand_mul() if isinstance(expr, Expr) else expr


@typechecked
def expand_pow(expr: ExprLike, /) -> ExprLike:
    """Distribute compact algebra in power."""
    return expr.expand_pow() if isinstance(expr, Expr) else expr


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
        expr = expand_mul(expr)
        if original == expr:
            break
    return expr


# TODO: expand to add, mul, pow, trig, log, func
# AKA factor
@typechecked
def reduce(expr: ExprLike, /) -> ExprLike:
    """Combine/factor flattened algebra in multiplication, power, and so on."""
    if isinstance(expr, Number):
        return expr
    return expr


# Reduce rational
@typechecked
def cancel(expr: ExprLike, /) -> ExprLike:
    if isinstance(expr, Number):
        return expr
    return expr


# Combine over common denominator
@typechecked
def together(expr: ExprLike, /) -> ExprLike:
    if isinstance(expr, Number):
        return expr
    return expr


# Partial fraction decomposition
@typechecked
def apart(expr: ExprLike, /) -> ExprLike:
    if isinstance(expr, Number):
        return expr
    return expr


# grouping the expression by variable
@typechecked
def collect(expr: ExprLike, var: Var, /) -> ExprLike:
    if isinstance(expr, Number):
        return expr
    return expr


@typechecked
def simplify(expr: ExprLike, /) -> ExprLike:
    if isinstance(expr, Number):
        return expr
    return expr


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


def main():
    x, y, z = symbols("xyz")
    expr = sum(x ** i for i in range(5))
    print(expr)
    print(expr.evalf())
    expr2 = x + y + z
    print(expr2)
    print(expr2.evalf())
    print(expr2.evalf({x: 1}))


if __name__ == "__main__":
    main()
