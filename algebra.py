from math import e as math_e, pi as math_pi, sqrt, prod, sin, cos, tan, log, asin, acos, atan
from numbers import Number
from types import FunctionType
from typing import Union

from typeguard import typechecked


VARIABLE_LETTERS = (
    "abcdfghijklmnopqrstuvwxyz"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "αβγδεζηθικλμνξοπρστυφχψω"
    "ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ"
)

RESERVED_LETTERS = "eπφΣΠ"


@typechecked
class Expr:
    def __init__(self, *args, **kwargs):
        raise Exception(
            "Cannot create an instance of base Expr class. "
            "To use a subclass, redefine the init code."
        )

    def apply(self, func: FunctionType, *args, **kwargs) -> any:
        return func(self, *args, **kwargs)

    def doit(self) -> Union["Expr", Number]:
        return self

    def subs(self, expr_map: dict["Var", Union["Expr", Number]] | None = None, /) -> Union["Expr", Number]:
        return self

    def evalf(self, value_map: dict["Var", Number] | None = None, /) -> Union["Expr", Number]:
        return self

    def __repr__(self) -> str:
        raise NotImplementedError

    def __str__(self) -> str:
        raise NotImplementedError

    def __hash__(self):
        raise NotImplementedError

    def __add__(self, other: Union["Expr", Number]) -> "Expr":
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

    def __radd__(self, other: Union["Expr", Number]) -> "Expr":
        if other == 0:
            return self
        elif isinstance(self, Add):
            return Add(other, self.const, *self.terms)
        else:
            return Add(other, self)

    def __sub__(self, other: Union["Expr", Number]) -> "Expr":
        if other == 0:
            return self
        elif isinstance(self, Add):
            return Add(Mul(-1, other), self.const, *self.terms)
        else:
            return Add(self, Mul(-1, other))

    def __rsub__(self, other: Union["Expr", Number]) -> "Expr":
        if other == 0:
            return Mul(-1, self)
        else:
            return Add(other, Mul(-1, self))

    def __mul__(self, other: Union["Expr", Number]) -> "Expr":
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

    def __rmul__(self, other: Union["Expr", Number]) -> "Expr":
        if other == 1:
            return self
        elif isinstance(self, Mul):
            return Mul(other, self.coef, *self.factors)
        else:
            return Mul(other, self)

    def __div__(self, other: Union["Expr", Number]) -> "Expr":
        if other == 1:
            return self
        elif isinstance(self, Mul):
            return Mul(self.const, *self.terms, Pow(other, -1))
        else:
            return Mul(self, Pow(other, -1))

    def __rdiv__(self, other: Union["Expr", Number]) -> "Expr":
        if other == 1:
            return Pow(self, -1)
        elif isinstance(other, Mul):
            return Mul(other.const, *other.terms, Pow(self, -1))
        else:
            return Mul(other, Pow(self, -1))

    def __pow__(self, other: Union["Expr", Number]) -> "Expr":
        if other == 0:
            # TODO: check for 0^0
            return 1
        elif other == 1:
            return self
        elif isinstance(self, Pow):
            return Pow(self.base, self.exp * other)
        else:
            return Pow(self, other)

    def __rpow__(self, other: Union["Expr", Number]) -> "Expr":
        # TODO: check for 0^0
        if other == 0 or other == 1:
            return other
        else:
            return Pow(other, self)


ExprLike = Expr | Number


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

    def __hash__(self):
        return hash(self.name)

    def subs(self, expr_map: dict["Var", ExprLike] | None = None, /) -> ExprLike:
        return (expr_map or {}).get(self, self)

    def evalf(self, value_map: dict["Var", Number] | None = None, /) -> ExprLike:
        return (value_map or {}).get(self, self)


ExprMap = dict[Var, ExprLike]
ValueMap = dict[Var, Number]


def doit(expr: ExprLike, /) -> ExprLike:
    """
    Recursively apply Expr.doit() method with ExprLike object.
    This function executes symbolic operations, like derivatives and
    integrals, and evaluates closed-form expressions, like Sin(π).

    :param expr: The expression to execute symbolic operations on.
    :type expr: ExprLike (Expr | Number)
    :rtype: ExprLike (Expr | Number)
    :raises TypeError: If any argument doesn't match expected type.
    """
    if isinstance(expr, Number):
        return expr
    elif isinstance(expr, Expr):
        return expr.doit()
    else:
        raise TypeError(f"must be expression or number,"
                        f" not {type(expr).__name__}")


def subs(expr: ExprLike, expr_map: ExprMap | None = None, /) -> ExprLike:
    """
    Recursively apply Expr.subs(mapping) method with ExprLike object.
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
    else:
        raise TypeError(f"must be expression or number,"
                        f" not {type(expr).__name__}")


def evalf(expr: ExprLike, value_map: ValueMap | None = None, /) -> ExprLike:
    """
    Recursively apply Expr.evalf(mapping) method with ExprLike object.
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
    else:
        raise TypeError(f"must be expression or number,"
                        f" not {type(expr).__name__}")


def symbols(letters, /):
    return [Var(letter) for letter in letters]


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

    def __hash__(self):
        return hash(self.name)

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return self

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        return self.value


e = Constant('e', math_e)
pi = π = Constant('π', math_pi)
phi = φ = Constant('φ', (1 + sqrt(5)) / 2)


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

    def __repr__(self):
        return (
            "Add(" + ", ".join(f"{term!r}" for term in self.terms) +
            (f", {self.const}" if self.const != 0 else "") + ")"
        )

    def __str__(self):
        return (
            " + ".join(f"{term}" for term in self.terms) +
            (f" + {self.const}" if self.const != 0 else "")
        )

    def doit(self) -> ExprLike:
        return sum(doit(term) for term in self.terms) + self.const

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return sum(subs(term, expr_map) for term in self.terms) + self.const

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        return sum(evalf(term, value_map) for term in self.terms) + self.const


class Mul(Expr):
    factors: list[Expr]
    coef: Number = 1

    def __init__(self, *factors: list[ExprLike]):
        self.factors = []
        self.coef = 1
        for factor in factors:
            if isinstance(factor, Number):
                self.coef *= factor
            else:
                self.factors.append(factor)

    def __repr__(self):
        return (
            "Mul(" + ", ".join(f"{factor!r}" for factor in self.factors) +
            (f", {self.coef}" if self.coef != 1 else "") + ")"
        )

    def __str__(self):
        return (
            " * ".join(
                f"({factor})" if isinstance(factor, Add)
                else f"{factor}"
                for factor in self.factors) +
            (f" * {self.coef}" if self.coef != 1 else "")
        )

    def doit(self) -> ExprLike:
        return prod(doit(term) for term in self.factors) * self.coef

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return prod(subs(term, expr_map) for term in self.factors) * self.coef

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        return prod(evalf(term, value_map) for term in self.factors) * self.coef


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
        if isinstance(self.exp, Add | Mul):
            exp_str = f"({exp_str})"
        return f"{base_str}^{exp_str}"

    def doit(self) -> ExprLike:
        return doit(self.base) ** doit(self.exp)

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return subs(self.base, expr_map) ** subs(self.exp, expr_map)

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        return evalf(self.base, value_map) ** evalf(self.exp, value_map)


class Function(Expr):
    name: str


class UnaryFunction(Function):
    name: str
    arg: ExprLike

    def __init__(self, arg: ExprLike, /):
        self.arg = arg

    def __repr__(self):
        return f"{self.name}({self.arg!r})"

    def __str__(self):
        return f"{self.name}({self.arg})"


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


class Sin(UnaryFunction):
    name = "Sin"

    def doit(self) -> ExprLike:
        arg = doit(self.arg)
        # TODO

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return Sin(subs(self.arg, expr_map))

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return sin(arg_value) if isinstance(arg_value, Number) else Sin(arg_value)


class Cos(UnaryFunction):
    name = "Cos"

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return Cos(subs(self.arg, expr_map))

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return cos(arg_value) if isinstance(arg_value, Number) else Cos(arg_value)


class Tan(UnaryFunction):
    name = "Tan"

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return Tan(subs(self.arg, expr_map))

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return tan(arg_value) if isinstance(arg_value, Number) else Tan(arg_value)


class Sec(UnaryFunction):
    name = "Sec"

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return Sec(subs(self.arg, expr_map))

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return (1 / cos(arg_value)) if isinstance(arg_value, Number) else Sec(arg_value)


class Csc(UnaryFunction):
    name = "Csc"

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return Csc(subs(self.arg, expr_map))

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return (1 / sin(arg_value)) if isinstance(arg_value, Number) else Csc(arg_value)


class Cot(UnaryFunction):
    name = "Cot"

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return Cot(subs(self.arg, expr_map))

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return (1 / tan(arg_value)) if isinstance(arg_value, Number) else Cot(arg_value)


class Ln(UnaryFunction):
    name = "Ln"

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return Ln(subs(self.arg, expr_map))

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return log(arg_value) if isinstance(arg_value, Number) else Ln(arg_value)


class Arcsin(UnaryFunction):
    name = "Arcsin"

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return Arcsin(subs(self.arg, expr_map))

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return asin(arg_value) if isinstance(arg_value, Number) else Arcsin(arg_value)


class Arccos(UnaryFunction):
    name = "Arccos"

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return Arccos(subs(self.arg, expr_map))

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return acos(arg_value) if isinstance(arg_value, Number) else Arccos(arg_value)


class Arctan(UnaryFunction):
    name = "Arctan"

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return Arctan(subs(self.arg, expr_map))

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg_value = evalf(self.arg, value_map)
        return atan(arg_value) if isinstance(arg_value, Number) else Arctan(arg_value)


class Log(BinaryFunction):
    name = "Log"

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        return Log(subs(self.arg1, expr_map), subs(self.arg2, expr_map))

    def evalf(self, value_map: ValueMap, /) -> ExprLike:
        arg1_value = evalf(self.arg1, value_map)
        arg2_value = evalf(self.arg2, value_map)
        if isinstance(arg1_value, Number) and isinstance(arg2_value, Number):
            return log(arg1_value, arg2_value)
        else:
            return Log(arg1_value, arg2_value)


class Limit(Expr):
    expr: ExprLike
    var: Var
    point: ExprLike
    direction: int = 0

    def __init__(self, expr: ExprLike, var: Var,
                 point: ExprLike, direction: int = 0):
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

    def subs(self, expr_map: ExprMap, /) -> ExprLike:
        new_var = subs(self.var, expr_map)
        if not isinstance(new_var, Var):
            pass
        return Limit(
            subs(self.expr, expr_map),
            subs(self.var, expr_map),
            subs(self.point, expr_map),
            self.direction,
        )


class Derivative(Expr):
    expr: ExprLike
    var: Var
    order: int = 1

    def __str__(self):
        order_str = '' if self.order == 1 else f"^{self.order}"
        return f"d{order_str}/d{self.var}{order_str} ({self.expr})"

    def __repr__(self):
        order_str = '' if self.order == 1 else f"^{self.order}"
        return f"d{order_str}/d{self.var}{order_str} ({self.expr!r})"


class Integral(Expr):
    expr: ExprLike
    var: Var
    a: ExprLike | None = None
    b: ExprLike | None = None

    def __str__(self):
        bound_str = '' if self.a is None else f"[{self.a}, {self.b}]"
        return f"∫{bound_str} {self.expr} d{self.var}"

    def __repr__(self):
        bound_str = '' if self.a is None else f"[{self.a!r}, {self.b!r}]"
        return f"∫{bound_str} {self.expr!r} d{self.var}"


# TODO: sorting and hashing to allow value comparison
# TODO: implement the "doit" for existing functions, namely trigs

# TODO: abs/piecewise

# TODO: (in)equality

# expand_complex is not planned

# TODO: expand to mul, pow, trig, log, func
def expand(expr: ExprLike, /) -> ExprLike:
    """Distribute compact algebra in multiplication, power, and so on."""
    if isinstance(expr, Number):
        return expr
    return expr


# TODO: expand to add, mul, pow, trig, log, func
def combine(expr: ExprLike, /) -> ExprLike:
    """Combine flattened algebra in multiplication, power, and so on."""
    if isinstance(expr, Number):
        return expr
    return expr


def factor(expr: ExprLike, /) -> ExprLike:
    if isinstance(expr, Number):
        return expr
    return expr


# Reduce rational
def cancel(expr: ExprLike, /) -> ExprLike:
    if isinstance(expr, Number):
        return expr
    return expr


# Combine over common denominator
def together(expr: ExprLike, /) -> ExprLike:
    if isinstance(expr, Number):
        return expr
    return expr


# Partial fraction decomposition
def apart(expr: ExprLike, /) -> ExprLike:
    if isinstance(expr, Number):
        return expr
    return expr


def simplify(expr: ExprLike, /) -> ExprLike:
    if isinstance(expr, Number):
        return expr
    return expr


def collect(expr: ExprLike, var: Var, /) -> ExprLike:
    if isinstance(expr, Number):
        return expr
    return expr


def main():
    x, y, z = symbols("xyz")
    expr = x + 2 * y
    values = {x: y, y: x}
    print(expr)
    print(expr.subs(values))
    print(expr.evalf({x: 1, y: 2}))


if __name__ == "__main__":
    main()
