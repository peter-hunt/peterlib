from functools import reduce as func_reduce
from numbers import Number, Integral, Rational, Real, Complex
from types import FunctionType
from typing import Iterable, Union as TypingUnion
from operator import and_, or_

from typeguard import typechecked

from expr import (
    CanonicalKey, Expr, Constant, Var, exprhash, is_infinite, is_finite, get_value,
    inf, undefined,
)
from tags import *
from utils import inherit_docstrings


__all__ = [
    "Set", "SetLike", "SetMap",
    "sethash", "setsorted", "subs", "reduce"
    "EmptySet", "UniversalSet", "FiniteSet",
    "FiniteSetType", "finite_set",
    "Union", "Intersection", "Difference", "Complement",
    "NumberDomain", "NaturalSet", "IntegerSet", "RationalSet", "RealSet", "ComplexSet",
    "Empty", "Universal", "NSet", "ZSet", "QSet", "RSet", "CSet",
    "ConstInterval",
]


@typechecked
class Set:
    def __init__(self, *args, **kwargs):
        raise Exception(
            "Cannot create an instance of base Set class. "
            "To use a subclass, redefine the init code."
        )

    def sethash(self) -> CanonicalKey:
        """
        Hash tuple generator for expressions and sets.
        Use the global function for supported numbers and expressions.
        This is used for sorting and hashing expressions and sets.

        :param self: The set instance.
        :return: The tuple of relevant information.
        :rtype: tuple[int | tuple, ...]
        """
        raise NotImplementedError

    def __hash__(self):
        return hash(self.sethash())

    def __eq__(self, other: TypingUnion["Set", Expr, Number]) -> bool:
        if isinstance(other, (Set, Expr, Number)):
            return sethash(self) == sethash(other)
        else:
            return NotImplemented

    def __sub__(self, other: "Set") -> "Set":
        if not isinstance(other, Set):
            return NotImplemented
        if isinstance(self, FiniteSetType) and isinstance(other, FiniteSetType):
            return Difference(self, other).reduce()
        elif isinstance(self, EmptySet) or isinstance(other, UniversalSet):
            return EmptySet()
        return Difference(self, other)

    def __and__(self, other: "Set") -> "Set":
        if not isinstance(other, Set):
            return NotImplemented
        if isinstance(self, FiniteSetType) and isinstance(other, FiniteSetType):
            return Intersection(self, other).reduce()
        elif isinstance(self, EmptySet) or isinstance(other, EmptySet):
            return EmptySet()
        return Intersection(self, other)

    def __or__(self, other: "Set") -> "Set":
        if not isinstance(other, Set):
            return NotImplemented
        if isinstance(self, FiniteSetType) and isinstance(other, FiniteSetType):
            return Union(self, other).reduce()
        elif isinstance(self, UniversalSet) or isinstance(other, UniversalSet):
            return UniversalSet()
        return Union(self, other)

    def __contains__(self, item: any) -> bool:
        return NotImplemented

    def apply(self, func: FunctionType, *args) -> "Set":
        """
        Recursively apply some functionality to the set instance.
        The function usually refers to some internally implemented
        method different for each subclass.
        Used for substitution, reduction, and so on.

        :param self: The set.
        :param func: The function to recursively apply.
        :type func: FunctionType
        :param args: Related arguments to pass to the function if any.
        :return: The result of the recursive function call.
        :rtype: Set
        """
        return func(self, *args)

    def subs(self, set_map: TypingUnion["Set", Expr, Number] | None = None) -> "Set":
        return self.apply(lambda x: x if isinstance(x, Number) else x._subs(set_map))

    def _subs(self, set_map: TypingUnion["Set", Expr, Number] | None = None) -> "Set":
        return self

    def reduce(self) -> "Set":
        return self.apply(lambda x: x if isinstance(x, Number) else x._reduce())

    def _reduce(self) -> "Set":
        return self


SetLike = Set | Expr | Number
SetMap = dict[Var, SetLike]


@typechecked
def sethash(item: SetLike) -> CanonicalKey:
    """
    Hash tuple generator for expressions and sets.
    Use the global function for supported numbers and expressions.
    This is used for sorting and hashing expressions and sets.

    :param item: The set, expression, or nnumber to hash.
    :type item: Set
    :return: The tuple of relevant information.
    :rtype: CanonicalKey
    """
    return item.sethash() if isinstance(item, Set) else exprhash(item)


@typechecked
def setsorted(iterable: Iterable[SetLike], /) -> list[SetLike]:
    """
    Sort iterable of expressions and sets by descending sethash.

    :param iterable: The iterable of expressions and sets to sort.
    :type iterable: Iterable[Set]
    :return: The list of sorted expressions and sets in descending order.
    :rtype: list[Set]
    """
    return sorted(iterable, key=sethash, reverse=True)


@typechecked
def subs(set_, set_map: SetMap, /) -> Set:
    return set_.subs(set_map) if isinstance(set_, Set | Expr) else set_


@typechecked
def reduce(set_, /) -> Set:
    return set_.reduce() if isinstance(set_, Set) else set_


@typechecked
@inherit_docstrings
class UniversalSet(Set):
    def __init__(self):
        pass

    def __repr__(self):
        return 'U'

    def sethash(self) -> CanonicalKey:
        return (KEY_SET_UNVSL,)

    def __contains__(self, item: SetLike) -> bool:
        return item != undefined


@typechecked
@inherit_docstrings
class EmptySet(Set):
    elements = []

    def __init__(self):
        pass

    def __repr__(self):
        return '∅'

    def sethash(self) -> CanonicalKey:
        return (KEY_SET_EMPTY,)

    def __iter__(self):
        return iter(self.elements)

    def __contains__(self, item: SetLike) -> bool:
        return False


@typechecked
@inherit_docstrings
class FiniteSet(Set):
    elements: list[SetLike]

    def __init__(self, *elements: SetLike):
        if len(elements) == 0:
            return ValueError("empty sets should be instantiated"
                              " with EmptySet type instead.")
        if inf in elements or -inf in elements:
            return ValueError("FiniteSet cannot have infinity values.")
        if undefined in elements:
            return ValueError("FiniteSet cannot have undefined values.")
        self.elements = setsorted(elements)

    def __repr__(self):
        return "FiniteSet(" + ", ".join(f"{element!r}" for element in self.elements) + ")"

    def __str__(self):
        return "{" + ", ".join(f"{element}" for element in self.elements) + "}"

    def sethash(self) -> CanonicalKey:
        return (KEY_SET_UNVSL,) + tuple(sethash(element) for element in self.elements)

    def __iter__(self):
        return iter(self.elements)

    def __contains__(self, item: SetLike) -> bool:
        return item in self.elements

    def apply(self, func: FunctionType, *args) -> "Set":
        return func(FiniteSet(func(element, *args) for element in self.elements), *args)


FiniteSetType = EmptySet | FiniteSet


def finite_set(*elements: SetLike) -> EmptySet | FiniteSet:
    """
    Construct either an EmptySet or a FiniteSet based on elements given.

    :param elements: The elements of the list
    :type elements: list[SetLike]
    :return: Description
    :rtype: EmptySet | FiniteSet
    """
    return EmptySet() if len(elements) == 0 else FiniteSet(*elements)


@typechecked
@inherit_docstrings
class Union(Set):
    sets: list[Set]

    def __init__(self, *sets: Set):
        if len(sets) == 0:
            raise ValueError("Union instance must have"
                             " at least one set")
        self.sets = setsorted({*sets})

    def __repr__(self):
        return "Union(" + ", ".join(f"{set_!r}" for set_ in self.sets) + ")"

    def __str__(self):
        return " | ".join(f"{set_}" for set_ in self.sets)

    def __contains__(self, item: SetLike) -> bool:
        return any(item in set_ for set_ in self.sets)

    def sethash(self) -> CanonicalKey:
        return (KEY_SET_UNION,) + tuple(sethash(set_) for set_ in self.sets)

    def apply(self, func: FunctionType, *args) -> "Set":
        return func(Union(*(func(set_, *args) for set_ in self.sets)), *args)

    def _reduce(self) -> "Set":
        if len(self.sets) == 1:
            return self.sets[0]
        elif Universal in self.sets:
            return Universal
        elif all(isinstance(set_, FiniteSetType) for set_ in self.sets):
            return finite_set(*func_reduce(or_, ({*set_.elements} for set_ in self.sets), *{}))
        else:
            return self


@typechecked
@inherit_docstrings
class Intersection(Set):
    sets: list[Set]

    def __init__(self, *sets: Set):
        if len(sets) == 0:
            raise ValueError("Intersection instance must have"
                             " at least one set")
        self.sets = setsorted({*sets})

    def __repr__(self):
        return "Intersection(" + ", ".join(f"{set_!r}" for set_ in self.sets) + ")"

    def __str__(self):
        return " & ".join(f"{set_}" for set_ in self.sets)

    def __contains__(self, item: SetLike) -> bool:
        return all(item in set_ for set_ in self.sets)

    def sethash(self) -> CanonicalKey:
        return (KEY_SET_INTSC,) + tuple(sethash(set_) for set_ in self.sets)

    def apply(self, func: FunctionType, *args) -> "Set":
        return func(Intersection(*(func(set_, *args) for set_ in self.sets)), *args)

    def _reduce(self) -> "Set":
        if len(self.sets) == 1:
            return self.sets[0]
        elif Empty in self.sets:
            return Empty
        elif all(isinstance(set_, FiniteSetType) for set_ in self.sets):
            return finite_set(*func_reduce(and_, ({*set_.elements} for set_ in self.sets), *{}))
        else:
            return self


@typechecked
@inherit_docstrings
class Difference(Set):
    left: Set
    right: Set

    def __init__(self, left: Set, right: Set, /):
        self.left = left
        self.right = right

    def __repr__(self):
        return f"Difference({self.left!r}, {self.right!r})"

    def __str__(self):
        return f"{self.left} - {self.right}"

    def __contains__(self, item: SetLike) -> bool:
        return item in self.left and item not in self.right

    def sethash(self) -> CanonicalKey:
        return (KEY_SET_DFRNC, sethash(self.left), sethash(self.right))

    def apply(self, func: FunctionType, *args) -> "Set":
        return func(Difference(func(self.left, *args), func(self.right, *args)), *args)

    def _reduce(self) -> "Set":
        if self.left == Empty or self.right == Universal:
            return Empty
        elif isinstance(self.left, FiniteSetType) and isinstance(self.right, FiniteSetType):
            return finite_set(*(element for element in self.left if element not in self.right))
        else:
            return self


@typechecked
@inherit_docstrings
class Complement(Set):
    set_: Set

    def __init__(self, set_: Set, /):
        self.set_ = set_

    def __repr__(self):
        return f"Complement({self.set_!r})"

    def __str__(self):
        return f"{self.set_}'"

    def __contains__(self, item: SetLike) -> bool:
        return item not in self.set_

    def sethash(self) -> CanonicalKey:
        return (KEY_SET_DFRNC, sethash(self.left), sethash(self.right))

    def apply(self, func: FunctionType, *args) -> "Set":
        return func(Complement(func(self.set_, *args)), *args)

    def _reduce(self) -> "Set":
        if self.set_ == Empty:
            return Universal
        elif self.set_ == Universal:
            return Empty
        else:
            return self


@typechecked
@inherit_docstrings
class NumberDomain(Set):
    letter: str

    def __init__(self, /):
        raise Exception(
            "Cannot create an instance of base NumberDomain class. "
            "To use a subclass, redefine the init code."
        )

    def __repr__(self):
        return f"{type(self).__name__}()"

    def __str__(self):
        return self.letter

    def sethash(self) -> CanonicalKey:
        return (KEY_SET_DOMAIN, self.domain_key)


@typechecked
@inherit_docstrings
class NaturalSet(NumberDomain):
    letter: str = 'ℕ'
    domain_key = DOMAIN_NATURAL

    def __init__(self, /):
        pass

    def __contains__(self, item: SetLike) -> bool:
        return isinstance(item, Integral) and item >= 0


@typechecked
@inherit_docstrings
class IntegerSet(NumberDomain):
    letter: str = 'ℤ'
    domain_key = DOMAIN_INTEGER

    def __init__(self, /):
        pass

    def __contains__(self, item: SetLike) -> bool:
        return isinstance(item, Integral)


@typechecked
@inherit_docstrings
class RationalSet(NumberDomain):
    letter: str = 'ℚ'
    domain_key = DOMAIN_RATIONAL

    def __init__(self, /):
        pass

    def __contains__(self, item: SetLike) -> bool:
        return isinstance(item, Rational)


@typechecked
@inherit_docstrings
class RealSet(NumberDomain):
    letter: str = 'ℝ'
    domain_key = DOMAIN_REAL

    def __init__(self, /):
        pass

    def __contains__(self, item: SetLike) -> bool:
        return isinstance(item, Real | Constant) and is_finite(item)


@typechecked
@inherit_docstrings
class ComplexSet(NumberDomain):
    letter: str = 'ℂ'
    domain_key = DOMAIN_COMPLEX

    def __init__(self, /):
        pass

    def __contains__(self, item: SetLike) -> bool:
        return isinstance(item, Complex | Constant) and is_finite(item)


Empty = EmptySet()
Universal = UniversalSet()
NSet = NaturalSet()
ZSet = IntegerSet()
QSet = RationalSet()
RSet = RealSet()
CSet = ComplexSet()


@typechecked
@inherit_docstrings
class ConstInterval(Set):
    a: Real | Constant  # lower bound
    b: Real | Constant  # upper bound
    a_inc: bool = True  # being false doesn't exclude infinity
    b_inc: bool = True  # being false doesn't exclude infinity
    domain: NumberDomain = CSet

    def __init__(self, a: Real | Constant, b: Real | Constant,
                 a_inc: bool = True, b_inc: bool = True,
                 domain: NumberDomain | None = None, /):
        if a == undefined or b == undefined:
            raise ValueError("cannot have undefined as bounds")
        if is_infinite(a) and is_infinite(b):
            if a == inf and b == -inf:
                return ValueError("a cannot be greater than b: Infinity > -Infinity")
            elif a == b:
                return ValueError(f"a and b cannot be the same infinity: {a} = {b}")
        elif is_infinite(a) and is_finite(b):
            if a == inf:
                return ValueError(f"a cannot be greater than b: Infinity > {b}")
            a_inc = True
        elif is_finite(a) and is_infinite(b):
            if b == -inf:
                return ValueError(f"a cannot be greater than b: {a} > -Infinity")
            b_inc = True
        else:
            if get_value(a) > get_value(b):
                return ValueError(f"a cannot be greater than b: {a} > {b}")
            elif get_value(a) == get_value(b) and (not a_inc or not b_inc):
                raise ValueError(f"single point interval must have"
                                 f" both boundaries inclusive: {a_inc}, {b_inc}")

        self.a = a
        self.b = b
        self.a_inc = a_inc
        self.b_inc = b_inc
        self.domain = CSet if domain is None else domain

    def __repr__(self):
        return f"ConstInterval({self.a!r}, {self.b!r}, {self.a_inc}, {self.b_inc}, {self.domain!r})"

    def __str__(self):
        left_bracket = '[' if self.a_inc else '('
        right_bracket = ']' if self.b_inc else ')'
        return f"Interval({left_bracket}{self.a}, {self.b}{right_bracket})"

    def sethash(self) -> CanonicalKey:
        return (KEY_SET_INTVL, exprhash(self.a), exprhash(self.b),
                int(self.a_inc), int(self.b_inc))

    def __contains__(self, item: SetLike) -> bool:
        if is_finite(item) and item not in self.domain:
            return False
        elif is_infinite(self.a) and is_infinite(self.b):
            return True
        elif is_infinite(self.a) and is_finite(self.b):
            return item == -inf or get_value(item) < get_value(self.b) or (get_value(item) == get_value(self.b) and self.b_inc)
        elif is_finite(self.a) and is_infinite(self.b):
            return item == inf or get_value(item) > get_value(self.a) or (get_value(item) == get_value(self.a) and self.a_inc)
        else:
            return ((get_value(item) < get_value(self.b) or (get_value(item) == get_value(self.b) and self.b_inc)) and
                    (get_value(item) > get_value(self.a) or (get_value(item) == get_value(self.a) and self.a_inc)))


def main():
    a = finite_set(1, 2, 3)
    b = finite_set(3, 4, 5)
    print(a - b)
    print(a & b)
    print(a | b)
    print(ConstInterval(1, 2))
    print(1.5 in ConstInterval(1, 2))
    print(ConstInterval(1, inf))
    print(100 in ConstInterval(1, inf))
    print(inf in ConstInterval(1, inf))


if __name__ == "__main__":
    main()
