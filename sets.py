from functools import reduce as func_reduce
from numbers import Number, Integral, Rational, Real, Complex
from types import FunctionType
from typing import Iterable, Union as TypingUnion
from operator import and_, or_

from typeguard import typechecked

from expr import (
    CanonicalKey, Expr, Constant, Var, exprhash, is_infinite, is_finite,
    inf, undefined,
)
from expr_utils import (
    const_cmp,
    const_lbound_min, const_lbound_max, const_rbound_min, const_rbound_max,
)
from tags import *
from utils import inherit_docstrings


__all__ = [
    "Set", "SetLike", "SetMap",
    "sethash", "setsorted", "subs", "reduce"
    "EmptySet", "UniversalSet", "FiniteSet",
    "FiniteSetType", "finite_set",
    "ConstInterval", "interval",
    "Union", "Intersection", "Difference", "Complement",
    "NumberDomain", "NaturalSet", "IntegerSet", "RationalSet", "RealSet", "ComplexSet",
    "Empty", "Universal", "NumberDomainSets", "NSet", "ZSet", "QSet", "RSet", "CSet",
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
        if isinstance(self, EmptySet) or isinstance(other, UniversalSet):
            return EmptySet()
        if isinstance(self, FiniteSetType) and isinstance(other, FiniteSetType):
            return Difference(self, other).reduce()
        return Difference(self, other)

    def __and__(self, other: "Set") -> "Set":
        if not isinstance(other, Set):
            return NotImplemented
        if isinstance(self, EmptySet) or isinstance(other, EmptySet):
            return EmptySet()
        for type_ in (FiniteSet, ConstInterval, NumberDomain):
            if isinstance(self, type_) and isinstance(other, type_):
                return Intersection(self, other).reduce()
        return Intersection(self, other)

    def __or__(self, other: "Set") -> "Set":
        if not isinstance(other, Set):
            return NotImplemented
        if isinstance(self, UniversalSet) or isinstance(other, UniversalSet):
            return UniversalSet()
        for type_ in (FiniteSet, ConstInterval, NumberDomain):
            if isinstance(self, type_) and isinstance(other, type_):
                return Intersection(self, other).reduce()
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


type SetLike = Set | Expr | Number
type SetMap = dict[Var, SetLike]


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


type FiniteSetType = EmptySet | FiniteSet


def finite_set(*elements: SetLike) -> FiniteSetType:
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
class ConstInterval(Set):
    a: Real | Constant  # lower bound
    b: Real | Constant  # upper bound
    a_inc: bool = True  # being false doesn't exclude infinity
    b_inc: bool = True  # being false doesn't exclude infinity

    def __init__(self, a: Real | Constant, b: Real | Constant,
                 a_inc: bool = True, b_inc: bool = True, /):
        if a == undefined or b == undefined:
            raise ValueError("cannot have undefined as bounds")
        if const_cmp(a, b) == 1:
            return ValueError(f"a cannot be greater than b: {a} > {b}")
        if a == b and is_infinite(a):
            return ValueError(f"a and b cannot be the same infinity: {a} = {b}")
        if const_cmp(a, b) == 0 and (not a_inc or not b_inc):
            raise ValueError(f"single point interval must have"
                             f" both boundaries inclusive: {a_inc}, {b_inc}")

        self.a = a
        self.b = b
        self.a_inc = a_inc or is_infinite(a)
        self.b_inc = b_inc or is_infinite(b)

    def __repr__(self):
        return f"ConstInterval({self.a!r}, {self.b!r}, {self.a_inc}, {self.b_inc})"

    def __str__(self):
        left_bracket = '[' if self.a_inc else '('
        right_bracket = ']' if self.b_inc else ')'
        return f"Interval({left_bracket}{self.a}, {self.b}{right_bracket})"

    def sethash(self) -> CanonicalKey:
        return (KEY_SET_INTVL, exprhash(self.a), exprhash(self.b),
                int(self.a_inc), int(self.b_inc))

    def __contains__(self, item: SetLike) -> bool:
        return (const_cmp(self.a, item) in (-1, 0) if self.a_inc else (-1,) and
                const_cmp(item, self.b) in (-1, 0) if self.a_inc else (-1,))


def interval(a: Real | Constant, b: Real | Constant, /,
             a_inc: bool = True, b_inc: bool = True) -> ConstInterval | EmptySet:
    """
    Creates a ConstInterval with the given bounds and domain.
    If the bound is empty, safety returns an EmptySet.

    :param a: The lower bound of the interval.
    :type a: Real | Constant
    :param b: The upper bound of the interval.
    :type b: Real | Constant
    :param a_inc: Whether if the lower bound is inclusive.
    :type a_inc: bool
    :param b_inc: Whether if the upper bound is inclusive.
    :type b_inc: bool
    :return: The ConstInterval or EmptySet instance.
    :rtype: ConstInterval | EmptySet
    """
    if const_cmp(a, b) == 1 or (const_cmp(a, b) == 0 and (not a_inc or not b_inc)):
        return EmptySet()
    elif const_cmp(a, b) == 0 and is_finite(a):
        return FiniteSet(a)
    else:
        return ConstInterval(a, b, a_inc, b_inc)


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
        # TODO: make this reduce multiple finite sets or intervals even if not all
        if len(self.sets) == 1:
            return self.sets[0]
        elif Universal in self.sets:
            return Universal
        elif any(isinstance(set_, Union) for set_ in self.sets):
            return Union(*(set_.sets if isinstance(set_, Union) else [set_] for set_ in self.sets)).reduce()
        elif all(isinstance(set_, FiniteSetType) for set_ in self.sets):
            return finite_set(*func_reduce(or_, ({*set_.elements} for set_ in self.sets), *{}))
        elif all(isinstance(set_, ConstInterval) for set_ in self.sets):
            # ! logic error of discontinuous regions that should be kept as union
            # min_bound = const_lbound_min(*((n.a, n.a_inc) for n in self.sets))
            # max_bound = const_rbound_max(*((n.b, n.b_inc) for n in self.sets))
            # return interval(min_bound[0], max_bound[0], min_bound[1], max_bound[1])
            return self
        elif all(isinstance(set_, NumberDomain) for set_ in self.sets):
            return NumberDomainSets[max(set_.domain_rank for set_ in self.sets)]
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
        # TODO: make this reduce multiple finite sets or intervals even if not all
        if len(self.sets) == 1:
            return self.sets[0]
        elif Empty in self.sets:
            return Empty
        elif all(isinstance(set_, FiniteSetType) for set_ in self.sets):
            return finite_set(*func_reduce(and_, ({*set_.elements} for set_ in self.sets), *{}))
        elif all(isinstance(set_, ConstInterval) for set_ in self.sets):
            min_bound = const_lbound_max(*((n.a, n.a_inc) for n in self.sets))
            max_bound = const_rbound_min(*((n.b, n.b_inc) for n in self.sets))
            return interval(min_bound[0], max_bound[0], min_bound[1], max_bound[1])
        elif all(isinstance(set_, NumberDomain) for set_ in self.sets):
            return NumberDomainSets[min(set_.domain_rank for set_ in self.sets)]
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
        # TODO: add this for intervals
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
        # TODO: add this for intervals
        if self.set_ == Empty:
            return Universal
        elif self.set_ == Universal:
            return Empty
        elif self.set_ == ConstInterval:
            if is_infinite(self.set_.a) and is_infinite(self.set_.b):
                return Empty
            elif is_infinite(self.set_.a):
                return interval(self.set_.b, inf, a_inc=not self.set_.b_inc)
            elif is_infinite(self.set_.b):
                return interval(-inf, self.set_.a, b_inc=not self.set_.a_inc)
            elif self.set_.a == self.set_.b:
                return Complement(FiniteSet(self.set_.a))
            else:
                return interval(-inf, self.set_.a, b_inc=not self.set_.a_inc) | interval(self.set_.b, inf, a_inc=not self.set_.b_inc)
        else:
            return self


@typechecked
@inherit_docstrings
class NumberDomain(Set):
    letter: str
    domain_rank: int

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
        return (KEY_SET_DOMAIN, self.domain_rank)


@typechecked
@inherit_docstrings
class NaturalSet(NumberDomain):
    letter: str = 'ℕ'
    domain_rank: int = DOMAIN_NATURAL

    def __init__(self, /):
        pass

    def __contains__(self, item: SetLike) -> bool:
        return isinstance(item, Integral) and item >= 0


@typechecked
@inherit_docstrings
class IntegerSet(NumberDomain):
    letter: str = 'ℤ'
    domain_rank: int = DOMAIN_INTEGER

    def __init__(self, /):
        pass

    def __contains__(self, item: SetLike) -> bool:
        return isinstance(item, Integral)


@typechecked
@inherit_docstrings
class RationalSet(NumberDomain):
    letter: str = 'ℚ'
    domain_rank: int = DOMAIN_RATIONAL

    def __init__(self, /):
        pass

    def __contains__(self, item: SetLike) -> bool:
        return isinstance(item, Rational)


@typechecked
@inherit_docstrings
class RealSet(NumberDomain):
    letter: str = 'ℝ'
    domain_rank: int = DOMAIN_REAL

    def __init__(self, /):
        pass

    def __contains__(self, item: SetLike) -> bool:
        return isinstance(item, Real | Constant) and is_finite(item)


@typechecked
@inherit_docstrings
class ComplexSet(NumberDomain):
    letter: str = 'ℂ'
    domain_rank: int = DOMAIN_COMPLEX

    def __init__(self, /):
        pass

    def __contains__(self, item: SetLike) -> bool:
        return isinstance(item, Complex | Constant) and is_finite(item)


Empty = EmptySet()
Universal = UniversalSet()
NumberDomainSets = (
    NSet := NaturalSet(),
    ZSet := IntegerSet(),
    QSet := RationalSet(),
    RSet := RealSet(),
    CSet := ComplexSet(),
)


def main():
    a = interval(1, 3)
    b = interval(2, 4)
    print(a | b)
    print(a & b)
    print(NSet | ZSet)
    print(QSet & CSet)


if __name__ == "__main__":
    main()
