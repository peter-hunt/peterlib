from numbers import Number, Real

from typeguard import typechecked

from expr import *


__all__ = [
    "get_value",
    "const_cmp", "const_min", "const_max",
    "const_lbound_cmp", "const_lbound_min", "const_lbound_max",
    "const_rbound_cmp", "const_rbound_min", "const_rbound_max",
]


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
    return expr.value if isinstance(expr, Constant) else expr


@typechecked
@number_casted
def const_cmp(a: Real | Constant, b: Real | Constant, /) -> int:
    """
    Compare two constant values and return the comparison state.
    Infinities of the same sign are considered equal.

    :param a: The first value.
    :type a: Real | Constant
    :param b: The second value.
    :type b: Real | Constant
    :return: A comparison state: -1 if a < b, 0 if a = b, and 1 if a > b.
    :rtype: int
    """
    if a == undefined or b == undefined:
        raise ValueError("cannot compare undefined values")
    av = get_value(a)
    bv = get_value(b)
    if isinstance(av, complex):
        if isinstance(bv, complex):
            raise ValueError(f"cannot compare complex values: a={av}, b={bv}")
        raise ValueError(f"cannot compare complex values: a={av}")
    if isinstance(bv, complex):
        raise ValueError(f"cannot compare complex values: b={bv}")
    return -1 if av < bv else 0 if av == bv else 1


@typechecked
@number_casted
def const_min(*values: Real | Constant) -> Real | Constant:
    """
    Find the minimum value of the constant values with const_cmp.

    :param values: The values to find the minimum from.
    :type values: tuple[Real | Constant]
    :return: The minimum value from the ones given.
    :rtype: Real | Constant
    """
    if len(values) == 0:
        raise TypeError("const_min expected at least 1 argument, got 0")
    result = values[0]
    for value in values[1:]:
        if const_cmp(value, result) == -1:
            result = value
    return result


@typechecked
@number_casted
def const_max(*values: Real | Constant) -> Real | Constant:
    """
    Find the maximum value of the constant values with const_cmp.

    :param values: The values to find the maximum from.
    :type values: tuple[Real | Constant]
    :return: The maximum value from the ones given.
    :rtype: Real | Constant
    """
    if len(values) == 0:
        raise TypeError("const_max expected at least 1 argument, got 0")
    result = values[0]
    for value in values[1:]:
        if const_cmp(value, result) == 1:
            result = value
    return result


@typechecked
@number_casted
def const_lbound_cmp(a: Real | Constant, a_inc: bool,
                     b: Real | Constant, b_inc: bool, /) -> int:
    """
    Compare two constant values and inclusivity as left bound
    and return the comparison state.
    Infinities of the same sign are considered equal.

    :param a: The first left bound value.
    :type a: Real | Constant
    :param a_inc: Whether the first left bound is inclusive.
    :type a_inc: bool
    :param b: The second left bound value.
    :type b: Real | Constant
    :param b_inc: Whether the second left bound is inclusive.
    :type b_inc: bool
    :return: A comparison state: -1 if a < b, 0 if a = b, and 1 if a > b.
    :rtype: int
    """
    cmp_result = const_cmp(a, b)
    if cmp_result != 0:
        return cmp_result
    else:
        return -1 if a_inc and not b_inc else 1 if b_inc and not a_inc else 0


@typechecked
def const_lbound_min(*pairs: tuple[Real | Constant, bool]) -> tuple[Real | Constant, bool]:
    """
    Find the minimum value of the constant values with const_lbound_cmp.

    :param values: The values to find the minimum from.
    :type values: tuple[Real | Constant]
    :return: The minimum value from the ones given.
    :rtype: Real | Constant
    """
    if len(pairs) == 0:
        raise TypeError("const_lbound_min expected at least 1 argument, got 0")
    result = pairs[0]
    for pair in pairs[1:]:
        if const_lbound_cmp(*pair, *result) == -1:
            result = pair
    return result


@typechecked
def const_lbound_max(*pairs: tuple[Real | Constant, bool]) -> tuple[Real | Constant, bool]:
    """
    Find the maximum value of the constant values with const_lbound_cmp.

    :param values: The values to find the maximum from.
    :type values: tuple[Real | Constant]
    :return: The maximum value from the ones given.
    :rtype: Real | Constant
    """
    if len(pairs) == 0:
        raise TypeError("const_lbound_max expected at least 1 argument, got 0")
    result = pairs[0]
    for pair in pairs[1:]:
        if const_lbound_cmp(*pair, *result) == 1:
            result = pair
    return result


@typechecked
@number_casted
def const_rbound_cmp(a: Real | Constant, a_inc: bool,
                     b: Real | Constant, b_inc: bool, /) -> int:
    """
    Compare two constant values and inclusivity as right bound
    and return the comparison state.
    Infinities of the same sign are considered equal.

    :param a: The first right bound value.
    :type a: Real | Constant
    :param a_inc: Whether the first right bound is inclusive.
    :type a_inc: bool
    :param b: The second right bound value.
    :type b: Real | Constant
    :param b_inc: Whether the second right bound is inclusive.
    :type b_inc: bool
    :return: A comparison state: -1 if a < b, 0 if a = b, and 1 if a > b.
    :rtype: int
    """
    cmp_result = const_cmp(a, b)
    if cmp_result != 0:
        return cmp_result
    else:
        return -1 if b_inc and not a_inc else 1 if a_inc and not b_inc else 0


@typechecked
def const_rbound_min(*pairs: tuple[Real | Constant, bool]) -> tuple[Real | Constant, bool]:
    """
    Find the minimum value of the constant values with const_rbound_cmp.

    :param values: The values to find the minimum from.
    :type values: tuple[Real | Constant]
    :return: The minimum value from the ones given.
    :rtype: Real | Constant
    """
    if len(pairs) == 0:
        raise TypeError("const_rbound_min expected at least 1 argument, got 0")
    result = pairs[0]
    for pair in pairs[1:]:
        if const_rbound_cmp(*pair, *result) == -1:
            result = pair
    return result


@typechecked
def const_rbound_max(*pairs: tuple[Real | Constant, bool]) -> tuple[Real | Constant, bool]:
    """
    Find the maximum value of the constant values with const_rbound_cmp.

    :param values: The values to find the maximum from.
    :type values: tuple[Real | Constant]
    :return: The maximum value from the ones given.
    :rtype: Real | Constant
    """
    if len(pairs) == 0:
        raise TypeError("const_rbound_max expected at least 1 argument, got 0")
    result = pairs[0]
    for pair in pairs[1:]:
        if const_rbound_cmp(*pair, *result) == 1:
            result = pair
    return result
