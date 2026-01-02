from expr import ExprLike
from sets import Set


__all__ = [
    "Statement",
    "Equality",
    "Inequality", "GreaterThan", "GreaterOrEqual", "LessThan", "LessOrEqual",
    "Membership",
    "LogicNeg", "LogicOr", "LogicAnd",
]


class Statement:
    def __init__(self, *args, **kwargs):
        raise Exception(
            "Cannot create an instance of base Statement class. "
            "To use a subclass, redefine the init code."
        )


class Equality(Statement):
    lhs: ExprLike
    rhs: ExprLike


class Inequality(Statement):
    lhs: ExprLike
    rhs: ExprLike


class GreaterThan(Inequality):
    lhs: ExprLike
    rhs: ExprLike


class GreaterOrEqual(Inequality):
    lhs: ExprLike
    rhs: ExprLike


class LessThan(Inequality):
    lhs: ExprLike
    rhs: ExprLike


class LessOrEqual(Inequality):
    lhs: ExprLike
    rhs: ExprLike


class Membership(Statement):
    element: ExprLike
    set_: Set


class LogicNeg(Statement):
    operand: Statement


class LogicOr(Statement):
    operands: list[Statement]


class LogicAnd(Statement):
    operands: list[Statement]
