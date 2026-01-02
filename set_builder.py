from sets import Set
from statement import Statement


__all__ = [
    "SetBuilder",
]


class SetBuilder(Set):
    base_set: Set
    condition: Statement
