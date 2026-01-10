__all__ = [
    "cmp",
    "inherit_docstrings",
]


def cmp(a, b):
    return -1 if a < b else 0 if a == b else 1


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
