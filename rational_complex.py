from fractions import Fraction
from math import sqrt
from numbers import Rational


__all__ = ["RationalComplex"]


def rdiv(a: Rational, b: Rational) -> Rational:
    if isinstance(a, int) and isinstance(b, int):
        return Fraction(a, b)
    else:
        return a / b


class RationalComplex:
    real: Rational
    imag: Rational

    def __init__(self, real: Rational = 0, imag: Rational = 0, /):
        if not isinstance(real, Rational):
            raise TypeError("RationalComplex() argument 'real'"
                            " must be a rational number")
        if not isinstance(imag, Rational):
            raise TypeError("RationalComplex() argument 'imag'"
                            " must be a rational number")
        self.real = Fraction(real)
        self.imag = Fraction(imag)

    def __repr__(self):
        if self.real == 0:
            if self.imag % 1 == 0:
                return f"{self.imag}j"
            else:
                return f"({self.imag})j"
        result = f"{self.real}"
        result += '+' if self.imag >= 0 else '-'
        if self.imag % 1 == 0:
            result += f"{abs(self.imag)}j"
        else:
            result += f"({abs(self.imag)})j"
        return result

    def __add__(self, other: Rational | RationalComplex):
        if isinstance(other, Rational):
            return self + RationalComplex(other)
        elif isinstance(other, RationalComplex):
            return RationalComplex(self.real + other.real, self.imag + other.imag)
        else:
            return NotImplemented

    def __radd__(self, other: Rational):
        if isinstance(other, Rational):
            return RationalComplex(other) + self
        else:
            return NotImplemented

    def __sub__(self, other: Rational | RationalComplex):
        if isinstance(other, Rational):
            return self - RationalComplex(other)
        elif isinstance(other, RationalComplex):
            return RationalComplex(self.real - other.real, self.imag - other.imag)
        else:
            return NotImplemented

    def __rsub__(self, other: Rational):
        if isinstance(other, Rational):
            return RationalComplex(other) - self
        else:
            return NotImplemented

    def __mul__(self, other: Rational | RationalComplex):
        if isinstance(other, Rational):
            return self * RationalComplex(other)
        elif isinstance(other, RationalComplex):
            return RationalComplex(self.real * other.real - self.imag * other.imag,
                                   self.real * other.imag + self.imag * other.real)
        else:
            return NotImplemented

    def __rmul__(self, other: Rational):
        if isinstance(other, Rational):
            return RationalComplex(other) * self
        else:
            return NotImplemented

    def __truediv__(self, other: Rational | RationalComplex):
        if isinstance(other, Rational):
            return self / RationalComplex(other)
        elif isinstance(other, RationalComplex):
            return RationalComplex(rdiv((self.real * other.real + self.imag * other.imag),
                                        (other.real ** 2 + other.imag ** 2)),
                                   rdiv((self.imag * other.real - self.real * other.imag),
                                        (other.real ** 2 + other.imag ** 2)))
        else:
            return NotImplemented

    def __rtruediv__(self, other: Rational):
        if isinstance(other, Rational):
            return RationalComplex(other) / self
        else:
            return NotImplemented

    def __pow__(self, other: int):
        if isinstance(other, int):
            if other < 0:
                raise ValueError(
                    "cannot raise RationalComplex to a negative power"
                    " without converting to float complex first")
            elif other == 0:
                if self.real == 0 and self.imag == 0:
                    return RationalComplex(1)
                else:
                    return RationalComplex()
            else:
                result = RationalComplex(1, 0)
                for _ in range(other):
                    result *= self
                return result
        else:
            return NotImplemented

    def __neg__(self):
        return RationalComplex(-self.real, -self.imag)

    def __pos__(self):
        return self

    def __abs__(self):
        return sqrt(self.real ** 2 + self.imag ** 2)

    def __complex__(self):
        return complex(self.real, self.imag)


def main():
    a = RationalComplex(3, 4)
    b = RationalComplex(1, 2)
    print(f"{a=}")
    print(f"{b=}")
    print(f"{+a=}")
    print(f"{-b=}")
    print(f"{abs(a)=}")
    print(f"{a + b=}")
    print(f"{a - b=}")
    print(f"{a * b=}")
    print(f"{a / b=}")
    print(f"{a / b * b=}")
    print(f"{b ** 2=}")


if __name__ == "__main__":
    main()
