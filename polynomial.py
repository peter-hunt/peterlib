from fractions import Fraction
from itertools import zip_longest
from math import inf
from numbers import Number

from .utils import cmp


__all__ = ["Polynomial", "X"]


DIVIDE_TO_FRACTION = True
RENDER_AS_FRACTION = False


if RENDER_AS_FRACTION:
    # may render 1.6666666666666667 as 7505999378950827/4503599627370496
    # if fractions are not used through calculations
    def format_float(number: Number, /) -> str:
        return f"{number:.0f}" if number % 1 == 0 else f"{Fraction(number)}"
else:
    def format_float(number: Number, /) -> str:
        return f"{number:.0f}" if number % 1 == 0 else f"{number}"
if DIVIDE_TO_FRACTION:
    def div(a: Number, b: Number) -> Number:
        a = int(a) if a % 1 == 0 else a
        b = int(b) if b % 1 == 0 else b
        if isinstance(a, int) and isinstance(b, int):
            return Fraction(a, b)
        else:
            return a / b
else:
    def div(a: Number, b: Number) -> Number:
        return a / b


class Polynomial:
    coefs: tuple[Number]

    def __init__(self, *coefs: Number):
        # constant and coeffcients of natural number power polynomial,
        # starting with constant and sorted in increasing power:
        # for example, x ** 2 - 3 * x + 2 is written as (2, -3, 1).
        # the trailing zeroes are always hidden.
        for coef in coefs:
            if not isinstance(coef, Number):
                raise TypeError("coefficients must be numbers")
        self.coefs = (0,) if len(coefs) == 0 else coefs
        while self.coefs[-1] == 0 and len(self.coefs) > 1:
            self.coefs = self.coefs[:-1]

    def __repr__(self):
        if len(self) == 0:
            return '0'
        elif len(self) == 1:
            return format_float(self[0])
        elif len(self) == 2:
            if abs(self[1]) == 1:
                result = 'x'
            else:
                result = format_float(self[1])
                if not self.is_int:
                    result += " * "
                result += 'x'
            if self[0] != 0:
                result += " + " if self[0] > 0 else " - "
                result += format_float(abs(self[0]))
            return result
        result = '-' if self[-1] < 0 else ''
        if abs(self[-1]) == 1:
            result += f"x^{len(self) - 1}"
        else:
            result += format_float(abs(self[-1]))
            if not self.is_int:
                result += " * "
            result += f"x^{len(self) - 1}"
        for expo, coef in reversed([*enumerate(self[:-1])]):
            if coef == 0:
                continue
            result += " + " if coef > 0 else " - "
            if expo == 0:
                result += format_float(abs(coef))
            elif expo == 1:
                if abs(coef) != 1:
                    result += format_float(abs(coef))
                    if not self.is_int:
                        result += " * "
                result += 'x'
            else:
                if abs(coef) != 1:
                    result += format_float(abs(coef))
                    if not self.is_int:
                        result += " * "
                result += f"x^{expo}"
        return result

    def __call__(self, value: Number, /):
        return sum(coef * value ** expo for expo, coef in enumerate(self))

    def __len__(self):
        return 0 if self.coefs == (0,) else len(self.coefs)

    def __getitem__(self, key: int) -> Number:
        return self.coefs[key]

    def __iter__(self):
        return iter(self.coefs)

    def __add__(self, other: Polynomial | Number) -> Polynomial:
        if isinstance(other, Polynomial):
            return Polynomial(*(scoef + ocoef for scoef, ocoef in zip_longest(self, other, fillvalue=0)))
        elif isinstance(other, Number):
            return Polynomial(*((self[0] + other,) + self[1:]))
        else:
            return NotImplemented

    def __radd__(self, other: Number) -> Polynomial:
        if isinstance(other, Number):
            return Polynomial(*((self[0] + other,) + self[1:]))
        else:
            return NotImplemented

    def __sub__(self, other: Polynomial | Number) -> Polynomial:
        if isinstance(other, Polynomial):
            return Polynomial(*(scoef - ocoef for scoef, ocoef in zip_longest(self, other, fillvalue=0)))
        elif isinstance(other, Number):
            return Polynomial(self[0] - other, *self[1:])
        else:
            return NotImplemented

    def __rsub__(self, other: Number) -> Polynomial:
        if isinstance(other, Number):
            return Polynomial(other - self[0], *(-coef for coef in self[1:]))
        else:
            return NotImplemented

    def __mul__(self, other: Polynomial | Number) -> Polynomial:
        if isinstance(other, Polynomial):
            return Polynomial(*(
                sum(
                    scoef * other[rexpo - sexpo]
                    for sexpo, scoef in [*enumerate(self)][max(rexpo - len(other) + 1, 0):rexpo + 1]
                )
                for rexpo in range(len(self) + len(other) - 1)
            ))
        elif isinstance(other, Number):
            return Polynomial(*(coef * other for coef in self))
        else:
            return NotImplemented

    def __rmul__(self, other: Number) -> Polynomial:
        if isinstance(other, Number):
            return Polynomial(*(coef * other for coef in self))
        else:
            return NotImplemented

    def __floordiv__(self, other: Polynomial | Number) -> Polynomial:
        if isinstance(other, Polynomial):
            return self.divmod(other)[0]
        elif isinstance(other, Number):
            return self.divmod(Polynomial(other))[0]
        else:
            return NotImplemented

    def __rfloordiv__(self, other: Number) -> Polynomial:
        if isinstance(other, Number):
            return Polynomial(other).divmod(self)[0]
        else:
            return NotImplemented

    def __mod__(self, other: Polynomial | Number) -> Polynomial:
        if isinstance(other, Polynomial):
            return self.divmod(other)[1]
        elif isinstance(other, Number):
            return self.divmod(Polynomial(other))[1]
        else:
            return NotImplemented

    def __rmod__(self, other: Number) -> Polynomial:
        if isinstance(other, Number):
            return Polynomial(other).divmod(self)[1]
        else:
            return NotImplemented

    def __divmod__(self, divisor: Polynomial | Number, /) -> tuple[Polynomial, Polynomial]:
        if isinstance(divisor, Number):
            divisor = Polynomial(divisor)
        elif not isinstance(divisor, Polynomial):
            return NotImplemented
        if len(divisor) == 0:
            raise ZeroDivisionError(f"division by zero: {self} /% {divisor}")
        dend_len = len(self)
        sor_len = len(divisor)
        rounds = max(dend_len - sor_len + 1, 0)
        remainder = [*self]
        quotient = [None for _ in range(rounds)]
        for scan_expo in range(rounds - 1, -1, -1):
            coef = div(remainder.pop(), divisor[-1])
            quotient[scan_expo] = coef
            for sor_expo in range(sor_len - 1):
                dend_expo = scan_expo + sor_expo
                remainder[dend_expo] -= coef * divisor[sor_expo]
        return (Polynomial(*quotient), Polynomial() if len(remainder) == 0 else Polynomial(*remainder))

    def __rdivmod__(self, divident: Number) -> tuple[Polynomial, Polynomial]:
        return divmod(Polynomial(divident), self)

    def __pow__(self, other: Polynomial | Number) -> Polynomial:
        if isinstance(other, Polynomial):
            if len(other) > 1:
                raise ValueError(
                    "cannot raise polynomial to power of polynomial")
            expo = other[0]
        elif isinstance(other, Number):
            expo = other
        else:
            return NotImplemented
        if expo < 0:
            raise ValueError("cannot raise polynomial to negative power")
        result = Polynomial(1)
        for _ in range(expo):
            result *= self
        return result

    def __neg__(self) -> Polynomial:
        return Polynomial(*(-coef for coef in self))

    def __pos__(self) -> Polynomial:
        return self

    @property
    def is_int(self) -> bool:
        return all(coef % 1 == 0 for coef in self)

    @property
    def derivative(self) -> Polynomial:
        return Polynomial(
            *(coef * expo
              for expo, coef in enumerate(self[1:], 1))
        )

    @property
    def integral(self) -> Polynomial:
        return Polynomial(
            0,
            *(div(coef, (expo + 1))
              for expo, coef in enumerate(self))
        )

    def get_sign_at(self, value: Number) -> int:
        if len(self) == 0:
            return 0
        elif value == -inf:
            return 1 if (self[-1] > 0) ^ (len(self) % 2 == 0) else -1
        elif value == inf:
            return 1 if self[-1] > 0 else -1
        else:
            return cmp(self(value), 0)

    def inc_mono_real_solve(self, lbound: Number, rbound: Number, /) -> Number | None:
        while True:
            mbound = (lbound + rbound) / 2
            if rbound == mbound or mbound == lbound:
                if self(lbound) == 0:
                    return lbound
                elif self(rbound) == 0:
                    return rbound
                elif self.get_sign_at(lbound) < 0 and self.get_sign_at(rbound) > 0:
                    return lbound
                else:
                    return
            msign = self.get_sign_at(mbound)
            if msign > 0:
                rbound = mbound
            elif msign < 0:
                lbound = mbound
            else:
                return mbound

    # solving within monotonic region guaranteed by either df/dx=0 bound or infinity bound
    def bound_real_solve(self, lbound: Number, rbound: Number, /) -> Number | None:
        lbound = float(lbound)
        rbound = float(rbound)
        lsign = self.get_sign_at(lbound)
        rsign = self.get_sign_at(rbound)
        if lsign == rsign and lsign in (1, -1) or lsign == 0 or rsign == 0:
            return
        if lbound == -inf and rbound == inf:
            lbound = -1
            try:
                while lbound != -inf and self.get_sign_at(lbound) != lsign:
                    lbound *= 2
            except OverflowError:
                return
            # ! return none if final sign swap takes farther than python number limit
            if lbound == -inf:
                return
        elif lbound == -inf:
            gap = 1
            try:
                while (rbound - gap) != -inf and self.get_sign_at(rbound - gap) != lsign:
                    gap *= 2
            except OverflowError:
                return
            lbound = rbound - gap
            if lbound == -inf:
                return
            if gap > 1:
                rbound = lbound + gap / 2
                if self.get_sign_at(rbound) == 0:
                    return rbound
        if rbound == inf:
            gap = 1
            try:
                while (lbound + gap) != -inf and self.get_sign_at(lbound + gap) != rsign:
                    gap *= 2
            except OverflowError:
                return
            rbound = lbound + gap
            if rbound == inf:
                return
            if gap > 1:
                lbound = rbound - gap / 2
                if self.get_sign_at(lbound) == 0:
                    return lbound
        if lsign < 0:
            return self.inc_mono_real_solve(lbound, rbound)
        else:
            return (-self).inc_mono_real_solve(lbound, rbound)

    @property
    def real_solutions(self) -> tuple[Number, ...]:
        if len(self) <= 1:
            return ()  # return none if polynomial is constant, either 0 or not
        elif len(self) == 2:
            return (-div(self[0], self[1]),)
        else:
            der_zeroes = self.derivative.real_solutions
            if len(der_zeroes) == 0:
                result = self.bound_real_solve(-inf, inf)
                return () if result is None else (result,)

            zeroes = []
            result = self.bound_real_solve(-inf, der_zeroes[0])
            if result is not None:
                zeroes.append(result)
            for i, lbound in enumerate(der_zeroes[:-1]):
                if self.get_sign_at(lbound) == 0:
                    zeroes.append(lbound)
                rbound = der_zeroes[i + 1]
                result = self.bound_real_solve(lbound, rbound)
                if result is not None:
                    zeroes.append(result)
            if self.get_sign_at(der_zeroes[-1]) == 0:
                zeroes.append(der_zeroes[-1])
            result = self.bound_real_solve(der_zeroes[-1], inf)
            if result is not None:
                zeroes.append(result)

            return tuple(zeroes)


X = Polynomial(0, 1)


def main():
    a = 5 * X ** 2 + 4 * X + 3
    b = 2 * X + 1
    print(f"{a=}")
    print(f"{b=}")
    print(f"{a + 2=}")
    print(f"{a + b=}")
    print(f"{a - 2=}")
    print(f"{3 - a=}")
    print(f"{a - b=}")
    print(f"{a * 5=}")
    print(f"{a * b=}")
    print(f"{a ** 2=}")
    print(f"{a(1)=}")
    print(f"{-a=}")
    print(f"{a.derivative=}")
    print(f"{a.integral=}")
    print(f"{divmod(a, b)=}")
    print(a - 6)
    print((a - 6).real_solutions)
    c = Polynomial(0, 1, 2, 3, 4, 5, 6, 7, 8)
    print(c)
    print(c.real_solutions)
    for x in c.real_solutions:
        print(divmod(c, (X - x)))


if __name__ == "__main__":
    main()
