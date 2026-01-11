from polynomial import *


def basic_test():
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
    c = RealPoly(0, 1, 2, 3, 4, 5, 6, 7, 8)
    print(c)
    print(c.real_solutions)
    for x in c.real_solutions:
        print(divmod(c, (X - x)))


def main():
    basic_test()


if __name__ == "__main__":
    main()
