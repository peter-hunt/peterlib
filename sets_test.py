from set_theory import *


def basics_test():
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


def main():
    basics_test()


if __name__ == "__main__":
    main()
