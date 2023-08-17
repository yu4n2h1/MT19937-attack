from sage.all import PolynomialRing, Zmod, Mod
from tqdm import tqdm


class Mt:
    """保证第0位对应数值的最低位，第32位对应数值的"""

    def __init__(self, num) -> None:
        self.length = 32
        if isinstance(num, int):
            self.mt = [Mod(1, 2) if (1 << i) & num else Mod(0, 2) for i in range(32)]

        elif isinstance(num, (list, tuple)):
            self.mt = [num[x] for x in range(32)]

        elif isinstance(num, Mt):
            self.mt = num.mt

    def __and__(self, other):
        res = [0 for i in range(624)]
        for i in range(self.length):
            res[i] = self.mt[i] * other.mt[i]
        return Mt(res)

    def __xor__(self, other):
        res = [0 for i in range(624)]
        for i in range(self.length):
            res[i] = self.mt[i] + other.mt[i]
        return Mt(res)

    def __rmul__(self, other):
        res = [0 for i in range(624)]
        for i in range(self.length):
            res[i] = other * self.mt[i]
        return Mt(res)

    def __lshift__(self, other):
        """为向高位移动"""
        res = [0 for i in range(624)]
        for i in range(self.length - 1, other - 1, -1):
            res[i] = self.mt[i - other]
        return Mt(res)

    def __rshift__(self, other):
        """向低位移动"""
        res = [0 for i in range(624)]
        for i in range(self.length - other):
            res[i] = self.mt[i + other]
        return Mt(res)

    def to_num(self):
        if isinstance(
            self.mt[0], (int, sage.rings.finite_rings.integer_mod.IntegerMod_int)
        ):
            res = 0
            for i in range(32):
                # if self.mt[i]:
                res += self.mt[i] << i
            return res
        return None

    def __repr__(self) -> str:
        return "Mt" + str(self.mt)

    def __str__(self) -> str:
        return "Mt" + str(self.mt)


class VactorMT19937:
    def __init__(self, state):
        self.mt = [Mt(x) for x in state]
        self.mti = 0

    def extract_number(self):
        if self.mti == 0:
            self.twist()

        y = self.mt[self.mti]
        # print(y >> 11)
        y = y ^ y >> 11
        y = y ^ y << 7 & Mt(2636928640)
        y = y ^ y << 15 & Mt(4022730752)
        y = y ^ y >> 18
        self.mti = (self.mti + 1) % 624
        return y

    def twist(self):
        for i in range(0, 624):
            y = (self.mt[i] & Mt(0x80000000)) ^ (
                self.mt[(i + 1) % 624] & Mt(0x7FFFFFFF)
            )
            self.mt[i] = (y >> 1) ^ self.mt[(i + 397) % 624]

            """
            等价于
            if y % 2 != 0:
                self.mt[i] = self.mt[i] ^ Mt(0x9908b0df)
            """
            self.mt[i] = self.mt[i] ^ (y.mt[0] * Mt(0x9908B0DF))


if __name__ == "__main__":
    x = PolynomialRing(Zmod(2), names=["x%05d" % i for i in range(624 * 32)]).gens()
    x = x[::-1]
    state = [list(x[32 * i : 32 * (i + 1)]) for i in range(624)]
    vr = VactorMT19937(state)
    p = open("output.txt", "a+")
    # print(x[-1])

    for i in tqdm(range(624)):
        tmp = vr.extract_number()
        for j in range(32):
            que_x = [str(x[1])[1:] for x in list(tmp.mt[j])]
            # vec = ['0' if x_args not in que_x else '1' for x_args in x]

            p.write(",".join(que_x) + "\n")
    p.close()
