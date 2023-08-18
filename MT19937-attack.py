from sage.all import *
from tqdm import tqdm
from tools import read_Matrix

E = read_Matrix("simulate/extract_Matix.txt")
T = read_Matrix("simulate/twist_Matirx.txt")

x = PolynomialRing(Zmod(2), names=["x%05d" % i for i in range(624 * 32)]).gens()
x = vector(x)
# print(x)g
print((E * T * x)[0])
