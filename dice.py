from sympy import symbols, factor, Poly
from scipy.signal import fftconvolve
import itertools
import sys
import numpy as np



def sicherman(sides):
    x = symbols("x")
    poly = sum([x**i for i in range(sides)])
    factors = list(factor(poly).args)*2
    factors = [Poly(factor,x).all_coeffs() for factor in factors]
    fs = [sum(factor) for factor in factors]
    lf = len(factors)
    indices = list(range(lf))
    print(lf)
    fcom = factorCombinations(sides)
    fcom = [com for com in fcom if all([i in fs for i in com])]
    maxiter = min(
        max([len(i) for i in fcom])
        +
        sum([1 for i in fs if i == 1])
        ,
        lf-min([len(i) for i in fcom]))
    results = []
    for r in range(len(fs)//2, maxiter+1):
        for a in itertools.combinations(indices, r):
            product = np.prod([fs[i] for i in a])
            if product == sides:
                b = [i for i in range(lf) if i not in a]
                ac = factors[a[0]]
                bc = factors[b[0]]
                if len(a) > 1:
                    for i in a[1:]:
                        ac = fftconvolve(ac,factors[i])
                    if ac.dtype == float:
                        ac = np.rint(ac)
                if len(b) > 1:
                    for i in b[1:]:
                        bc = fftconvolve(bc,factors[i])
                    if bc.dtype == float:
                        bc = np.rint(bc)
                if min(ac) >= 0 and min(bc) >= 0 and sum(ac) == sides and sum(bc) == sides:
                    found1 = True
                    found2 = True
                    a4 = []
                    b4 = []
                    for i,c in enumerate(ac):
                        if c:
                            a4.extend([i+1]*int(c))
                    for i,c in enumerate(bc):
                        if c:
                            b4.extend([i+1]*int(c))
                    result = (a4,b4) if a4 > b4 else (b4,a4)
                    if result not in results:
                        results.append(result)
        if found1 and not found2:
            break
    if len(results) > 1:
        results.sort()
        contents = "\n".join([str(line)[1:-1] for line in results])

        with open(f"results/sicherman-d{sides:03}.txt","w") as f:
            f.write(contents)
            f.close()
if __name__ == "__main__":
    try:
        sides = int(sys.argv[1])
    except IndexError:
        sides = input("Enter the number of sides: ") or 8
        sides = int(sides)
    sicherman(sides)