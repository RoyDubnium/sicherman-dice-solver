from sympy import symbols, factor, Poly
from scipy.signal import fftconvolve
import itertools
import sys
import numpy as np



def sicherman(sides):
    x = symbols("x")
    poly = sum([x**i for i in range(1,sides+1)])
    factors = list(factor(poly).args)*2
    factors = [Poly(factor,x).all_coeffs() for factor in factors]

    lf = len(factors)
    indices = list(range(lf))
    print(lf)
    combosidx = [itertools.combinations(indices, i) for i in range(len(factors)//2,len(factors))]
    
    results = []
    found1 = False
    for combs in combosidx:
        combs = list(combs)
        altidx = [tuple([i for i in range(lf) if i not in elem]) for elem in combs]
        found2 = False
        for a,b in zip(combs,altidx):
            a1 = [factors[i] for i in a]
            b1 = [factors[i] for i in b]
            ac = [1]
            bc = [1]
            for i in a:
                ac = fftconvolve(ac,factors[i])
            for i in b:
                bc = fftconvolve(bc,factors[i])
            if ac.dtype == float:
                ac = np.rint(ac)
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
    results.sort()
    contents = "\n".join([str(line)[1:-1] for line in results])

    with open(f"results/sicherman-d{sides:03}-test.txt","w") as f:
        f.write(contents)
        f.close()
if __name__ == "__main__":
    try:
        sides = int(sys.argv[1])
    except IndexError:
        sides = input("Enter the number of sides: ") or 8
        sides = int(sides)
    sicherman(sides)