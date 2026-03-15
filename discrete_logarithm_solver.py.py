import math
import random
import itertools
from primefac import primefac
from collections import Counter
import sympy as sp

# Define the parameters
p = 18443
g = 37
h = 211
B = 5

def is_Bsmooth(b, n):
    """Determine if n is B-smooth (uses fast factoring)."""
    P = list(primefac(n))
    if len(P) != 0 and P[-1] <= b: 
        return True, P
    return False, P

def factorlist_to_explist(L):
    """Convert a list of prime factors to a dictionary of bases and exponents."""
    return dict(Counter(int(n) for n in L))

def find_congruences():
    """Generate random smooth relations."""
    congruences = []
    bases_set = set()
    
    print('Searching for congruences...')
    while True:                                                 
        k = random.randint(2, p - 1)                                                
        smooth, factors = is_Bsmooth(B, pow(g, k, p))                                     
        if smooth:                                                                         
            congruences.append((factorlist_to_explist(factors), k))                    
            
            # Update known bases with any unique prime factors found
            for base in congruences[-1][0].keys():
                bases_set.add(base)
            
            # We need at least as many equations as bases to solve the system.
            # Generating a few extra increases the chance of finding an invertible matrix.
            if len(bases_set) > 0 and len(congruences) >= len(bases_set) + 3: 
                break
                
    bases = list(bases_set)
    print(f'Congruences found: {len(congruences)}\nBases: {len(bases)}')
    return bases, congruences

def solve_system(bases, congruences):
    """Convert the linear system to matrices and solve modulo p-1 using SymPy."""
    n = len(bases)
    
    # Iterate through subsets of congruences to find a square, invertible submatrix
    for subset in itertools.combinations(congruences, n):
        M_list = [[c[0].get(base, 0) for base in bases] for c in subset]
        b_list = [[c[1]] for c in subset]
        
        M_sp = sp.Matrix(M_list)
        
        # A matrix is only invertible modulo (p-1) if gcd(det(M), p-1) == 1
        if math.gcd(int(M_sp.det()), p - 1) == 1:
            # Invert matrix modulo p-1
            M_inv = M_sp.inv_mod(p - 1)
            b_sp = sp.Matrix(b_list)
            
            # Solve for the discrete log exponents: x = (M^-1 * b) mod (p-1)
            x = (M_inv * b_sp) % (p - 1)
            return [int(val) for val in x]
            
    print("Could not find an invertible submatrix modulo p-1. Try running again.")
    exit()

def evaluate(eq_dict, dlogs):
    """Evaluate a linear equation using known discrete logs."""
    # Renamed 'exp' to 'exponent' to avoid shadowing math.exp (if imported later)
    return sum([dlogs[term] * exponent for term, exponent in eq_dict.items()]) % (p - 1)

def check_congruences(congruences, dlogs):
    """Verify that the calculated discrete logs satisfy the original congruences."""
    print('Checking congruences:', end=" ")
    for c in congruences:
        if evaluate(c[0], dlogs) != c[1]:
            print('Failed! Congruence mismatch.')
            return False
    print('Passed!\n')
    return True

def check_dlogs(exponents, bases):
    """Verify that g^x = base (mod p)."""
    print('Checking dlog exponents:')
    for exponent, base in zip(exponents, bases):
        if pow(g, exponent, p) != base:
            print(f'Failed for base {base}.')
            return False
        else:
            print(f'{g}^{exponent} = {base} (mod {p})')
    print('Passed!\n')
    return True

def main():
    print(f'p: {p}, g: {g}, h: {h}, B: {B}')
    
    # 1. Generate and extract bases/congruences
    bases, congruences = find_congruences()
    
    # 2. Solve the linear system
    print('Solving linear system with SymPy...')
    exponents = solve_system(bases, congruences)
    print('SymPy done.\n')

    # Map bases to their solved discrete logarithms
    dlogs = {b: exp for b, exp in zip(bases, exponents)}

    # 3. Verify the intermediate results
    check_congruences(congruences, dlogs)
    check_dlogs(exponents, bases)

    # 4. Find k such that h * g^-k is B-smooth
    print('Searching for k such that h * g^-k is B-smooth...')
    for _ in range(10**6):
        k = random.randint(2, p - 1)
        # Using pow(g, -k, p) is standard in Python 3.8+ to find inverse exponents
        val = (h * pow(g, -k, p)) % p
        smooth, factors = is_Bsmooth(B, val)
        if smooth:
            print(f'Found k = {k}')
            break
    else:
        print("Failed to find a smooth value after 1,000,000 iterations.")
        return

    # 5. Solve the main DLP
    print('\nSolving the main dlog problem:')
    factor_dict = factorlist_to_explist(factors)
    soln = (evaluate(factor_dict, dlogs) + k) % (p - 1)
    
    if pow(g, soln, p) == h:
        print(f'{g}^{soln} = {h} (mod {p}) holds!')
        print(f'DLP solution: {soln}')
    else:
        print('Failed to verify final solution.')

if __name__ == '__main__':
    main()
