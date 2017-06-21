def alt(n):
        return (-1)^n

def I_even(n):
    return (1 + alt(n))/2

def I_odd(n):
    return (1 - alt(n))/2

@cached_function
def H(n):
    if n <= 0:
        return 0
    return H(n-1)+1/n

@cached_function
def H_odd(n):
    if n <= 0:
        return 0
    return H_odd(n-1) + bool(n%2 == 1)/n

@cached_function
def H_alt(n):
    if n <= 0:
        return 0
    return H_alt(n-1) + alt(n)/n

def X_NE(n):
    if n <= 1:
        return 0
    return H_odd(n-2)/2 - 1/8 + alt(n) / 8 / (n - I_even(n))
def P_CT(n):
    if n <= 1:
        return 0
    return 3*n/2 - 9/4 + 1/4/(n - I_even(n)) + X_NE(n)

def comparisons_quicksort_dual(P):
    @cached_function
    def C(n):
        if n == 0:
            return 0
        return P(n) + (3 / binomial(n, 2) *
                       sum((n-1-k) * C(k) for k in srange(1, n-2+1))
                       if n >= 2 else 0)
    return C
comparisons_quicksort_dual_count = comparisons_quicksort_dual(P_CT)

data = [(n, factorial(n) * comparisons_quicksort_dual_count(n))
        for n in srange(100+1)]

#print(', '.join(str(c) for n, c in data))

for n, c in data:
    print('{} {}'.format(n, c))
