@cached_function
def H(n):
    if n <= 0:
        return 0
    return H(n-1)+1/n

def comparisons_quicksort_classic(n):
    return 2 * (n+1) * H(n) - 4*n

data = [(n, factorial(n) * comparisons_quicksort_classic(n))
        for n in srange(100+1)]

#print(', '.join(str(c) for n, c in data))

for n, c in data:
    print('{} {}'.format(n, c))
