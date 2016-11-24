r"""

Test-Suites
===========

Preparation
-----------

    sage: def alt(n):
    ....:         return (-1)^n

    sage: def I_even(n):
    ....:     return (1 + alt(n))/2

    sage: def I_odd(n):
    ....:     return (1 - alt(n))/2

    sage: @cached_function
    ....: def H(n):
    ....:     if n <= 0:
    ....:         return 0
    ....:     return H(n-1)+1/n

    sage: @cached_function
    ....: def H_odd(n):
    ....:     if n <= 0:
    ....:         return 0
    ....:     return H_odd(n-1) + bool(n%2 == 1)/n

    sage: @cached_function
    ....: def H_alt(n):
    ....:     if n <= 0:
    ....:         return 0
    ....:     return H_alt(n-1) + alt(n)/n

Classical Quicksort
-------------------

::

    sage: def comparisons_quicksort_classic(n):
    ....:     return 2 * (n+1) * H(n) - 4*n

    sage: all(QuickStar().avg_comparisons_quicksort(
    ....:         n, 'partitioned_classic', verbose=True) ==
    ....:     comparisons_quicksort_classic(n)
    ....:     for n in srange(9))
    n=0, cmp=0, avg=0
    n=1, cmp=0, avg=0
    n=2, cmp=2, avg=1
    n=3, cmp=16, avg=8/3
    n=4, cmp=116, avg=29/6
    n=5, cmp=888, avg=37/5
    n=6, cmp=7416, avg=103/10
    n=7, cmp=67968, avg=472/35
    n=8, cmp=682272, avg=2369/140
    True

Dual-pivot Quicksort "Count"
----------------------------

    sage: def comparisons_quicksort_dual(P):
    ....:     @cached_function
    ....:     def C(n):
    ....:         if n == 0:
    ....:             return 0
    ....:         return P(n) + (3 / binomial(n, 2) *
    ....:                        sum((n-1-k) * C(k) for k in srange(1, n-2+1))
    ....:                        if n >= 2 else 0)
    ....:     return C

    sage: def X_NE(n):
    ....:     if n <= 1:
    ....:         return 0
    ....:     return H_odd(n-2)/2 - 1/8 + alt(n) / 8 / (n - I_even(n))
    sage: def P_CT(n):
    ....:     if n <= 1:
    ....:         return 0
    ....:     return 3*n/2 - 9/4 + 1/4/(n - I_even(n)) + X_NE(n)
    sage: comparisons_quicksort_dual_count = comparisons_quicksort_dual(P_CT)

    sage: all(QuickStar().avg_comparisons_quicksort(
    ....:         n, 'partitioned_dual_count', verbose=True) ==
    ....:     comparisons_quicksort_dual_count(n)
    ....:     for n in srange(9))
    n=0, cmp=0, avg=0
    n=1, cmp=0, avg=0
    n=2, cmp=2, avg=1
    n=3, cmp=16, avg=8/3
    n=4, cmp=114, avg=19/4
    n=5, cmp=866, avg=433/60
    n=6, cmp=7188, avg=599/60
    n=7, cmp=65580, avg=1093/84
    n=8, cmp=655872, avg=244/15
    True

"""
class QuickStar(object):
    r"""

    EXAMPLES:

    Classical quicksort::

        sage: qs = QuickStar()
        sage: qs.quicksorted([2,1,3,5]), qs.comparisons
        ([1, 2, 3, 5], 4)

        sage: qs = QuickStar()
        sage: n = 4
        sage: E = srange(n)
        sage: pis = Permutations(E)
        sage: assert all(qs.quicksorted(pi) == E for pi in pis)
        sage: qs.comparisons
        116

    Dual-pivot quicksort "count"::

        sage: qs = QuickStar()
        sage: qs.quicksorted([2,1,3,5], 'partitioned_dual_count'), qs.comparisons
        ([1, 2, 3, 5], 5)

        sage: qs = QuickStar()
        sage: n = 4
        sage: E = srange(n)
        sage: pis = Permutations(E)
        sage: assert all(qs.quicksorted(pi, 'partitioned_dual_count') == E for pi in pis)
        sage: qs.comparisons
        114

    Classical quickselect::

        sage: qs = QuickStar()
        sage: qs.quickselected([2,1,3,0], 0), qs.comparisons
        (0, 4)

    Dual-pivot quickselect "count"::

        sage: qs = QuickStar()
        sage: n = 4
        sage: E = srange(n)
        sage: pis = Permutations(E)
        sage: assert all(qs.quickselected(pi, j, 'partitioned_dual_count') == j for pi in pis for j in srange(n))
        sage: qs.comparisons
        432
    """

    comparisons = 0
    

    def partitioned_classic(self, L):
        small = list()
        pivot = list()
        large = list()
        
        iterL = iter(L)
        try:
            p = next(iterL)
        except StopIteration:
            return small, list(), large
        else:
            pivot.append(p)
            
        for element in iterL:
            if element == p:
                pivot.append(element)
            else:
                self.comparisons += 1
                if element < p:
                    small.append(element)
                else:
                    large.append(element)
        return small, pivot, large

    
    def partitioned_dual_count(self, L):
        small = list()
        pivot = list()
        medium = list()
        qivot = list()
        large = list()

        iterL = iter(L)
        try:
            p = next(iterL)
        except StopIteration:
            return small, list(), medium, list(), large
        else:
            pivot.append(p)
        try:
            q = next(iterL)
        except StopIteration:
            return small, pivot, medium, list(), large
        else:
            qivot.append(q)

        if q < p:
            p, q = (q, p)
            pivot, qivot = (qivot, pivot)
        self.comparisons += 1

        s = 0
        ell = 0
        for element in iterL:
            if element == p:
                pivot.append(element)
            elif element == q:
                qivot.append(element)
            else:
                if s >= ell:
                    self.comparisons += 1
                    if element < p:
                        small.append(element)
                        s += 1
                    else:
                        self.comparisons += 1
                        if element > q:
                            large.append(element)
                            ell += 1
                        else:
                            medium.append(element)
                else:
                    self.comparisons += 1
                    if element > q:
                        large.append(element)
                        ell += 1
                    else:
                        self.comparisons += 1
                        if element < p:
                            small.append(element)
                            s += 1
                        else:
                            medium.append(element)
        return small, pivot, medium, qivot, large

        
    def partitioned_dual_clairvoyant(self, L):
        small = list()
        pivot = list()
        medium = list()
        qivot = list()
        large = list()

        iterL = iter(L)
        try:
            p = next(iterL)
        except StopIteration:
            return small, list(), medium, list(), large
        else:
            pivot.append(p)
        try:
            q = next(iterL)
        except StopIteration:
            return small, pivot, medium, list(), large
        else:
            qivot.append(q)

        if q < p:
            p, q = (q, p)
            pivot, qivot = (qivot, pivot)
        self.comparisons += 1

        s = 0
        ell = 0
        for element in L:
            if element < p:
                s += 1
            if element > q:
                ell += 1

        for element in iterL:
            if element == p:
                pivot.append(element)
            elif element == q:
                qivot.append(element)
            else:
                if s >= ell:
                    self.comparisons += 1
                    if element < p:
                        small.append(element)
                        s -= 1
                    else:
                        self.comparisons += 1
                        if element > q:
                            large.append(element)
                            ell -= 1
                        else:
                            medium.append(element)
                else:
                    self.comparisons += 1
                    if element > q:
                        large.append(element)
                        ell -= 1
                    else:
                        self.comparisons += 1
                        if element < p:
                            small.append(element)
                            s -= 1
                        else:
                            medium.append(element)
        assert s == 0 and ell == 0
        return small, pivot, medium, qivot, large

        
    def partitioned_dual_p_first(self, L):
        small = list()
        pivot = list()
        medium = list()
        qivot = list()
        large = list()

        iterL = iter(L)
        try:
            p = next(iterL)
        except StopIteration:
            return small, list(), medium, list(), large
        else:
            pivot.append(p)
        try:
            q = next(iterL)
        except StopIteration:
            return small, pivot, medium, list(), large
        else:
            qivot.append(q)

        if q < p:
            p, q = (q, p)
            pivot, qivot = (qivot, pivot)
        self.comparisons += 1

        for element in iterL:
            if element == p:
                pivot.append(element)
            elif element == q:
                qivot.append(element)
            else:
                self.comparisons += 1
                if element < p:
                    small.append(element)
                else:
                    self.comparisons += 1
                    if element > q:
                        large.append(element)
                    else:
                        medium.append(element)
        return small, pivot, medium, qivot, large

    
    def quicksorted(self, L, partitioning_strategy='partitioned_classic'):
        partitioned = getattr(self, partitioning_strategy)
        if len(L) <= 1:
            return L
        return sum((self.quicksorted(partition, partitioning_strategy=partitioning_strategy)
                    if is_even(i) else partition
                    for i, partition in enumerate(partitioned(L))),
                   list())


    def quickselected(self, L, j, partitioning_strategy='partitioned_classic'):
        partitioned = getattr(self, partitioning_strategy)
        if len(L) == 0:
            raise ValueError('empty list')
        if len(L) == 1:
            return L[0]
        partitions = partitioned(L)
        m = 0
        for partition in partitions:
            m += len(partition)
            if j < m:
                return self.quickselected(partition, j - m + len(partition),
                                          partitioning_strategy=partitioning_strategy)
        raise ValueError('index out of range')


    def avg_comparisons_quicksort(self, n, partitioning_strategy, verbose=False):
        r"""
        """
        E = srange(n)
        pis = Permutations(E)
        assert all(self.quicksorted(list(pi), partitioning_strategy) == E
                   for pi in pis)
        if verbose:
            print('n={}, cmp={}, avg={}'.format(
                n, self.comparisons, self.comparisons / n.factorial()))
        return self.comparisons / n.factorial()

