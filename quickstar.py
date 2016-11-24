r"""

Test-Suites
===========

Preparation
-----------

::

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

::

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

Dual-pivot Quicksort "Clairvoyant"
----------------------------------

::

    sage: def X_SE(n):
    ....:     if n <= 1:
    ....:         return 0
    ....:     return X_NE(n) - 1/2 + 1/2/(n-I_even(n))
    sage: def P_CV(n):
    ....:     if n <= 1:
    ....:         return 0
    ....:     return 3*n/2 - 9/4 + 1/4/(n - I_even(n)) - X_SE(n)
    sage: comparisons_quicksort_dual_clairvoyant = comparisons_quicksort_dual(P_CV) 
    sage: all(QuickStar().avg_comparisons_quicksort(
    ....:         n, 'partitioned_dual_clairvoyant', verbose=True) ==
    ....:     comparisons_quicksort_dual_clairvoyant(n)
    ....:     for n in srange(9))
    n=0, cmp=0, avg=0
    n=1, cmp=0, avg=0
    n=2, cmp=2, avg=1
    n=3, cmp=14, avg=7/3
    n=4, cmp=102, avg=17/4
    n=5, cmp=778, avg=389/60
    n=6, cmp=6492, avg=541/60
    n=7, cmp=59484, avg=4957/420
    n=8, cmp=597216, avg=6221/420
    True

Dual-pivot Quicksort "p first"
------------------------------

::

    sage: def P_pF(n):
    ....:     if n <= 1:
    ....:         return 0
    ....:     return 5/3*(n-2) + 1
    sage: comparisons_quicksort_dual_p_first = comparisons_quicksort_dual(P_pF) 
    sage: all(QuickStar().avg_comparisons_quicksort(
    ....:         n, 'partitioned_dual_p_first', verbose=True) ==
    ....:     comparisons_quicksort_dual_p_first(n)
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

Classical Quickselect
---------------------

::

    sage: def comparisons_quickselect_classic(n, j):
    ....:     return 2 * (n + 3 + (n+1)*H(n) - (j+2)*H(j) - (n+3-j)*H(n+1-j))
    sage: all(QuickStar().avg_comparisons_quickselect(
    ....:         n, j, 'partitioned_classic', verbose=True) ==
    ....:     comparisons_quickselect_classic(n, j+1)
    ....:     for n in srange(7) for j in srange(n))
    n=1, j=0, cmp=0, avg=0
    n=2, j=0, cmp=2, avg=1
    n=2, j=1, cmp=2, avg=1
    n=3, j=0, cmp=14, avg=7/3
    n=3, j=1, cmp=16, avg=8/3
    n=3, j=2, cmp=14, avg=7/3
    n=4, j=0, cmp=92, avg=23/6
    n=4, j=1, cmp=108, avg=9/2
    n=4, j=2, cmp=108, avg=9/2
    n=4, j=3, cmp=92, avg=23/6
    n=5, j=0, cmp=652, avg=163/30
    n=5, j=1, cmp=768, avg=32/5
    n=5, j=2, cmp=808, avg=101/15
    n=5, j=3, cmp=768, avg=32/5
    n=5, j=4, cmp=652, avg=163/30
    n=6, j=0, cmp=5112, avg=71/10
    n=6, j=1, cmp=6000, avg=25/3
    n=6, j=2, cmp=6456, avg=269/30
    n=6, j=3, cmp=6456, avg=269/30
    n=6, j=4, cmp=6000, avg=25/3
    n=6, j=5, cmp=5112, avg=71/10
    True

Dual-pivot Quickselect "Count"
------------------------------

::

    sage: def comparisons_quickselect_dual(P):
    ....:     @cached_function
    ....:     def C(n, j):
    ....:         if n <= 0 or j <= 0 or n < j:
    ....:             return 0
    ....:         return P(n) + (1 / binomial(n, 2) * \
    ....:             (sum((n-1-s) * C(s, j)
    ....:                  for s in srange(j, n-2+1)) +
    ....:              sum(C(n-2-s-ell, j-(s+1))
    ....:                  for s in srange(0, j-2+1)
    ....:                  for ell in srange(0, n-j-1+1)) +
    ....:              sum((n-1-ell) * C(ell, j-(n-ell))
    ....:                  for ell in srange(n-j+1, n-2+1))
    ....:             ) if n >= 2 else 0)
    ....:     return C
    sage: comparisons_quickselect_dual_count = comparisons_quickselect_dual(P_CT)

    sage: all(QuickStar().avg_comparisons_quickselect(
    ....:         n, j, 'partitioned_dual_count', verbose=True) ==
    ....:     comparisons_quickselect_dual_count(n, j+1)
    ....:     for n in srange(7) for j in srange(n))
    n=1, j=0, cmp=0, avg=0
    n=2, j=0, cmp=2, avg=1
    n=2, j=1, cmp=2, avg=1
    n=3, j=0, cmp=16, avg=8/3
    n=3, j=1, cmp=16, avg=8/3
    n=3, j=2, cmp=16, avg=8/3
    n=4, j=0, cmp=106, avg=53/12
    n=4, j=1, cmp=110, avg=55/12
    n=4, j=2, cmp=110, avg=55/12
    n=4, j=3, cmp=106, avg=53/12
    n=5, j=0, cmp=754, avg=377/60
    n=5, j=1, cmp=798, avg=133/20
    n=5, j=2, cmp=818, avg=409/60
    n=5, j=3, cmp=798, avg=133/20
    n=5, j=4, cmp=754, avg=377/60
    n=6, j=0, cmp=5916, avg=493/60
    n=6, j=1, cmp=6312, avg=263/30
    n=6, j=2, cmp=6564, avg=547/60
    n=6, j=3, cmp=6564, avg=547/60
    n=6, j=4, cmp=6312, avg=263/30
    n=6, j=5, cmp=5916, avg=493/60
    True

Dual-pivot Quickselect "Clairvoyant"
------------------------------------

::

    sage: comparisons_quickselect_dual_clairvoyant = comparisons_quickselect_dual(P_CV)

    sage: all(QuickStar().avg_comparisons_quickselect(
    ....:         n, j, 'partitioned_dual_clairvoyant', verbose=True) ==
    ....:     comparisons_quickselect_dual_clairvoyant(n, j+1)
    ....:     for n in srange(7) for j in srange(n))
    n=1, j=0, cmp=0, avg=0
    n=2, j=0, cmp=2, avg=1
    n=2, j=1, cmp=2, avg=1
    n=3, j=0, cmp=14, avg=7/3
    n=3, j=1, cmp=14, avg=7/3
    n=3, j=2, cmp=14, avg=7/3
    n=4, j=0, cmp=94, avg=47/12
    n=4, j=1, cmp=98, avg=49/12
    n=4, j=2, cmp=98, avg=49/12
    n=4, j=3, cmp=94, avg=47/12
    n=5, j=0, cmp=674, avg=337/60
    n=5, j=1, cmp=714, avg=119/20
    n=5, j=2, cmp=730, avg=73/12
    n=5, j=3, cmp=714, avg=119/20
    n=5, j=4, cmp=674, avg=337/60
    n=6, j=0, cmp=5332, avg=1333/180
    n=6, j=1, cmp=5688, avg=79/10
    n=6, j=2, cmp=5900, avg=295/36
    n=6, j=3, cmp=5900, avg=295/36
    n=6, j=4, cmp=5688, avg=79/10
    n=6, j=5, cmp=5332, avg=1333/180
    True

Dual-pivot Quickselect "p first"
--------------------------------

::

    sage: comparisons_quickselect_dual_p_first = comparisons_quickselect_dual(P_pF)

    sage: all(QuickStar().avg_comparisons_quickselect(
    ....:         n, j, 'partitioned_dual_p_first', verbose=True) ==
    ....:     comparisons_quickselect_dual_p_first(n, j+1)
    ....:     for n in srange(7) for j in srange(n))
    n=1, j=0, cmp=0, avg=0
    n=2, j=0, cmp=2, avg=1
    n=2, j=1, cmp=2, avg=1
    n=3, j=0, cmp=16, avg=8/3
    n=3, j=1, cmp=16, avg=8/3
    n=3, j=2, cmp=16, avg=8/3
    n=4, j=0, cmp=108, avg=9/2
    n=4, j=1, cmp=112, avg=14/3
    n=4, j=2, cmp=112, avg=14/3
    n=4, j=3, cmp=108, avg=9/2
    n=5, j=0, cmp=776, avg=97/15
    n=5, j=1, cmp=820, avg=41/6
    n=5, j=2, cmp=840, avg=7
    n=5, j=3, cmp=820, avg=41/6
    n=5, j=4, cmp=776, avg=97/15
    n=6, j=0, cmp=6136, avg=767/90
    n=6, j=1, cmp=6536, avg=817/90
    n=6, j=2, cmp=6792, avg=283/30
    n=6, j=3, cmp=6792, avg=283/30
    n=6, j=4, cmp=6536, avg=817/90
    n=6, j=5, cmp=6136, avg=767/90
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
        from sage.arith.srange import srange
        from sage.combinat.permutation import Permutations
        E = srange(n)
        pis = Permutations(E)
        assert all(self.quicksorted(list(pi), partitioning_strategy) == E
                   for pi in pis)
        if verbose:
            print('n={}, cmp={}, avg={}'.format(
                n, self.comparisons, self.comparisons / n.factorial()))
        return self.comparisons / n.factorial()


    def avg_comparisons_quickselect(self, n, j, partitioning_strategy, verbose=False):
        r"""
        """
        from sage.arith.srange import srange
        from sage.combinat.permutation import Permutations
        E = srange(n)
        pis = Permutations(E)
        assert all(self.quickselected(list(pi), j, partitioning_strategy) == j
                   for pi in pis)
        if verbose:
            print('n={}, j={}, cmp={}, avg={}'.format(
                n, j, self.comparisons, self.comparisons / n.factorial()))
        return self.comparisons / n.factorial()
