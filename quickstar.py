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

    Dual-pivot quicksort::

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

    Dual-pivot quickselect::

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
