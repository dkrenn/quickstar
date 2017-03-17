from __future__ import absolute_import
from __future__ import print_function

from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import SageObject


class ClassificationTree(SageObject):

    PREFIX = 's'

    def __init__(self, top, left, right):
        self.top = top
        elements = []
        if left is not None:
            elements += list(left)
        elements += [top]
        if right is not None:
            elements += list(right)
        self.minimum = min(elements)
        self.maximum = max(elements)
        if not elements == range(self.minimum, self.maximum + 1):
            raise ValueError("elements {} do not form a consecutive range")
        self.left = left
        self.right = right
        self.polyhedron = None


    def __iter__(self):
        if self.left is not None:
            for x in self.left:
                yield x
        yield self.top
        if self.right is not None:
            for x in self.right:
                yield x


    def __repr__(self):
        def no_None(e):
            if e is None:
                return ""
            else:
                return e
        return "{} ({}) ({})".format(self.top, no_None(self.left), no_None(self.right))


    def __hash__(self):
        return hash(repr(self))


    def height(self, i):
        if i < self.minimum - 1 or i > self.maximum:
            raise ValueError("{} must be in [{}, {}]".format(i, self.minimum-1, self.maximum))
        if i < self.top:
            if self.left is None:
                return 1
            else: 
                return 1 + self.left.height(i)
        else:
            if self.right is None:
                return 1
            else:
                return 1 + self.right.height(i)


    def set_polyhedron(self, polyhedron):
        self.polyhedron = polyhedron
        

    def _coefficients_inequalities_(self):
        for i in self.polyhedron.Hrepresentation():
            if not i.is_inequality():
                raise NotImplementedError("{} is not an inequality".format(i))
            yield tuple(i.vector())


    def indices(self):
        return tuple(range(self.minimum - 1, self.maximum + 1))


    def pretty_inequalities(self):
        from sage.geometry.polyhedron.representation import repr_pretty
        for coeffs in self._coefficients_inequalities_():
            print(repr_pretty(coeffs, 0, indices=self.indices(), prefix=self.PREFIX))


    def repr_pretty_polyhedron(self, strict_inequality=False):
        return repr_pretty_Hrepresentation(
            self.polyhedron,
            strict_inequality=strict_inequality,
            indices=self.indices(), prefix=self.PREFIX, separator='\n')


    def generating_function(self):
        return generating_function_of_partitioners((self,))


    @cached_method        
    def partition_cost(self):
        from sage.symbolic.ring import SR
        return sum(self.height(i) * (SR(self.PREFIX + "{}".format(i)) + 1)
                  for i in range(self.minimum - 1, self.maximum + 1))
    

    def classify_element(self, element, pivots):
        r"""
        EXAMPLES::

            sage: from partitioner import classification_strategy

            sage: p = classification_strategy(2, verbose=False)[0]
            sage: p.classify_element(5, [10, 20])
            (0, 1)
            sage: p.classify_element(15, [10, 20])
            (1, 2)
            sage: p.classify_element(25, [10, 20])
            (2, 2)

            sage: p = classification_strategy(3, verbose=False)[1]
            sage: p.classify_element(5, [10, 20, 30])
            (0, 1)
            sage: p.classify_element(15, [10, 20, 30])
            (1, 3)
            sage: p.classify_element(25, [10, 20, 30])
            (2, 3)
            sage: p.classify_element(35, [10, 20, 30])
            (3, 2)
        """
        comparisons = 0
        current = self
        while True:
            comparisons += 1
            if element < pivots[current.top-1]:
                if current.left is None:
                    return current.top-1, comparisons
                current = current.left
            else:
                if current.right is None:
                    return current.top, comparisons
                current = current.right


def classification_trees(r):
    if not r:
        yield None
        return
    for i, top in enumerate(r):
        for left in classification_trees(r[:i]):
            for right in classification_trees(r[i+1:]):
                yield ClassificationTree(top, left, right)


def break_tie(values, undo=False):
    r"""
    EXAMPLES::

        sage: from partitioner import classification_strategy, break_tie

        sage: for p in classification_strategy(2, verbose=False):
        ....:     print(Polyhedron(
        ....:         ieqs=[break_tie(tuple(ieq))
        ....:               for ieq in p.polyhedron.inequalities()]).repr_pretty_Hrepresentation())
        x2 >= 0, x1 >= 0, x0 >= x2
        x2 >= x0 + 1, x1 >= 0, x0 >= 0
        sage: for p in classification_strategy(3, verbose=False):
        ....:     print(Polyhedron(
        ....:         ieqs=[break_tie(tuple(ieq))
        ....:               for ieq in p.polyhedron.inequalities()]).repr_pretty_Hrepresentation())
        x3 >= 0, x2 >= 0, x1 >= x3, x0 >= x2 + x3 + 2
        x3 >= x1 + 1, x2 >= 0, x1 >= 0, x0 >= x1 + x2 + 2, x0 >= x3
        x2 + x3 + 1 >= x0, x1 + x2 + 1 >= x0, x3 >= 0, x2 >= 0, x1 >= 0, x1 + x2 + 1 >= x3, x0 >= 0, x0 + x1 + 1 >= x3
        x3 >= x0 + 1, x3 >= x1 + x2 + 2, x2 >= 0, x1 >= 0, x0 >= x2 + 1
        x3 >= x0 + x1 + 2, x2 >= x0, x1 >= 0, x0 >= 0
    """
    if _is_strict_(values):
        return _make_strict_(values, undo=undo)
    else:
        return tuple(values)


def _is_strict_(values):
    nc_values = values[1:]
    k = len(nc_values)
    left = tuple(i for i, t in enumerate(nc_values) if t > 0)
    right = tuple(i for i, t in enumerate(nc_values) if t < 0)
    return strict_inequality_symmetric_choice(k, left, right)


def _make_strict_(values, undo=False):
    from sage.arith.misc import lcm
    from sage.rings.rational_field import QQ

    nc_values = values[1:]
    den = lcm([QQ(t).denominator() for t in nc_values])
    if undo:
        den = -den
    return (values[0] - den,) + tuple(nc_values)


def strict_inequality_symmetric_choice(k, left, right):
    r"""
    TESTS::

        sage: from partitioner import strict_inequality_symmetric_choice
        sage: S = Set(srange(5))
        sage: for s in S.subsets():
        ....:     L = tuple(s)
        ....:     R = tuple(S.difference(s))
        ....:     a = strict_inequality_symmetric_choice(5, L, R)
        ....:     b = strict_inequality_symmetric_choice(5, R, L)
        ....:     assert a != b

        sage: strict_inequality_symmetric_choice(1, [], [0])
        True
    """
    from sage.rings.rational_field import QQ

    assert not set(left) & set(right)
    center = QQ(k - 1) / QQ(2)
    def weight(t):
        return -(center + 1 - (t-center).abs())
    def total(T):
        return sum(weight(t) for t in T)
    def rule(T):
        return total(T), sorted((weight(t), t) for t in T)
    return rule(left) > rule(right)


def polyhedron_break_tie(polyhedron, undo=False):
    from sage.geometry.polyhedron.constructor import Polyhedron

    if polyhedron.equations():
        raise NotImplementedError
    return Polyhedron(ieqs=[break_tie(tuple(ieq), undo=undo)
                            for ieq in polyhedron.inequalities()])


def split_polyhedra(dim):
    r"""
    ::

        sage: from partitioner import split_polyhedra, repr_pretty_Hrepresentation

        sage: for P in split_polyhedra(2):
        ....:     print(repr_pretty_Hrepresentation(P, strict_inequality=True,
        ....:                                       prefix='s'))
        s1 > s0
        s0 >= s1

        sage: for P in split_polyhedra(3):
        ....:     print(repr_pretty_Hrepresentation(P, strict_inequality=True,
        ....:                                       prefix='s'))
        s1 >= s0, s2 > s1
        s2 > s0, s1 >= s2
        s2 > s0, s0 > s1
        s2 > s1, s0 >= s2
        s1 >= s0, s0 >= s2
        s1 >= s2, s0 > s1

        sage: for P in split_polyhedra(4):
        ....:     print(repr_pretty_Hrepresentation(P, strict_inequality=True,
        ....:                                       prefix='s'))
        s1 >= s0, s2 > s1, s3 > s2
        s1 >= s0, s3 > s1, s2 >= s3
        s2 >= s0, s3 > s1, s1 >= s2
        s2 >= s0, s3 > s2, s1 >= s3
        s3 > s0, s2 > s1, s1 >= s3
        s3 > s0, s2 >= s3, s1 >= s2
        s2 >= s0, s3 > s2, s0 > s1
        s3 > s0, s2 >= s3, s0 > s1
        s3 > s0, s2 > s1, s0 > s2
        s2 > s1, s3 > s2, s0 >= s3
        s2 >= s0, s3 > s1, s0 >= s3
        s3 > s1, s2 >= s3, s0 > s2
        s1 >= s0, s3 > s1, s0 > s2
        s3 > s0, s1 >= s3, s0 > s2
        s3 > s0, s1 >= s2, s0 > s1
        s3 > s1, s1 >= s2, s0 >= s3
        s1 >= s0, s3 > s2, s0 >= s3
        s3 > s2, s1 >= s3, s0 > s1
        s1 >= s0, s2 > s1, s0 >= s3
        s2 >= s0, s1 >= s2, s0 >= s3
        s2 >= s0, s1 >= s3, s0 > s1
        s2 > s1, s1 >= s3, s0 > s2
        s1 >= s0, s2 >= s3, s0 > s2
        s2 >= s3, s1 >= s2, s0 > s1
    """
    from sage.combinat.permutation import Permutations
    from sage.geometry.polyhedron.constructor import Polyhedron

    return iter(
        polyhedron_break_tie(
            Polyhedron(
                ieqs=[tuple(1 if i==b else (-1 if i==a else 0)
                            for i in range(dim+1))
                      for a, b in zip(pi[:-1], pi[1:])]))
        for pi in Permutations(dim))


def repr_pretty_Hrepresentation(self, separator=', ',
                                strict_inequality=False,
                                **kwds):
    if not strict_inequality:
        return self.repr_pretty_Hrepresentation(separator=separator,
                                                **kwds)

    P = polyhedron_break_tie(self, undo=True)
    def repr_ieq(ieq):
        rp = ieq.repr_pretty(**kwds)
        if _is_strict_(tuple(ieq)):
            rp = rp.replace('>=', '>')
        return rp

    return separator.join(repr_ieq(h) for h in P.Hrepresentation())


def new_ClassificationStrategy(arg, disjoint=False):
    cs = ClassificationStrategy(arg)
    cs._disjoint_ = disjoint
    return cs


class ClassificationStrategy(tuple):

    def __repr__(self):
        return '\n'.join(
            '*********************************' + '\n' +
            'classification tree {}'.format(tree) + '\n' +
            tree.repr_pretty_polyhedron(
                strict_inequality=self._disjoint_)
            for tree in self)

    def indices(self):
        assert len(set(tree.indices() for tree in self)) <= 1
        return self[0].indices()

    def d(self):
        from sage.rings.integer_ring import ZZ
        d = len(self.indices()) - 1
        assert d == self[0].polyhedron.ambient_dim() - 1
        return ZZ(d)

    def H(self):
        from sage.modules.free_module_element import vector
        return {i: vector(tree.height(i) for tree in self)
                for i in self.indices()}

    def next_classification_tree_by_counts(self, counts):
        return next(tree for tree in self
                    if tree.polyhedron.contains(counts))

    def nonempty_subsets(self):
        from sage.misc.misc import subsets
        return iter(new_ClassificationStrategy(s, disjoint=self._disjoint_)
                    for s in subsets(self)
                    if s)


def classification_strategy(dimension,
                            make_disjoint=False,
                            verbose=True,
                            trees=None):
    r"""
    EXAMPLES::

        sage: from partitioner import classification_strategy

        sage: p2 = classification_strategy(2)
        *********************************
        classification tree 1 () (2 () ())
        s2 >= 0
        s1 >= 0
        s0 >= s2
        *********************************
        classification tree 2 (1 () ()) ()
        s2 >= s0
        s0 >= 0
        s1 >= 0

        sage: p3 = classification_strategy(3)
        *********************************
        classification tree 1 () (2 () (3 () ()))
        s3 >= 0
        s2 >= 0
        s1 >= s3
        s0 >= s2 + s3 + 1
        *********************************
        classification tree 1 () (3 (2 () ()) ())
        s3 >= s1
        s0 >= s3
        s2 >= 0
        s1 >= 0
        s0 >= s1 + s2 + 1
        *********************************
        classification tree 2 (1 () ()) (3 () ())
        s2 + s3 + 1 >= s0
        s1 + s2 + 1 >= s0
        s3 >= 0
        s2 >= 0
        s1 >= 0
        s1 + s2 + 1 >= s3
        s0 >= 0
        s0 + s1 + 1 >= s3
        *********************************
        classification tree 3 (1 () (2 () ())) ()
        s0 >= s2
        s3 >= s0
        s3 >= s1 + s2 + 1
        s1 >= 0
        s2 >= 0
        *********************************
        classification tree 3 (2 (1 () ()) ()) ()
        s0 >= 0
        s1 >= 0
        s3 >= s0 + s1 + 1
        s2 >= s0

    TESTS::

        sage: p2 = classification_strategy(2, make_disjoint=True)
        *********************************
        classification tree 1 () (2 () ())
        s2 >= 0
        s1 >= 0
        s0 >= s2
        *********************************
        classification tree 2 (1 () ()) ()
        s2 > s0
        s1 >= 0
        s0 >= 0

        sage: p3 = classification_strategy(3, make_disjoint=True)
        *********************************
        classification tree 1 () (2 () (3 () ()))
        s3 >= 0
        s2 >= 0
        s1 >= s3
        s0 > s2 + s3 + 1
        *********************************
        classification tree 1 () (3 (2 () ()) ())
        s3 > s1
        s2 >= 0
        s1 >= 0
        s0 > s1 + s2 + 1
        s0 >= s3
        *********************************
        classification tree 2 (1 () ()) (3 () ())
        s2 + s3 + 1 >= s0
        s1 + s2 + 1 >= s0
        s3 >= 0
        s2 >= 0
        s1 >= 0
        s1 + s2 + 1 >= s3
        s0 >= 0
        s0 + s1 + 1 >= s3
        *********************************
        classification tree 3 (1 () (2 () ())) ()
        s3 > s0
        s3 > s1 + s2 + 1
        s2 >= 0
        s1 >= 0
        s0 > s2
        *********************************
        classification tree 3 (2 (1 () ()) ()) ()
        s3 > s0 + s1 + 1
        s2 >= s0
        s1 >= 0
        s0 >= 0
    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    from sage.symbolic.ring import SR

    def get_vector(inequality, vars):
        coefficients = list(inequality.coefficient(var) for var in vars)
        constant = inequality - sum(c*v for c, v in zip(coefficients, vars))
        return [constant] + coefficients
    
    if trees is None:
        trees = iter(p for p in classification_trees(range(1, dimension + 1)))
    strategy = new_ClassificationStrategy(trees, disjoint=make_disjoint)

    prefix = strategy[0].PREFIX
    vars = list(SR(prefix + "{}".format(j)) for j in range(dimension+1))
    for this in strategy:
        others = list(strategy)
        others.remove(this)
        ineqs = [other.partition_cost() - this.partition_cost() for other in others] + vars
        ineq_matrix = [get_vector(ineq, vars) for ineq in ineqs]
        P = Polyhedron(ieqs=ineq_matrix)
        if make_disjoint:
            P = polyhedron_break_tie(P)
        this.polyhedron = P

    dim = strategy[0].polyhedron.ambient_dim()
    nonnegative_orthant = Polyhedron(ieqs=[dd*(0,) + (1,) + (dim-dd)*(0,)
                                           for dd in range(1, dim+1)])
    assert all(A.polyhedron & nonnegative_orthant == A.polyhedron
               for A in strategy)
    if make_disjoint:
        assert all((A.polyhedron & B.polyhedron).is_empty()
                   for A in strategy for B in strategy if A != B)
    if verbose:
        print(repr(strategy))
    return strategy
