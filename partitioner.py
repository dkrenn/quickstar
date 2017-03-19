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

            sage: from partitioner import ClassificationStrategy

            sage: T = ClassificationStrategy(2)[0]
            sage: T.classify_element(5, [10, 20])
            (0, 1)
            sage: T.classify_element(15, [10, 20])
            (1, 2)
            sage: T.classify_element(25, [10, 20])
            (2, 2)

            sage: T = ClassificationStrategy(3)[1]
            sage: T.classify_element(5, [10, 20, 30])
            (0, 1)
            sage: T.classify_element(15, [10, 20, 30])
            (1, 3)
            sage: T.classify_element(25, [10, 20, 30])
            (2, 3)
            sage: T.classify_element(35, [10, 20, 30])
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

    def tree_equals(self, other, mirrored=False):
        r"""
        EXAMPLES:

            sage: from partitioner import ClassificationStrategy

            sage: cs2 = ClassificationStrategy(2)
            sage: Matrix([[a.tree_equals(b) for a in cs2.trees]
            ....:         for b in cs2.trees])
            [1 0]
            [0 1]
            sage: Matrix([[a.tree_equals(b, mirrored=True) for a in cs2.trees]
            ....:         for b in cs2.trees])
            [0 1]
            [1 0]

            sage: cs3 = ClassificationStrategy(3)
            sage: Matrix([[a.tree_equals(b) for a in cs3.trees]
            ....:         for b in cs3.trees])
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            sage: Matrix([[a.tree_equals(b, mirrored=True) for a in cs3.trees]
            ....:         for b in cs3.trees])
            [0 0 0 0 1]
            [0 0 0 1 0]
            [0 0 1 0 0]
            [0 1 0 0 0]
            [1 0 0 0 0]
        """
        indices = self.indices()
        if indices != other.indices():
            return False
        indices = indices[1:]
        mapping = {j: i for i, j
                   in zip(indices, reversed(indices) if mirrored else indices)}

        def eq(s, o):
            if s is None and o is None:
                return True
            if (s is None) != (o is None):
                return False
            if s.top != mapping[o.top]:
                return False
            if mirrored:
                return eq(s.left, o.right) and eq(s.right, o.left)
            else:
                return eq(s.left, o.left) and eq(s.right, o.right)

        return eq(self, other)

    def polyhedron_equals(self, other, mirrored=False):
        if mirrored:
            return self.polyhedron == mirror_polyhedron(other.polyhedron)
        else:
            return self.polyhedron == other.polyhedron

    def __eq__(self, other):
        return self.tree_equals(other) and self.polyhedon == other.polyhedon

    def is_mirroring_tree_of(self, other):
        return self.tree_equals(other, mirrored=True)

    def is_mirroring_polyhedron_of(self, other):
        return self.polyhedron_equals(other, mirrored=True)

    def is_mirroring_of(self, other, polyhedra=True):
        return (self.is_mirroring_tree_of(other) and
                (not polyhedra or self.is_mirroring_polyhedron_of(other)))


def classification_trees(r):
    if not r:
        yield None
        return
    for i, top in enumerate(r):
        for left in classification_trees(r[:i]):
            for right in classification_trees(r[i+1:]):
                yield ClassificationTree(top, left, right)


def mirror_polyhedron(polyhedron):
    r"""
    EXAMPLES::

        sage: from partitioner import mirror_polyhedron

        sage: P = Polyhedron(
        ....:         ieqs=[(1,2,3,0), (0,1,1,0)])
        sage: P.Hrepresentation()
        (An inequality (1, 1, 0) x + 0 >= 0, An inequality (2, 3, 0) x + 1 >= 0)
        sage: mirror_polyhedron(P).Hrepresentation()
        (An inequality (0, 1, 1) x + 0 >= 0, An inequality (0, 3, 2) x + 1 >= 0)

        sage: P = Polyhedron(
        ....:         ieqs=[(1,2,3,0), (0,1,1,0)],
        ....:         eqns=[(2,1,0,1)])
        sage: P.Hrepresentation()
        (An equation (1, 0, 1) x + 2 == 0,
         An inequality (1, 1, 0) x + 0 >= 0,
         An inequality (2, 3, 0) x + 1 >= 0)
        sage: mirror_polyhedron(P).Hrepresentation()
        (An equation (1, 0, 1) x + 2 == 0,
         An inequality (-1, 1, 0) x - 2 >= 0,
         An inequality (-2, 3, 0) x - 3 >= 0)
    """
    from sage.geometry.polyhedron.constructor import Polyhedron

    if polyhedron is None:
        return polyhedron

    def mirror(vec):
        return vec[:1] + tuple(reversed(vec[1:]))

    return Polyhedron(
        eqns=[mirror(tuple(eqn)) for eqn in polyhedron.equations()],
        ieqs=[mirror(tuple(ieq)) for ieq in polyhedron.inequalities()])


def break_tie(values, undo=False):
    r"""
    EXAMPLES::

        sage: from partitioner import ClassificationStrategy, break_tie

        sage: for cs in ClassificationStrategy(2):
        ....:     print(Polyhedron(
        ....:         ieqs=[break_tie(tuple(ieq))
        ....:               for ieq in cs.polyhedron.inequalities()]).repr_pretty_Hrepresentation())
        x2 >= 0, x1 >= 0, x0 >= x2
        x2 >= x0 + 1, x1 >= 0, x0 >= 0
        sage: for cs in ClassificationStrategy(3):
        ....:     print(Polyhedron(
        ....:         ieqs=[break_tie(tuple(ieq))
        ....:               for ieq in cs.polyhedron.inequalities()]).repr_pretty_Hrepresentation())
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


class ClassificationStrategy(SageObject):
    r"""
    EXAMPLES::

        sage: from partitioner import ClassificationStrategy

        sage: cs2 = ClassificationStrategy(2)
        sage: cs2
        ====================================
        classification strategy with 2 trees
        ------------------------------------
        classification tree 1 () (2 () ())
        s2 >= 0
        s1 >= 0
        s0 >= s2
        ------------------------------------
        classification tree 2 (1 () ()) ()
        s2 >= s0
        s0 >= 0
        s1 >= 0

        sage: cs3 = ClassificationStrategy(3)
        sage: cs3
        ====================================
        classification strategy with 5 trees
        ------------------------------------
        classification tree 1 () (2 () (3 () ()))
        s3 >= 0
        s2 >= 0
        s1 >= s3
        s0 >= s2 + s3 + 1
        ------------------------------------
        classification tree 1 () (3 (2 () ()) ())
        s3 >= s1
        s0 >= s3
        s2 >= 0
        s1 >= 0
        s0 >= s1 + s2 + 1
        ------------------------------------
        classification tree 2 (1 () ()) (3 () ())
        s2 + s3 + 1 >= s0
        s1 + s2 + 1 >= s0
        s3 >= 0
        s2 >= 0
        s1 >= 0
        s1 + s2 + 1 >= s3
        s0 >= 0
        s0 + s1 + 1 >= s3
        ------------------------------------
        classification tree 3 (1 () (2 () ())) ()
        s0 >= s2
        s3 >= s0
        s3 >= s1 + s2 + 1
        s1 >= 0
        s2 >= 0
        ------------------------------------
        classification tree 3 (2 (1 () ()) ()) ()
        s0 >= 0
        s1 >= 0
        s3 >= s0 + s1 + 1
        s2 >= s0

    TESTS::

        sage: cs2 = ClassificationStrategy(2, make_disjoint=True)
        sage: cs2
        ====================================
        classification strategy with 2 trees
        ------------------------------------
        classification tree 1 () (2 () ())
        s2 >= 0
        s1 >= 0
        s0 >= s2
        ------------------------------------
        classification tree 2 (1 () ()) ()
        s2 > s0
        s1 >= 0
        s0 >= 0

        sage: cs3 = ClassificationStrategy(3, make_disjoint=True)
        sage: cs3
        ====================================
        classification strategy with 5 trees
        ------------------------------------
        classification tree 1 () (2 () (3 () ()))
        s3 >= 0
        s2 >= 0
        s1 >= s3
        s0 > s2 + s3 + 1
        ------------------------------------
        classification tree 1 () (3 (2 () ()) ())
        s3 > s1
        s2 >= 0
        s1 >= 0
        s0 > s1 + s2 + 1
        s0 >= s3
        ------------------------------------
        classification tree 2 (1 () ()) (3 () ())
        s2 + s3 + 1 >= s0
        s1 + s2 + 1 >= s0
        s3 >= 0
        s2 >= 0
        s1 >= 0
        s1 + s2 + 1 >= s3
        s0 >= 0
        s0 + s1 + 1 >= s3
        ------------------------------------
        classification tree 3 (1 () (2 () ())) ()
        s3 > s0
        s3 > s1 + s2 + 1
        s2 >= 0
        s1 >= 0
        s0 > s2
        ------------------------------------
        classification tree 3 (2 (1 () ()) ()) ()
        s3 > s0 + s1 + 1
        s2 >= s0
        s1 >= 0
        s0 >= 0
    """

    def __init__(self, dimension,
                       make_disjoint=False,
                       trees=None):
        from sage.rings.integer_ring import ZZ

        super(ClassificationStrategy, self).__init__()

        self._disjoint_ = make_disjoint
        self._dimension_ = ZZ(dimension)

        if trees is None:
            trees = iter(p for p in classification_trees(range(1, dimension + 1)))
        self.trees = tuple(trees)
        assert self.dimension() == len(self.indices()) - 1

        self.populate_polyhedra()
        assert self.dimension() == self.trees[0].polyhedron.ambient_dim() - 1

    def populate_polyhedra(self):
        from sage.geometry.polyhedron.constructor import Polyhedron
        from sage.symbolic.ring import SR

        def get_vector(inequality, vars):
            coefficients = list(inequality.coefficient(var) for var in vars)
            constant = inequality - sum(c*v for c, v in zip(coefficients, vars))
            return [constant] + coefficients

        d = self.dimension()
        prefix = self.trees[0].PREFIX
        vars = list(SR(prefix + "{}".format(j)) for j in range(d+1))
        for tree in self.trees:
            others = list(self.trees)
            others.remove(tree)
            ineqs = [other.partition_cost() - tree.partition_cost()
                     for other in others] + vars
            ineq_matrix = [get_vector(ineq, vars) for ineq in ineqs]
            P = Polyhedron(ieqs=ineq_matrix)
            if self._disjoint_:
                P = polyhedron_break_tie(P)
            tree.polyhedron = P

        nonnegative_orthant = Polyhedron(ieqs=[dd*(0,) + (1,) + (d+1-dd)*(0,)
                                               for dd in range(1, d+1+1)])
        assert all(A.polyhedron & nonnegative_orthant == A.polyhedron
                   for A in self.trees)
        if self._disjoint_:
            assert all((A.polyhedron & B.polyhedron).is_empty()
                       for A in self.trees for B in self.trees if A != B)

    def __iter__(self):
        return iter(self.trees)

    def __getitem__(self, i):
        return self.trees[i]

    def _repr_(self):
        return (
            '====================================' + '\n' +
            'classification strategy with {} trees'.format(len(self.trees)) +
            '\n' +
            '\n'.join(
            '------------------------------------' + '\n' +
            'classification tree {}'.format(tree) + '\n' +
            tree.repr_pretty_polyhedron(
                strict_inequality=self._disjoint_)
            for tree in self.trees))

    def indices(self):
        assert len(set(tree.indices() for tree in self.trees)) <= 1
        return self.trees[0].indices()

    def dimension(self):
        return self._dimension_

    def H(self):
        from sage.modules.free_module_element import vector
        return {i: vector(tree.height(i) for tree in self.trees)
                for i in self.indices()}

    def next_classification_tree_by_counts(self, counts):
        return next(tree for tree in self.trees
                    if tree.polyhedron.contains(counts))

    def is_symmetric(self, polyhedra=True):
        r"""
        EXAMPLES::

            sage: from partitioner import ClassificationStrategy
            sage: ClassificationStrategy(2).is_symmetric()
            True
            sage: ClassificationStrategy(3).is_symmetric()
            True
            sage: ClassificationStrategy(4).is_symmetric()
            True
        """
        c = (len(self.trees) + 1) // 2
        return all(a.is_mirroring_tree_of(b)
                   for a, b in zip(self.trees[:c],
                                   tuple(reversed(self.trees))[:c]))

    def nonempty_subsets(self):
        r"""
        EXAMPLES::

            sage: from partitioner import ClassificationStrategy
            sage: cs2 = ClassificationStrategy(2, make_disjoint=True)
            sage: for cs in cs2.nonempty_subsets():
            ....:     print(cs)
            ====================================
            classification strategy with 1 trees
            ------------------------------------
            classification tree 1 () (2 () ())
            s2 >= 0
            s1 >= 0
            s0 >= 0
            ====================================
            classification strategy with 1 trees
            ------------------------------------
            classification tree 2 (1 () ()) ()
            s2 >= 0
            s1 >= 0
            s0 >= 0
            ====================================
            classification strategy with 2 trees
            ------------------------------------
            classification tree 1 () (2 () ())
            s2 >= 0
            s1 >= 0
            s0 >= s2
            ------------------------------------
            classification tree 2 (1 () ()) ()
            s2 > s0
            s1 >= 0
            s0 >= 0
        """
        from sage.misc.misc import subsets
        d = self.dimension()
        return iter(ClassificationStrategy(d,
                                           make_disjoint=self._disjoint_,
                                           trees=srees)
                    for srees in subsets(self.trees)
                    if srees)

    def symmetric_subsets(self):
        r"""
        EXAMPLES::

            sage: from partitioner import ClassificationStrategy

            sage: cs2 = ClassificationStrategy(2, make_disjoint=True)
            sage: for cs in cs2.symmetric_subsets():
            ....:     print(cs)
            ====================================
            classification strategy with 2 trees
            ------------------------------------
            classification tree 1 () (2 () ())
            s2 >= 0
            s1 >= 0
            s0 >= s2
            ------------------------------------
            classification tree 2 (1 () ()) ()
            s2 > s0
            s1 >= 0
            s0 >= 0

            sage: cs3 = ClassificationStrategy(3, make_disjoint=True)
            sage: for cs in cs3.symmetric_subsets():
            ....:     print(cs)
            ====================================
            classification strategy with 1 trees
            ------------------------------------
            classification tree 2 (1 () ()) (3 () ())
            s3 >= 0
            s2 >= 0
            s1 >= 0
            s0 >= 0
            ====================================
            classification strategy with 2 trees
            ------------------------------------
            classification tree 1 () (3 (2 () ()) ())
            s3 >= 0
            s2 >= 0
            s1 >= 0
            s0 >= s3
            ------------------------------------
            classification tree 3 (1 () (2 () ())) ()
            s3 > s0
            s2 >= 0
            s1 >= 0
            s0 >= 0
            ====================================
            classification strategy with 3 trees
            ------------------------------------
            classification tree 1 () (3 (2 () ()) ())
            s3 >= 0
            s2 >= 0
            s1 >= 0
            s0 > s1 + s2 + 1
            s0 >= s3
            ------------------------------------
            classification tree 2 (1 () ()) (3 () ())
            s1 + s2 + 1 >= s0
            s3 >= 0
            s2 >= 0
            s1 >= 0
            s1 + s2 + 1 >= s3
            s0 >= 0
            ------------------------------------
            classification tree 3 (1 () (2 () ())) ()
            s3 > s0
            s3 > s1 + s2 + 1
            s2 >= 0
            s1 >= 0
            s0 >= 0
            ====================================
            classification strategy with 2 trees
            ------------------------------------
            classification tree 1 () (2 () (3 () ()))
            s3 >= 0
            s2 >= 0
            s1 >= 0
            s0 >= 0
            2*s0 + s1 >= s2 + 2*s3
            ------------------------------------
            classification tree 3 (2 (1 () ()) ()) ()
            s2 + 2*s3 > 2*s0 + s1
            s3 >= 0
            s2 >= 0
            s1 >= 0
            s0 >= 0
            ====================================
            classification strategy with 3 trees
            ------------------------------------
            classification tree 1 () (2 () (3 () ()))
            s3 >= 0
            s2 >= 0
            s1 >= 0
            s0 > s2 + s3 + 1
            ------------------------------------
            classification tree 2 (1 () ()) (3 () ())
            s2 + s3 + 1 >= s0
            s3 >= 0
            s2 >= 0
            s1 >= 0
            s0 >= 0
            s0 + s1 + 1 >= s3
            ------------------------------------
            classification tree 3 (2 (1 () ()) ()) ()
            s3 > s0 + s1 + 1
            s2 >= 0
            s1 >= 0
            s0 >= 0
            ====================================
            classification strategy with 4 trees
            ------------------------------------
            classification tree 1 () (2 () (3 () ()))
            s3 >= 0
            s2 >= 0
            s1 >= s3
            s0 >= 0
            s0 + s1 >= 2*s3
            2*s0 + s1 >= s2 + 2*s3
            ------------------------------------
            classification tree 1 () (3 (2 () ()) ())
            s3 > s1
            s2 >= 0
            s1 >= 0
            s0 >= s3
            2*s0 > s2 + s3
            ------------------------------------
            classification tree 3 (1 () (2 () ())) ()
            2*s3 > s0 + s1
            s3 > s0
            s2 >= 0
            s1 >= 0
            s0 > s2
            ------------------------------------
            classification tree 3 (2 (1 () ()) ()) ()
            s2 + 2*s3 > 2*s0 + s1
            s2 + s3 >= 2*s0
            s2 >= s0
            s3 >= 0
            s1 >= 0
            s0 >= 0
            ====================================
            classification strategy with 5 trees
            ------------------------------------
            classification tree 1 () (2 () (3 () ()))
            s3 >= 0
            s2 >= 0
            s1 >= s3
            s0 > s2 + s3 + 1
            ------------------------------------
            classification tree 1 () (3 (2 () ()) ())
            s3 > s1
            s2 >= 0
            s1 >= 0
            s0 > s1 + s2 + 1
            s0 >= s3
            ------------------------------------
            classification tree 2 (1 () ()) (3 () ())
            s2 + s3 + 1 >= s0
            s1 + s2 + 1 >= s0
            s3 >= 0
            s2 >= 0
            s1 >= 0
            s1 + s2 + 1 >= s3
            s0 >= 0
            s0 + s1 + 1 >= s3
            ------------------------------------
            classification tree 3 (1 () (2 () ())) ()
            s3 > s0
            s3 > s1 + s2 + 1
            s2 >= 0
            s1 >= 0
            s0 > s2
            ------------------------------------
            classification tree 3 (2 (1 () ()) ()) ()
            s3 > s0 + s1 + 1
            s2 >= s0
            s1 >= 0
            s0 >= 0
        """
        if not self.is_symmetric(polyhedra=True):
            raise ValueError('classification strategy is not symmetric')
        return iter(cs
                    for cs in self.nonempty_subsets()
                    if cs.is_symmetric(polyhedra=False))
