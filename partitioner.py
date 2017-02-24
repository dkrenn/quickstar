from __future__ import absolute_import
from __future__ import print_function

from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import SageObject


class Partitioner(SageObject):

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
        self.polytope = None


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


    def set_polytope(self, polytope):
        self.polytope = polytope
        

    def _coefficients_inequalities_(self):
        for i in self.polytope.Hrepresentation():
            if not i.is_inequality():
                raise NotImplementedError("{} is not an inequality".format(i))
            yield tuple(i.vector())


    def _indices_(self):
        return tuple(range(self.minimum - 1, self.maximum + 1))


    def pretty_inequalities(self):
        from sage.geometry.polyhedron.representation import repr_pretty
        for coeffs in self._coefficients_inequalities_():
            print(repr_pretty(coeffs, 0, indices=self._indices_(), prefix='b'))


    def repr_pretty_polytope(self):
        return self.polytope.repr_pretty_Hrepresentation(
            indices=self._indices_(), prefix='b', separator='\n')


    def generating_function(self):
        return generating_function_of_partitioners((self,))


    @cached_method        
    def partition_cost(self):
        from sage.symbolic.ring import SR
        return sum(self.height(i) * (SR("b{}".format(i)) + 1)
                  for i in range(self.minimum - 1, self.maximum + 1))
    

    def classify_element(self, element, pivots):
        r"""
        EXAMPLES::

            sage: p = get_polytopes(2, verbose=False)[0]
            sage: p.classify_element(5, [10, 20])
            (0, 1)
            sage: p.classify_element(15, [10, 20])
            (1, 2)
            sage: p.classify_element(25, [10, 20])
            (2, 2)

            sage: p = get_polytopes(3, verbose=False)[1]
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


def get_partitioners(r):
    if not r:
        yield None
        return
    for i, top in enumerate(r):
        for left in get_partitioners(r[:i]):
            for right in get_partitioners(r[i+1:]):
                yield Partitioner(top, left, right)


def break_tie(values):
    r"""
    EXAMPLES::

        sage: for p in get_polytopes(2, verbose=False):
        ....:     print(Polyhedron(
        ....:         ieqs=[break_tie(tuple(ieq))
        ....:               for ieq in p.polytope.inequalities()]).repr_pretty_Hrepresentation())
        x2 >= 0, x1 >= 0, x0 >= x2
        x2 >= x0 + 1, x1 >= 0, x0 >= 0
        sage: for p in get_polytopes(3, verbose=False):
        ....:     print(Polyhedron(
        ....:         ieqs=[break_tie(tuple(ieq))
        ....:               for ieq in p.polytope.inequalities()]).repr_pretty_Hrepresentation())
        x3 >= 0, x2 >= 0, x1 >= x3 + 1, x0 >= x2 + x3 + 1
        x3 >= x1, x2 >= 0, x1 >= 0, x0 >= x1 + x2 + 1, x0 >= x3
        x2 + x3 >= x0, x1 + x2 >= x0, x3 >= 0, x2 >= 0, x1 >= 0, x1 + x2 >= x3, x0 >= 0, x0 + x1 >= x3
        x3 >= x0 + 1, x3 >= x1 + x2 + 1, x2 >= 0, x1 >= 0, x0 >= x2
        x3 >= x0 + x1 + 1, x2 >= x0 + 1, x1 >= 0, x0 >= 0
    """
    nc_values = values[1:]
    k = len(nc_values)
    left = tuple(i for i, t in enumerate(nc_values) if t > 0)
    right = tuple(i for i, t in enumerate(nc_values) if t < 0)
    if strict_inequality_symmetric_choice(k, left, right):
        den = lcm([QQ(t).denominator() for t in nc_values])
        return (values[0] - den,) + tuple(nc_values)
    else:
        return tuple(values)


def strict_inequality_symmetric_choice(k, left, right):
    r"""
    TESTS::

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
    assert not set(left) & set(right)
    k = QQ(k - 1) / QQ(2)
    def distances_to_center(T):
        return sorted((((t-k).abs().ceil(), -t) for t in T))
    return distances_to_center(left) < distances_to_center(right)


def polyhedron_break_tie(polyhedron):
    return Polyhedron(ieqs=[break_tie(tuple(ieq))
                            for ieq in polyhedron.inequalities()])


def get_polytopes(dimension, make_disjoint=False, verbose=True):
    r"""
    EXAMPLES::

        sage: p2 = get_polytopes(2)
        *********************************
        Partitioner 1 () (2 () ())
        b2 >= 0
        b1 >= 0
        b0 >= b2
        *********************************
        Partitioner 2 (1 () ()) ()
        b2 >= b0
        b0 >= 0
        b1 >= 0
        sage: p3 = get_polytopes(3)
        *********************************
        Partitioner 1 () (2 () (3 () ()))
        b3 >= 0
        b2 >= 0
        b1 >= b3
        b0 >= b2 + b3 + 1
        *********************************
        Partitioner 1 () (3 (2 () ()) ())
        b3 >= b1
        b0 >= b3
        b2 >= 0
        b1 >= 0
        b0 >= b1 + b2 + 1
        *********************************
        Partitioner 2 (1 () ()) (3 () ())
        b2 + b3 + 1 >= b0
        b1 + b2 + 1 >= b0
        b3 >= 0
        b2 >= 0
        b1 >= 0
        b1 + b2 + 1 >= b3
        b0 >= 0
        b0 + b1 + 1 >= b3
        *********************************
        Partitioner 3 (1 () (2 () ())) ()
        b0 >= b2
        b3 >= b0
        b3 >= b1 + b2 + 1
        b1 >= 0
        b2 >= 0
        *********************************
        Partitioner 3 (2 (1 () ()) ()) ()
        b0 >= 0
        b1 >= 0
        b3 >= b0 + b1 + 1
        b2 >= b0

    TESTS::

        sage: p2 = get_polytopes(2, make_disjoint=True)
        *********************************
        Partitioner 1 () (2 () ())
        b2 >= 0
        b1 >= 0
        b0 >= b2
        *********************************
        Partitioner 2 (1 () ()) ()
        b2 >= b0 + 1
        b1 >= 0
        b0 >= 0
        sage: p3 = get_polytopes(3, make_disjoint=True)
        *********************************
        Partitioner 1 () (2 () (3 () ()))
        b3 >= 0
        b2 >= 0
        b1 >= b3 + 1
        b0 >= b2 + b3 + 1
        *********************************
        Partitioner 1 () (3 (2 () ()) ())
        b3 >= b1
        b2 >= 0
        b1 >= 0
        b0 >= b1 + b2 + 1
        b0 >= b3
        *********************************
        Partitioner 2 (1 () ()) (3 () ())
        b2 + b3 >= b0
        b1 + b2 >= b0
        b3 >= 0
        b2 >= 0
        b1 >= 0
        b1 + b2 >= b3
        b0 >= 0
        b0 + b1 >= b3
        *********************************
        Partitioner 3 (1 () (2 () ())) ()
        b3 >= b0 + 1
        b3 >= b1 + b2 + 1
        b2 >= 0
        b1 >= 0
        b0 >= b2
        *********************************
        Partitioner 3 (2 (1 () ()) ()) ()
        b3 >= b0 + b1 + 1
        b2 >= b0 + 1
        b1 >= 0
        b0 >= 0
    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    from sage.symbolic.ring import SR

    def get_vector(inequality, vars):
        coefficients = list(inequality.coefficient(var) for var in vars)
        constant = inequality - sum(c*v for c, v in zip(coefficients, vars))
        return [constant] + coefficients
    
    partitioners = list(p for p in get_partitioners(range(1, dimension + 1)))
    vars = list(SR("b{}".format(j)) for j in range(dimension+1))
    for this in partitioners:
        others = list(partitioners)
        others.remove(this)
        ineqs = [other.partition_cost() - this.partition_cost() for other in others] + vars
        ineq_matrix = [get_vector(ineq, vars) for ineq in ineqs]
        P = Polyhedron(ieqs=ineq_matrix)
        if make_disjoint:
            P = polyhedron_break_tie(P)
        this.polytope = P
        if verbose:
            print("*********************************")
            print("Partitioner {}".format(this))
            print(this.repr_pretty_polytope())

    dim = partitioners[0].polytope.ambient_dim()
    nonnegative_orthant = Polyhedron(ieqs=[dd*(0,) + (1,) + (dim-dd)*(0,)
                                           for dd in range(1, dim+1)])
    assert all(A.polytope & nonnegative_orthant == A.polytope
               for A in partitioners)
    if make_disjoint:
        assert all((A.polytope & B.polytope).is_empty()
                   for A in partitioners for B in partitioners if A != B)
    return partitioners
