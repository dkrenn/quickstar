from partitioner import ClassificationStrategy
from quickstar import QuickStar

cs = ClassificationStrategy(3)
Q = [(n, QuickStar().avg_comparisons_quicksort(n,
                                               'partitioned_trees',
                                               strategy=cs))
     for n in srange(9)]

print(', '.join(str(q * n.factorial()) for n, q in Q))
