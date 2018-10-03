from structopt.tools import disjoint_set_merge

elements = [1,2,3,4,5]
equivalent_pairs = [(1,2), (2,3), (4,5)]
super_sets = disjoint_set_merge(elements, equivalent_pairs)
assert {1,2,3} in super_sets and {4,5} in super_sets and len(super_sets) == 2

equivalent_pairs = [(1,2), (4,5), (2,5)]
super_sets = disjoint_set_merge(elements, equivalent_pairs)
assert {1,2,4,5} in super_sets and {3} in super_sets and len(super_sets) == 2

elements = [1,2,3,4,5,6,7,8,9,0]
equivalent_pairs = [(1,3), (5,6), (8,0), (0, 8), (8, 0), (1, 3), (4, 5), (5,3)]
super_sets = disjoint_set_merge(elements, equivalent_pairs)
assert super_sets == [{8, 0}, {1, 3, 4, 5, 6}, {2}, {7}, {9}]
