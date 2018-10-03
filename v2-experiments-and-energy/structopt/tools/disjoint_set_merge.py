def disjoint_set_merge(elements, equivalent_pairs):
    # See https://www.topcoder.com/community/data-science/data-science-tutorials/disjoint-set-data-structures/
    def _merge_sets(sets, set1, set2):
        if set1 is not set2:
            for elem in set2:
                set1.add(elem)
            sets.remove(set2)

    def _find_set(sets, elem):
        for s in sets:
            if elem in s:
                return s

    sets = []
    for x in elements:
        sets.append({x})
    for x, y in equivalent_pairs:
        setx = _find_set(sets, x)
        sety = _find_set(sets, y)
        if setx != sety:
            _merge_sets(sets, setx, sety)
    sets.sort(key=lambda s: min(s))
    return sets

