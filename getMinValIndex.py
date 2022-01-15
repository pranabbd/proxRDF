def minValIndex(l):
    minVal = min(l)
    if l.count(minVal) > 1:
        return [i for i, x in enumerate(l) if x == min(l)]
    else:
        return l.index(min(l))
