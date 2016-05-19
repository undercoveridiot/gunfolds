import sys

sys.path.append('./tools/')


def osumnum(s, num):
    return set(num + x for x in s)


def osumset(s1, s2):
    s = set()
    for e in s1:
        s.update(osumnum(s2, e))
    return s


class PathTree:
    def __init__(self, lset, pre=None):
        self.loopset = lset
        if pre is None:
            self.preset = {0}
        else:
            assert (type(pre) == set)
            self.preset = pre

    def __add__(self, other):
        if type(other) is int:
            return PathTree(self.loopset, pre=osumnum(self.preset, other))
        if type(other) is set:
            return PathTree(self.loopset, pre=osumset(self.preset, other))
        return PathTree(self.loopset.union(other.loopset), pre=osumset(self.preset, other.preset))

    def __radd__(self, other):
        return self.__add__(other)

    def __str__(self):
        s = '('
        comma = False
        if not self.preset == {0}:
            s += '{' + ', '.join(map(str, self.preset)) + '} '
            comma = True
        if not self.loopset == {0}:
            if comma:
                s += ', '
            s += '<{ ' + ', '.join(map(str, self.loopset)) + ' }>'
        s += ')'
        return s
