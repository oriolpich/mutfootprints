
import collections
import functools


def leveled_default_dict(inner, level=1):
    level = level - 1
    if level == 0:
        return collections.defaultdict(inner)
    else:
        f = functools.partial(leveled_default_dict, inner=inner, level=level)
        return collections.defaultdict(f)
