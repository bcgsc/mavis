import glob
import os


def glob_exists(*pos, strict=True, n=1):
    globexpr = os.path.join(*pos)
    l = glob.glob(globexpr)
    if strict and len(l) == n:
        return l[0] if len(l) == 1 else l
    elif not strict and len(l) > 0:
        return l
    else:
        print(globexpr)
        print(l)
        return False
