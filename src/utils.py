import sys

import argparse
import configparser as cp



def print2(msg="", end='\n'):
    print(msg, file=sys.stderr, end=end)
    if end=='':
        sys.stderr.flush()

def withP(cur, all, add = None, fl = 0, fm="%.1f"):
    add = "" if (add is None) else (" " + add)
    return int2str(cur,fl=fl) + " (" + fm % (cur * 100.0 / all) + "%" + add + ")"







def int2str(v, fl = 0, delim="'"):
    s = str(abs(v))
    ans = ""
    while len(s) > 0:
        if len(ans) > 0:
            ans = delim + ans
        ans = s[-3:] + ans
        s = s[:-3]
    if v < 0:
        ans = "-" + ans
    return fl_func(ans, fl)
def fl(v, fl):
    ans = str(v)
    while len(ans) < fl:
        ans = " " + ans
    return ans
fl_func = fl



def str2int(s):
    if s[-1] == "k":
        value = int(s[:-1]) * 1000
    elif s[-1] == "M":
        value = int(s[:-1]) * 1000000
    elif s[-1] == "G":
        value = int(s[:-1]) * 1000000000
    else:
        value = int(s)
    return value

def compressInt(v, toAdd_format=" {suf}b"):
    if v < 1000:
        return str(v)
    # elif v < 100000:
    #     return str(v // 1000) + "k"
    else:
        return int2str(v // 1000, delim=",") + toAdd_format.format(suf="K")


def count_lines(file):
    lines = 0
    with open(file, "r") as f:
        for line in f:
            lines += 1
    return lines






def parseEdgeNumberWithDirection(str):
    forward = (str[-1] == '+')
    v = int(str[:-1])
    return v, forward

def parseEdgeNumber(name, ID="EDGE"):
    info = name.split('_')
    if info[0] != ID:
        raise ValueError("Unknown edge name's structure = '" + name + "'")
    return int(info[1])

def parseContigNumber(name):
    info = name.split('_')
    if info[0] not in ["contig", "NODE"]:
        raise ValueError("Unknown contig name's structure = '" + name + "'")
    return int(info[1])


def PERC(a, b):
    if b == 0:
        if a == 0:
            return "0.0%"
        else:
            return " ?? "
    else:
        return "%.1f%%" % (a * 100.0 / b)



class ExpandRanges(argparse.Action):
    """Translate a str like 1,2,3-5,40 to [1,2,3,4,5,40]"""
    def __call__(self, parser, namespace, values, option_string=None):
        import re
        elts = []
        for item in values.replace(' ', '').split(','):
            mo = re.search(r'(\d+)-(\d+)', item)
            if mo is not None:
                rng = [int(x) for x in mo.groups()]
                elts.extend(list(range(rng[0], rng[1] + 1)))
            else:
                elts.append(int(item))
        setattr(namespace, self.dest, elts)


def load_config(file, verbosePrint=True):
    if verbosePrint:
        print("Loading config... ", end='', flush=True)
    config = cp.ConfigParser()
    loaded = config.read(file)
    if verbosePrint:
        print("Loaded = " + str(loaded))
    assert loaded, "Config file not loaded!"

    return config

