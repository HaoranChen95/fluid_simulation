import re


class file_name:
    def __init__(self, fn: "str") -> None:
        self.Nm = None
        self.kT = None
        self.phi = None
        self.gamma = None
        self.aV = None
        self.str = fn
        self.__set_file_param(fn)
        self.pattern = self.__pattern()

    def __str__(self) -> str:
        return self.str

    def __set_file_param(self, fn: "str"):
        print("[setting] global parameter from file name")
        print(fn)
        if fn.find("aV") != -1:
            self.aV = find_param_int("aV", fn)
        if fn.find("phi") != -1:
            self.phi = find_param_float("phi", fn)
        if fn.find("Nm") != -1:
            self.Nm = find_param_int("Nm", fn)
        if fn.find("kT") != -1:
            self.kT = find_param_float("kT", fn)
        if fn.find("gam") != -1:
            self.gamma = find_param_float("gam", fn)

    def __pattern(self) -> "str":
        return f"phi_{self.phi}_Nm_{self.Nm}_kT_{self.kT}_gam_{self.gamma}"


def find_param_int(ps, s):
    """find the int parameter in the string s after the pattern ps_

    Args:
        ps (string): pattern for search
        s (string): string to analyse

    Returns:
        int: the first value after pattern
    """
    match = re.search(ps + "_(-?\d+)", s)
    return int(match.group(1))


def find_param_float(ps, s):
    """find the float parameter in the string s after the pattern ps_

    Args:
        ps (string): pattern for search
        s (string): string to analyse

    Returns:
        float: the first value after pattern
    """
    match = re.search(ps + "_(-?\d+.?\d*[e]?-?\d+)", s)
    res = 0
    try:
        res = float(match.group(1))
    except AttributeError:
        res = find_param_int(ps, s)
    return float(res)
