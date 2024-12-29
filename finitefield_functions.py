from configparser import ConfigParser
from math import sqrt

config = ConfigParser()
config.read("config.cfg")

IRREDUCIBLE_POLYS_PATH = str(config.get("General", "IRREDUCIBLE_POLYS_PATH"))
ALWAYS_REDUCE = (config.get("General", "ALWAYS_REDUCE") == "True")

def isPrime(n: int) -> bool:
    if n == 2:
        return True
    if (not isinstance(n, int)) or n<0 or n%2 == 0:
        return False
    for i in range(3, int(sqrt(n))+1, 2):
        if n%i == 0:
            return False
    return True

def largest_index(lst: list):
    for i in range(len(lst)-1, 0, -1):
        if lst[i] != 0:
            return i
    return 0