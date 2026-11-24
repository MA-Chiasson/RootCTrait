"""Natural order of sample names: S2 before S10, and B2T3G1S2 before B2T3G1S10.
Every number in the name is compared as a number, from left to right."""
import re


def natural_key(name):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(name))]
