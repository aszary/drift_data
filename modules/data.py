from astropy.table import Table
import pandas as pa


def read_p0_edot(filename):
    df = pa.read_csv(filename, sep=";")
    p0s = df["P0"][1:].to_list() # strings
    p0s = [float(p0) for p0 in p0s] # floats
    edots = df["EDOT"][1:].to_list() # strings
    edots = [float(ed) for ed in edots] # floats
    return [p0s, edots]
