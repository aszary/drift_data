from astropy.io import ascii
from astropy.table import Table, vstack
import pandas as pd


def latexify(table):
    return ascii.write(table, format='latex')


def read_p0_edot(filename):
    df = pd.read_csv(filename, sep=";")
    p0s = df["P0"][1:].to_list() # strings
    p0s = [float(p0) for p0 in p0s] # floats
    edots = df["EDOT"][1:].to_list() # strings
    edots = [float(ed) for ed in edots] # floats
    return [p0s, edots]


def read_highedot(filename):
    """ using pandas - obsolete? """
    df = pd.read_csv(filename, sep=",", low_memory=False) #  have mixed types.Specify dtype option on import or set low_memory=False.
    # prints columns information
    """
    descriptions = df.iloc[0]
    ##pd.set_option('display.max_columns', None)
    #pd.set_option('display.max_colwidth', -1) # uncomment for full width
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', None)
    #print(descriptions)
    """

    # pulsars in census
    ce = df.loc[df["Census"]=="YES"]

    # high edot only
    ce["Edot [ergs/s]"] = pd.to_numeric(ce["Edot [ergs/s]"]) # convert Edot to numbers
    ce = ce.loc[ce["Edot [ergs/s]"] > 5e32]

    # drifting only # checking all components
    ces = []
    cols = []
    for i in range(1,3):
        for j in range(1, 6):
            ces.append(ce.loc[ce["MP_C{}_F{}".format(i, j)] == "drift"])
            cols.append("MP_C{}_F{}".format(i, j))
    # interpulses
    for i in range(1,3):
        for j in range(1, 6):
            ces.append(ce.loc[ce["IP_C{}_F{}".format(i, j)] == "drift"])
            cols.append("MP_C{}_F{}".format(i, j))

    ce = pd.concat(ces)
    ce.rename(columns={"Edot [ergs/s]":"dot E", "JName_paper":"PSRJ", "BName":"PSRB"}, inplace=True)
    res = ce.loc[:, ["PSRJ", "PSRB", "dot E"]]
    print(res.to_latex(index=False)) # prints latex table on screen
    return res


def read_highedot2(filename):
    st = Table.read(filename, format='ascii', header_start=0, data_start=2) #read in the table
    ce = st[st['Census']=='YES'] # select the census observation only
    ce.rename_column("Edot [ergs/s]", r"$\dot E$")
    ce.rename_column("JName_paper", "PSRJ")
    ce.rename_column("BName", "PSRB")
    ce[r"$\dot E$"] = ce[r"$\dot E$"].astype(float)
    # drifting only # checking all features for MP
    ces = []
    cols = []
    for i in range(1,3):
        for j in range(1, 6):
            try:
                mask = ce["MP_C{}_F{}".format(i, j)] == "drift"
                ces.append(ce[mask])
                cols.append("MP_C{}_F{}".format(i, j))
            except ValueError:
                print("Error for ", "MP_C{}_F{}".format(i, j), "keyword")
    ce = vstack(ces)
    res = ce["PSRJ", "PSRB", r"$\dot E$"]
    res = res[res[r"$\dot E$"] > 5e32]
    return res
