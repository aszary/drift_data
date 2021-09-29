from astropy.io import ascii
from astropy.table import Table, vstack, setdiff
import pandas as pd
import numpy as np


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


def high_edot2(filename):
    """
    st0 = Table.read(filename, format='ascii', header_start=0, data_start=1) # read in the table
    a = list(st0.colnames)
    b = list(st0[0])
    for i in range(len(a)):
        print('"',a[i],'" :: ', b[i])
    return
    """
    st = Table.read(filename, format='ascii', header_start=0, data_start=2) # read in the table
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
    #print(ce.info)
    # create latex columns
    c1 = ce["C1 Power"]
    c1err = ce["C1 PowerErr"]
    rec1 = ["${:.1f} \\pm {:.1f}$".format(float(c1[i]), float(c1err[i])) for i in range(len(c1))]
    ce["C1 Power-Latex"] = rec1

    c1 = ce["MP C1 F1: P3_value"]
    c1err = ce["MP C1 F1: P3_error"]
    rec1 = ["${:.1f} \\pm {:.1f}$".format(float(c1[i]), float(c1err[i])) for i in range(len(c1))]
    ce["$P_3$"] = rec1
    ce["P_3"] = c1

    ce["C2 Power"][ce["C2 Power"]=="???"] = "$--$"  # changes value in table
    ce["C3 Power"][ce["C3 Power"]=="???"] = "$--$"  # changes value in table
    ce["C4 Power"][ce["C4 Power"]=="???"] = "$--$"  # changes value in table

    res = ce["PSRJ", "PSRB", r"$\dot E$", "P_3", "$P_3$", "SNR", "C1 Power", "C1 Power-Latex", "C2 Power", "C3 Power", "C4 Power"]
    res = res[res[r"$\dot E$"] > 5e32]
    return res


def drifting_p3only(filename="data/stats.csv"):
    st = Table.read(filename, format='ascii', header_start=0, data_start=2) # read in the table
    ce = st[st['Census']=='YES'] # select the census observation only
    # drifting only # checking all features for MP no IP included
    ces = []
    for i in range(0, 2):
        for j in range(0, 5):
            try:
                mask = ce["MP_C{}_F{}".format(i+1, j+1)] == "drift"
                ces.append(ce[mask])
            except ValueError:
                print("Error for ", "MP_C{}_F{}".format(i+1, j+1), "keyword...")
    drifting = vstack(ces)
    # P3only # checking all features for MP no IP included
    ces2 = []
    for i in range(0, 2):
        for j in range(0, 5):
            try:
                mask = ce["MP_C{}_F{}".format(i+1, j+1)] == "P3only"
                ces2.append(ce[mask])
            except ValueError:
                print("Error for ", "MP_C{}_F{}".format(i+1, j+1), "keyword...")
    p3only = vstack(ces2)

    # unique only
    for i in range(0, 2):
        for j in range(0, 5):
            try:
                mask1 = drifting["MP_C{}_F{}".format(i+1, j+1)] == "P3only"
                mask2 = p3only["MP_C{}_F{}".format(i+1, j+1)] == "drift"
                drifting = drifting[~mask1]
                p3only = p3only[~mask2]
            except ValueError:
                print("Error for ", "MP_C{}_F{}".format(i+1, j+1), "keyword...")

    drifting["Edot [ergs/s]"] = drifting["Edot [ergs/s]"].astype(float)
    p3only["Edot [ergs/s]"] = p3only["Edot [ergs/s]"].astype(float)

    return drifting, p3only


def positive_negative_mixed(filename="data/stats.csv"):
    st = Table.read(filename, format='ascii', header_start=0, data_start=2) # read in the table
    ce = st[st['Census']=='YES'] # select the census observation only

    # drifting only # checking all features for MP no IP included
    ces = []
    for i in range(0, 2):
        for j in range(0, 5):
            try:
                mask = ce["MP_C{}_F{}".format(i+1, j+1)] == "drift"
                ces.append(ce[mask])
            except ValueError:
                print("Error for ", "MP_C{}_F{}".format(i+1, j+1), "keyword...")
    dr = vstack(ces)

    # checking P2 / positive, negative, mixed
    pos = []
    neg = []
    for i in range(0, 2):
        for j in range(0, 5):
            #print("{} {} {}".format(i, j, dr["MP C{} F{}: P2_value".format(i+1, j+1)]))
            try:
                mask1 = ce["MP C{} F{}: P2_value".format(i+1, j+1)] > 0.
                pos.append(ce[mask1])
                mask2 = ce["MP C{} F{}: P2_value".format(i+1, j+1)] < 0.
                neg.append(ce[mask2])
            except ValueError:
                print("Error for ", "MP C{} F{}: P2_value".format(i+1, j+1), "keyword...")
    positive = vstack(pos)
    negative = vstack(neg)
    #print(len(positive))

    mix = []
    for i in range(0, 2):
        for j in range(0, 5):
            mask1 = negative["MP C{} F{}: P2_value".format(i+1, j+1)] > 0.
            mask2 = positive["MP C{} F{}: P2_value".format(i+1, j+1)] < 0.
            #print(positive[mask2])
            mix.append(positive[mask2])
            mix.append(negative[mask1])
            # unique records
            positive = positive[~mask2]
            negative = negative[~mask1]

    print(len(positive))
    print(len(negative))
    mixed = vstack(mix)

    positive["Edot [ergs/s]"] = positive["Edot [ergs/s]"].astype(float)
    negative["Edot [ergs/s]"] = negative["Edot [ergs/s]"].astype(float)
    mixed["Edot [ergs/s]"] = mixed["Edot [ergs/s]"].astype(float)

    #print(positive.info())
    #print(list(negative["JName_paper"]))

    return positive, negative, mixed


def positive_negative_mixed3(filename="data/stats.csv"):
    st = Table.read(filename, format='ascii', header_start=0, data_start=2) # read in the table
    ce = st[st['Census']=='YES'] # select the census observation only
    ce["Edot [ergs/s]"] = ce["Edot [ergs/s]"].astype(float)

    jnames_pos = []
    p3s_pos = []
    ep3s_pos = []
    edots_pos = []
    dps_pos = []

    jnames_neg = []
    p3s_neg = []
    ep3s_neg = []
    edots_neg = []
    dps_neg = []

    for row in ce:
        # get MP component numbers
        mcs = [int(x) for x in row["MP 2dfs nrs"].split(",")]
        for i,c in enumerate(mcs):
            try:
                # drift feature
                #print("ip1 ", i+1, " c ", c)
                f = int(row["MPdominantDriftFeature_C{}".format(i+1)]) # not c only two records
                # P3
                if row["MP_C{}_F{}".format(i+1, f)] == "drift":
                    p3 = float(row["MP C{} F{}: P3_value".format(i+1, f)])
                    p3error = float(row["MP C{} F{}: P3_error".format(i+1, f)])
                    p2 = float(row["MP C{} F{}: P2_value".format(i+1, f)])
                    edot = float(row["Edot [ergs/s]"])
                    jname = row["JName_paper"]
                    driftpower = float(row["C{} Power".format(c)]) # here c is fine..
                    if p2 > 0:
                        jnames_pos.append(jname)
                        p3s_pos.append(p3)
                        ep3s_pos.append(p3error)
                        edots_pos.append(edot)
                        dps_pos.append(driftpower)
                    elif p2 < 0:
                        jnames_neg.append(jname)
                        p3s_neg.append(p3)
                        ep3s_neg.append(p3error)
                        edots_neg.append(edot)
                        dps_neg.append(driftpower)
            except:
                pass
        # get IP components # adds only one record?!
        ics = []
        if row["IP 2dfs nrs"] != "???":
            ics = [int(x) for x in row["IP 2dfs nrs"].split(",")]
        for i,c in enumerate(ics):
            try:
                # drift feature
                f = int(row["IPdominantDriftFeature_C{}".format(i+1)])
                # P3
                if row["IP_C{}_F{}".format(i+1, f)] == "drift":
                    p3 = float(row["IP C{} F{}: P3_value".format(i+1, f)])
                    p3error = float(row["IP C{} F{}: P3_error".format(i+1, f)])
                    p2 = float(row["IP C{} F{}: P2_value".format(i+1, f)])
                    edot = float(row["Edot [ergs/s]"])
                    jname = row["JName_paper"]
                    driftpower = float(row["C{} Power".format(c)])
                    if p2 > 0:
                        jnames_pos.append(jname)
                        p3s_pos.append(p3)
                        ep3s_pos.append(p3error)
                        edots_pos.append(edot)
                        dps_pos.append(driftpower)
                    elif p2 < 0:
                        jnames_neg.append(jname)
                        p3s_neg.append(p3)
                        ep3s_neg.append(p3error)
                        edots_neg.append(edot)
                        dps_neg.append(driftpower)
            except np.ma.core.MaskError:
                pass

    jnames_mix = []
    p3s_mix = []
    ep3s_mix = []
    edots_mix = []
    dps_mix = []

    mixed = []
    for name in jnames_pos:
        if name in jnames_neg:
           mixed.append(name)

    # wow so lame...
    for name in mixed:
        i = jnames_pos.index(name)
        jnames_mix.append(name)
        p3s_mix.append(p3s_pos[i])
        ep3s_mix.append(ep3s_pos[i])
        edots_mix.append(edots_pos[i])
        dps_mix.append(dps_pos[i])
        jnames_pos.pop(i)
        p3s_pos.pop(i)
        ep3s_pos.pop(i)
        edots_pos.pop(i)
        dps_pos.pop(i)
        j = jnames_neg.index(name)
        jnames_mix.append(name)
        p3s_mix.append(p3s_neg[j])
        ep3s_mix.append(ep3s_neg[j])
        edots_mix.append(edots_neg[j])
        dps_mix.append(dps_neg[j])
        jnames_neg.pop(j)
        p3s_neg.pop(j)
        ep3s_neg.pop(j)
        edots_neg.pop(j)
        dps_neg.pop(j)

    print("Positive: ", len(jnames_pos))
    print("Negative: ", len(jnames_neg))
    print("Mixed: ", len(jnames_mix))

    return [jnames_pos, p3s_pos, ep3s_pos, edots_pos, dps_pos], [jnames_neg, p3s_neg, ep3s_neg, edots_neg, dps_neg], [jnames_mix, p3s_mix, ep3s_mix, edots_mix, dps_mix]



def drifting_p3only2(filename="data/stats.csv"):
    st = Table.read(filename, format='ascii', header_start=0, data_start=2) # read in the table
    ce = st[st['Census']=='YES'] # select the census observation only
    #ce["Edot [ergs/s]"] = ce["Edot [ergs/s]"].astype(float)

    jnames_dr = []
    p3s_dr = []
    ep3s_dr = []
    ps_dr = []
    dps_dr = []
    pdots_dr = []

    jnames_p3o = []
    p3s_p3o = []
    ep3s_p3o = []
    ps_p3o = []
    dps_p3o = []
    pdots_p3o = []

    for row in ce:
        # get MP component numbers
        mcs = [int(x) for x in row["MP 2dfs nrs"].split(",")]
        for i,c in enumerate(mcs):
            try:
                fd = int(row["MPdominantDriftFeature_C{}".format(i+1)]) # not c only two records
            except:
                fd = None
            try:
                f3 = int(row["MPdominantP3Feature_C{}".format(i+1)]) # not c only two records
            except:
                f3 = None
            # drift feature
            p = float(row["Period [s]"])
            pdot = float(row["Pdot [s/s]"])
            jname = row["JName_paper"]
            driftpower = float(row["C{} Power".format(c)]) # here c is fine..
            if (fd is not None) and (row["MP_C{}_F{}".format(i+1, fd)] == "drift"):
                p3 = float(row["MP C{} F{}: P3_value".format(i+1, fd)])
                p3error = float(row["MP C{} F{}: P3_error".format(i+1, fd)])
                jnames_dr.append(jname)
                p3s_dr.append(p3)
                ep3s_dr.append(p3error)
                ps_dr.append(p)
                dps_dr.append(driftpower)
                pdots_dr.append(pdot)
            elif (f3 is not None) and (row["MP_C{}_F{}".format(i+1, f3)] == "P3only"):
                p3 = float(row["MP C{} F{}: P3_value".format(i+1, f3)])
                p3error = float(row["MP C{} F{}: P3_error".format(i+1, f3)])
                jnames_p3o.append(jname)
                p3s_p3o.append(p3)
                ep3s_p3o.append(p3error)
                ps_p3o.append(p)
                dps_p3o.append(driftpower)
                pdots_p3o.append(pdot)
            #except np.ma.core.MaskError:
            #    pass
        # get IP components #
        ics = []
        if row["IP 2dfs nrs"] != "???":
            ics = [int(x) for x in row["IP 2dfs nrs"].split(",")]
        for i,c in enumerate(ics):
            try:
                fd = int(row["IPdominantDriftFeature_C{}".format(i+1)]) # not c only two records
            except:
                fd = None
            try:
                f3 = int(row["IPdominantP3Feature_C{}".format(i+1)]) # not c only two records
            except:
                f3 = None
            # drift feature
            p = float(row["Period [s]"])
            pdot = float(row["Pdot [s/s]"])
            jname = row["JName_paper"]
            driftpower = float(row["C{} Power".format(c)]) # here c is fine..
            if (fd is not None) and (row["MP_C{}_F{}".format(i+1, fd)] == "drift"):
                p3 = float(row["IP C{} F{}: P3_value".format(i+1, fd)])
                p3error = float(row["IP C{} F{}: P3_error".format(i+1, fd)])
                jnames_dr.append(jname)
                p3s_dr.append(p3)
                ep3s_dr.append(p3error)
                ps_dr.append(p)
                dps_dr.append(driftpower)
                pdots_dr.append(pdot)
            elif (f3 is not None) and (row["IP_C{}_F{}".format(i+1, f3)] == "P3only"):
                p3 = float(row["IP C{} F{}: P3_value".format(i+1, f3)])
                p3error = float(row["IP C{} F{}: P3_error".format(i+1, f3)])
                jnames_p3o.append(jname)
                p3s_p3o.append(p3)
                ep3s_p3o.append(p3error)
                ps_p3o.append(p)
                dps_p3o.append(driftpower)
                pdots_p3o.append(pdot)
    print("Drift: ", len(jnames_dr))
    print("P3only: ", len(jnames_p3o))

    return [jnames_dr, p3s_dr, ep3s_dr, ps_dr, dps_dr, pdots_dr], [jnames_p3o, p3s_p3o, ep3s_p3o, ps_p3o, dps_p3o, pdots_p3o]
