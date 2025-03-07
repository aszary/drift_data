from astropy.io import ascii
from astropy.table import Table, vstack, setdiff, Column, join
from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
import pandas as pd
import numpy as np
from uncertainties import ufloat
import matplotlib.pyplot as pl
from scipy.stats import ttest_ind
#import psrqpy

from modules.functions import *
import modules.plot as plo


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
    ps_pos = []

    jnames_neg = []
    p3s_neg = []
    ep3s_neg = []
    edots_neg = []
    dps_neg = []
    ps_neg = []

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
                    p = float(row["Period [s]"])
                    edot = float(row["Edot [ergs/s]"])
                    jname = row["JName_paper"]
                    driftpower = float(row["C{} Power".format(c)]) # here c is fine..
                    if p2 > 0:
                        jnames_pos.append(jname)
                        p3s_pos.append(p3)
                        ep3s_pos.append(p3error)
                        edots_pos.append(edot)
                        dps_pos.append(driftpower)
                        ps_pos.append(p)
                    elif p2 < 0:
                        jnames_neg.append(jname)
                        p3s_neg.append(p3)
                        ep3s_neg.append(p3error)
                        edots_neg.append(edot)
                        dps_neg.append(driftpower)
                        ps_neg.append(p)
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
                    p = float(row["Period [s]"])
                    jname = row["JName_paper"]
                    driftpower = float(row["C{} Power".format(c)])
                    if p2 > 0:
                        jnames_pos.append(jname)
                        p3s_pos.append(p3)
                        ep3s_pos.append(p3error)
                        edots_pos.append(edot)
                        dps_pos.append(driftpower)
                        ps_pos.append(p)
                    elif p2 < 0:
                        jnames_neg.append(jname)
                        p3s_neg.append(p3)
                        ep3s_neg.append(p3error)
                        edots_neg.append(edot)
                        dps_neg.append(driftpower)
                        ps_neg.append(p)
            except np.ma.core.MaskError:
                pass

    jnames_mix = []
    p3s_mix = []
    ep3s_mix = []
    edots_mix = []
    dps_mix = []
    ps_mix = []

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
        ps_mix.append(ps_pos[i])
        jnames_pos.pop(i)
        p3s_pos.pop(i)
        ep3s_pos.pop(i)
        edots_pos.pop(i)
        dps_pos.pop(i)
        ps_pos.pop(i)
        j = jnames_neg.index(name)
        jnames_mix.append(name)
        p3s_mix.append(p3s_neg[j])
        ep3s_mix.append(ep3s_neg[j])
        edots_mix.append(edots_neg[j])
        dps_mix.append(dps_neg[j])
        ps_mix.append(ps_neg[j])
        jnames_neg.pop(j)
        p3s_neg.pop(j)
        ep3s_neg.pop(j)
        edots_neg.pop(j)
        dps_neg.pop(j)
        ps_neg.pop(j)

    print("Positive: ", len(jnames_pos))
    print("Negative: ", len(jnames_neg))
    print("Mixed: ", len(jnames_mix))

    return [jnames_pos, p3s_pos, ep3s_pos, edots_pos, dps_pos, ps_pos], [jnames_neg, p3s_neg, ep3s_neg, edots_neg, dps_neg, ps_neg], [jnames_mix, p3s_mix, ep3s_mix, edots_mix, dps_mix, ps_mix]


def positive_negative_mixed4(filename="data/stats.csv", edot_max=2e32):
    edot_max = 0 # take all!
    st = Table.read(filename, format='ascii', header_start=0, data_start=2) # read in the table
    ce = st[st['Census']=='YES'] # select the census observation only
    ce["Edot [ergs/s]"] = ce["Edot [ergs/s]"].astype(float)

    jnames_pos = []
    p3s_pos = []
    ep3s_pos = []
    edots_pos = []
    dps_pos = []
    ps_pos = []

    jnames_neg = []
    p3s_neg = []
    ep3s_neg = []
    edots_neg = []
    dps_neg = []
    ps_neg = []

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
                    p = float(row["Period [s]"])
                    edot = float(row["Edot [ergs/s]"])
                    jname = row["JName_paper"]
                    driftpower = float(row["C{} Power".format(c)]) # here c is fine..
                    if p2 > 0:
                        if edot >= edot_max:
                            jnames_pos.append(jname)
                            p3s_pos.append(p3)
                            ep3s_pos.append(p3error)
                            edots_pos.append(edot)
                            dps_pos.append(driftpower)
                            ps_pos.append(p)
                    elif p2 < 0:
                        if edot >= edot_max:
                            jnames_neg.append(jname)
                            p3s_neg.append(p3)
                            ep3s_neg.append(p3error)
                            edots_neg.append(edot)
                            dps_neg.append(driftpower)
                            ps_neg.append(p)
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
                    p = float(row["Period [s]"])
                    jname = row["JName_paper"]
                    driftpower = float(row["C{} Power".format(c)])
                    if p2 > 0:
                        if edot >= edot_max:
                            jnames_pos.append(jname)
                            p3s_pos.append(p3)
                            ep3s_pos.append(p3error)
                            edots_pos.append(edot)
                            dps_pos.append(driftpower)
                            ps_pos.append(p)
                    elif p2 < 0:
                        if edot >= edot_max:
                            jnames_neg.append(jname)
                            p3s_neg.append(p3)
                            ep3s_neg.append(p3error)
                            edots_neg.append(edot)
                            dps_neg.append(driftpower)
                            ps_neg.append(p)
            except np.ma.core.MaskError:
                pass

    jnames_mix = []
    p3s_mix = []
    ep3s_mix = []
    edots_mix = []
    dps_mix = []
    ps_mix = []

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
        ps_mix.append(ps_pos[i])
        jnames_pos.pop(i)
        p3s_pos.pop(i)
        ep3s_pos.pop(i)
        edots_pos.pop(i)
        dps_pos.pop(i)
        ps_pos.pop(i)
        j = jnames_neg.index(name)
        jnames_mix.append(name)
        p3s_mix.append(p3s_neg[j])
        ep3s_mix.append(ep3s_neg[j])
        edots_mix.append(edots_neg[j])
        dps_mix.append(dps_neg[j])
        ps_mix.append(ps_neg[j])
        jnames_neg.pop(j)
        p3s_neg.pop(j)
        ep3s_neg.pop(j)
        edots_neg.pop(j)
        dps_neg.pop(j)
        ps_neg.pop(j)

    print("Positive: ", len(jnames_pos))
    print("Negative: ", len(jnames_neg))
    print("Mixed: ", len(jnames_mix))

    return [jnames_pos, p3s_pos, ep3s_pos, edots_pos, dps_pos, ps_pos], [jnames_neg, p3s_neg, ep3s_neg, edots_neg, dps_neg, ps_neg], [jnames_mix, p3s_mix, ep3s_mix, edots_mix, dps_mix, ps_mix]


def positive_negative_manual(filename="data/stats.csv", posf="data/positive.txt", negf="data/negative.txt"):

    st = Table.read(filename, format='ascii', header_start=0, data_start=2) # read in the table
    ce = st[st['Census']=='YES'] # select the census observation only
    ce["Edot [ergs/s]"] = ce["Edot [ergs/s]"].astype(float)

    # positive drift
    jnames_pos = []
    p3s_pos = []
    ep3s_pos = []
    edots_pos = []

    for line in open(posf).readlines():
        jname = line.split()[0]
        rec = ce[ce['JName']==jname]
        fe = rec['PulsarDominantP3Feature'].value[0].split()
        #print(fe)
        feature = "{}dominantDriftFeature_{}".format(fe[0], fe[1])
        f = int(rec[feature].value[0])
        p3 =  float(rec["{} {} F{}: P3_value".format(fe[0], fe[1], f)])
        ep3 =  float(rec["{} {} F{}: P3_error".format(fe[0], fe[1], f)])
        edot = float(rec["Edot [ergs/s]"])
        jnames_pos.append(jname)
        p3s_pos.append(p3)
        ep3s_pos.append(ep3)
        edots_pos.append(edot)
    

    # negative drift
    jnames_neg = []
    p3s_neg = []
    ep3s_neg = []
    edots_neg = []



    for line in open(negf).readlines():
        jname = line.split()[0]
        rec = ce[ce['JName']==jname]
        fe = rec['PulsarDominantP3Feature'].value[0].split()
        #print(fe)
        feature = "{}dominantDriftFeature_{}".format(fe[0], fe[1])
        f = int(rec[feature].value[0])
        p3 =  float(rec["{} {} F{}: P3_value".format(fe[0], fe[1], f)])
        ep3 =  float(rec["{} {} F{}: P3_error".format(fe[0], fe[1], f)])
        edot = float(rec["Edot [ergs/s]"])
        jnames_neg.append(jname)
        p3s_neg.append(p3)
        ep3s_neg.append(ep3)
        edots_neg.append(edot)
 
    return [jnames_pos, p3s_pos, ep3s_pos, edots_pos], [jnames_neg, p3s_neg, ep3s_neg, edots_neg]


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


def allfeaturesp3_xiaoxi(filename="data/allfeaturesp3.csv"):
    st = Table.read(filename, format='ascii', header_start=0, data_start=1) # read in the table
    #print(st.info)
    jnames = list(st["JName"])
    p3s = list(st["P3"])
    ep3s = list(st["P3_error"])
    ps = list(st["Period [s]"])
    edots = list(st["Edot [ergs/s]"])
    pdots = list(st["Pdot [s/s]"])

    return [jnames, p3s, ep3s, ps, pdots, edots]


def extractdomp3s(filename="data/stats.csv"):
    """ Written by Xiaoxi modified by me """
    """
    # HEADER header
    st0 = Table.read(filename, format='ascii', header_start=0, data_start=1) # read in the table
    a = list(st0.colnames)
    b = list(st0[0])
    for i in range(len(a)):
        print('"',a[i],'" :: ', b[i])
    return
    """
    statstab = Table.read(filename, format='ascii', header_start=0, data_start=0)
    #input: stats.csv, census observations and the ones with features only
    statsndomfeature = statstab[(statstab['Census']=='YES')&(statstab['PulsarDominantP3Feature'].mask==False)]
    statsndomfeature['DomFeature'] = 'XXXXXXXX'
    for i in range(len(statsndomfeature)):
        compn = statsndomfeature['PulsarDominantP3Feature'][i]
        mip,comps = compn.split(' ')
        dname = mip+'dominantDriftFeature_' + comps
        pname = mip+'dominantP3Feature_' + comps
        # if DriftFeature and P3feature both defined, choose the P3 feature; if only DriftFeature or P3feature, take this.
        if (float(statsndomfeature[dname][i]) > 0):
            if (float(statsndomfeature[pname][i]) > 0):
                fn = 'F'+str(statsndomfeature[pname][i])
            else:
                fn = 'F'+str(statsndomfeature[dname][i])
            statsndomfeature['DomFeature'][i] = mip+' '+comps+' '+fn
        elif (float(statsndomfeature[pname][i]) > 0):
            fn = 'F'+str(statsndomfeature[pname][i])
            statsndomfeature['DomFeature'][i] = mip+' '+comps+' '+fn
    p3s = np.zeros(len(statsndomfeature))
    p3serr = np.zeros(len(statsndomfeature))
    for j in range(len(statsndomfeature)):
        p3s[j] = statsndomfeature[statsndomfeature['DomFeature'][j]+': P3_value'][j]
        #print(statsndomfeature['DomFeature'][j])
        p3serr[j] = statsndomfeature[statsndomfeature['DomFeature'][j]+': P3_error'][j]
    edots = statsndomfeature['Edot [ergs/s]'].astype(float) #can extract p, pdot, age etc as necessary
    jnames = statsndomfeature["JName_paper"]
    return p3s, p3serr, edots, jnames, statsndomfeature #return the table too



def p3plusp2(filename="data/stats.csv"):
    p3s = []
    p3serr = []
    p3ws = []
    p3wserr = []
    p2s = []
    p2serr_p = []
    p2serr_m = []
    edots = []
    jnames = []
    p2as = []
    p2aserr = []


    st = Table.read(filename, format='ascii', header_start=0, data_start=0)

    # get P2only pulsars first - no PulsarDominantP3Feature but P2only components (2 pulsars?)
    st0 = st[(st['Census']=='YES') & (st["IgnoreBecauseOfScattering"]!="ignore") & (st['PulsarDominantP3Feature'].mask==True)]
    # checking all components
    for i in range(1,3):
        for j in range(1, 6):
            fename = "MP_C{}_F{}".format(i, j)
            st0 = st0[st0[fename] == "P2only"]
            for row in st0:
                feature = fename.replace("_", " ")
                p3s.append("--")
                p3serr.append("--")
                # dominant (and only one probably)
                p2s.append(float(row[feature + ': P2_value']))
                p2serr_p.append(float(row[feature + ': P2_error_plus']))
                p2serr_m.append(float(row[feature + ': P2_error_minus']))
                p3ws.append("--")
                p3wserr.append("--")
                edots.append(row['Edot [ergs/s]'])
                jnames.append(row["JName_paper"])
                # P2 assymetry
                c = feature[4]
                # the dominant feature
                p2as.append(float(row["C{} Power".format(c)]))
                p2aserr.append(float(row["C{} PowerErr".format(c)]))
    # interpulses # not implemented - not really needed
    for i in range(1,3):
        for j in range(1, 6):
            fename = "IP_C{}_F{}".format(i, j)

    # get rest of the pulsars
    st = st[(st['Census']=='YES')&(st['PulsarDominantP3Feature'].mask==False) & (st["IgnoreBecauseOfScattering"]!="ignore")]
    st["Edot [ergs/s]"] = st["Edot [ergs/s]"].astype(float)
    for row in st:
        compn = row['PulsarDominantP3Feature']
        mip, comp = compn.split()
        dname = mip+'dominantDriftFeature_' + comp
        pname = mip+'dominantP3Feature_' + comp
        # if DriftFeature and P3feature both defined, choose the drift feature; if only DriftFeature or P3feature, take that one.
        # kind of strange, but works
        if float(row[dname]) > 0 and float(row[pname]) > 0:
            fn = 'F'+str(row[dname])
        elif float(row[dname]) > 0:
            fn = 'F'+str(row[dname])
        elif float(row[pname]) > 0:
            fn = 'F'+str(row[pname])
        feature = mip + ' ' + comp + ' ' + fn
        p3s.append(float(row[feature + ': P3_value']))
        p3serr.append(float(row[feature + ': P3_error']))
        # P2 the lowest value
        """
        if row[feature.replace(" ", "_")] == "P3only":
            p2s.append("--")
            p2serr_p.append("--")
            p2serr_m.append("--")
        else:
            p2, p2ep, p2em = get_smallest_p2(row)
            p2s.append(p2)
            p2serr_p.append(p2ep)
            p2serr_m.append(p2em)
        """
        # p2 for the dominant drift or p3 feature (not the lowest value?!)
        if row[feature.replace(" ", "_")] == "P3only": # removes very P2 values! ASK ABOUT IT!!!
            p2s.append("--")
            p2serr_p.append("--")
            p2serr_m.append("--")
        else:
            p2s.append(float(row[feature + ': P2_value']))
            p2serr_p.append(float(row[feature + ': P2_error_plus']))
            p2serr_m.append(float(row[feature + ': P2_error_minus']))
        p3ws.append(float(row[feature + ': P3_stddev_value']))
        p3wserr.append(float(row[feature + ': P3_stddev_error']))
        edots.append(row['Edot [ergs/s]'])
        jnames.append(row["JName_paper"])
        # P2 assymetry
        comps = row["{} 2dfs nrs".format(mip)].split(",") # get component numbers
        p2a_all = []
        p2ae_all = []
        # TODO change to the dominant feature
        # maximum value of P2assymetry
        """
        for c in comps:
            p2a_all.append(float(row["C{} Power".format(c)]))
            p2ae_all.append(float(row["C{} PowerErr".format(c)]))
        ind = np.argmax(p2a_all)
        p2as.append(p2a_all[ind])
        p2aserr.append(p2ae_all[ind])
        """
        # dominant component value
        p2as.append(float(row["{} Power".format(comp)]))
        p2aserr.append(float(row["{} PowerErr".format(comp)]))
        """
        for ii in range(4):
            print(mip, comp)
            print(row["MP 2dfs nrs"])
            print(ii, row["JName_paper"], row["C{} phase left".format(ii+1)])
        """
    #print(st)
    #print(len(p3s), len(p3serr), len(p3ws), len(p3wserr), len(p2s), len(p2serr_p), len(p2serr_m), len(p2as), len(p2aserr), len(edots), len(jnames))
    #exit()
    ind = np.argsort(jnames)
    p3s = np.array(p3s)[ind]
    p3serr = np.array(p3serr)[ind]
    p3ws = np.array(p3ws)[ind]
    p3wserr = np.array(p3wserr)[ind]
    p2s = np.array(p2s)[ind]
    p2serr_p = np.array(p2serr_p)[ind]
    p2serr_m = np.array(p2serr_m)[ind]
    p2as = np.array(p2as)[ind]
    p2aserr = np.array(p2aserr)[ind]
    edots = np.array(edots)[ind]
    jnames = np.array(jnames)[ind]
    return p3s, p3serr, p3ws, p3wserr, p2s, p2serr_p, p2serr_m, p2as, p2aserr, edots, jnames


def p3plusp2_full(filename="data/stats.csv"):

    jnames = []
    ps = []
    pdots = []
    edots = []
    bs = []
    blcs = []
    ages = []
    s1400s = []
    w10s = []
    w50s = []
    ip = []

    obsnames = []
    nants = []
    npulses = []
    snrs = []
    fftlengths = []
    avmods = []
    avmodeserr = []
    #avmods_ip = []
    #avmodserr_ip = []
    mp_data = []
    ip_data = []

    # P2only pulsars should be included
    st = Table.read(filename, format='ascii', header_start=0, data_start=0)
    st = st[(st['Census']=='YES')& (st["IgnoreBecauseOfScattering"]!="ignore")]
    st["Edot [ergs/s]"] = st["Edot [ergs/s]"].astype(float)
    for row in st:
        data1 = []
        comps = row["MP 2dfs nrs"].split(",")
        #mask = ce["MP_C{}_F{}".format(i, j)] == "drift"
        for c in comps: # using 2dfs numbers as component numbers
            df = row["MPdominantDriftFeature_C{}".format(c)]
            try:
                df = int(df)
            except:
                df = 0
            p3f = row["MPdominantP3Feature_C{}".format(c)]
            try:
                p3f = int(p3f)
            except:
                p3f = 0
            for f in range(1,6):
                if row["MP_C{}_F{}".format(c, f)] == "P3only":
                    p3 = row["MP C{} F{}: P3_value".format(c, f)]
                    p3err = row["MP C{} F{}: P3_error".format(c, f)]
                    type_ = "P3only"
                    if f == p3f:
                        data1.insert(0, [type_, p3, p3err, None, None, None])
                    else:
                        data1.append([type_, p3, p3err, None, None, None])
                elif row["MP_C{}_F{}".format(c, f)] == "drift":
                    p3 = row["MP C{} F{}: P3_value".format(c, f)]
                    p3err = row["MP C{} F{}: P3_error".format(c, f)]
                    p2 = row["MP C{} F{}: P2_value".format(c, f)]
                    p2pl = row["MP C{} F{}: P2_error_plus".format(c, f)]
                    p2mi = row["MP C{} F{}: P2_error_minus".format(c, f)]
                    type_ = "drift"
                    if f == df:
                        data1.insert(0, [type_, p3, p3err, p2, p2pl, p2mi])
                    else:
                        data1.append([type_, p3, p3err, p2, p2pl, p2mi])
                elif row["MP_C{}_F{}".format(c, f)] == "P2only":
                    p2 = row["MP C{} F{}: P2_value".format(c, f)]
                    p2pl = row["MP C{} F{}: P2_error_plus".format(c, f)]
                    p2mi = row["MP C{} F{}: P2_error_minus".format(c, f)]
                    type_ = "P2only"
                    data1.append([type_, None, None, p2, p2pl, p2mi])
        data2 = []
        comps = row["IP 2dfs nrs"].split(",")
        #mask = ce["MP_C{}_F{}".format(i, j)] == "drift"
        for c, cc in enumerate(comps): # using index of 2dfs numbers as components numbers
            try:
                cc = int(cc)
            except:
                break
            c += 1
            df = row["IPdominantDriftFeature_C{}".format(c)]
            p3f = row["IPdominantP3Feature_C{}".format(c)]
            try:
                df = int(df)
            except:
                df = 0
            try:
                p3f = int(p3f)
            except:
                p3f = 0
            for f in range(1,6):
                if row["IP_C{}_F{}".format(c, f)] == "P3only":
                    p3 = row["IP C{} F{}: P3_value".format(c, f)]
                    p3err = row["IP C{} F{}: P3_error".format(c, f)]
                    type_ = "P3only"
                    if f == p3f:
                        data2.insert(0, [type_, p3, p3err, None, None, None])
                    else:
                        data2.append([type_, p3, p3err, None, None, None])
                elif row["IP_C{}_F{}".format(c, f)] == "drift":
                    p3 = row["IP C{} F{}: P3_value".format(c, f)]
                    p3err = row["IP C{} F{}: P3_error".format(c, f)]
                    p2 = row["IP C{} F{}: P2_value".format(c, f)]
                    p2pl = row["IP C{} F{}: P2_error_plus".format(c, f)]
                    p2mi = row["IP C{} F{}: P2_error_minus".format(c, f)]
                    type_ = "drift"
                    if f == df:
                        data2.insert(0, [type_, p3, p3err, p2, p2pl, p2mi])
                    else:
                        data2.append([type_, p3, p3err, p2, p2pl, p2mi])
                elif row["IP_C{}_F{}".format(c, f)] == "P2only":
                    p2 = row["IP C{} F{}: P2_value".format(c, f)]
                    p2pl = row["IP C{} F{}: P2_error_plus".format(c, f)]
                    p2mi = row["IP C{} F{}: P2_error_minus".format(c, f)]
                    type_ = "P2only"
                    data2.append([type_, None, None, p2, p2pl, p2mi])
        if len(data1) > 0:
            mp_data.append(data1)
            if len(data2) > 0:
                ip_data.append(data2)
            else:
                ip_data.append([])
            jnames.append(row["JName_paper"])
            ps.append(float(row["Period [s]"]))
            pdots.append(float(row["Pdot [s/s]"]))
            edots.append(float(row['Edot [ergs/s]']))
            bs.append(float(row["Bsurf [G]"]))
            blcs.append(float(row["Blc [G]"]))
            #if float(row["Blc [G]"]) < 1e-3:
            #    print(row["JName_paper"], row["Blc [G]"])
            ages.append(float(row["Age [yr]"]))
            s1400s.append(row["S1400 [mJy]"])
            w10s.append(row["W10 [ms]"])
            w50s.append(row["W50 [ms]"])
            obsnames.append(row["Obsname"])
            nants.append(int(row["NrAnts"]))
            npulses.append(int(row["NrPulsesUsuable"]))
            snrs.append(float(row["SNRclean"]))
            fftlengths.append(int(row["FFT length (used)"]))
            avmods.append(row["AvMod"])
            avmodeserr.append(row["AvModErr"])


    # TODO add sorting here (if needed)
    #ind = np.argsort(jnames)


    return jnames, ps, pdots, edots, bs, blcs, ages, s1400s, w10s, w50s, obsnames, nants, npulses, snrs, fftlengths, avmods, avmodeserr, mp_data, ip_data



def p3dominant(filename="data/stats.csv"):
    p3s = []
    p3serr = []
    p3ws = []
    p3wserr = []
    p2s = []
    p2serr_p = []
    p2serr_m = []
    edots = []
    jnames = []
    p2as = []
    p2aserr = []


    st = Table.read(filename, format='ascii', header_start=0, data_start=0)
    st = st[(st['Census']=='YES')&(st['PulsarDominantP3Feature'].mask==False)]
    st["Edot [ergs/s]"] = st["Edot [ergs/s]"].astype(float)
    for row in st:
        compn = row['PulsarDominantP3Feature']
        mip, comp = compn.split()
        dname = mip+'dominantDriftFeature_' + comp
        pname = mip+'dominantP3Feature_' + comp
        # if DriftFeature and P3feature both defined, choose the P3 feature; if only DriftFeature or P3feature, take this.
        # kind of strange, but works
        if float(row[dname]) > 0 and float(row[pname])>0:
            fn = 'F'+str(row[pname])
        elif float(row[dname]) > 0:
            fn = 'F'+str(row[dname])
        elif float(row[pname]) > 0:
            fn = 'F'+str(row[pname])
        feature = mip+' '+comp+' '+fn
        p3s.append(float(row[feature + ': P3_value']))
        p3serr.append(float(row[feature + ': P3_error']))
        # P2 the lowest value
        if row[feature.replace(" ", "_")] == "P3only":
            p2s.append("--")
            p2serr_p.append("--")
            p2serr_m.append("--")
        else:
            p2, p2ep, p2em = get_smallest_p2(row)
            p2s.append(p2)
            p2serr_p.append(p2ep)
            p2serr_m.append(p2em)
        # p2 for the dominant feature (not the lowest value)
        #p2s.append(float(row[feature + ': P2_value']))
        #p2serr_p.append(float(row[feature + ': P2_error_plus']))
        #p2serr_m.append(float(row[feature + ': P2_error_minus']))
        p3ws.append(float(row[feature + ': P3_stddev_value']))
        p3wserr.append(float(row[feature + ': P3_stddev_error']))
        edots.append(row['Edot [ergs/s]'])
        jnames.append(row["JName_paper"])
        # P2 assymetry
        comps = row["{} 2dfs nrs".format(mip)].split(",") # get component numbers
        p2a_all = []
        p2ae_all = []
        # TODO change to the dominant feature
        for c in comps:
            p2a_all.append(float(row["C{} Power".format(c)]))
            p2ae_all.append(float(row["C{} PowerErr".format(c)]))
        ind = np.argmax(p2a_all)
        p2as.append(p2a_all[ind])
        p2aserr.append(p2ae_all[ind])
        """
        for ii in range(4):
            print(mip, comp)
            print(row["MP 2dfs nrs"])
            print(ii, row["JName_paper"], row["C{} phase left".format(ii+1)])
        """
    #print(st)
    return p3s, p3serr, p3ws, p3wserr, p2s, p2serr_p, p2serr_m, p2as, p2aserr, edots, jnames


def p3dominant_driftonly(filename="data/stats.csv"):
    p3s = []
    p3serr = []
    p3ws = []
    p3wserr = []
    p2s = []
    p2serr_p = []
    p2serr_m = []
    edots = []
    jnames = []
    p2as = []
    p2aserr = []
    age = []
    bsurf = []
    blc = []
    s1400 = []
    w50 = []
    w10 = []
    periods = []
    pdots = []

    st = Table.read(filename, format='ascii', header_start=0, data_start=0)
    st = st[(st['Census']=='YES')&(st['PulsarDominantP3Feature'].mask==False)]
    st["Edot [ergs/s]"] = st["Edot [ergs/s]"].astype(float)
    for row in st:
        compn = row['PulsarDominantP3Feature']
        mip, comp = compn.split()
        dname = mip+'dominantDriftFeature_' + comp
        pname = mip+'dominantP3Feature_' + comp
        # if DriftFeature and P3feature both defined, choose the P3 feature; if only DriftFeature or P3feature, take this.
        # kind of strange, but works
        if float(row[dname]) > 0 and float(row[pname])>0:
            fn = 'F'+str(row[pname])
        elif float(row[dname]) > 0:
            fn = 'F'+str(row[dname])
        elif float(row[pname]) > 0:
            fn = 'F'+str(row[pname])
        feature = mip+' '+comp+' '+fn
        if row[feature.replace(" ", "_")] != "P3only":
            p3s.append(float(row[feature + ': P3_value']))
            p3serr.append(float(row[feature + ': P3_error']))
            # P2 the lowest value
            p2, p2ep, p2em = get_smallest_p2(row)
            p2s.append(p2)
            p2serr_p.append(p2ep)
            p2serr_m.append(p2em)
            # p2 for the dominant feature (not the lowest value)
            #p2s.append(float(row[feature + ': P2_value']))
            #p2serr_p.append(float(row[feature + ': P2_error_plus']))
            #p2serr_m.append(float(row[feature + ': P2_error_minus']))
            p3ws.append(float(row[feature + ': P3_stddev_value']))
            p3wserr.append(float(row[feature + ': P3_stddev_error']))
            edots.append(row['Edot [ergs/s]'])
            jnames.append(row["JName_paper"])
            age.append(float(row["Age [yr]"]))
            bsurf.append(float(row["Bsurf [G]"])) # float important!
            blc.append(float(row["Blc [G]"])) # float important!
            periods.append(float(row["Period [s]"]))
            pdots.append(float(row["Pdot [s/s]"]))
            try:
                s1400.append(float(row["S1400 [mJy]"])) # float important!
            except:
                s1400.append(1)
            try:
                w50.append(float(row["W50 [ms]"])) # float important!
            except:
                w50.append(1)
            try:
                w10.append(float(row["W10 [ms]"])) # float important!
            except:
                w10.append(1) # float important!

            # P2 assymetry
            comps = row["{} 2dfs nrs".format(mip)].split(",") # get component numbers
            p2a_all = []
            p2ae_all = []
            # TODO change to the dominant feature
            for c in comps:
                p2a_all.append(float(row["C{} Power".format(c)]))
                p2ae_all.append(float(row["C{} PowerErr".format(c)]))
            ind = np.argmax(p2a_all)
            p2as.append(p2a_all[ind])
            p2aserr.append(p2ae_all[ind])

    #print(st)
    return p3s, p3serr, p3ws, p3wserr, p2s, p2serr_p, p2serr_m, p2as, p2aserr, edots, jnames, age, bsurf, blc, s1400, w50, w10, periods, pdots


def p3dominant_driftonly_best(filename="data/stats.csv"):
    p3s = []
    p3serr = []
    p3ws = []
    p3wserr = []
    p2s = []
    p2serr_p = []
    p2serr_m = []
    edots = []
    jnames = []
    p2as = []
    p2aserr = []
    age = []
    bsurf = []
    blc = []
    s1400 = []
    w50 = []
    w10 = []
    periods = []
    pdots = []
    snrs = []

    st = Table.read(filename, format='ascii', header_start=0, data_start=0)
    st = st[(st['Census']=='YES')&(st['PulsarDominantP3Feature'].mask==False)]
    st["Edot [ergs/s]"] = st["Edot [ergs/s]"].astype(float)
    for row in st:
        compn = row['PulsarDominantP3Feature']
        mip, comp = compn.split()
        dname = mip+'dominantDriftFeature_' + comp
        pname = mip+'dominantP3Feature_' + comp
        # if DriftFeature and P3feature both defined, choose the P3 feature; if only DriftFeature or P3feature, take this.
        # kind of strange, but works
        if float(row[dname]) > 0 and float(row[pname])>0:
            fn = 'F'+str(row[pname])
        elif float(row[dname]) > 0:
            fn = 'F'+str(row[dname])
        elif float(row[pname]) > 0:
            fn = 'F'+str(row[pname])
        feature = mip+' '+comp+' '+fn
        if row[feature.replace(" ", "_")] != "P3only":
            p3s.append(float(row[feature + ': P3_value']))
            p3serr.append(float(row[feature + ': P3_error']))
            # P2 the lowest value
            p2, p2ep, p2em = get_smallest_p2(row)
            p2s.append(p2)
            p2serr_p.append(p2ep)
            p2serr_m.append(p2em)
            # p2 for the dominant feature (not the lowest value)
            #p2s.append(float(row[feature + ': P2_value']))
            #p2serr_p.append(float(row[feature + ': P2_error_plus']))
            #p2serr_m.append(float(row[feature + ': P2_error_minus']))
            p3ws.append(float(row[feature + ': P3_stddev_value']))
            p3wserr.append(float(row[feature + ': P3_stddev_error']))
            edots.append(row['Edot [ergs/s]'])
            jnames.append(row["JName_paper"])
            age.append(float(row["Age [yr]"]))
            bsurf.append(float(row["Bsurf [G]"])) # float important!
            blc.append(float(row["Blc [G]"])) # float important!
            periods.append(float(row["Period [s]"]))
            pdots.append(float(row["Pdot [s/s]"]))
            snrs.append(float(row["SNRclean"]))
            try:
                s1400.append(float(row["S1400 [mJy]"])) # float important!
            except:
                s1400.append(1)
            try:
                w50.append(float(row["W50 [ms]"])) # float important!
            except:
                w50.append(1)
            try:
                w10.append(float(row["W10 [ms]"])) # float important!
            except:
                w10.append(1) # float important!

            # P2 assymetry
            comps = row["{} 2dfs nrs".format(mip)].split(",") # get component numbers
            p2a_all = []
            p2ae_all = []
            # TODO change to the dominant feature
            for c in comps:
                p2a_all.append(float(row["C{} Power".format(c)]))
                p2ae_all.append(float(row["C{} PowerErr".format(c)]))
            ind = np.argmax(p2a_all)
            p2as.append(p2a_all[ind])
            p2aserr.append(p2ae_all[ind])

    #print(st)
    #return p3s, p3serr, p3ws, p3wserr, p2s, p2serr_p, p2serr_m, p2as, p2aserr, edots, jnames, age, bsurf, blc, s1400, w50, w10, periods, pdots
    return jnames, edots, p3ws, p3wserr, p2as, p2aserr, snrs



def p3dominant_p3only(filename="data/stats.csv"):
    p3s = []
    p3serr = []
    p3ws = []
    p3wserr = []
    p2s = []
    p2serr_p = []
    p2serr_m = []
    edots = []
    jnames = []
    p2as = []
    p2aserr = []
    age = []
    bsurf = []
    blc = []
    s1400 = []
    w50 = []
    w10 = []


    st = Table.read(filename, format='ascii', header_start=0, data_start=0)
    st = st[(st['Census']=='YES')&(st['PulsarDominantP3Feature'].mask==False)]
    st["Edot [ergs/s]"] = st["Edot [ergs/s]"].astype(float)
    for row in st:
        compn = row['PulsarDominantP3Feature']
        mip, comp = compn.split()
        dname = mip+'dominantDriftFeature_' + comp
        pname = mip+'dominantP3Feature_' + comp
        # if DriftFeature and P3feature both defined, choose the P3 feature; if only DriftFeature or P3feature, take this.
        # kind of strange, but works
        if float(row[dname]) > 0 and float(row[pname])>0:
            fn = 'F'+str(row[pname])
        elif float(row[dname]) > 0:
            fn = 'F'+str(row[dname])
        elif float(row[pname]) > 0:
            fn = 'F'+str(row[pname])
        feature = mip+' '+comp+' '+fn
        if row[feature.replace(" ", "_")] == "P3only":
            p3s.append(float(row[feature + ': P3_value']))
            p3serr.append(float(row[feature + ': P3_error']))
            # P2 the lowest value
            p2, p2ep, p2em = get_smallest_p2(row)
            p2s.append(p2)
            p2serr_p.append(p2ep)
            p2serr_m.append(p2em)
            # p2 for the dominant feature (not the lowest value)
            #p2s.append(float(row[feature + ': P2_value']))
            #p2serr_p.append(float(row[feature + ': P2_error_plus']))
            #p2serr_m.append(float(row[feature + ': P2_error_minus']))
            p3ws.append(float(row[feature + ': P3_stddev_value']))
            p3wserr.append(float(row[feature + ': P3_stddev_error']))
            edots.append(row['Edot [ergs/s]'])
            jnames.append(row["JName_paper"])
            age.append(float(row["Age [yr]"])) # float important!
            bsurf.append(float(row["Bsurf [G]"])) # float important!
            blc.append(float(row["Blc [G]"])) # float important!
            try:
                s1400.append(float(row["S1400 [mJy]"])) # float important!
            except:
                s1400.append(1) # float important!
            try:
                w50.append(float(row["W50 [ms]"])) # float important!
            except:
                w50.append(1)
            try:
                w10.append(float(row["W10 [ms]"])) # float important!
            except:
                w10.append(1)
            #print(row["JName_paper"], row['Edot [ergs/s]'], row["Age [yr]"])
            # P2 assymetry
            comps = row["{} 2dfs nrs".format(mip)].split(",") # get component numbers
            p2a_all = []
            p2ae_all = []
            # TODO change to the dominant feature
            for c in comps:
                p2a_all.append(float(row["C{} Power".format(c)]))
                p2ae_all.append(float(row["C{} PowerErr".format(c)]))
            ind = np.argmax(p2a_all)
            p2as.append(p2a_all[ind])
            p2aserr.append(p2ae_all[ind])


    #print(st)
    return p3s, p3serr, p3ws, p3wserr, p2s, p2serr_p, p2serr_m, p2as, p2aserr, edots, jnames, age, bsurf, blc, s1400, w50, w10




def get_smallest_p2(row):
    p2s = []
    p2serr_p = []
    p2serr_m = []
    for i in range(0, 2):
        for j in range(0, 5):
            # MP - main pulse value
            p2 = row["MP C{} F{}: P2_value".format(i+1, j+1)]
            if not np.isnan(float(p2)):
                #if row["MP C{} F{}: P2_value".format(i+1, j+1)].mask is False:
                p2err_p = row["MP C{} F{}: P2_error_plus".format(i+1, j+1)]
                p2err_m = row["MP C{} F{}: P2_error_minus".format(i+1, j+1)]
                p2s.append(float(p2))
                p2serr_p.append(float(p2err_p))
                p2serr_m.append(float(p2err_m))
            # IP - inter pulse value
            p2 = row["IP C{} F{}: P2_value".format(i+1, j+1)]
            if not np.isnan(float(p2)):
                #if row["MP C{} F{}: P2_value".format(i+1, j+1)].mask is False:
                p2err_p = row["IP C{} F{}: P2_error_plus".format(i+1, j+1)]
                p2err_m = row["IP C{} F{}: P2_error_minus".format(i+1, j+1)]
                p2s.append(float(p2))
                p2serr_p.append(float(p2err_p))
                p2serr_m.append(float(p2err_m))
    p2s = np.array(p2s)
    p2serr_p = np.array(p2serr_p)
    p2serr_m = np.array(p2serr_m)
    ind = np.argmin(np.fabs(p2s))
    return float(p2s[ind]), float(p2serr_p[ind]), float(p2serr_m[ind])


def table_1(filename="data/stats.csv"):
    p3s, p3serr, edots, jnames, tabb = extractdomp3s(filename)
    # manual table generation

    p3s = np.array_split(np.array(p3s), 3)
    p3serr = np.array_split(np.array(p3serr), 3)
    jnames = np.array_split(np.array(jnames), 3)

    tab_st = r"""
\begin{table*}
\caption{List of Pulsars with detectable $P_3$}
\centering
\begin{tabular}{cccc|cccc|cccc}
\hline
  No. & Pulsar & $\dot{E}$ & $P_3$ &  No.  & Pulsar & $\dot{E}$ & $P_3$ & No. & Pulsar & $\dot{E}$ & $P_3$  \\
   &   & (10$^{30}$erg~s$^{-1}$) & (in $P$) &   &   & (10$^{30}$erg~s$^{-1}$) & (in $P$)  &   &   & (10$^{30}$erg~s$^{-1}$) & (in $P$) \\
\hline
    """

    tab_end = r"""\hline
\end{tabular}
\label{tab:p3s}
\medskip
\end{table*}
    """

    tab = """ """
    size = len(p3s[0])

    for i in range(size):
        try:
            tab += "{} & {} & {} & {:.1f} & {} & {} & {} & {:.1f} & {} & {} & {} & {:.1f}\\\\ \n ".format(i, jnames[0][i], -1, p3s[0][i], i+size, jnames[1][i], -1, p3s[1][i], i+2*size, jnames[2][i], -1, p3s[2][i])
        except IndexError:
            #tab += "{} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {}\\".format(i, jnames[0][i], -1, -1, i+size, jnames[1][i], -1, -1, "", "", "", "")
            tab += "{} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {}\\\\ \n".format(i, jnames[0][i], -1, -1, "", "", "", "", "", "", "", "")

    #print(len(p3s[2]))
    print(tab_st + tab+ tab_end)


def table_2(filename="data/stats.csv"):
    #p3s, p3serr, edots, jnames, tabb = extractdomp3s(filename)
    p3s, p3serr, p3ws, p3wserr, p2s, p2serr_p, p2serr_m, p2as, p2aserr, edots, jnames = p3dominant(filename)

    # manual table generation
    tab_st = r"""
\begin{table*}
\caption{List of Pulsars with detectable $P_3$}
\centering
\begin{tabular}{ccccc|ccccc|ccccc|ccccc}
\hline
  Pulsar & $P_3$ & $\sigma_{P/P_3}$ &  $P_2$ &  $P_{\rm asym}$ & Pulsar & $P_3$ & $\sigma_{P/P_3}$ &  $P_2$ &  $P_{\rm asym}$ & Pulsar & $P_3$ & $\sigma_{P/P_3}$ &  $P_2$ &  $P_{\rm asym}$ & Pulsar & $P_3$ & $\sigma_{P/P_3}$ &  $P_2$ &  $P_{\rm asym}$ \\
   & (in $P$) & ($10^{-3}$cpp) &  $(^{\circ})$ &  $(\%)$ &  &  (in $P$) & ($10^{-3}$cpp) & $(^{\circ})$ & $(\%)$ &  & (in $P$) & ($10^{-3}$ cpp) &  $(^{\circ})$ &  $(\%)$ &  & (in $P$) & ($10^{-3}$cpp) &  $(^{\circ})$ &  $(\%)$ \\
   % & & $10^{-2}$ & & & & & $10^{-2}$ & & & & & $10^{-2}$ & & & & & $10^{-2}$ & & \\
\hline
    """

    tab_end = r"""\hline
\end{tabular}
\label{tab:p3s}
\medskip
\end{table*}
    """

    pages = 3

    p3s = np.array_split(np.array(p3s), pages)
    p3serr = np.array_split(np.array(p3serr), pages)
    p3ws = np.array_split(np.array(p3ws), pages)
    p3wserr = np.array_split(np.array(p3wserr), pages)
    p2s = np.array_split(np.array(p2s), pages)
    p2serr_p = np.array_split(np.array(p2serr_p), pages)
    p2serr_m = np.array_split(np.array(p2serr_m), pages)
    p2as = np.array_split(np.array(p2as), pages)
    p2aserr = np.array_split(np.array(p2aserr), pages)
    jnames = np.array_split(np.array(jnames), pages)



    f = open("output/table_2.tex", "w")

    for i in range(pages):

        colnum = 4

        p3s_ = np.array_split(np.array(p3s[i]), colnum)
        p3serr_ = np.array_split(np.array(p3serr[i]), colnum)
        p3ws_ = np.array_split(np.array(p3ws[i]), colnum)
        p3wserr_ = np.array_split(np.array(p3wserr[i]), colnum)
        p2s_ = np.array_split(np.array(p2s[i]), colnum)
        p2serr_p_ = np.array_split(np.array(p2serr_p[i]), colnum)
        p2serr_m_ = np.array_split(np.array(p2serr_m[i]), colnum)
        p2as_ = np.array_split(np.array(p2as[i]), colnum)
        p2aserr_ = np.array_split(np.array(p2aserr[i]), colnum)
        jnames_ = np.array_split(np.array(jnames[i]), colnum)

        tab = """ """
        size = len(p3s_[0])

        for j in range(size):
            for k in range(colnum):
                try:
                    tab += "{} & {} & {} & {} & {}".format(jnames_[k][j], P3_error(p3s_[k][j], p3serr_[k][j]), P3w_error(p3ws_[k][j], p3wserr_[k][j]), P2_e(p2s_[k][j], p2serr_p_[k][j], p2serr_m_[k][j]), P2a_error(p2as_[k][j], p2aserr_[k][j]))
                except IndexError:
                    tab += " & & & & "
                if k < colnum-1:
                     tab += " & "
                else:
                     tab += " \\\\ \n "

                """
            try:
                tab += "{} & {} & {} & {} & {:.1f} & {} & {} & {} & {:.1f} & {:.1f} & {} & {} & {} & {:.1f} & {:.1f} & {} & {} & {} & {:.1f} & {:.1f} \\\\ \n ".format(
                    ,
                    jnames_[1][j], P3_error(p3s_[1][j], p3serr_[1][j]), P3w_error(p3ws_[1][j], p3wserr_[1][j]), P2_e(p2s_[1][j], p2serr_p_[1][j], p2serr_m_[1][j]), p2as_[1][j],
                    jnames_[2][j], P3_error(p3s_[2][j], p3serr_[2][j]), P3w_error(p3ws_[2][j], p3wserr_[2][j]), p2s_[2][j], p2as_[2][j],
                    jnames_[3][j], P3_error(p3s_[3][j], p3serr_[3][j]), P3w_error(p3ws_[3][j], p3wserr_[3][j]), p2s_[3][j], p2as_[3][j])
            except IndexError:
                #tab += "{} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {}\\".format(i, jnames[0][i], -1, -1, i+size, jnames[1][i], -1, -1, "", "", "", "")
                tab += ""
                # TODO some records are missing!
                """
        #print(len(p3s[2]))
        f.write(79*"%" + "\n")
        f.write("% PAGE {}".format(i+1) + "\n")
        f.write(79*"%" + "\n")
        f.write(tab_st + tab + tab_end + "\n")


    f.close()


def table_3(filename="data/stats.csv"):
    #p3s, p3serr, edots, jnames, tabb = extractdomp3s(filename)
    p3s, p3serr, p3ws, p3wserr, p2s, p2serr_p, p2serr_m, p2as, p2aserr, edots, jnames = p3plusp2(filename)

    # manual table generation
    tab_st = r"""
\begin{table*}
\caption{List of Pulsars with detectable $P_3$ or $P_2$}
\centering
\begin{tabular}{ccccc|ccccc|ccccc|ccccc}
\hline
  Pulsar & $P_3$ & $\sigma_{P/P_3}$ &  $P_2$ &  $P_{\rm asym}$ & Pulsar & $P_3$ & $\sigma_{P/P_3}$ &  $P_2$ &  $P_{\rm asym}$ & Pulsar & $P_3$ & $\sigma_{P/P_3}$ &  $P_2$ &  $P_{\rm asym}$ & Pulsar & $P_3$ & $\sigma_{P/P_3}$ &  $P_2$ &  $P_{\rm asym}$ \\
   & (in $P$) & ($10^{-3}$cpp) &  $(^{\circ})$ &  $(\%)$ &  &  (in $P$) & ($10^{-3}$cpp) & $(^{\circ})$ & $(\%)$ &  & (in $P$) & ($10^{-3}$ cpp) &  $(^{\circ})$ &  $(\%)$ &  & (in $P$) & ($10^{-3}$cpp) &  $(^{\circ})$ &  $(\%)$ \\
   % & & $10^{-2}$ & & & & & $10^{-2}$ & & & & & $10^{-2}$ & & & & & $10^{-2}$ & & \\
\hline
    """

    tab_end = r"""\hline
\end{tabular}
\label{tab:p3s}
\medskip
\end{table*}
    """

    pages = 4

    p3s = np.array_split(np.array(p3s), pages)
    p3serr = np.array_split(np.array(p3serr), pages)
    p3ws = np.array_split(np.array(p3ws), pages)
    p3wserr = np.array_split(np.array(p3wserr), pages)
    p2s = np.array_split(np.array(p2s), pages)
    p2serr_p = np.array_split(np.array(p2serr_p), pages)
    p2serr_m = np.array_split(np.array(p2serr_m), pages)
    p2as = np.array_split(np.array(p2as), pages)
    p2aserr = np.array_split(np.array(p2aserr), pages)
    jnames = np.array_split(np.array(jnames), pages)


    f = open("output/table_3.tex", "w")

    for i in range(pages):

        colnum = 4

        p3s_ = np.array_split(np.array(p3s[i]), colnum)
        p3serr_ = np.array_split(np.array(p3serr[i]), colnum)
        p3ws_ = np.array_split(np.array(p3ws[i]), colnum)
        p3wserr_ = np.array_split(np.array(p3wserr[i]), colnum)
        p2s_ = np.array_split(np.array(p2s[i]), colnum)
        p2serr_p_ = np.array_split(np.array(p2serr_p[i]), colnum)
        p2serr_m_ = np.array_split(np.array(p2serr_m[i]), colnum)
        p2as_ = np.array_split(np.array(p2as[i]), colnum)
        p2aserr_ = np.array_split(np.array(p2aserr[i]), colnum)
        jnames_ = np.array_split(np.array(jnames[i]), colnum)

        tab = """ """
        size = len(p3s_[0])

        for j in range(size):
            for k in range(colnum):
                #print(type(p3ws_[k][j])) # why string?
                #print(type(p3s_[k][j])) # why strings? array_split?
                try:
                    tab += "{} & {} & {} & {} & {}".format(jnames_[k][j], P3_error(p3s_[k][j], p3serr_[k][j]), P3w_error(p3ws_[k][j], p3wserr_[k][j]), P2_e(p2s_[k][j], p2serr_p_[k][j], p2serr_m_[k][j]), P2a_error(p2as_[k][j], p2aserr_[k][j]))
                except IndexError:
                    tab += " & & & & "
                if k < colnum-1:
                     tab += " & "
                else:
                     tab += " \\\\ \n "

                """
            try:
                tab += "{} & {} & {} & {} & {:.1f} & {} & {} & {} & {:.1f} & {:.1f} & {} & {} & {} & {:.1f} & {:.1f} & {} & {} & {} & {:.1f} & {:.1f} \\\\ \n ".format(
                    ,
                    jnames_[1][j], P3_error(p3s_[1][j], p3serr_[1][j]), P3w_error(p3ws_[1][j], p3wserr_[1][j]), P2_e(p2s_[1][j], p2serr_p_[1][j], p2serr_m_[1][j]), p2as_[1][j],
                    jnames_[2][j], P3_error(p3s_[2][j], p3serr_[2][j]), P3w_error(p3ws_[2][j], p3wserr_[2][j]), p2s_[2][j], p2as_[2][j],
                    jnames_[3][j], P3_error(p3s_[3][j], p3serr_[3][j]), P3w_error(p3ws_[3][j], p3wserr_[3][j]), p2s_[3][j], p2as_[3][j])
            except IndexError:
                #tab += "{} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {}\\".format(i, jnames[0][i], -1, -1, i+size, jnames[1][i], -1, -1, "", "", "", "")
                tab += ""
                # TODO some records are missing!
                """
        #print(len(p3s[2]))
        f.write(79*"%" + "\n")
        f.write("% PAGE {}".format(i+1) + "\n")
        f.write(79*"%" + "\n")
        f.write(tab_st + tab + tab_end + "\n")

    f.close()



def table_4(filename="data/stats.csv"):
    jnames, ps, pdots, edots, bs, blcs, ages, s1400s, w10s, w50s, obsnames, nants, npulses, snrs, fftlengths, avmods, avmodeserr, mp_data, ip_data = p3plusp2_full(filename)

    # manual table generation
    tab_st = r"""
\begin{table*}
\caption{TODO add caption...}
\centering
\begin{tabular}{cccccccccc}
\toprule
Pulsar & $P$ & $\dot{P}$ & $\dot{E}$ & $B_{\rm s}$ & $B_{\rm LC}$ & Age & $S_{1400}$ & $W_{10}$ & $W_{50}$ \\
 & (s)  &  ($10^{-15}$ s/s)  &  ($10^{33}$ ergs/s) & ($10^{12}$G) & ? & (Myrs) & (mJy) & (ms) & (ms) \\ \cmidrule[0.1pt]{3-10}
\multicolumn{3}{r}{Obsname} & $ N_{\rm ants}$ & $N_{\rm pulses}$ & $S/N$ & FFT & $\bar{m}$ \\
\toprule
    """

    tab_end = r"""\hline
\end{tabular}
\label{tab:all_data}
\medskip
\end{table*}
    """

    tab = " "

    for i, name in enumerate(jnames):
        tab += "\specialrule{0.1pt}{1pt}{1pt} \n"
        # \cmidrule[0.1pt]{{3-10}}
        # \specialrule{0.1pt}{1pt}{1pt}
        tab +=  "{} & {} & {} & {} & {} & {} & {} & {} & {} & {} \\\\ \\cmidrule[0.1pt]{{3-10}}  \n".format(name, ps[i], pdots[i]/1e-15, edots[i]/1e33, bs[i]/1e12, blcs[i], ages[i]/1e6, tof(s1400s[i]), tof(w10s[i]), tof(w50s[i]))
        #nants, npulses, snrs, fftlengths, avmods, avmodeserr
        #tab += "\multicolumn{3}{r}{Obsname} & $ N_{\\rm ants}$ & $N_{\\rm pulses}$ & $S/N$ & FFT & $\\bar{m}$ \\\\ \n "
        tab += "\multicolumn{{3}}{{r}}{{ {} }} & {} &  {} & {} & {} & {} \\\\ \n".format(obsnames[i], nants[i], npulses[i], snrs[i], fftlengths[i], avmods[i])
        # main pulse components
        for mp in mp_data[i]:
            if mp[0] == "drift":
                tab += "\multicolumn{{4}}{{r}}{{ drift }} & {} &  {} & {} & {} & {} \\\\ \n".format(mp[1], mp[2], mp[3], mp[4], mp[5])
            elif mp[0] == "P3only":
                tab += "\multicolumn{{4}}{{r}}{{ P3only }} & {} &  {} & -- & -- & -- \\\\ \n".format(mp[1], mp[2])
            elif mp[0] == "P2only":
                tab += "\multicolumn{{4}}{{r}}{{ P2only }} & -- &  -- & {} & {} & {} \\\\ \n".format(mp[3], mp[4], mp[5])
        for ip in ip_data[i]:
            if ip[0] == "drift":
                tab += "\multicolumn{{4}}{{r}}{{ IP drift }} & {} &  {} & {} & {} & {} \\\\ \n".format(ip[1], ip[2], ip[3], ip[4], ip[5])
            elif ip[0] == "P3only":
                tab += "\multicolumn{{4}}{{r}}{{ IP P3only }} & {} &  {} & -- & -- & -- \\\\ \n".format(ip[1], ip[2])
            elif ip[0] == "P2only":
                tab += "\multicolumn{{4}}{{r}}{{ IP P2only }} & -- &  -- & {} & {} & {} \\\\ \n".format(ip[3], ip[4], ip[5])
        # limit the table records (for test)
        if i == 10:
            break

    f = open("output/table_4.tex", "w")
    f.write(tab_st + tab + tab_end + "\n")
    f.close()

def tof(s14):
    try:
        return float(s14)
    except ValueError:
        return "--"


def P3_error(p3, p3err):
    if p3 == "--":
        return "--"
    p3u = ufloat(float(p3), float(p3err))
    #res = "${:.1uf}$".format(p3u).replace("+/-", "\pm") # OK
    #res = "${:.1ufL}$".format(p3u) # # also OK
    res = "${:.1ufS}$".format(p3u)
    return res


def P2a_error(p2a, p2aerr):
    if p2a == "--":
        return "--"
    p2au = ufloat(float(p2a), float(p2aerr))
    res = "${:.1ufS}$".format(p2au)
    return res

def P3_error_old(p3, p3err):
    if p3err < 1:
        err = float("{:.1g}".format(p3err))
        p3 = "{:.1f}".format(p3)
    else:
        err = int(p3err)
        p3 = int(p3)
    return "${} \pm {}$".format(p3, err)

def P3w_error(p3w, p3werr):
    if p3w == "--":
        return "--"
    times = 1000
    p3w = float(p3w) * times
    p3werr = float(p3werr) * times
    p3wu = ufloat(float(p3w), float(p3werr))
    return "${:.1ufS}$".format(p3wu)
    """
    if p3w < 1:
        return "${:.1g}$".format(p3w)
    else:
        return "${:d}$".format(int(p3w))
    """

def P2_e(p2, p2err_p, p2err_m):
    """ # shorter
    if p2 < 0:
        p2u = ufloat(float(p2), float(p2err_m))
    else:
        p2u = ufloat(float(p2), float(p2err_p))
    res = "${:.1ufS}$".format(p2u)
    return res
    """
    if p2 == "--":
        return "--"
    else:
        p2 = float(p2)
        p2err_p = float(p2err_p)
        p2err_m = float(p2err_m)

    if (p2err_p < 1):
        err_p = float("{:.1g}".format(p2err_p))
    else:
        err_p = int(p2err_p)
    if (p2err_m < 1):
        err_m = float("{:.1g}".format(p2err_m))
    else:
        err_m = int(p2err_m)
    if np.fabs(p2) < 1:
        p2 = "{:.1g}".format(p2)
    else:
        p2 = int(p2)
    return r"${}^{{ +{}}}_{{ -{}}}$".format(p2, err_p, err_m)


def generate_p3_edot(p3s, ep3s, edots, edot_min, edot_max):

    # to fit linear dependence
    p3s_fit = []
    ep3s_fit = []
    edots_fit = []

    for i,edot in enumerate(edots):
        if (edot >= edot_min) and (edot <= edot_max):
            p3s_fit.append(p3s[i])
            ep3s_fit.append(ep3s[i])
            edots_fit.append(edots[i])

    # fitting linear dependence
    fun = lambda v, x: v[0] * x + v[1]
    v0 = [-0.6, 20]
    x, y, v = least_sq(np.log10(np.array(edots_fit)), np.log10(np.array(p3s_fit)), fun, v0, xmax=None)
    sigma = np.std(np.log10(np.array(p3s_fit))) #/ 2 # TODO this is wrong? you lazy bastard
    #print("Fitted parameters: ", v)
    x = 10 ** x
    y = 10 ** y

    p3pred = lambda x: x ** v[0] * 10 ** v[1]

    # generating data
    p3s_model = []
    p3s_model_obs = []
    edots_model = []

    for i,edot in enumerate(edots):
        edot_lin = np.log10(edot)
        p3_lin = np.random.normal(fun(v, edot_lin), sigma)
        p3 = 10 ** p3_lin
        edots_model.append(edots[i])
        p3s_model.append(p3)
        if p3 > 2:
            p3s_model_obs.append(p3)
        else:
            for n in range(1000):
                p3obs = np.abs(p3 / (1 - n * p3))
                if p3obs > 2:
                    break
            p3s_model_obs.append(p3obs)

    divs = np.abs((np.log10(p3s) - np.log10(p3s_model)) / np.log10(ep3s)) # divided by errors
    divergence = np.sum(divs)
    return p3s_model, p3s_model_obs, edots_model, v, divergence


def calculate_divergance_old(data1, data2, xdata, err):
    # NOPE
    #divs = np.abs((data1 - data2)) / err # divided by errors
    #divergence = np.sum(divs)
    #return divergence
    data2 = np.copy(data2)
    #print(err)

    sum = 0
    for i, y1 in enumerate(data1):
        min = 1e50
        x1 = xdata[i]
        for j,y2 in enumerate(data2):
            x2 = xdata[j]
            div = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) #/ err[i] Need to normalize err first (some smaller than one)
            if div < min:
                min = div
                ind = j
        sum += min
        data2[ind] = 1e50 # do not include this point anymore
    return sum


def calculate_divergance(data1, data2, xdata, err, sample=30):
    idx = np.argsort(xdata)

    d1 = data1[idx]
    d2 = data2[idx]
    xd = xdata[idx]
    er = err[idx]

    rngs = np.array(range(0, len(xd)+1, sample))
    if rngs[-1] < len(xd):
        rngs[-1] = len(xd)

    pvals = np.empty(len(rngs)-1)

    for i in range(len(rngs)-1):
        dd1 = d1[rngs[i]:rngs[i+1]]
        dd2 = d2[rngs[i]:rngs[i+1]]
        p = ttest_ind(dd1, dd2).pvalue
        #p = ttest_ind(10**np.array(dd1), 10**np.array(dd2)).pvalue
        pvals[i] = p

    xds = (rngs + int(sample/2))[:-1]
    bins = len(pvals) * 2
    pvals_bin = np.empty(bins)
    xpvals_bin = np.empty(bins)
    for i in range(0, bins, 2):
        pvals_bin[i] = pvals[i//2]
        pvals_bin[i+1] = pvals[i//2]

    xpvals_bin[0] = xd[0]
    xpvals_bin[-1] = xd[-1]
    j = 1
    for i in range(1, bins - 1, 2):
        xpvals_bin[i] = xd[xds[j]]
        xpvals_bin[i+1] = xd[xds[j]]
        j += 1

    return pvals, xds, pvals_bin, xpvals_bin


def find_bestp3edot(data, y=(-1.5, 1.5), a=(-1.0, 1.0), size=20, repeat=20, sig_thresh_proc=10):
    # TODO repeat not implemented yet
    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])

    """
    # YES this is the JiuJitsu way :D
    f = open("data/p3_edot.txt", "w")
    for i in range(len(p3s)):
        f.write("{} {} {}\n".format(p3s[i], ep3s[i], edots[i]))
    f.close()
    return
    #"""

    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))

    sig_thresh = int(len(edots) * sig_thresh_proc / 100)

    a_ = np.linspace(a[0], a[1], num=size)
    ys = np.linspace(y[0], y[1], num=size)
    X = np.log10(1e32)

    edots_model = np.empty(len(edots))
    p3s_model = np.empty(len(p3s))
    edots_model.fill(np.nan)
    p3s_model.fill(np.nan)
    p3_sigma = np.empty(len(edots))
    sigs = np.empty(len(edots))
    sig_pos = 0
    sig_neg = 0
    SIGMA = np.std(p3s10)

    two10 = np.log10(2)
    t_tests = np.empty([size, size, repeat])
    t_test = np.empty([size, size])
    t_test.fill(np.nan)
    acc_min = 1e50
    acc_ind = None

    plot = False

    for i in range(size):
        for j in range(size):
            for ii in range(repeat):
                a = a_[j]
                b = ys[i] - a * X
                p3fun = lambda x: a * x + b
                xline = np.linspace(28, 37, num=100)
                yline = p3fun(xline)
                """
                # estimate sigma? not needed?
                for k in range(len(edots10)):
                    p3_sigma[k] = p3fun(edots10[k])
                    sigs[k] = p3_sigma[k] - p3s10[k]
                    if sigs[k] > 0:
                        sig_pos += 1
                    else:
                        sig_neg += 1
                if sig_pos < sig_thresh or sig_neg < sig_thresh:
                    sigma = SIGMA
                else:
                    sigma = np.min([np.abs(np.min(sigs)), np.max(sigs)]) / 3 # TODO improve that!
                sig_pos, sig_neg = 0, 0
                """
                sigma = SIGMA
                # model data
                p3s_notobs = []
                edots_notobs = []
                skip_a = False
                for k in range(len(edots10)):
                    p3 = np.random.normal(p3fun(edots10[k]), sigma)
                    if p3 >= two10:
                        p3s_model[k] = p3
                        edots_model[k] = edots10[k]
                    else:
                        p3s_notobs.append(p3)
                        edots_notobs.append(edots10[k])
                        p3obs = None
                        for n in range(1000):
                            p3nolog = 10 ** p3
                            # is it working p3obs_n removed
                            p3obss = np.abs(p3nolog / (1 - n * p3nolog))
                            if p3obss > 2:
                                p3obs = p3obss
                                break
                        if p3obs != None:
                            p3s_model[k] = np.log10(p3obs) # it is ok
                            edots_model[k] = edots10[k]
                        else:
                            skip_a = True
                            t_tests[i, j, ii] = 11
                            #print("skip", i, j, ii, "p3obss", p3obss, "n", n, edots10[k])
                            pass
                """
                if a > 0:
                    skip_a = True # TODO HERE!!!
                    t_tests[i, j, ii] = 10
                """
                if skip_a is False:
                    # randomize observations # something is wrong here
                    #for zz in range(len(p3s10)):
                    #    p3s10_rand[zz] = np.log10(np.random.normal(p3s[zz], np.abs(ep3s[zz])))

                    #pvals, xpvals, pvals_bin, xpvals_bin = calculate_divergance(p3s10, p3s_model, edots10, ep3s) # Note ep3s not ep3s10 (this is fine..?)
                    pvals, xpvals, pvals_bin, xpvals_bin = calculate_divergance(p3s10_rand, p3s_model, edots10, ep3s) # Note ep3s not ep3s10 (this is fine..?)
                    acc = len([pv for pv in pvals if pv < 0.05])
                    t_tests[i, j, ii] = acc
                    if acc < acc_min:
                        acc_min = acc
                        acc_ind = (i, j)
                        p3s_best = np.copy(p3s_model)
                        edots_best = np.copy(edots_model)
                        p3s_notobs_best = np.copy(p3s_notobs)
                        edots_notobs_best = np.copy(edots_notobs)
                        xline_best = xline
                        yline_best = yline
                        pvals_best = pvals
                        xpvals_best = xpvals
                        pvals_binbest = pvals_bin
                        xpvals_binbest = xpvals_bin
                plot = False
                if plot is True:
                    pl.figure(figsize=(13.149606, 13.946563744))
                    pl.subplot(2,1,1)
                    pl.minorticks_on()
                    pl.scatter(edots10, p3s10)
                    pl.scatter(edots_model, p3s_model)
                    pl.scatter(edots_notobs, p3s_notobs)
                    pl.plot(xline, yline)
                    xlims = pl.xlim()
                    pl.subplot(2,1,2)
                    pl.minorticks_on()
                    #pl.scatter(xpvals, pvals)
                    try:
                        pl.plot(xpvals_bin, pvals_bin, lw=2, c="C1")
                    except:
                        pass
                    pl.axhline(y=0.05, c="C2", ls="--")
                    pl.xlim(xlims[0], xlims[1])
                    #pl.scatter(edots10, sigs, c="green")
                    pl.show()
                print(i, j, ii)
            #print(t_tests[i, j])
            t_test[i, j] = np.mean(t_tests[i, j]) # OK?

    #s = 0.7

    """
    print(acc_ind, )
    print("y= ", ys[acc_ind[0]], "a=", a_[acc_ind[1]])

    pl.figure(figsize=(13.149606, 13.946563744))
    pl.subplot(2,1,1)
    pl.minorticks_on()
    pl.scatter(edots10, p3s10)
    pl.scatter(edots_best, p3s_best)
    pl.scatter(edots_notobs_best, p3s_notobs_best, c="tab:grey")
    pl.plot(xline_best, yline_best)
    pl.scatter([X for i in range(len(ys))], ys, c="black")
    pl.scatter(X, ys[acc_ind[0]], c='tab:red')
    xlims = pl.xlim()
    pl.subplot(2,1,2)
    pl.minorticks_on()
    #pl.scatter(xpvals_best, pvals_best)
    pl.plot(xpvals_binbest, pvals_binbest, lw=2, c="C1")
    pl.axhline(y=0.05, c="C2", ls="--")
    pl.xlim(xlims[0], xlims[1])
    #pl.scatter(xg, yg, s=0.1, c="red")
    pl.show()
    #"""

    pl.figure(figsize=(13.149606, 13.946563744))
    pl.minorticks_on()
    #pl.imshow(t_test.transpose(), origin="lower", extent=[ys[0], ys[-1], a_[0], a_[-1]])
    pl.imshow(t_test, origin="lower", extent=[a_[0], a_[-1], ys[0], ys[-1]])
    pl.xlabel("a")
    pl.ylabel("y")
    ticks = np.linspace(0, 2, num=10)
    #pl.colorbar(ticks=ticks)
    pl.colorbar()
    pl.show()


def check_p3edot(data, y=-0.1, a=-0.6, sig_thresh_proc=10):
    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])

    """
    # YES this is the JiuJitsu way :D
    f = open("data/p3_edot.txt", "w")
    for i in range(len(p3s)):
        f.write("{} {} {}\n".format(p3s[i], ep3s[i], edots[i]))
    f.close()
    return
    #"""

    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))

    #sig_thresh = int(len(edots) * sig_thresh_proc / 100)

    X = np.log10(1e32)

    edots_model = np.empty(len(edots))
    p3s_model = np.empty(len(p3s))
    edots_model.fill(np.nan)
    p3s_model.fill(np.nan)
    #p3_sigma = np.empty(len(edots))
    #sigs = np.empty(len(edots))
    #sig_pos = 0
    #sig_neg = 0

    sigma = np.std(p3s10)
    two10 = np.log10(2)

    # adding nsp dependence
    nsp_fun = lambda v, x: v[0] * x + v[1] + x ** 2 * v[2]
    v0 = [1, 1, 1]
    xpoints = np.array([29, 31, 33, 35, 37]) # np.array([28, 30, 33, 38])
    #ypoints = np.array([45, 30, 15, 10, 5])# np.array([60, 30, 15, 2])
    #ypoints = np.array([5, 10, 15, 30, 45])# np.array([60, 30, 15, 2])
    #ypoints = np.array([10, 10, 10, 10, 10])# np.array([60, 30, 15, 2])
    ypoints = np.array([20, 20, 20, 20, 20])# np.array([60, 30, 15, 2])
    x_, y_, v_nsp = least_sq(xpoints, ypoints, nsp_fun, v0, xmax=None)

    b = y - a * X

    p3fun = lambda x: a * x + b
    xline = np.linspace(28, 37, num=100)
    yline = p3fun(xline)

    cpps = []
    ecpps = []

    # model data
    p3s_notobs = []
    edots_notobs = []
    skip_a = False
    for k in range(len(edots10)):
        p3 = np.random.normal(p3fun(edots10[k]), sigma)
        nsp = int(nsp_fun(v_nsp, edots10[k]))
        print(nsp)
        if p3 >= two10:
            p3s_model[k] = p3
            edots_model[k] = edots10[k]
            cpp = 1 / 10 ** p3
            cpps.append(cpp)
            ecpps.append(edots10[k])
        elif 10 ** p3 < 1 / nsp:
            p3 = np.random.normal(1 / nsp, sigma)
            while 10 ** p3 < 1 / nsp:
                p3 = np.random.normal(1 / nsp, sigma)
            if p3 >= two10:
                p3s_model[k] = p3
                edots_model[k] = edots10[k]
                cpp = 1 / 10 ** p3
                cpps.append(cpp)
                ecpps.append(edots10[k])
            else:
                p3s_notobs.append(p3)
                edots_notobs.append(edots10[k])
                p3obs = None # not in log scale
                for n in range(1000):
                    p3nolog = 10 ** p3
                    p3obss = np.abs(p3nolog / (1 - n * p3nolog))
                    if p3obss > 2:
                        p3obs = p3obss
                        cpp = 1 / p3obs
                        cpps.append(cpp)
                        ecpps.append(edots10[k])
                        if edots10[k] > 34:
                            print(" n ", n, " ", cpp, " edot ", edots10[k])
                        break
                if p3obs != None:
                    p3s_model[k] = np.log10(p3obs) # it is ok
                    edots_model[k] = edots10[k]
        else:
            p3s_notobs.append(p3)
            edots_notobs.append(edots10[k])
            p3obs = None # not in log scale
            for n in range(1000):
                p3nolog = 10 ** p3
                p3obss = np.abs(p3nolog / (1 - n * p3nolog))
                if p3obss > 2:
                    p3obs = p3obss
                    cpp = 1 / p3obs
                    cpps.append(cpp)
                    ecpps.append(edots10[k])
                    if edots10[k] > 34:
                        print(" n ", n, " ", cpp, " edot ", edots10[k])
                    break
            if p3obs != None:
                p3s_model[k] = np.log10(p3obs) # it is ok
                edots_model[k] = edots10[k]

    """
    hi1, hi1b, me1, err1, xs1, med1 = plo.create_scatter_xy(ecpps, cpps, 15)
    pl.figure()
    #pl.scatter(ecpps, cpps)
    pl.plot(hi1b, hi1, c="tab:red", lw=3, ls="--", label=r"$P_3$ width")
    pl.errorbar(xs1, me1, yerr=err1, c="tab:red", lw=0, marker="o", ms=5, elinewidth=1)
    pl.plot(xs1, med1, c="tab:red", lw=0, marker="o", ms=3)
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"(cpp)")
    pl.minorticks_on()
    pl.show()
    pl.close()
    """

    # randomize observations # turned off... but do not comment..
    #"""
    for zz in range(len(p3s10)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s10_rand[zz] = np.log10(ran)
        #print("p3 ", p3s[zz], "new p3", ran, "error", ep3s[zz])
        #p3s10_rand[zz] = p3s10[zz] # not random
    #"""

    # saving synthetic sample to a file
    """
    f = open("data/p3_edot-synthetic.txt", "w")
    for i in range(len(p3s)):
        f.write("{} {} {}\n".format(10**p3s_model[i], 0.1*10**p3s_model[i], edots[i]))
    f.close()
    #return
    #"""


    pvals, xpvals, pvals_bin, xpvals_bin = calculate_divergance(p3s10_rand, p3s_model, edots10, ep3s, sample=50)
    #pvals, xpvals, pvals_bin, xpvals_bin = calculate_divergance(p3s10, p3s_model, edots10, ep3s) # Note ep3s not ep3s10 (is it fine..?)

    divergence = len([pv for pv in pvals if pv < 0.05])
    print("Divergence: ", divergence)

    pl.figure(figsize=(13.149606, 13.946563744))
    pl.subplot(2,1,1)
    pl.minorticks_on()
    pl.scatter(edots10, p3s10_rand, alpha=0.5, ec="None")
    pl.scatter(edots_model, p3s_model, alpha=0.5, ec="None")
    pl.scatter(edots_notobs, p3s_notobs)
    pl.plot(xline, yline)
    pl.axhline(y=two10, ls=":", c="black")
    xlims = pl.xlim()
    pl.subplot(2,1,2)
    pl.minorticks_on()
    #pl.scatter(xpvals, pvals)
    try:
        pl.plot(xpvals_bin, pvals_bin, lw=2, c="C1")
    except:
        pass
    pl.axhline(y=0.05, c="C2", ls="--")
    pl.xlim(xlims[0], xlims[1])
    #pl.scatter(edots10, sigs, c="green")
    filename = "output/check_p3edot.pdf"
    print(filename)
    pl.savefig(filename)

    pl.show()


def check_p3edot_fun(data, sig_thresh_proc=10):
    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])

    """
    # YES this is the JiuJitsu way :D
    f = open("data/p3_edot.txt", "w")
    for i in range(len(p3s)):
        f.write("{} {} {}\n".format(p3s[i], ep3s[i], edots[i]))
    f.close()
    return
    #"""

    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))

    #sig_thresh = int(len(edots) * sig_thresh_proc / 100)

    X = np.log10(1e32)

    edots_model = np.empty(len(edots))
    p3s_model = np.empty(len(p3s))
    edots_model.fill(np.nan)
    p3s_model.fill(np.nan)
    #p3_sigma = np.empty(len(edots))
    #sigs = np.empty(len(edots))
    #sig_pos = 0
    #sig_neg = 0

    sigma = np.std(p3s10)
    two10 = np.log10(2)

    # adding nsp dependence
    nsp_fun = lambda v, x: v[0] * x + v[1] + x ** 2 * v[2]
    v0 = [1, 1, 1]
    xpoints = np.array([28, 31, 33, 35, 37]) # np.array([28, 30, 33, 38])
    #ypoints = np.array([45, 30, 15, 10, 5])# np.array([60, 30, 15, 2])
    #ypoints = np.array([5, 10, 15, 30, 45])# np.array([60, 30, 15, 2])
    ypoints = np.array([20, 20, 20, 20, 20])# np.array([60, 30, 15, 2])
    x, y, v = least_sq(xpoints, ypoints, nsp_fun, v0, xmax=None)

    """
    xline = np.linspace(28, 38, num=100)
    yline = []
    for x in xline:
        yline.append(nsp_fun(v, x))
    pl.plot(xline, yline)
    pl.scatter(xpoints, ypoints)
    pl.xlabel(r"$\dot{E}$")
    pl.ylabel(r"$n_{\rm sp}$")
    pl.show()
    #return
    # """

    #p3fun = p3_rs_log10n1
    p3fun = p3_log10n2
    xline = np.linspace(28, 37, num=100)
    yline = []
    for x in xline:
        yline.append(p3fun(x, nsp_fun, v))

    cpps = []
    ecpps = []

    # model data
    p3s_notobs = []
    edots_notobs = []
    skip_a = False
    for k in range(len(edots10)):
        nsp = int(nsp_fun(v, edots10[k]))
        #print(edots10[k], "nsp=", nsp)
        p3 = np.random.normal(p3fun(edots10[k], nsp_fun, v), sigma)
        if p3 >= two10:
            p3s_model[k] = p3
            edots_model[k] = edots10[k]
            cpp = 1 / 10 ** p3
            cpps.append(cpp)
            ecpps.append(edots10[k])
        elif p3 < 1 / nsp:
            p3 = np.random.normal(1 / nsp, sigma)
            while p3 < 1 / nsp:
                p3 = np.random.normal(1 / nsp, sigma)
            if p3 >= two10:
                p3s_model[k] = p3
                edots_model[k] = edots10[k]
                cpp = 1 / 10 ** p3
                cpps.append(cpp)
                ecpps.append(edots10[k])
            else:
                p3s_notobs.append(p3)
                edots_notobs.append(edots10[k])
                p3obs = None # not in log scale
                for n in range(1000):
                    p3nolog = 10 ** p3
                    p3obss = np.abs(p3nolog / (1 - n * p3nolog))
                    if p3obss > 2:
                        p3obs = p3obss
                        cpp = 1 / p3obs
                        cpps.append(cpp)
                        ecpps.append(edots10[k])
                        if edots10[k] > 34:
                            print(" n ", n, " ", cpp, " edot ", edots10[k])
                        break
                if p3obs != None:
                    p3s_model[k] = np.log10(p3obs) # it is ok
                    edots_model[k] = edots10[k]
        else:
            p3s_notobs.append(p3)
            edots_notobs.append(edots10[k])
            p3obs = None # not in log scale
            for n in range(1000):
                p3nolog = 10 ** p3
                p3obss = np.abs(p3nolog / (1 - n * p3nolog))
                if p3obss > 2:
                    p3obs = p3obss
                    cpp = 1 / p3obs
                    cpps.append(cpp)
                    ecpps.append(edots10[k])
                    if edots10[k] > 34:
                        print(" n ", n, " ", cpp, " edot ", edots10[k])
                    break
            if p3obs != None:
                p3s_model[k] = np.log10(p3obs) # it is ok
                edots_model[k] = edots10[k]


    """
    hi1, hi1b, me1, err1, xs1, med1 = plo.create_scatter_xy(ecpps, cpps, 15)
    pl.figure()
    #pl.scatter(ecpps, cpps)
    pl.plot(hi1b, hi1, c="tab:red", lw=3, ls="--", label=r"$P_3$ width")
    pl.errorbar(xs1, me1, yerr=err1, c="tab:red", lw=0, marker="o", ms=5, elinewidth=1)
    pl.plot(xs1, med1, c="tab:red", lw=0, marker="o", ms=3)
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"(cpp)")
    pl.minorticks_on()
    pl.show()
    pl.close()
    """

    # randomize observations # turned off?... but do not comment..
    #"""
    for zz in range(len(p3s10)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s10_rand[zz] = np.log10(ran)
        #print("p3 ", p3s[zz], "new p3", ran, "error", ep3s[zz])
        #p3s10_rand[zz] = p3s10[zz] # not random
    #"""

    # saving synthetic sample to a file
    """
    f = open("data/p3_edot-synthetic.txt", "w")
    for i in range(len(p3s)):
        f.write("{} {} {}\n".format(10**p3s_model[i], 0.1*10**p3s_model[i], edots[i]))
    f.close()
    #return
    #"""


    pvals, xpvals, pvals_bin, xpvals_bin = calculate_divergance(p3s10_rand, p3s_model, edots10, ep3s, sample=50)
    #pvals, xpvals, pvals_bin, xpvals_bin = calculate_divergance(p3s10, p3s_model, edots10, ep3s) # Note ep3s not ep3s10 (is it fine..?)

    divergence = len([pv for pv in pvals if pv < 0.05])
    print("Divergence: ", divergence)

    pl.figure(figsize=(13.149606, 13.946563744))
    pl.subplot(2,1,1)
    pl.minorticks_on()
    pl.scatter(edots10, p3s10_rand, alpha=0.5, ec="None")
    pl.scatter(edots_model, p3s_model, alpha=0.5, ec="None")
    pl.scatter(edots_notobs, p3s_notobs)
    pl.plot(xline, yline)
    pl.axhline(y=two10, ls=":", c="black")
    xlims = pl.xlim()
    pl.subplot(2,1,2)
    pl.minorticks_on()
    #pl.scatter(xpvals, pvals)
    try:
        pl.plot(xpvals_bin, pvals_bin, lw=2, c="C1")
    except:
        pass
    pl.axhline(y=0.05, c="C2", ls="--")
    pl.xlim(xlims[0], xlims[1])
    #pl.scatter(edots10, sigs, c="green")
    filename = "output/check_p3edot.pdf"
    print(filename)
    pl.savefig(filename)

    pl.show()



def drift_nodrift(filename="data/stats.csv"):
    # drift includes p3only
    st = Table.read(filename, format='ascii', header_start=0, data_start=0)
    dr = st[(st['Census']=='YES')&(st['PulsarDominantP3Feature'].mask==False)]
    nodr = st[(st['Census']=='YES')&(st['PulsarDominantP3Feature'].mask==True)]
    dr["Edot [ergs/s]"] = dr["Edot [ergs/s]"].astype(float)
    nodr["Edot [ergs/s]"] = nodr["Edot [ergs/s]"].astype(float)
    dr["Age [yr]"] = dr["Age [yr]"].astype(float)
    nodr["Age [yr]"] = nodr["Age [yr]"].astype(float)
    return dr, nodr


def p3dominant_driftonly_edotinfo(filename="data/stats.csv", thresh=5e32):

    names = []
    p2as = []
    edots = []


    st = Table.read(filename, format='ascii', header_start=0, data_start=0)
    st = st[(st['Census']=='YES')&(st['PulsarDominantP3Feature'].mask==False)]
    st["Edot [ergs/s]"] = st["Edot [ergs/s]"].astype(float)
    for row in st:
        compn = row['PulsarDominantP3Feature']
        mip, comp = compn.split()
        dname = mip+'dominantDriftFeature_' + comp
        pname = mip+'dominantP3Feature_' + comp
        # if DriftFeature and P3feature both defined, choose the P3 feature; if only DriftFeature or P3feature, take this.
        # kind of strange, but works
        if float(row[dname]) > 0 and float(row[pname])>0:
            fn = 'F'+str(row[pname])
        elif float(row[dname]) > 0:
            fn = 'F'+str(row[dname])
        elif float(row[pname]) > 0:
            fn = 'F'+str(row[pname])
        feature = mip+' '+comp+' '+fn
        if row[feature.replace(" ", "_")] != "P3only":

            edot = float(row['Edot [ergs/s]'])
            if edot <= thresh:
                # P2 assymetry
                comps = row["{} 2dfs nrs".format(mip)].split(",") # get component numbers
                p2a_all = []
                p2ae_all = []
                # TODO change to the dominant feature
                for c in comps:
                    p2a_all.append(float(row["C{} Power".format(c)]))
                    p2ae_all.append(float(row["C{} PowerErr".format(c)]))
                    ind = np.argmax(p2a_all)
                p2as.append(p2a_all[ind])
                names.append(row["JName_paper"])
                edots.append(row['Edot [ergs/s]'])

                #p2aserr.append(p2ae_all[ind])
            """
            p3s.append(float(row[feature + ': P3_value']))
            p3serr.append(float(row[feature + ': P3_error']))
            # P2 the lowest value
            p2, p2ep, p2em = get_smallest_p2(row)
            p2s.append(p2)
            p2serr_p.append(p2ep)
            p2serr_m.append(p2em)
            # p2 for the dominant feature (not the lowest value)
            #p2s.append(float(row[feature + ': P2_value']))
            #p2serr_p.append(float(row[feature + ': P2_error_plus']))
            #p2serr_m.append(float(row[feature + ': P2_error_minus']))
            p3ws.append(float(row[feature + ': P3_stddev_value']))
            p3wserr.append(float(row[feature + ': P3_stddev_error']))
            edots.append(row['Edot [ergs/s]'])
            jnames.append(row["JName_paper"])
            age.append(float(row["Age [yr]"]))
            bsurf.append(float(row["Bsurf [G]"])) # float important!
            blc.append(float(row["Blc [G]"])) # float important!
            # P2 assymetry
            comps = row["{} 2dfs nrs".format(mip)].split(",") # get component numbers
            p2a_all = []
            p2ae_all = []
            # TODO change to the dominant feature
            for c in comps:
                p2a_all.append(float(row["C{} Power".format(c)]))
                p2ae_all.append(float(row["C{} PowerErr".format(c)]))
            ind = np.argmax(p2a_all)
            p2as.append(p2a_all[ind])
            p2aserr.append(p2ae_all[ind])
            """
    import subprocess

    idx = np.argsort(p2as)
    names = np.array(names)[idx]
    edots = np.array(edots)[idx]
    p2as = np.array(p2as)[idx]
    ht = "https://www.atnf.csiro.au/people/Simon.Johnston/meerkat/"
    for i in range(len(names)):

        print("{}_{}".format(i+1, names[i]))
        cmd = "firefox {}{}{}".format(ht, names[i], ".html")
        subprocess.call(cmd, shell = True)
    #print(edots)
    #print(st)
    return

""" fits line to output/line.txt (result of best_p3edot_geoff in bestp3_edot.jl)"""
def fit_line():
    x = []
    y = []
    f = open("output/line.txt")
    lines = f.readlines()
    for line in lines:
        res = line.split()
        x.append(float(res[0]))
        y.append(float(res[1]))
    x = np.array(x)
    y = np.array(y)
    fun = lambda v, x: v[0] * x + v[1]
    v0 = [-1, 1]
    xn, yn, v = least_sq(x, y, fun, v0, xmax=None)
    print("Fitted parameters:", v)
    f.close()
    return fun, v


def best_drifters(d, filename="data/drifters_list.txt", csvname="data/stats.csv" ):
    pd.set_option("display.max_rows", None)
    """
    jnames, edots, p3ws, p3wserr, p2as, p2aserr, snrs = d
    dataset = {"jnames":jnames, "edots":edots, "p3ws":p3ws, "p3wserr":p3wserr, "p2as":p2as, "p2aserr":p2aserr, "snrs":snrs}
    df = pd.DataFrame(dataset)
    df.sort_values(by=["p3ws", "snrs"], inplace=True)
    print(df)
    """

    df = pd.read_csv(filename)
    #print(df.iloc[1]["jname"])
    #return
    st = Table.read(csvname, format='ascii', header_start=0, data_start=0)
    edots = []
    w50s = []
    snrs = []

    for i in range(len(df)):

        rec = st[st['JName'] == df.iloc[i]["jname"]]
        edots.append(float(rec["Edot [ergs/s]"]))
        w50s.append(float(rec["bp_w50"]))
        snrs.append(float(rec["SNRclean"]))

    df["edots"] = edots
    df["w50"] = w50s
    df["SNR"] = snrs

    df = df.sort_values(by=["SNR"], ascending=False)
    #df = df.sort_values(by=["w50"], ascending=False)
    print(df)
    print("Number of sources: ", len(df))

def get_coordinates(st, filename="data/pulsars_ra_dec.csv"):
    # get RAJ. DECJ. from ATNF database and save to file
    query = psrqpy.QueryATNF()
    psrs = query.get_pulsars()
    tab = "NAME, RAJ, DECJ\n"
    for name in st["JName"]:
        try:
            tab += "{}, {}, {}\n".format(name,  psrs[name].RAJ, psrs[name].DECJ)
            print(name)
        except:
            print("No pulsar {} found...".format(name))
    f = open(filename, "w")
    f.write(tab)
    f.close()

def read_coordinates(filename="data/pulsars_ra_dec.csv"):
    tab = Table.read(filename)
    return tab

def all_drifters(filename="data/stats.csv"):
    #pd.set_option("display.max_rows", None)
    
    # MeerKAT drifters!
    st = Table.read(filename, format='ascii', header_start=0, data_start=0)
    st = st[(st['Census']=='YES')&(st['PulsarDominantP3Feature'].mask==False)] # get all drifters
    st["Edot [ergs/s]"] = st["Edot [ergs/s]"].astype(float)

    # PRESS observations
    press_psrs = np.loadtxt("data/press.txt", dtype=str)
    press_table = Table({"JName":press_psrs})


    cross_match = join(st, press_table, keys="JName")


    print("MeerKAT drifters: ", len(st))
    print("PRESS pulsars: ", len(press_psrs))
    #print(cross_match)
    print("MeerKAT/Press cross match: ", len(cross_match))


    # LOFAR potential obs
    #get_coordinates(st) # get coordinates of st pulsars
    locs = read_coordinates()
    #locs["DECJ (float)"] = locs["DECJ"].astype(float)
    #col = Column(range(len(locs)), name="DECJ (float)")
    #locs.add_column(col)
    ras = []
    decs =[]
    for i in range(len(locs["DECJ"])):
        dec = locs["DECJ"][i]
        ra = locs["RAJ"][i]
        coord = SkyCoord(ra=ra, dec=dec, frame="fk5", unit=(u.hourangle, u.deg))
        ras.append(coord.ra.hour)
        decs.append(coord.dec.degree)
    locs["RAJ (float)"] = ras
    locs["DECJ (float)"] = decs
    #print(st)
    ranges = [-10, 0, 10, 20]
    nums = []
    for r in ranges:
        num = len(locs[locs["DECJ (float)"] > r])
        nums.append(num)
        print("Dec. > {} Num:{}".format(r, num))

    #pl.plot(ranges, nums)
    #pl.show()
    

def check_drifters(filename="data/stats.csv"):
    import paramiko
    import os
    from getpass import getpass
    
    # MeerKAT drifters!
    st = Table.read(filename, format='ascii', header_start=0, data_start=0)
    st = st[(st['Census']=='YES')&(st['PulsarDominantP3Feature'].mask==False)] # get all drifters
    
    # Create paths for each pulsar
    base_path = "/fred/oz005/search_processed"
    paths = []
    for row in st:
        path = f"{base_path}/{row['JName']}/{row['Obsname']}/{row['Freqname']}/single"
        paths.append(path)
    
    # Setup SFTP connection
    hostname = "ozstar.swin.edu.au"
    username = "aszary"
    
    # Create SSH client
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.load_system_host_keys()
    
    try:
        # Connect using SSH key authentication
        ssh.connect(hostname, username=username)
        sftp = ssh.open_sftp()
        
        # Check each path
        existing_paths = []
        path_sizes = {}
        for path in paths:
            try:
                # Try to check if the directory exists
                sftp.stat(path)
                print(f"EXISTS: {path}")
                existing_paths.append(path)
                
                # Get disk usage for the directory
                stdin, stdout, stderr = ssh.exec_command(f"du -sh {path}")
                size_output = stdout.read().decode().strip()
                if size_output:
                    size = size_output.split()[0]
                    path_sizes[path] = size
                    print(f"Size: {size}")
                
            except FileNotFoundError:
                print(f"NOT FOUND: {path}")
        
        # Calculate total size if there are existing paths
        if existing_paths:
            paths_str = " ".join(existing_paths)
            stdin, stdout, stderr = ssh.exec_command(f"du -sch {paths_str} | tail -n1")
            total_output = stdout.read().decode().strip()
            if total_output:
                total_size = total_output.split()[0]
                print(f"\nTotal size: {total_size}")
                path_sizes['total'] = total_size
        
        return existing_paths, path_sizes
        
    except Exception as e:
        print(f"Error: {str(e)}")
        return [], {}
        
    finally:
        if 'sftp' in locals():
            sftp.close()
        if 'ssh' in locals():
            ssh.close()


def download_drifters(filename="data/stats.csv", outdir="output"):
    import paramiko
    import os
    from astropy.table import Table
    
    print("\n=== Starting download_drifters ===")
    
    # MeerKAT drifters!
    print(f"Reading data from: {filename}")
    st = Table.read(filename, format='ascii', header_start=0, data_start=0)
    st = st[(st['Census']=='YES')&(st['PulsarDominantP3Feature'].mask==False)] # get all drifters
    print(f"Found {len(st)} drifters in the table")
    
    # Create paths for each pulsar
    base_path = "/fred/oz005/search_processed"
    paths = []
    for row in st:
        path = f"{base_path}/{row['JName']}/{row['Obsname']}/{row['Freqname']}/single"
        paths.append(path)
    
    print(f"Generated {len(paths)} paths to check")
    
    if not paths:
        print("No paths found to download from")
        return []
    
    # Create output directory
    outdir_abs = outdir # outdir
    print(f"Files will be downloaded to: {outdir_abs}")
    #os.makedirs(outdir_abs, exist_ok=True)
    
    # Setup SFTP connection
    hostname = "ozstar.swin.edu.au"
    username = "aszary"
    print(f"Connecting to {hostname} as {username}...")
    
    # Create SSH client
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.load_system_host_keys()
    
    downloaded_files = []
    
    try:
        # Connect using SSH key authentication
        ssh.connect(hostname, username=username)
        print("SSH connection established")
        sftp = ssh.open_sftp()
        print("SFTP session opened")
        
        # Download files from each path
        for path in paths:
            try:
                # Check if directory exists
                sftp.stat(path)
                print(f"\nFound directory: {path}")
                
                # List and download files
                remote_files = sftp.listdir(path)
                print(f"Found {len(remote_files)} files")
                
                for remote_file in remote_files:
                    remote_path = f"{path}/{remote_file}"
                    local_path = os.path.join(outdir_abs, remote_file)
                    
                    # Download the file
                    print(f"Downloading: {remote_file}")
                    sftp.get(remote_path, local_path)
                    print(f"Saved to: {local_path}")
                    downloaded_files.append(local_path)
                    
            except FileNotFoundError:
                print(f"Directory not found: {path}")
                continue
            except Exception as e:
                print(f"Error downloading from {path}: {str(e)}")
                continue
        
        print(f"\nTotal files downloaded: {len(downloaded_files)}")
        return downloaded_files
        
    except Exception as e:
        print(f"Connection error: {str(e)}")
        return []
        
    finally:
        if 'sftp' in locals():
            print("Closing SFTP session")
            sftp.close()
        if 'ssh' in locals():
            print("Closing SSH connection")
            ssh.close()
        print("=== Finished download_drifters ===\n")
