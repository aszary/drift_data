import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Slider
from scipy.stats import ttest_1samp, ks_2samp
from scipy.optimize import leastsq
import scipy.stats as stats

from modules.functions import least_sq, least_sq1D, odr, least_sq_samex
import modules.functions as fun
import modules.data as da

def test_plot(data):
    pl.figure()
    pl.scatter(data[0], data[1])
    pl.semilogy()
    pl.xlabel("P (s)")
    pl.ylabel("$\dot{E}$ (ergs/s)")
    pl.savefig("output/p0_edot.pdf")
    pl.show()


def p3_edot(da1, da2):
    #print(da.info)
    # all drifting
    p3s = []
    ep3s = []
    edots = []
    # all main components / features
    for i in range(0,2):
        p3s.append([])
        ep3s.append([])
        edots.append([])
        for j in range(0, 5):
            # get valid P3 value
            da_tmp = da1[da1["MP C{} F{}: P3_value".format(i+1, j+1)]>0.0] # has a proper P3
            p3s[i].append(da_tmp["MP C{} F{}: P3_value".format(i+1, j+1)])
            ep3s[i].append(da_tmp["MP C{} F{}: P3_error".format(i+1, j+1)])
            edots[i].append(da_tmp["Edot [ergs/s]"])

    # P3 only
    p3s2 = []
    ep3s2 = []
    edots2 = []
    # all main components / features
    for i in range(0,2):
        p3s2.append([])
        ep3s2.append([])
        edots2.append([])
        for j in range(0, 5):
            # get valid P3 value
            da_tmp = da2[da2["MP C{} F{}: P3_value".format(i+1, j+1)] > 0.0] # has a proper P3
            p3s2[i].append(da_tmp["MP C{} F{}: P3_value".format(i+1, j+1)])
            ep3s2[i].append(da_tmp["MP C{} F{}: P3_error".format(i+1, j+1)])
            edots2[i].append(da_tmp["Edot [ergs/s]"])

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    """
    # all components / features
    for i in range(0,2):
        for j in range(0, 5):
            pl.scatter(edots[i][j], p3s[i][j], color="C{}".format(1+i+j))
    """
    pl.scatter(edots[0][0], p3s[0][0], color="tab:blue", label="drift", s=5, zorder=1)
    pl.errorbar(edots[0][0], p3s[0][0], fmt='none', yerr=np.array(ep3s[0][0]), color="tab:blue", zorder=2)

    pl.scatter(edots2[0][0], p3s2[0][0], color="tab:red", label="P3only", s=5, zorder=3)
    pl.errorbar(edots2[0][0], p3s2[0][0], fmt='none', yerr=np.array(ep3s2[0][0]), color="tab:red", zorder=4)
    #pl.scatter(edots[1][0], p3s[1][0], color="tab:blue")
    #pl.errorbar(edots[1][0], p3s[1][0], fmt='none', yerr=np.array(ep3s[1][0]), color="tab:blue")
    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([1, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    print("output/p3_edot.pdf")
    pl.savefig("output/p3_edot.pdf")
    #pl.show()



def get_p3s_ps(da):

    p3s = []
    ep3s = []
    ps = []

    # all main components / features
    for i in range(0,2):
        p3s.append([])
        ep3s.append([])
        ps.append([])
        for j in range(0, 5):
            # get valid P3 value
            da_tmp = da[da["MP C{} F{}: P3_value".format(i+1, j+1)] > 0.0] # has a proper P3
            p3s[i].append(da_tmp["MP C{} F{}: P3_value".format(i+1, j+1)])
            ep3s[i].append(da_tmp["MP C{} F{}: P3_error".format(i+1, j+1)])
            ps[i].append(da_tmp["Period [s]"])
    return p3s, ep3s, ps



def get_p3s_edots(da):

    p3s = []
    ep3s = []
    edots = []

    # all main components / features
    for i in range(0,2):
        p3s.append([])
        ep3s.append([])
        edots.append([])
        for j in range(0, 5):
            # get valid P3 value
            da_tmp = da[da["MP C{} F{}: P3_value".format(i+1, j+1)] > 0.0] # has a proper P3
            p3s[i].append(da_tmp["MP C{} F{}: P3_value".format(i+1, j+1)])
            ep3s[i].append(da_tmp["MP C{} F{}: P3_error".format(i+1, j+1)])
            edots[i].append(da_tmp["Edot [ergs/s]"])
    return p3s, ep3s, edots


def get_p3s_edots_driftpowers(da):

    p3s = []
    ep3s = []
    edots = []
    dps = []

    # all main components / features
    for i in range(0,2):
        p3s.append([])
        ep3s.append([])
        edots.append([])
        dps.append([])
        for j in range(0, 5):
            # get valid P3 value
            da_tmp = da[da["MP C{} F{}: P3_value".format(i+1, j+1)] > 0.0] # has a proper P3
            p3s[i].append(da_tmp["MP C{} F{}: P3_value".format(i+1, j+1)])
            ep3s[i].append(da_tmp["MP C{} F{}: P3_error".format(i+1, j+1)])
            edots[i].append(da_tmp["Edot [ergs/s]"])
            dps[i].append(da_tmp["C{} Power".format(i+1)]) # drift power here
            print(da_tmp["JName_paper"], da_tmp["MP 2dfs nrs"])
    return p3s, ep3s, edots, dps


def p3_edot2(datas, labels):
    #print(da.info)

    p3s = []
    ep3s = []
    edots = []
    dps = []

    for da in datas:
        p3_, ep3_, edots_, dps_ = get_p3s_edots_driftpowers(da)
        p3s.append(p3_)
        ep3s.append(ep3_)
        edots.append(edots_)
        dps.append(dps_)

    sets = len(datas)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)


    pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    cm = pl.cm.get_cmap('RdYlBu')
    for i in range(sets):
        # first component and feature only
        for j in range(1):
            #print("le ", len(edots[i][0][j])) # TODO different numbers for unique but plots are the same!
            sc = pl.scatter(edots[i][0][j], p3s[i][0][j], c=dps[i][0][j].astype(float), label=labels[i], s=5, zorder=1)
            pl.errorbar(edots[i][0][j], p3s[i][0][j], fmt='none', yerr=np.array(ep3s[i][0][j]), color="C{}".format(i+1), zorder=2)
        pl.colorbar(sc)
    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([1, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    print("output/p3_edot2.pdf")
    pl.savefig("output/p3_edot2.pdf")
    #pl.show()


def p3_edot3(datas, labels):
    p3s = []
    ep3s = []
    edots = []
    dps = []

    for d in datas:
        p3s.append(d[1])
        ep3s.append(d[2])
        edots.append(d[3])
        dps.append(d[4])

    sets = len(datas)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)

    for i in range(sets):
        #sc = pl.scatter(edots[i], p3s[i], c=dps[i], s=5, cmap=cmaps[i], zorder=1)
        sc = pl.scatter(edots[i], p3s[i], c=colors[i], s=5, cmap=cmaps[i], zorder=1)
        pl.errorbar(edots[i], p3s[i], fmt='none', yerr=ep3s[i], color=colors[i], zorder=2, label=labels[i])
        #co = pl.colorbar(sc, shrink=0.5)
        #pl.clim([0,70])
    pl.legend()
    pl.loglog()
    #pl.semilogx()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    #pl.ylim([0.7, 50])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_edot3.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3_edot_sec(datas, labels):
    p3s = []
    ep3s = []
    edots = []
    dps = []
    ps = []

    for d in datas:
        p3s.append(d[1])
        ep3s.append(d[2])
        edots.append(d[3])
        dps.append(d[4])
        ps.append(d[5])

    sets = len(datas)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)

    for i in range(sets):
        #sc = pl.scatter(edots[i], np.array(p3s[i]) * np.array(ps[i]), c=dps[i], s=5, cmap=cmaps[i], zorder=1)
        sc = pl.scatter(edots[i], np.array(p3s[i]) * np.array(ps[i]), c=colors[i], s=5, cmap=cmaps[i], zorder=1)
        pl.errorbar(edots[i], np.array(p3s[i]) * np.array(ps[i]), fmt='none', yerr=np.array(ep3s[i]) * np.array(ps[i]), color=colors[i], zorder=2, label=labels[i])
        #co = pl.colorbar(sc, shrink=0.5)
        #pl.clim([0,70])
    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    #pl.ylim([0.1, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ ($s$)")
    filename = "output/p3_edot_sec.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()



def p3_edot_rahul(datas, labels):

    p3s = []
    ep3s = []
    edots = []

    for da in datas:
        p3_, ep3_, edots_ = get_p3s_edots(da)
        p3s.append(p3_)
        ep3s.append(ep3_)
        edots.append(edots_)

    # assuming aliasing for negative drift
    n = 1
    p3s[1] = p3s_rahul(p3s[1], n=n)
    ep3s[1] = p3s_rahul(ep3s[1], n=n)

    #print(ep3s[1][0][0])

    sets = len(datas)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    for i in range(sets):
        # first component and feature only
        for j in range(1):
            #print("le ", len(edots[i][0][j])) # TODO different numbers for unique but plots are the same!
            pl.scatter(edots[i][0][j], p3s[i][0][j], color="C{}".format(i+1), label=labels[i], s=5, zorder=1)
            pl.errorbar(edots[i][0][j], p3s[i][0][j], fmt='none', yerr=np.array(ep3s[i][0][j]), color="C{}".format(i+1), zorder=2)
    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.1, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    print("output/p3_edot_rahul.pdf")
    pl.savefig("output/p3_edot_rahul.pdf")
    #pl.show()


def p3_edot_rahul2(datas, labels):
    p3s = []
    ep3s = []
    edots = []
    dps = []

    for d in datas:
        p3s.append(d[1])
        ep3s.append(d[2])
        edots.append(d[3])
        dps.append(d[4])

    # assuming aliasing for negative drift
    n = 1
    p3s[1] = p3s_rahul2(p3s[1], n=n)
    ep3s[1] = p3s_rahul2(ep3s[1], n=n)

    sets = len(datas)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)

    for i in range(sets):
        sc = pl.scatter(edots[i], p3s[i], c=dps[i], s=5, cmap=cmaps[i], zorder=1)
        pl.errorbar(edots[i], p3s[i], fmt='none', yerr=ep3s[i], color=colors[i], zorder=2, label=labels[i])
        #co = pl.colorbar(sc, shrink=0.5)
        #pl.clim([0,70])
    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.1, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_edot_rahul_2.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3_p(datas, labels):
    #print(da.info)

    p3s = []
    ep3s = []
    ps = []

    for da in datas:
        p3_, ep3_, p_ = get_p3s_ps(da)
        p3s.append(p3_)
        ep3s.append(ep3_)
        ps.append(p_)

    #print(ps[1])
    #return
    sets = len(datas)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    for i in range(sets):
        # first component and feature only
        pl.scatter(ps[i][0][0], p3s[i][0][0], color="C{}".format(i+1), label=labels[i], s=5, zorder=1)
        pl.errorbar(ps[i][0][0], p3s[i][0][0], fmt='none', yerr=np.array(ep3s[i][0][0]), color="C{}".format(i+1), zorder=2)
    pl.legend()
    #pl.semilogx()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([1, yl[1]])
    pl.xlabel("$P$ (s)")
    pl.ylabel(r"$P_3$ in $P$")
    print("output/p3_p.pdf")
    pl.savefig("output/p3_p.pdf")
    #pl.show()


def p3_pnew(datas, labels):
    #print(da.info)

    p3s = []
    ep3s = []
    ps = []
    dps = []
    pdots = []

    for d in datas:
        p3s.append(d[1])
        ep3s.append(d[2])
        ps.append(d[3])
        dps.append(d[4])
        pdots.append(d[5])

    #print(ps[1])
    #return
    sets = len(datas)

    #cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    #colors = ["tab:red", "tab:blue", "tab:green"]
    colors = ["C1", "C2", "C3"]

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    for i in range(sets):
        #pl.scatter(ps[i], p3s[i], c=dps[i], cmap=cmaps[i], s=5, zorder=1)
        #pl.scatter(np.array(ps[i])**2.1, p3s[i], c=colors[i], s=5, zorder=1, alpha=0.77)
        #pl.scatter(np.array(pdots[i])**(-0.015) * np.array(ps[i]), p3s[i], c=colors[i], s=5, zorder=1, alpha=0.77)
        pl.scatter(np.array(ps[i]), np.array(p3s[i]) , c=colors[i], s=5, zorder=1, alpha=0.77)
        pl.errorbar(ps[i], p3s[i], fmt='none', yerr=ep3s[i], label=labels[i], color=colors[i], zorder=2)
    #pl.legend()
    #pl.semilogx()
    #pl.semilogy()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([1, yl[1]])
    #pl.xlim([1e-20, 1e-11])
    pl.xlabel("$P$ (s)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_pnew.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3_psec(datas, labels):
    #print(da.info)

    p3s = []
    ep3s = []
    ps = []
    dps = []
    pdots = []

    for d in datas:
        p3s.append(d[1])
        ep3s.append(d[2])
        ps.append(d[3])
        dps.append(d[4])
        pdots.append(d[5])

    #print(ps[1])
    #return
    sets = len(datas)

    #cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    #colors = ["tab:red", "tab:blue", "tab:green"]
    colors = ["C1", "C2", "C3"]

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    for i in range(sets):
        pl.scatter(np.array(ps[i]), np.array(p3s[i]) * np.array(ps[i]) , c=colors[i], s=5, zorder=1, alpha=0.77)
    pl.plot(np.array(ps[0]), np.array(ps[0]) * 2., color="black", ls="--" , dashes=(10, 10))
    #pl.legend()
    #pl.semilogx()
    #pl.semilogy()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([7e-2, yl[1]])
    #pl.xlim([1e-20, 1e-11])
    pl.xlabel("$P$ (s)")
    pl.ylabel(r"$P_3$ (s)")
    filename = "output/p3_psec.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()



def p3s_rahul(p3s, n=1):
    # why so many indices i, j, k? does p3_edot_rahul() work?
    new_p3s = []
    for i in range(len(p3s)):
        new_p3s.append([])
        for j in range(len(p3s[i])):
            new_p3s[-1].append([])
            for k in range(len(p3s[i][j])):
                #new_p3s[-1][-1].append([])
                p3_obs = p3s[i][j][k]
                p3 = 1.0 / (n + 1.0 / (p3_obs)) # p3_obs in Periods (taken from Julia version)
                new_p3s[-1][-1].append(p3)
    """
    for p3_obs in p3s: # p3_obs in Periods (taken from Julia version)
        p3 = 1 / (n + 1 / (p3_obs))
        new_p3s.append(p3)
    """
    return new_p3s


def p3s_rahul2(p3s, n=1):
    new_p3s = []
    for i in range(len(p3s)):
        #new_p3s[-1][-1].append([])
        p3_obs = p3s[i]
        p3 = 1.0 / (n + 1.0 / (p3_obs)) # p3_obs in Periods (taken from Julia version)
        new_p3s.append(p3)

    """
    for p3_obs in p3s: # p3_obs in Periods (taken from Julia version)
        p3 = 1 / (n + 1 / (p3_obs))
        new_p3s.append(p3)
    """
    return new_p3s


def p3_p_xcheck(datas, xdata):
    #print(da.info)

    jnames = list(datas[0][0]) + list(datas[1][0])
    p3s = list(datas[0][1]) + list(datas[1][1])
    ep3s = list(datas[0][2]) + list(datas[1][2])
    ps = list(datas[0][3]) + list(datas[1][3])

    xp3s = xdata[1]
    xep3s = xdata[2]
    xps = xdata[3]

    # find blue points
    for i in range(len(p3s)):
        res_list = [j for j in range(len(xp3s)) if xp3s[j] == p3s[i]]
        for ind in res_list:
            #print(ps[i], xps[ind])
            if ps[i] != xps[ind]:
                print(jnames[i], ps[i], xps[ind])

    sets = len(datas)

    colors = ["C0", "C1", "C3"]

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    pl.scatter(np.array(ps), np.array(p3s) , c=colors[0], s=5, zorder=1, alpha=0.77)
    pl.scatter(np.array(xps), np.array(xp3s) , c=colors[1], s=5, zorder=2, alpha=0.5)
    #pl.legend()
    #pl.semilogx()
    #pl.semilogy()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([1, yl[1]])
    #pl.xlim([1e-20, 1e-11])
    pl.xlabel("$P$ (s)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_p_xcheck.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3mean_p(datas, xdata):
    #print(da.info)

    p3s = [list(datas[0][1]) + list(datas[1][1])]
    ep3s = [list(datas[0][2]) + list(datas[1][2])]
    ps = [list(datas[0][3]) + list(datas[1][3])]

    p3s.append(xdata[1])
    ep3s.append(xdata[2])
    ps.append(xdata[3])

    sets = len(p3s)

    mi = np.min([np.min(ps[0]), np.min(ps[1])])
    ma = np.max([np.max(ps[0]), np.max(ps[1])])
    #print(mi)
    #print(ma)

    bins = 5

    pbins = np.logspace(np.log10(mi), np.log10(ma), num=bins+1)

    p3_ = []
    for i in range(sets):
        p3_.append([])
        for j in range(bins):
            p3_[-1].append([])
            for k in range(len(p3s[i])):
                if ps[i][k] >= pbins[j] and ps[i][k] < pbins[j+1]:
                    p3_[-1][-1].append(p3s[i][k])
    p3means = []
    ep3means = []
    for i in range(sets):
        p3means.append([])
        ep3means.append([])
        for j in range(bins):
            p3means[-1].append(np.mean(p3_[i][j]))
            #ep3means[-1].append(np.std(p3_[i][j]) / len(p3_[i][j]))
            ep3means[-1].append(np.std(p3_[i][j]))
    # last step?
    p3means[0].append(p3means[0][-1])
    p3means[1].append(p3means[1][-1])

    # for errorbars
    dlog = (np.log10(pbins[1]) - np.log10(pbins[0])) / 2
    epbins = np.logspace(np.log10(pbins[0]) + dlog, np.log10(pbins[-2]) + dlog, num=bins)

    colors = ["C1", "C2", "C3"]

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    pl.scatter(ps[0], p3s[0] , c="C0", s=5, zorder=1, alpha=0.57, label="Andrzej")
    pl.scatter(ps[1], p3s[1] , c="C1", s=5, zorder=2, alpha=0.57, label="Xiaoxi")
    pl.step(pbins, p3means[0], where='post', c="C0", lw=3, alpha=0.7)
    pl.step(pbins, p3means[1], where='post', c="C1", lw=3, alpha=0.7)
    pl.errorbar(epbins, p3means[1][:-1],  yerr=ep3means[1], c="C1", fmt='none', lw=2)
    pl.legend()
    pl.semilogx()
    #pl.semilogy()
    #pl.loglog()
    yl = pl.ylim()
    #pl.ylim([1, yl[1]])
    #pl.ylim([0, 25])
    pl.ylim([-10, 50])
    #pl.xlim([1e-20, 1e-11])
    pl.xlabel("$P$ (s)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3mean_p.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3_edotp(datas, labels):
    p3s = []
    ep3s = []
    edots = []
    dps = []
    ps = []

    for d in datas:
        p3s.append(d[1])
        ep3s.append(d[2])
        edots.append(d[3])
        dps.append(d[4])
        ps.append(d[5])

    sets = len(datas)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)

    for i in range(sets):
        #sc = pl.scatter(edots[i], np.array(p3s[i]) * np.array(ps[i]), c=dps[i], s=5, cmap=cmaps[i], zorder=1)
        sc = pl.scatter(np.sqrt(np.array(edots[i]) * np.array(ps[i])), np.array(p3s[i]), c=colors[i], s=5, cmap=cmaps[i], zorder=1)
        pl.errorbar(np.sqrt(np.array(edots[i]) * np.array(ps[i])), np.array(p3s[i]), fmt='none', yerr=np.array(ep3s[i]) * np.array(ps[i]), color=colors[i], zorder=2, label=labels[i])
        #co = pl.colorbar(sc, shrink=0.5)
        #pl.clim([0,70])
    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    #pl.ylim([0.1, yl[1]])
    pl.xlabel("$\sqrt{\dot{E} P}$ ($\sqrt{\\rm ergs}$)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_edotp.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def generate_p3_edot(datas, labels, edot_min=5e30, edot_max=2e31, fmod=""): # 3e30-2e31, 5e30-2e31
    # data
    p3s = []
    ep3s = []
    edots = []
    dps = []
    ps = []

    for d in datas:
        p3s += list(d[1])
        ep3s += list(d[2])
        edots += list(d[3])
        dps += list(d[4])
        ps += list(d[5])

    # to fit linear dependence
    p3s_fit = []
    ep3s_fit = []
    edots_fit = []
    ps_fit = []

    for i,edot in enumerate(edots):
        if (edot >= edot_min) and (edot <= edot_max):
            p3s_fit.append(p3s[i])
            ep3s_fit.append(ep3s[i])
            edots_fit.append(edots[i])
            ps_fit.append(ps[i])

    # fitting linear dependence
    fun = lambda v, x: v[0] * x + v[1]
    v0 = [-0.6, 20]
    x, y, v = least_sq(np.log10(np.array(edots_fit)), np.log10(np.array(p3s_fit)), fun, v0, xmax=None)
    sigma = np.std(np.log10(np.array(p3s_fit))) #/ 2 # TODO this is wrong? you lazy bastard
    print(v)
    x = 10 ** x
    y = 10 ** y


    p3pred = lambda x: x ** v[0] * 10 ** v[1]
    xline = np.logspace(30, 35, num=100)
    yline = p3pred(xline)

    # to be modeled only high or all?
    all = True
    p3s_high = []
    ep3s_high = []
    edots_high = []
    ps_high = []
    for i,edot in enumerate(edots):
        if all is True:
            p3s_high.append(p3s[i])
            ep3s_high.append(ep3s[i])
            edots_high.append(edots[i])
            ps_high.append(ps[i])
        else:
            if (edot > edot_max):
                p3s_high.append(p3s[i])
                ep3s_high.append(ep3s[i])
                edots_high.append(edots[i])
                ps_high.append(ps[i])

    # generating data
    p3s_model = []
    p3s_model_obs = []
    edots_model = []

    for i,edot in enumerate(edots_high):
        edot_lin = np.log10(edot)
        p3_lin = np.random.normal(fun(v, edot_lin), sigma)
        p3 = 10 ** p3_lin
        edots_model.append(edots_high[i])
        p3s_model.append(p3)
        if p3 > 2:
            p3s_model_obs.append(p3)
        else:
            for n in range(1000):
                p3obs_p = np.abs(p3 / (1 - n * p3))
                p3obs_n = np.abs(p3 / (1 - (-n) * p3))
                if p3obs_p > 2:
                    p3obs = p3obs_p
                    #print(edots_high[i], n)
                    break
                if p3obs_n > 2:
                    p3obs = p3obs_n
                    #print(edots_high[i], n)
                    break
            p3s_model_obs.append(p3obs)

    pl.rc("font", size=8)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]
    size = 2
    al = 1.0

    #fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    fig = pl.figure(figsize=(3.149606, 1.946563744))  # 8 cm x  cm # golden ratio
    # pl.minorticks_on() # does not work
    pl.subplots_adjust(left=0.19, bottom=0.215, right=0.99, top=0.99)

    sc = pl.scatter(np.array(edots), np.array(p3s), c="tab:blue", s=size, alpha=al, label="obs. data", ec="none")
    #sc = pl.scatter(np.array(edots_fit), np.array(p3s_fit), c="tab:green", s=5, alpha=0.7)
    sc = pl.scatter(np.array(edots_model), np.array(p3s_model), c="tab:grey", s=size, alpha=al, ec="none")
    sc = pl.scatter(np.array(edots_model), np.array(p3s_model_obs), c="tab:red", s=size, alpha=al, label="modeled data", ec="none")
    #print(np.min(edots))
    #sc = pl.scatter(np.array(edots_high), np.array(p3s_high), c="tab:pink", s=5)
    pl.plot(x, y, lw=1, alpha=0.7, color="C1") # fitted dependence
    #pl.plot(xline, yline, lw=2, alpha=0.7, color="C0") # works
    pl.legend(fontsize=5,loc='lower left')
    pl.loglog()
    yl = pl.ylim()
    #pl.ylim([0.007, yl[1]])
    pl.axvline(x=edot_min, c="black", lw=0.3, ls="--")
    pl.axvline(x=edot_max, c="black", lw=0.3, ls="--")
    pl.ylim([0.001, 7e2])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/synthetic_p3_edot{}.pdf".format(fmod)
    print(filename)
    pl.savefig(filename)
    #pl.show()


def generate_p3_edot2(data, edot_min=5e30, edot_max=2e31, fmod=""): # 3e30-2e31, 5e30-2e31
    # data
    p3s = []
    ep3s = []
    edots = []

    p3s += list(data[0])
    ep3s += list(data[1])
    edots += list(data[9])

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
    sigma = np.std(np.log10(np.array(p3s_fit))) #/ 2 # TODO this is wrong or not.  You lazy bastard
    print("Fitted parameters: ", v)
    x = 10 ** x
    y = 10 ** y


    p3pred = lambda x: x ** v[0] * 10 ** v[1]
    xline = np.logspace(28, 37, num=100)
    yline = p3pred(xline)

    # to be modeled only high or all?
    all = True
    p3s_high = []
    ep3s_high = []
    edots_high = []
    for i,edot in enumerate(edots):
        if all is True:
            p3s_high.append(p3s[i])
            ep3s_high.append(ep3s[i])
            edots_high.append(edots[i])
        else:
            if (edot > edot_max):
                p3s_high.append(p3s[i])
                ep3s_high.append(ep3s[i])
                edots_high.append(edots[i])

    # generating data
    p3s_model = []
    p3s_model_obs = []
    edots_model = []

    for i,edot in enumerate(edots_high):
        edot_lin = np.log10(edot)
        p3_lin = np.random.normal(fun(v, edot_lin), sigma)
        p3 = 10 ** p3_lin
        edots_model.append(edots_high[i])
        p3s_model.append(p3)
        if p3 > 2:
            p3s_model_obs.append(p3)
        else:
            for n in range(1000):
                p3obs_p = np.abs(p3 / (1 - n * p3))
                p3obs_n = np.abs(p3 / (1 - (-n) * p3))
                if p3obs_p > 2:
                    p3obs = p3obs_p
                    #print(edots_high[i], n)
                    break
                if p3obs_n > 2:
                    p3obs = p3obs_n
                    #print(edots_high[i], n)
                    break
            p3s_model_obs.append(p3obs)

    #divs = np.abs(np.log10(p3s) - np.log10(p3s_model))
    #divergence = np.sum(divs)

    divs = []
    for i in range(100):
        ge = da.generate_p3_edot(p3s, ep3s, edots, edot_min, edot_max)
        div = ge[4]
        divs.append(div)

    print("Divergence {}: {} +/- {}".format(fmod, np.mean(divs), np.std(divs)))
    #return

    pl.rc("font", size=7)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]
    size = 2
    al = 1.0

    #fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    fig = pl.figure(figsize=(3.149606, 1.946563744))  # 8 cm x  cm # golden ratio
    # pl.minorticks_on() # does not work
    pl.subplots_adjust(left=0.19, bottom=0.215, right=0.99, top=0.99)

    sc = pl.scatter(np.array(edots), np.array(p3s), c="tab:blue", s=size, alpha=al, label="obs. data", ec="none")
    #sc = pl.scatter(np.array(edots_fit), np.array(p3s_fit), c="tab:green", s=5, alpha=0.7)
    sc = pl.scatter(np.array(edots_model), np.array(p3s_model), c="tab:grey", s=size, alpha=al, ec="none")
    sc = pl.scatter(np.array(edots_model), np.array(p3s_model_obs), c="tab:red", s=size, alpha=al, label="modeled data", ec="none")
    #print(np.min(edots))
    #sc = pl.scatter(np.array(edots_high), np.array(p3s_high), c="tab:pink", s=5)
    pl.plot(x, y, lw=1, alpha=0.7, color="C1") # fitted dependence
    pl.plot(xline, yline, lw=2, alpha=0.3, color="C0") # works
    pl.legend(fontsize=5,loc='lower left')
    pl.loglog()
    yl = pl.ylim()
    #pl.ylim([0.007, yl[1]])
    pl.axvline(x=edot_min, c="black", lw=0.3, ls="--")
    pl.axvline(x=edot_max, c="black", lw=0.3, ls="--")
    pl.ylim([7e-3, 7e2])
    #pl.ylim([0.7, 5e2])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/synthetic_p3_edot2{}.pdf".format(fmod)
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3_edot4(data):
    # data
    p3s = []
    ep3s = []
    edots = []

    p3s += list(data[0])
    ep3s += list(data[1])
    edots += list(data[9])


    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["C0", "tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    sc = pl.scatter(edots, p3s, c=colors[0], s=5, zorder=1)
    pl.errorbar(edots, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2, label="dominant feature (drift only)")
    pl.legend()
    pl.loglog()
    #pl.semilogx()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    #pl.ylim([0.7, 50])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_edot4.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3_edot_driftp3only(data, data2):
    # data
    p3s = []
    ep3s = []
    edots = []
    ages = []

    p3s += list(data[0])
    ep3s += list(data[1])
    edots += list(data[9])
    ages += list(data[11])

    p3s2 = []
    ep3s2 = []
    edots2 = []
    ages2 = []

    p3s2 += list(data2[0])
    ep3s2 += list(data2[1])
    edots2 += list(data2[9])
    ages2 += list(data2[11])


    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["C0", "tab:red", "tab:blue", "tab:green"]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    sc = pl.scatter(edots, p3s, c=colors[0], s=5, zorder=1, label="drift")
    pl.errorbar(edots, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2)
    sc = pl.scatter(edots2, p3s2, c="C1", s=5, zorder=1, label="P3only")
    pl.errorbar(edots2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2)

    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_edot_driftp3only.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()



def p3_age(data, data2):
    # data
    p3s = []
    ep3s = []
    edots = []
    ages = []

    p3s += list(data[0])
    ep3s += list(data[1])
    edots += list(data[9])
    ages += list(data[11])

    p3s2 = []
    ep3s2 = []
    edots2 = []
    ages2 = []

    p3s2 += list(data2[0])
    ep3s2 += list(data2[1])
    edots2 += list(data2[9])
    ages2 += list(data2[11])

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["C0", "tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    sc = pl.scatter(ages, p3s, c=colors[0], s=5, zorder=1, label="drift")
    pl.errorbar(ages, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2)
    sc = pl.scatter(ages2, p3s2, c="C1", s=5, zorder=1, label="P3only")
    pl.errorbar(ages2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2)

    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    pl.xlabel("Age (yr)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_age.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3_bsurf(data, data2):
    # data
    p3s = []
    ep3s = []
    edots = []
    ages = []
    bsurf = []

    p3s += list(data[0])
    ep3s += list(data[1])
    edots += list(data[9])
    ages += list(data[11])
    bsurf += list(data[12])

    p3s2 = []
    ep3s2 = []
    edots2 = []
    ages2 = []
    bsurf2 = []

    p3s2 += list(data2[0])
    ep3s2 += list(data2[1])
    edots2 += list(data2[9])
    ages2 += list(data2[11])
    bsurf2 += list(data2[12])

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["C0", "tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    sc = pl.scatter(bsurf, p3s, c=colors[0], s=5, zorder=1, label="drift")
    pl.errorbar(bsurf, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2)
    sc = pl.scatter(bsurf2, p3s2, c="C1", s=5, zorder=1, label="P3only")
    pl.errorbar(bsurf2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2)

    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    pl.xlabel(r"$B_{\rm s}$ (G)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_bsurf.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3_blc(data, data2):
    # data
    p3s = []
    ep3s = []
    blc = []

    p3s += list(data[0])
    ep3s += list(data[1])
    blc += list(data[13])

    p3s2 = []
    ep3s2 = []
    blc2 = []

    p3s2 += list(data2[0])
    ep3s2 += list(data2[1])
    blc2 += list(data2[13])

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["C0", "tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    sc = pl.scatter(blc, p3s, c=colors[0], s=5, zorder=1, label="drift")
    pl.errorbar(blc, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2)
    sc = pl.scatter(blc2, p3s2, c="C1", s=5, zorder=1, label="P3only")
    pl.errorbar(blc2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2)

    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    pl.xlabel(r"$B_{\rm lc}$ (G)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_blc.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3_s1400(data, data2):
    # data
    p3s = []
    ep3s = []
    s1400 = []

    p3s += list(data[0])
    ep3s += list(data[1])
    s1400 += list(data[14])

    p3s2 = []
    ep3s2 = []
    s14002 = []

    p3s2 += list(data2[0])
    ep3s2 += list(data2[1])
    s14002 += list(data2[14])

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["C0", "tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    sc = pl.scatter(s1400, p3s, c=colors[0], s=5, zorder=1, label="drift")
    pl.errorbar(s1400, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2)
    sc = pl.scatter(s14002, p3s2, c="C1", s=5, zorder=1, label="P3only")
    pl.errorbar(s14002, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2)

    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    pl.xlabel(r"$S1400$ (mJy)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_s1400.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3_w50(data, data2):
    # data
    p3s = []
    ep3s = []
    w50 = []

    p3s += list(data[0])
    ep3s += list(data[1])
    w50 += list(data[15])

    p3s2 = []
    ep3s2 = []
    w502 = []

    p3s2 += list(data2[0])
    ep3s2 += list(data2[1])
    w502 += list(data2[15])

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["C0", "tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    sc = pl.scatter(w50, p3s, c=colors[0], s=5, zorder=1, label="drift")
    pl.errorbar(w50, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2)
    sc = pl.scatter(w502, p3s2, c="C1", s=5, zorder=1, label="P3only")
    pl.errorbar(w502, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2)

    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    pl.xlabel(r"$W50$ (ms)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_w50.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3_tests(data, data2):
    # data
    p3s = []
    ep3s = []
    edots = []
    ages = []

    p3s += list(data[0])
    ep3s += list(data[1])
    edots += list(data[9])
    ages += list(data[11])

    p3s2 = []
    ep3s2 = []
    edots2 = []
    ages2 = []

    p3s2 += list(data2[0])
    ep3s2 += list(data2[1])
    edots2 += list(data2[9])
    ages2 += list(data2[11])


    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["C0", "tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    #"""
    sc = pl.scatter(edots, p3s, c=colors[0], s=5, zorder=1, label="drift")
    pl.errorbar(edots, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2)
    sc = pl.scatter(edots2, p3s2, c="C1", s=5, zorder=1, label="P3only")
    pl.errorbar(edots2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2)
    #"""

    #sc = pl.scatter(ages, p3s, c=colors[0], s=5, zorder=1, label="drift")
    #pl.errorbar(ages, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2)
    #sc = pl.scatter(ages2, p3s2, c="C1", s=5, zorder=1, label="P3only")
    #pl.errorbar(ages2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2)

    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    #pl.xlabel("Age (yr)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_tests.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()

def create_hist(dr, nodr, bins, xval="Edot [ergs/s]"):
    # sort tables
    dr.sort(xval)
    nodr.sort(xval)
    mi1 = dr[0][xval]
    ma1 = dr[-1][xval]
    mi2 = nodr[0][xval]
    ma2 = nodr[-1][xval]
    min = np.min([mi1, mi2])
    max = np.max([ma1, ma2])

    ebins = np.logspace(np.log10(min), np.log10(max), num=bins)
    hi1 = np.histogram(dr[xval], bins=ebins)[0]
    hi2 = np.histogram(nodr[xval], bins=ebins)[0]
    hisum = hi1 + hi2
    hifr = hi1 / hisum
    hifrac = []
    for hi in hifr:
        hifrac.append(hi)
        hifrac.append(hi)
    hibins = []
    for i in range(len(ebins) - 1):
        hibins.append(ebins[i])
        hibins.append(ebins[i+1])
    return hibins, hifrac


def create_scatter(dr, bins, xval="Edot [ergs/s]", val="Period [s]"):
    # sort tables
    dr.sort(xval)
    #nodr.sort("Edot [ergs/s]")
    min = dr[0][xval] # mi1
    max = dr[-1][xval] # ma2
    #mi2 = nodr[0]["Edot [ergs/s]"]
    #ma2 = nodr[-1]["Edot [ergs/s]"]
    #min = np.min([mi1, mi2])
    #max = np.max([ma1, ma2])

    ebins = np.logspace(np.log10(min), np.log10(max), num=bins)
    hibins = []
    for i in range(len(ebins) - 1):
        hibins.append(ebins[i])
        hibins.append(ebins[i+1])
    hi = np.zeros(len(hibins))

    vals = [[] for i in range(len(ebins) - 1)]

    for row in dr:
        try:
            v = float(row[val])
        except:
            continue
        if ~np.isnan(v):
            # HERE
            #"""
            if val == "W10 [ms]":# or val == "bp_w50":
                v = v / 1000 / float(row["Period [s]"])
            if val == val == "bp_w50":
                v = v / 360. # assuming in deg.
                #v = v / 1000 / float(row["Period [s]"])
            #"""
            edot = row[xval]
            for i in range(len(ebins) - 1):
                if edot >= ebins[i] and edot <= ebins[i+1]:
                    vals[i].append(v)
                    break

    means = [np.mean(v) for v in vals]
    medians = [np.median(v) for v in vals]
    stds = [np.std(v) for v in vals]
    xs = []
    # get xs for means
    for i in range(len(ebins) - 1):
        ed = 10 ** (np.log10(ebins[i]) + (np.log10(ebins[i+1]) - np.log10(ebins[i])) / 2)
        xs.append(ed)

    # calculate histogram
    for i in range(int(len(hi) / 2)):
        hi[2*i] = means[i]
        hi[2*i + 1] = means[i]
    errs = np.zeros(len(stds))
    # standard error
    for i in range(len(errs)):
        errs[i] = stds[i] / np.sqrt(len(vals[i]))

    return np.array(hi), np.array(hibins), np.array(means), np.array(errs), np.array(xs), np.array(medians)


def create_scatter_xy(edots, ys, bins):
    min = np.min(edots)
    max = np.max(edots)

    ebins = np.logspace(np.log10(min), np.log10(max), num=bins)
    hibins = []
    for i in range(len(ebins) - 1):
        hibins.append(ebins[i])
        hibins.append(ebins[i+1])
    hi = np.zeros(len(hibins))

    vals = [[] for i in range(len(ebins) - 1)]

    for (i,y) in enumerate(ys):
        x = edots[i]
        for j in range(len(ebins) - 1):
            if x >= ebins[j] and x <= ebins[j+1]:
                vals[j].append(y)
                break

    means = [np.mean(v) for v in vals]
    medians = [np.median(v) for v in vals]
    stds = [np.std(v) for v in vals]
    xs = []
    # get xs for means
    for i in range(len(ebins) - 1):
        ed = 10 ** (np.log10(ebins[i]) + (np.log10(ebins[i+1]) - np.log10(ebins[i])) / 2)
        xs.append(ed)

    # calculate histogram
    for i in range(int(len(hi) / 2)):
        hi[2*i] = means[i]
        hi[2*i + 1] = means[i]
    errs = np.zeros(len(stds))
    # standard error
    for i in range(len(errs)):
        errs[i] = stds[i] / np.sqrt(len(vals[i]))

    return np.array(hi), np.array(hibins), np.array(means), np.array(errs), np.array(xs), np.array(medians)



def create_npulsars(dr, bins, xval="Edot [ergs/s]"):
    # sort tables
    dr.sort(xval)
    min = dr[0][xval] # mi1
    max = dr[-1][xval] # ma2

    ebins = np.logspace(np.log10(min), np.log10(max), num=bins)
    hibins = []
    for i in range(len(ebins) - 1):
        hibins.append(ebins[i])
        hibins.append(ebins[i+1])
    hi = np.zeros(len(hibins))

    for row in dr:
        edot = row[xval]
        for i in range(len(hibins) - 1):
            if edot >= hibins[i] and edot <= hibins[i+1]:
                hi[i] += 1
                hi[i+1] += 1
                break
    return np.array(hi), np.array(hibins)


def p3_edot_fraction(data, data2, dr, nodr, bins=15):
    # data
    p3s = []
    ep3s = []
    edots = []
    ages = []

    p3s += list(data[0])
    ep3s += list(data[1])
    edots += list(data[9])
    ages += list(data[11])

    p3s2 = []
    ep3s2 = []
    edots2 = []
    ages2 = []

    p3s2 += list(data2[0])
    ep3s2 += list(data2[1])
    edots2 += list(data2[9])
    ages2 += list(data2[11])

    hibins, hifrac = create_hist(dr, nodr, bins)
    hi1, hib1, me1, err1, xs1, med1 = create_scatter(dr, bins, val="W10 [ms]")# val="W10 [ms]") #val="S1400 [mJy]")
    #hi1, hib1, me1, err1, xs1, med1 = create_scatter(dr, bins, val="bp_w50")# val="W10 [ms]") #val="S1400 [mJy]")
    #hi2, hib2, me2, err2, xs2, med2 = create_scatter(dr, bins, val="W50 [ms]") # P
    hi2, hib2, me2, err2, xs2, med2 = create_scatter(dr, bins, val="SNRclean")
    hi7, hib7 = create_npulsars(dr, bins)

    # changing period to polar cap radius
    #hi1 = 150 * np.power(hi1, -0.5)
    #me1 = 150 * np.power(me1, -0.5)
    #hi2 = 150 * np.power(hi2, -0.5)
    #me2 = 150 * np.power(me2, -0.5)

    #std1 = 150 * np.power(std1, -0.5)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["C0", "tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]
    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.91, top=0.99)
    sc = pl.scatter(edots, p3s, c=colors[0], s=5, zorder=1, label="drift", alpha=0.7)
    pl.errorbar(edots, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2, alpha=0.7)
    sc = pl.scatter(edots2, p3s2, c="C1", s=5, zorder=1, label="P3only", alpha=0.7)
    pl.errorbar(edots2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2, alpha=0.7)
    #sc = pl.scatter(ages, p3s, c=colors[0], s=5, zorder=1, label="drift")
    #pl.errorbar(ages, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2)
    #sc = pl.scatter(ages2, p3s2, c="C1", s=5, zorder=1, label="P3only")
    #pl.errorbar(ages2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2)

    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    #pl.xlabel("Age (yr)")
    pl.ylabel(r"$P_3$ in $P$")
    ax1 = fig.get_axes()[0]

    ax2 = ax1.twinx()
    ax2.plot(hibins, hifrac, c="C3", lw=2, ls="--", label="fraction")
    pl.semilogx()
    pl.ylabel("fraction of drifters")
    #ax2.set_yticks([])
    pl.legend(loc="upper center")

    ax3 = ax1.twinx()
    ax3.plot(hib1, hi1, c="tab:green", lw=2, ls="--", label="W10")
    pl.errorbar(xs1, me1, yerr=err1, c="tab:green", lw=0, marker="o", ms=5, elinewidth=1)
    pl.semilogx()
    #pl.ylabel("W10 [ms]")
    ax3.set_yticks([])
    pl.legend(loc="lower left")

    #ax4 = ax1.twinx()
    #ax4.plot(hib7, hi7, c="tab:grey", lw=2, ls="--", label="N-pulsars")
    #ax4.plot(hib2, hi2, c="tab:blue", lw=2, ls="--", label="SNR")
    #pl.errorbar(xs2, me2, yerr=err2, c="tab:blue", lw=0, marker="o", ms=5, elinewidth=1)
    #pl.loglog()
    #ax4.set_yticks([])
    #pl.legend(loc="upper left")

    filename = "output/p3_edot_frac.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3_age_fraction(data, data2, dr, nodr, bins=15):
    # data
    p3s = []
    ep3s = []
    edots = []
    ages = []

    p3s += list(data[0])
    ep3s += list(data[1])
    edots += list(data[9])
    ages += list(data[11])

    p3s2 = []
    ep3s2 = []
    edots2 = []
    ages2 = []

    p3s2 += list(data2[0])
    ep3s2 += list(data2[1])
    edots2 += list(data2[9])
    ages2 += list(data2[11])

    hibins, hifrac = create_hist(dr, nodr, bins, xval="Age [yr]")
    hi1, hib1, me1, err1, xs1, med1 = create_scatter(dr, bins, xval="Age [yr]", val="W10 [ms]")# val="W10 [ms]") #val="S1400 [mJy]")
    #hi2, hib2, me2, err2, xs2, med2 = create_scatter(dr, bins)# val="W50 [ms]") # P
    hi2, hib2, me2, err2, xs2, med2 = create_scatter(dr, bins, xval="Age [yr]", val="SNRclean")
    hi7, hib7 = create_npulsars(dr, bins, xval="Age [yr]")

    # changing period to polar cap radius
    #hi1 = 150 * np.power(hi1, -0.5)
    #me1 = 150 * np.power(me1, -0.5)
    #hi2 = 150 * np.power(hi2, -0.5)
    #me2 = 150 * np.power(me2, -0.5)

    #std1 = 150 * np.power(std1, -0.5)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["C0", "tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.91, top=0.99)
    sc = pl.scatter(ages, p3s, c=colors[0], s=5, zorder=1, label="drift", alpha=0.7)
    pl.errorbar(ages, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2, alpha=0.7)
    sc = pl.scatter(ages2, p3s2, c="C1", s=5, zorder=1, label="P3only", alpha=0.7)
    pl.errorbar(ages2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2, alpha=0.7)

    #sc = pl.scatter(ages, p3s, c=colors[0], s=5, zorder=1, label="drift")
    #pl.errorbar(ages, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2)
    #sc = pl.scatter(ages2, p3s2, c="C1", s=5, zorder=1, label="P3only")
    #pl.errorbar(ages2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2)

    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    pl.xlabel("Age (yr)")
    #pl.xlabel("Age (yr)")
    pl.ylabel(r"$P_3$ in $P$")
    ax = fig.get_axes()[0]

    ax2 = ax.twinx()
    ax2.plot(hibins, hifrac, c="C3", lw=2, ls="--", label="fraction")
    #pl.loglog()
    pl.ylabel("fraction of drifters")
    #ax2.set_yticks([])
    pl.legend(loc="upper left")

    ax3 = ax.twinx()
    ax3.plot(hib1, hi1, c="tab:green", lw=2, ls="--", label="W10")
    pl.errorbar(xs1, me1, yerr=err1, c="tab:green", lw=0, marker="o", ms=5, elinewidth=1)
    #pl.loglog()
    #pl.ylabel("W10 [ms]")
    ax3.set_yticks([])
    pl.legend(loc="lower left")

    ax4 = ax.twinx()
    #ax4.plot(hib7, hi7, c="tab:grey", lw=2, ls="--", label="N-pulsars")
    #ax4.plot(hib2, hi2, c="tab:blue", lw=2, ls="--", label="SNR")
    #pl.errorbar(xs2, me2, yerr=err2, c="tab:blue", lw=0, marker="o", ms=5, elinewidth=1)
    #pl.loglog()
    ax4.set_yticks([])
    #pl.legend(loc="upper left")

    filename = "output/p3_age_frac.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def edot_fraction(data, data2, dr, nodr, bins=15):
    # data
    p3s = []
    ep3s = []
    edots = []
    ages = []

    p3s += list(data[0])
    ep3s += list(data[1])
    edots += list(data[9])
    ages += list(data[11])

    p3s2 = []
    ep3s2 = []
    edots2 = []
    ages2 = []

    p3s2 += list(data2[0])
    ep3s2 += list(data2[1])
    edots2 += list(data2[9])
    ages2 += list(data2[11])

    hibins, hifrac = create_hist(dr, nodr, bins)
    hi3, hib3, me3, err3, xs3, med3 = create_scatter(dr, bins, val="W10 [ms]")# val="W10 [ms]") #val="S1400 [mJy]")
    #hi3, hib3, me3, err3, xs3, med3 = create_scatter(da.vstack([dr, nodr]), bins, val="W10 [ms]")# val="W10 [ms]") #val="S1400 [mJy]")
    #hi1, hib1, me1, err1, xs1, med1 = create_scatter(dr, bins, val="bp_w50")# val="W10 [ms]") #val="S1400 [mJy]")
    #hi5, hib5, me5, err5, xs5, med5 = create_scatter(dr, bins, val="W50 [ms]") # P
    #hi5, hib5, me5, err5, xs5, med5 = create_scatter(da.vstack([dr, nodr]), bins, val="W10 [ms]") # P
    hi5, hib5, me5, err5, xs5, med5 = create_scatter(dr, bins, val="bp_w50")# val="W10 [ms]") #val="S1400 [mJy]")
    #hi5, hib5, me5, err5, xs5, med5 = create_scatter(da.vstack([dr, nodr]), bins, val="bp_w50")# val="W10 [ms]") #val="S1400 [mJy]")

    #hi4, hib4 = create_npulsars(da.vstack([dr, nodr] ), bins)
    hi4, hib4, me4, err4, xs4, med4 = create_scatter(dr, bins, val="SNRclean")

    # changing period to polar cap radius
    #hi1 = 150 * np.power(hi1, -0.5)
    #me1 = 150 * np.power(me1, -0.5)
    #hi2 = 150 * np.power(hi2, -0.5)
    #me2 = 150 * np.power(me2, -0.5)

    #std1 = 150 * np.power(std1, -0.5)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["C0", "tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]
    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.10, bottom=0.15, right=0.89, top=0.99)
    """
    sc = pl.scatter(edots, p3s, c=colors[0], s=5, zorder=1, label="drift", alpha=0.7)
    pl.errorbar(edots, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2, alpha=0.7)
    sc = pl.scatter(edots2, p3s2, c="C1", s=5, zorder=1, label="P3only", alpha=0.7)
    pl.errorbar(edots2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2, alpha=0.7)
    #"""
    #sc = pl.scatter(ages, p3s, c=colors[0], s=5, zorder=1, label="drift")
    #pl.errorbar(ages, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2)
    #sc = pl.scatter(ages2, p3s2, c="C1", s=5, zorder=1, label="P3only")
    #pl.errorbar(ages2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2)

    #pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    #pl.xlabel("Age (yr)")
    pl.ylabel(r"$P_3$ in $P$")
    ax1 = fig.get_axes()[0]

    ax2 = ax1.twinx()
    ax2.plot(hibins, hifrac, c="C3", lw=3, ls="--", label="fraction")
    pl.ylabel("fraction of drifters")
    ax2.yaxis.label.set_color("C3")
    ax2.tick_params(colors="C3", which="both")
    #pl.loglog()
    #pl.ylim((0.07, 1.01))

    ax3 = ax1.twinx()
    ax3.plot(hib3, hi3, c="tab:green", lw=2, ls="--", label="W10")
    pl.errorbar(xs3, me3, yerr=err3, c="tab:green", lw=0, marker="o", ms=5, elinewidth=1)
    pl.plot(xs3, med3, c="tab:green", lw=0, marker="o", ms=3)
    pl.ylabel("W10 [phase]", labelpad=-32)
    ax3.yaxis.label.set_color("tab:green")
    ax3.tick_params(axis="y",direction="in", pad=-22, colors="tab:green", which="both")
    pl.loglog()
    pl.ylim((0.01, 0.6))
    pl.legend(loc="lower left")

    """
    ax4 = ax1.twinx()
    ax4.plot(hib4, hi4, c="tab:blue", lw=2, ls="--")
    pl.errorbar(xs4, me4, yerr=err4, c="tab:blue", lw=0, marker="o", ms=5, elinewidth=1)
    pl.loglog()
    ax4.set_yticks([])
    #pl.legend(loc="upper left")
    #"""

    ax5 = ax1.twinx()
    ax5.plot(hib5, hi5, c="tab:purple", lw=2, ls="--")
    pl.errorbar(xs5, me5, yerr=err5, c="tab:purple", lw=0, marker="o", ms=5, elinewidth=1)
    pl.plot(xs5, med5, c="tab:purple", lw=0, marker="o", ms=3)
    ax5.tick_params(axis="y",direction="in", pad=-25, colors="tab:purple", which="both", right=False, left=True, labelright=False, labelleft=True)
    pl.ylabel("W50 [phase] (B.P.)", labelpad=-35)
    ax5.yaxis.set_label_position("left")
    ax5.yaxis.label.set_color("tab:purple")
    pl.loglog()
    pl.ylim((0.01, 0.6))

    filename = "output/edot_fraction.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def edot_fraction_snr(data, data2, dr, nodr, bins=15):
    # data
    p3s = []
    ep3s = []
    edots = []
    ages = []

    p3s += list(data[0])
    ep3s += list(data[1])
    edots += list(data[9])
    ages += list(data[11])

    p3s2 = []
    ep3s2 = []
    edots2 = []
    ages2 = []

    p3s2 += list(data2[0])
    ep3s2 += list(data2[1])
    edots2 += list(data2[9])
    ages2 += list(data2[11])

    hibins, hifrac = create_hist(dr, nodr, bins)
    #hi3, hib3, me3, err3, xs3, med3 = create_scatter(dr, bins, val="W10 [ms]")# val="W10 [ms]") #val="S1400 [mJy]")
    #hi3, hib3, me3, err3, xs3, med3 = create_scatter(da.vstack([dr, nodr]), bins, val="W10 [ms]")# val="W10 [ms]") #val="S1400 [mJy]")
    #hi1, hib1, me1, err1, xs1, med1 = create_scatter(dr, bins, val="bp_w50")# val="W10 [ms]") #val="S1400 [mJy]")
    #hi5, hib5, me5, err5, xs5, med5 = create_scatter(dr, bins, val="W50 [ms]") # P
    #hi5, hib5, me5, err5, xs5, med5 = create_scatter(da.vstack([dr, nodr]), bins, val="W10 [ms]") # P
    #hi5, hib5, me5, err5, xs5, med5 = create_scatter(dr, bins, val="bp_w50")# val="W10 [ms]") #val="S1400 [mJy]")
    #hi5, hib5, me5, err5, xs5, med5 = create_scatter(da.vstack([dr, nodr]), bins, val="bp_w50")# val="W10 [ms]") #val="S1400 [mJy]")

    hi4, hib4, me4, err4, xs4, med4 = create_scatter(dr, bins, val="SNRclean")
    hi5, hib5 = create_npulsars(da.vstack([dr, nodr] ), bins)

    # changing period to polar cap radius
    #hi1 = 150 * np.power(hi1, -0.5)
    #me1 = 150 * np.power(me1, -0.5)
    #hi2 = 150 * np.power(hi2, -0.5)
    #me2 = 150 * np.power(me2, -0.5)

    #std1 = 150 * np.power(std1, -0.5)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    cmaps = [pl.cm.get_cmap('Reds'), pl.cm.get_cmap('Blues'), pl.cm.get_cmap('Greens')]
    colors = ["C0", "tab:red", "tab:blue", "tab:green"]
    #cmaps = [pl.cm.get_cmap('viridis'), pl.cm.get_cmap('inferno')]
    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.10, bottom=0.15, right=0.89, top=0.99)
    #"""
    sc = pl.scatter(edots, p3s, c=colors[0], s=5, zorder=1, label="drift", alpha=0.7)
    pl.errorbar(edots, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2, alpha=0.7)
    sc = pl.scatter(edots2, p3s2, c="C1", s=5, zorder=1, label="P3only", alpha=0.7)
    pl.errorbar(edots2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2, alpha=0.7)
    #"""
    #sc = pl.scatter(ages, p3s, c=colors[0], s=5, zorder=1, label="drift")
    #pl.errorbar(ages, p3s, fmt='none', yerr=ep3s, color=colors[0], zorder=2)
    #sc = pl.scatter(ages2, p3s2, c="C1", s=5, zorder=1, label="P3only")
    #pl.errorbar(ages2, p3s2, fmt='none', yerr=ep3s2, color="C1", zorder=2)
    #pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    #pl.xlabel("Age (yr)")
    pl.ylabel(r"$P_3$ in $P$")
    ax1 = fig.get_axes()[0]

    ax2 = ax1.twinx()
    ax2.plot(hibins, hifrac, c="C3", lw=3, ls="--", label="fraction")
    pl.ylabel("fraction of drifters")
    ax2.yaxis.label.set_color("C3")
    ax2.tick_params(colors="C3", which="both")
    pl.loglog()
    pl.ylim((0.07, 1.01))

    ax3 = ax1.twinx()
    ax3.plot(hib4, hi4, c="tab:green", lw=2, ls="--")
    pl.errorbar(xs4, me4, yerr=err4, c="tab:green", lw=0, marker="o", ms=5, elinewidth=1)
    pl.plot(xs4, med4, c="tab:green", lw=0, marker="o", ms=3)
    pl.axhline(y=2e3, ls=":", c="tab:green", lw=1)
    pl.ylabel("SNRclean", labelpad=-31)
    ax3.yaxis.label.set_color("tab:green")
    ax3.tick_params(axis="y",direction="in", pad=-22, colors="tab:green", which="both")
    pl.loglog()
    #pl.ylim((0.01, 0.6))
    #pl.legend(loc="lower left")


    ax5 = ax1.twinx()
    ax5.plot(hib5, hi5, c="tab:purple", lw=2, ls="--")
    ax5.tick_params(axis="y",direction="in", pad=-25, colors="tab:purple", which="both", right=False, left=True, labelright=False, labelleft=True)
    pl.ylabel("Number of pulsars", labelpad=-35)
    ax5.yaxis.set_label_position("left")
    ax5.yaxis.label.set_color("tab:purple")
    pl.loglog()
    #pl.ylim((0.01, 0.6))



    filename = "output/edot_fraction_snr.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def edot_cpp(data, data2, bins=15):
    # data
    p3s = []
    ep3s = []
    p3ws = []
    p3wserr = []
    edots = []
    ages = []

    p3s += list(data[0])
    ep3s += list(data[1])
    p3ws += list(data[2])
    p3wserr += list(data[3])
    edots += list(data[9])
    ages += list(data[11])

    p3s2 = []
    ep3s2 = []
    p3ws2 = []
    p3wserr2 = []
    edots2 = []
    ages2 = []

    p3s2 += list(data2[0])
    ep3s2 += list(data2[1])
    p3ws2 += list(data2[2])
    p3wserr2 += list(data2[3])
    edots2 += list(data2[9])
    ages2 += list(data2[11])

    #print(p3ws)
    #print(p3ws2)

    # create histograms
    hi1, hi1b, me1, err1, xs1, med1 = create_scatter_xy(edots, p3ws, bins)
    hi2, hi2b, me2, err2, xs2, med2 = create_scatter_xy(edots, p3wserr, bins)
    #hi2, hi2b, me2, err2, xs2, med2 = create_scatter_xy(edots2, p3wserr2, bins)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    fig = pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.13, bottom=0.15, right=0.99, top=0.99)
    #pl.loglog()
    pl.semilogx()
    #pl.errorbar(edots, p3ws, yerr=p3wserr, c="tab:red", lw=0, marker="o", ms=5, elinewidth=1)
    #pl.errorbar(edots2, p3ws2, yerr=p3wserr2, c="tab:blue", lw=0, marker="o", ms=5, elinewidth=1)
    pl.plot(hi1b, hi1, c="tab:red", lw=3, ls="--", label=r"$P_3$ width")
    pl.errorbar(xs1, me1, yerr=err1, c="tab:red", lw=0, marker="o", ms=5, elinewidth=1)
    pl.plot(xs1, med1, c="tab:red", lw=0, marker="o", ms=3)
    #"""
    pl.plot(hi2b, hi2, c="tab:blue", lw=3, ls="--", label=r"$P_3$ width error")
    pl.errorbar(xs2, me2, yerr=err2, c="tab:blue", lw=0, marker="o", ms=5, elinewidth=1)
    pl.plot(xs2, med2, c="tab:blue", lw=0, marker="o", ms=3)
    #"""
    #yl = pl.ylim()
    #pl.ylim([0.7, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    #pl.xlabel("Age (yr)")
    pl.ylabel(r"(cpp)")
    pl.legend()
    ax1 = fig.get_axes()[0]

    #ax2 = ax1.twinx()
    #ax2.plot(hibins, hifrac, c="C3", lw=3, ls="--", label="fraction")
    #pl.ylabel("fraction of drifters")
    #ax2.yaxis.label.set_color("C3")
    #ax2.tick_params(colors="C3", which="both")
    #pl.loglog()
    #pl.ylim((0.07, 1.01))

    """
    ax3 = ax1.twinx()
    ax3.plot(hib4, hi4, c="tab:green", lw=2, ls="--")
    pl.errorbar(xs4, me4, yerr=err4, c="tab:green", lw=0, marker="o", ms=5, elinewidth=1)
    pl.plot(xs4, med4, c="tab:green", lw=0, marker="o", ms=3)
    pl.axhline(y=2e3, ls=":", c="tab:green", lw=1)
    pl.ylabel("SNRclean", labelpad=-31)
    ax3.yaxis.label.set_color("tab:green")
    ax3.tick_params(axis="y",direction="in", pad=-22, colors="tab:green", which="both")
    pl.loglog()
    #pl.ylim((0.01, 0.6))
    #pl.legend(loc="lower left")

    ax5 = ax1.twinx()
    ax5.plot(hib5, hi5, c="tab:purple", lw=2, ls="--")
    ax5.tick_params(axis="y",direction="in", pad=-25, colors="tab:purple", which="both", right=False, left=True, labelright=False, labelleft=True)
    pl.ylabel("Number of pulsars", labelpad=-35)
    ax5.yaxis.set_label_position("left")
    ax5.yaxis.label.set_color("tab:purple")
    pl.loglog()
    #pl.ylim((0.01, 0.6))
    """

    filename = "output/edot_cpp.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()



def fit_p3_edot(data, a=-0.5, y0=0.7, edot_break=2.5e31, edot_min=2e27):
    """
        using two linear fits for low- and high-edot regions seperatly
    """
    # TODO repeat not implemented yet
    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])

    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))

    p3s_low = []
    ep3s_low = []
    edots_low = []

    p3s_high = []
    ep3s_high = []
    edots_high = []

    for i, edot in enumerate(edots):
        if edot <= edot_break:
            if edot > edot_min:
                p3s_low.append(p3s[i])
                ep3s_low.append(ep3s[i])
                edots_low.append(edots[i])
        else:
            p3s_high.append(p3s[i])
            ep3s_high.append(ep3s[i])
            edots_high.append(edots[i])

    p3s_low = np.array(p3s_low)
    ep3s_low = np.array(ep3s_low)
    edots_low = np.array(edots_low)
    p3s_high = np.array(p3s_high)
    ep3s_high = np.array(ep3s_high)
    edots_high = np.array(edots_high)


    # fitting linear dependence
    fun = lambda v, x: v[0] * x + v[1]
    v0 = [1, 1]
    xl, yl, vl = least_sq(np.log10(edots_low), np.log10(p3s_low), fun, v0, xmax=None)
    xh, yh, vh = least_sq(np.log10(edots_high), np.log10(p3s_high), fun, v0, xmax=None)
    # steeper # Patrick's idea
    #vh = [0.4, -12]
    #vh = [0.6, -18.5]
    #vh = [0.8, -24.5]
    #yh = fun(vh, xh)


    # continue line # Geoff's idea
    xf = np.array(list(xl) + list(xh))
    yh2 = fun(vh, xf)
    """
    for i,y in enumerate(yh2):
        if y < np.log10(2):
            yh2[i] = np.abs(1 / yh2[i])
            if yh2[i] > 1.5:
                yh2[i] = 1.5
    #"""

    #p3 = lambda p3obs, n: p3obs / (n* p3obs + 1)
    def p3(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            elif p3 == 1:
                return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                elif p3_ == 0:
                    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)


    yintr = np.empty(yh.size)
    n = 1
    for i, y in enumerate(yh):
        yintr[i] = np.log10(p3(10**y, n))

    para = lambda v, x: v[0] * x ** 2 + v[1] * x + v[2] + v[3] ** x **3
    v0 = [1,1,1, 1]
    xp, yp, vp = least_sq(xh, yintr, para, v0, xmax=None)
    #exit()

    b = y0 - a * np.log10(1e31)
    true = lambda x: a * x + b
    ytrue = true(np.log10(edots))

    p3ob_al10 = lambda p3true10, n: np.log10(np.abs(10 ** p3true10 / (1 - n * 10 ** p3true10)))


    # adding nsp dependence
    nsp_fun = lambda v, x: v[0] * x + v[1] + x ** 2 * v[2]
    v0 = [1, 1, 1]
    xpoints = np.array([29, 31, 33, 35, 36.4])
    #ypoints = np.array([60, 30, 15, 10, 2])
    #ypoints = np.array([5, 10, 15, 30, 45])
    #ypoints = np.array([10, 10, 10, 10, 10])
    #ypoints = np.array([20, 20, 20, 20, 20])
    ypoints = np.array([15, 10, 5, 3, 2])
    #ypoints = np.array([2, 2, 2, 2, 2])
    xnsp, ynsp, vnsp = least_sq(xpoints, ypoints, nsp_fun, v0, xmax=None)
    p3minimum = 1. / ynsp


    def genp3(edots, true, para, vp, nsp_fun, vnsp, sigma, edot_min):
        #sigma_old = sigma
        p3s_obs = []
        edots_obs = []
        p3s_intr = []
        edots_intr = []
        for edot in edots:
            x = np.log10(edot)
            p3i = true(x)
            p3p = para(vp, x)
            if p3i < p3p:
                p3 = p3p
            else:
                p3 = p3i
            p3m = p3

            """ # standard proposition
            sig = sigma
            p3 = np.random.normal(p3m, sig)
            while 10 ** p3 < 1 / nsp_fun(vnsp, x):
                p3 = np.random.normal(p3m, sig)
            #"""

            """ # Geoff's proposition
            if p3m < 1:
                sig = np.abs(p3) # sigma
            else:
                sig = sigma
            p3 = np.random.normal(p3m, sig)
            while 10 ** p3 < 1 / nsp_fun(vnsp, x):
                p3 = np.random.normal(p3m, sig)
            #"""

            #""" #  p3_lin just for tests, but should be fine
            if x < 32.5:
                sig = 1.0 * 10 ** p3
            else:
                sig = 0.1 * 10 ** p3
            p3_lin = np.abs(np.random.normal(10 ** p3, sig))
            while p3_lin == 0:
                p3_lin = np.abs(np.random.normal(10 ** p3, sig))
            p3 = np.log10(p3_lin)

            while 10 ** p3 < 1 / nsp_fun(vnsp, x):
                p3_lin = np.abs(np.random.normal(10 ** p3, sig))
                while p3_lin == 0:
                    p3_lin = np.abs(np.random.normal(10 ** p3, sig))
                p3 =  np.log10(p3_lin)
            #"""
            alias = False
            n = 1
            p3obs = p3
            while (10 ** p3obs < 2):
                p3obs = np.log10(np.abs(10 ** p3 / (1 - n * 10 ** p3)))
                """
                print("\t", 10 ** p3obs, n)
                if np.isnan(p3obs):
                    #p3 = np.random.normal(p3, sigma)
                    #n = 0
                    #p3obs = p3
                    print(10 ** p3)
                    print(10 ** p3obs, n)
                #print("p3", 10 ** p3)
                #print((1 - n * 10 ** p3))
                #"""
                n += 1
                alias = True
                #print(n)
            if alias is True:
                p3s_intr.append(10 ** p3)
                edots_intr.append(edot)

            p3s_obs.append(10 ** p3obs)
            edots_obs.append(edot)
        return np.array(p3s_obs), np.array(edots_obs), np.array(p3s_intr), np.array(edots_intr)

    sigma = np.std(np.log10(p3s_low))
    #sigma = sigma / (sigma + 1) # alised sigma?
    print("SIGMA: ", sigma)
    p3s_model, edots_model, p3s_true, edots_true = genp3(edots, true, para, vp, nsp_fun, vnsp, sigma, edot_min)


    # randomize observations
    #"""
    for zz in range(len(p3s10)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s10_rand[zz] = np.log10(ran)
        #p3s10_rand[zz] = p3s10[zz] # not random
    #"""

    pvals, xpvals, pvals_bin, xpvals_bin = da.calculate_divergance(p3s10_rand, np.log10(p3s_model), edots10, ep3s, sample=30)
    #pvals, xpvals, pvals_bin, xpvals_bin = calculate_divergance(p3s10, p3s_model, edots10, ep3s) # Note ep3s not ep3s10 (is it fine..?)

    divergence = len([pv for pv in pvals if pv < 0.05])
    print("Divergence: ", divergence)


    print("LOW: ", vl)
    print("HIGH: ", vh)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    fig = pl.figure(figsize=(7.086614, 1.3*4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    pl.subplot(2,1,1)
    pl.minorticks_on()
    sc = pl.scatter(edots_low, p3s_low, c="tab:green", s=5, zorder=1)
    pl.plot(10**xl, 10**yl, c="tab:green")
    sc = pl.scatter(edots_high, p3s_high, c="tab:blue", s=5, zorder=1)
    pl.plot(10**xh, 10**yh, c="tab:blue")
    pl.plot(10**xf, 10**yh2, c="tab:blue", lw=2, alpha=0.4)
    pl.plot(10**xh, 10**yintr, c="tab:red")
    pl.plot(10**xp, 10**yp, c="tab:orange", alpha=0.3, lw=3)
    pl.plot(10**xnsp, p3minimum, c="black", lw=1, ls="--")
    sc = pl.plot(edots, 10**ytrue, c="tab:pink")
    sc = pl.scatter(edots_model, p3s_model, c="tab:orange", s=5)
    sc = pl.scatter(edots_true, p3s_true, c="tab:grey", s=5)
    #pl.legend()
    pl.loglog()
    #pl.axhline(y=1)
    yl = pl.ylim()
    #pl.ylim([0.01, yl[1]])
    #pl.ylim([0.7, 50])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    xlims = pl.xlim()

    pl.subplot(2,1,2)
    pl.minorticks_on()
    #pl.scatter(xpvals, pvals)
    try:
        pl.plot(xpvals_bin, pvals_bin, lw=2, c="C1")
    except:
        pass
    pl.axhline(y=0.05, c="C2", ls="--")
    pl.xlim(np.log10(xlims[0]), np.log10(xlims[1]))

    filename = "output/fit_p3_edot.pdf"
    print(filename)
    pl.savefig(filename)
    pl.show()


def fit_p3_edot2(data, a=-0.7, y0=0.5):
#def fit_p3_edot2(data, a=0.3, y0=0.3):
    """
        using parabolic fit for the whole edot range
    """
    # TODO repeat not implemented yet
    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])

    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))

    # dependence in question
    b = y0 - a * np.log10(1e31)
    true = lambda x: a * x + b
    yt = true(np.log10(edots))

    # poarabolic fit to whole data
    para = lambda v, x: v[0] * x ** 2 + v[1] * x + v[2]
    xp, yp, vp = least_sq(np.log10(edots), np.log10(p3s), para, [1,1,1], xmax=None)

    # unalised values assuming all observations are aliased
    p3 = lambda p3obs, n: p3obs / (n * p3obs + 1)
    yintr = np.empty(yp.size)
    n = 1
    for i, y in enumerate(yp):
        yintr[i] = np.log10(p3(10**y, n))

    # polynominal fit to intrinsic p3
    poly = lambda v, x: v[0] * x **3 +  v[1] * x ** 2 + v[2] * x + v[3]
    xp2, yp2, vp2 = least_sq(xp, yintr, poly, [1,1,1,1], xmax=None)
    #print(vp2)

    # adding nsp dependence
    nsp_fun = lambda v, x: v[0] * x + v[1] + x ** 2 * v[2]
    v0 = [1, 1, 1]
    xpoints = np.array([29, 31, 33, 35, 36.4])
    #ypoints = np.array([60, 30, 15, 10, 2])
    #ypoints = np.array([5, 10, 15, 30, 45])
    #ypoints = np.array([10, 10, 10, 10, 10])
    #ypoints = np.array([20, 20, 20, 20, 20])
    ypoints = np.array([15, 10, 5, 3, 2])
    #ypoints = np.array([2, 2, 2, 2, 2])
    xnsp, ynsp, vnsp = least_sq(xpoints, ypoints, nsp_fun, v0, xmax=None)
    p3minimum = 1. / ynsp


    def genp3(edots, true, poly, vp, nsp_fun, vnsp, sigma):
        #sigma_old = sigma
        p3s_obs = []
        edots_obs = []
        p3s_intr = []
        edots_intr = []
        for edot in edots:
            x = np.log10(edot)
            p3i = true(x)
            p3p = poly(vp, x)
            p3 = np.max([p3i, p3p])
            p3m = p3
            sig = sigma
            #sig = np.abs(p3)
            p3 = np.random.normal(p3m, sig)
            while 10 ** p3 < 1 / nsp_fun(vnsp, x):
                p3 = np.random.normal(p3m, sig)

            #p3 = np.abs(np.random.normal(p3, sig))
            #if 10 ** p3 > 2:
            #    p3 = np.abs(np.random.normal(p3, sig))
            #else:
            #    p3 = np.random.normal(p3, sig)

            #if 10 ** p3 < 1:
            #    print(p3)
            #while 10 ** p3 == 1:
            #    if sig == 0:
            #        sig = 1
            #    p3 = np.abs(np.random.normal(p3, sig)) # HACK

            alias = False
            n = 1
            p3obs = p3
            while (10 ** p3obs < 2):
                p3obs = np.log10(np.abs(10 ** p3 / (1 - n * 10 ** p3)))
                n += 1
                alias = True
                #print(n)
            if alias is True:
                p3s_intr.append(10 ** p3)
                edots_intr.append(edot)

            p3s_obs.append(10 ** p3obs)
            edots_obs.append(edot)
        return np.array(p3s_obs), np.array(edots_obs), np.array(p3s_intr), np.array(edots_intr)

    sigma = np.std(np.log10(p3s))
    #sigma = sigma / (sigma + 1) # alised sigma?
    print("SIGMA: ", sigma)
    p3s_model, edots_model, p3s_true, edots_true = genp3(edots, true, poly, vp2, nsp_fun, vnsp, sigma)


    # randomize observations
    for zz in range(len(p3s10)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s10_rand[zz] = np.log10(ran)
        #p3s10_rand[zz] = p3s10[zz] # not random

    pvals, xpvals, pvals_bin, xpvals_bin = da.calculate_divergance(p3s10_rand, np.log10(p3s_model), edots10, ep3s, sample=30)
    #pvals, xpvals, pvals_bin, xpvals_bin = calculate_divergance(p3s10, p3s_model, edots10, ep3s) # Note ep3s not ep3s10 (is it fine..?)

    divergence = len([pv for pv in pvals if pv < 0.05])
    print("Divergence: ", divergence)



    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    fig = pl.figure(figsize=(7.086614, 1.3*4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    pl.subplot(2,1,1)
    pl.minorticks_on()
    sc = pl.scatter(edots, p3s, c="tab:blue", s=5, zorder=1)
    #pl.plot(10**xh, 10**yh, c="tab:blue")
    #pl.plot(10**xh, 10**yintr, c="tab:red")
    pl.plot(10**xp, 10**yp, c="tab:orange", alpha=0.3, lw=3)
    #pl.plot(10**xp, 10**yintr, c="tab:orange", alpha=0.3, lw=3)
    pl.plot(10**xp2, 10**yp2, c="tab:red", alpha=0.3, lw=2)
    pl.plot(10**xnsp, p3minimum, c="black", lw=1, ls="--")

    pl.plot(edots, 10**yt, c="tab:green")
    sc = pl.scatter(edots_model, p3s_model, c="tab:orange", s=5)
    sc = pl.scatter(edots_true, p3s_true, c="tab:grey", s=5)
    #pl.legend()
    pl.loglog()
    #pl.axhline(y=1)
    yl = pl.ylim()
    pl.ylim([0.01, yl[1]])
    #pl.ylim([0.7, 50])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    xlims = pl.xlim()

    pl.subplot(2,1,2)
    pl.minorticks_on()
    #pl.scatter(xpvals, pvals)
    try:
        pl.plot(xpvals_bin, pvals_bin, lw=2, c="C1")
    except:
        pass
    pl.axhline(y=0.05, c="C2", ls="--")
    pl.xlim(np.log10(xlims[0]), np.log10(xlims[1]))

    filename = "output/fit_p3_edot2.pdf"
    print(filename)
    pl.savefig(filename)
    pl.show()



def fit_p3_edot3(data, a=-0.77, y0=0.5):
#def fit_p3_edot2(data, a=0.3, y0=0.3):
    """
        using only evolution of number of sparks
    """
    # TODO repeat not implemented yet
    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])

    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))

    # dependence in question
    b = y0 - a * np.log10(1e31)
    p3fun = lambda x: a * x + b
    yt = p3fun(np.log10(edots))

    # adding nsp dependence
    nsp_fun = lambda v, x: v[0] * x + v[1] + x ** 2 * v[2]
    v0 = [1, 1, 1]
    xpoints = np.array([29, 31, 33, 35, 36.4])
    #ypoints = np.array([30, 20, 6, 4, 2])
    #ypoints = np.array([5, 10, 15, 30, 45])
    #ypoints = np.array([10, 10, 10, 10, 10])
    #ypoints = np.array([20, 20, 20, 20, 20])
    #ypoints = np.array([25, 15, 10, 3, 2])
    #ypoints = np.array([1, 1, 1, 1, 1])
    #ypoints = np.array([2, 2, 2, 2, 2])
    ypoints = np.array([15, 10, 5, 3, 2])
    xnsp, ynsp, vnsp = least_sq(xpoints, ypoints, nsp_fun, v0, xmax=None)
    p3minimum = 1. / ynsp


    def genp3(edots, p3fun, nsp_fun, vnsp, sigma):
        #sigma_old = sigma
        p3s_obs = []
        edots_obs = []
        p3s_intr = []
        edots_intr = []
        for edot in edots:
            x = np.log10(edot)
            p3 = p3fun(x)
            p3min =  np.log10(1 / nsp_fun(vnsp, x))
            p3m = max([p3, p3min])
            if p3 > p3min:
                sig = 1.0 * sigma
            else:
                sig = 1.0 * sigma
            #sig = np.max([np.abs(p3m), 0.1])
            #sig = np.abs(p3m - np.log10(10**p3m * 0.7))
            #print(sig, "p3m ", p3m)
            p3 = np.random.normal(p3m, sig)
            while 10 ** p3 < 1 / nsp_fun(vnsp, x):
                p3 = np.random.normal(p3m, sig)
            alias = False
            n = 1
            p3obs = p3
            while (10 ** p3obs < 2):
                p3obs = np.log10(np.abs(10 ** p3 / (1 - n * 10 ** p3)))
                n += 1
                alias = True
                #print(n)
            if alias is True:
                p3s_intr.append(10 ** p3)
                edots_intr.append(edot)
            p3s_obs.append(10 ** p3obs)
            edots_obs.append(edot)
        return np.array(p3s_obs), np.array(edots_obs), np.array(p3s_intr), np.array(edots_intr)

    sigma = np.std(np.log10(p3s))
    #sigma = sigma / (sigma + 1) # alised sigma?
    print("SIGMA: ", sigma)
    p3s_model, edots_model, p3s_true, edots_true = genp3(edots, p3fun, nsp_fun, vnsp, sigma)


    # randomize observations
    for zz in range(len(p3s10)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s10_rand[zz] = np.log10(ran)
        #p3s10_rand[zz] = p3s10[zz] # not random

    pvals, xpvals, pvals_bin, xpvals_bin = da.calculate_divergance(p3s10_rand, np.log10(p3s_model), edots10, ep3s, sample=30)
    #pvals, xpvals, pvals_bin, xpvals_bin = calculate_divergance(p3s10, p3s_model, edots10, ep3s) # Note ep3s not ep3s10 (is it fine..?)

    divergence = len([pv for pv in pvals if pv < 0.05])
    print("Divergence: ", divergence)


    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    fig = pl.figure(figsize=(7.086614, 1.3*4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    pl.subplot(2,1,1)
    pl.minorticks_on()
    sc = pl.scatter(edots, p3s, c="tab:blue", s=5, zorder=1)
    #pl.plot(10**xh, 10**yh, c="tab:blue")
    #pl.plot(10**xh, 10**yintr, c="tab:red")
    #pl.plot(10**xp, 10**yintr, c="tab:orange", alpha=0.3, lw=3)
    pl.plot(10**xnsp, p3minimum, c="black", lw=1, ls="--")

    pl.plot(edots, 10**yt, c="tab:green")
    sc = pl.scatter(edots_model, p3s_model, c="tab:orange", s=5)
    sc = pl.scatter(edots_true, p3s_true, c="tab:grey", s=5)
    #pl.legend()
    pl.loglog()
    #pl.axhline(y=1)
    yl = pl.ylim()
    pl.ylim([0.01, yl[1]])
    #pl.ylim([0.7, 50])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    xlims = pl.xlim()

    pl.subplot(2,1,2)
    pl.minorticks_on()
    #pl.scatter(xpvals, pvals)
    try:
        pl.plot(xpvals_bin, pvals_bin, lw=2, c="C1")
    except:
        pass
    pl.axhline(y=0.05, c="C2", ls="--")
    pl.xlim(np.log10(xlims[0]), np.log10(xlims[1]))

    filename = "output/fit_p3_edot3.pdf"
    print(filename)
    pl.savefig(filename)
    pl.show()


def fit_p3_edot4(data):
    """
        using Geoff's 1 - 1/x dependence
    """
    # TODO repeat not implemented yet
    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])

    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))


    # alised values
    #p3aliased = lambda p3, n: p3 / (1  - n * p3)
    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            elif p3 == 1:
                return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                elif p3_ == 0:
                    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    """ #not needed! just for tests!
    def p3aliased(p3, n):
        # n is ignored here!
        if type(p3) == np.float64 or type(p3) == float:
            p3obs = 1 / (1 - 1 / p3)
            if p3obs < 2:
                p3obs = 2
            return p3obs
        else:
            res = []
            for p3_ in p3:
                p3obs = 1 / (1 - 1 / p3_)
                if p3obs < 2:
                    p3obs = 2
                    res.append(p3_obs)
            return np.array(res)
    #"""

    # dependence in question
    #p3fun = lambda x, w, th: np.fabs(1 / (1 - w / (x - th)))
    def p3fun(xs, w, th):
        if type(xs) == np.float64:
            if xs != th:
                return np.fabs(1 / (1 - w / (xs - th)))
            else:
                return 1.
        else:
            res = []
            for x in xs:
                if x != th:
                    res.append(np.fabs(1 / (1 - w / (x - th))))
                else:
                    res.append(1)
            return np.array(res)

    # observed values
    def p3_obs(v, edot10):
        w = v[0]
        th = v[1]
        ym = p3fun(edot10, w, th)
        y = ym
        #print(edot10, ym)
        n = 1
        if ym < 1e-3:
            return 1e-3
        while y < 2:
            #print(edot10, ym, p3aliased(ym, n), n)
            y = np.abs(p3aliased(ym, n))
            #print(y, n)
            n += 1
        return y

    # fitting
    def errfunc(v, xs, ys):
        size = len(xs)
        diff = np.empty(size)
        for i in range(size):
            diff[i] = p3_obs(v, xs[i]) - ys[i]
        return diff


    #""" # FITTING DOES NOT WORK check commented plot below
    #w = 1.8
    #th = 30.5 - w
    w = 0.7
    th = 30.06
    v0 = [w, th]
    #vobs = leastsq(errfunc, v0, args=(edots10, p3s), maxfev=10000, full_output=True)[0]
    vobs = v0
    print(vobs)
    xt = np.logspace(np.log10(28), np.log10(37), num=1000)
    yt = p3fun(xt, vobs[0], vobs[1])
    #"""

    yo = np.empty(xt.size)
    for i, x in enumerate(xt):
        yo[i] = p3_obs(vobs, x)

    """
    szi = 50
    szj = 50

    ws = np.linspace(0, 2, num=szi)
    ths = np.linspace(29, 32, num=szj)

    diffs = np.empty([szi, szj])

    for i in range(szi):
        for j in range(szj):
            w = ws[i]
            th = ths[j]
            v = [w, th]
            diff = errfunc(v, edots10, p3s)
            diffs[i, j] = np.log10(np.abs(np.sum(diff)))
            #diffs[i, j] = np.abs(np.sum(diff))
            #print("w ", w, " th ", th,  )

    pl.imshow(diffs, origin="lower", extent=[ths[0], ths[-1], ws[0], ws[-1]])
    pl.colorbar()
    pl.show()
    #"""



    """
    # test function
    #xx = np.logspace(np.log10(1e29), np.log10(1e39), num=10000)
    #xx = np.logspace(np.log10(1), np.log10(30), num=1000)
    xx = np.linspace(-4, 30, num=1000)
    yx = np.empty(xx.size)
    for i, x in enumerate(xx):
        yx[i] = p3_obs([1, 2], x)
        #yx[i] = p3fun(x, 1, -3)
    #print(xx)
    #print(yx)

    print(xx)
    print(yx)

    #xx = np.logspace(np.log10(1e28), np.log10(1e37), num=1000)

    #pl.plot(xx-1e44, yx)
    pl.plot(xx, yx)
    pl.loglog()
    #pl.semilogy()
    pl.ylim(-1, 100)
    pl.show()
    return
    #"""


    """
    # test function 2
    #xx = np.linspace(27, 38)
    xx = np.linspace(0.7, 30, num=100)
    yx = np.empty(xx.size)
    for i, x in enumerate(xx):
        yx[i] = p3_obs([1, 7], x)

    #pl.plot(xx-1e44, yx)
    pl.plot(xx, yx)
    #pl.scatter(edots10, p3s)
    pl.loglog()
    #pl.semilogy()
    #pl.ylim(-1, 100)
    pl.show()
    return
    #"""

    #"""
    # test function 3
    f = lambda x, u: 1 / (1 - u[0] / (x - u[1])) #p3
    g = lambda x, u: 1 / (1 - 1 / f(x, u)) # aliased value
    x = np.linspace(1, 10, num=100)
    u = [1, 1]
    fs = f(x, u)
    gs = g(x, u)
    pl.subplot(2,1,1)
    pl.plot(x, fs)
    pl.plot(x, gs)
    pl.subplot(2,1,2)
    pl.plot(x, fs)
    pl.plot(x, gs)
    pl.loglog()
    pl.show()
    return
    # LINE in log-log scale only for u[1]=0!
    #"""


    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    fig = pl.figure(figsize=(7.086614, 1.3*4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    pl.subplot(2,1,1)
    pl.minorticks_on()
    sc = pl.scatter(edots, p3s, c="tab:blue", s=5, zorder=1)
    #pl.plot(10**xh, 10**yh, c="tab:blue")
    #pl.plot(10**xh, 10**yintr, c="tab:red")
    #pl.plot(10**xp, 10**yintr, c="tab:orange", alpha=0.3, lw=3)
    lt, = pl.plot(10**xt, yt, c="black", lw=1, ls="-")
    lo, = pl.plot(10**xt, yo, c="red", lw=2, ls="-", alpha=0.4)
    #pl.legend()
    pl.loglog()
    #pl.axhline(y=1)
    yl = pl.ylim()
    #pl.ylim([0.01, yl[1]])
    pl.ylim([0.7, 300])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    xlims = pl.xlim()

    pl.figtext(0.5, 0.4, r"$f(x) = \frac{1}{1-\frac{w}{x-t}}$", size=15)

    # Sliders
    axcolor = 'lightgoldenrodyellow'
    wax = pl.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    thax = pl.axes([0.25, 0.3, 0.65, 0.03], facecolor=axcolor)

    sw = Slider(wax, 'w', 0.0, 3, valinit=w, valstep=0.01)
    sth = Slider(thax, 't', 27, 32, valinit=th, valstep=0.01)
    def update(val):
        w = sw.val
        th = sth.val
        yt = p3fun(xt, w, th)
        yo = np.empty(xt.size)
        for i, x in enumerate(xt):
            yo[i] = p3_obs([w, th], x)
        lt.set_ydata(yt)
        lo.set_ydata(yo)
        fig.canvas.draw_idle()
        #print(w, th)

    sw.on_changed(update)
    sth.on_changed(update)

    filename = "output/fit_p3_edot4.pdf"
    print(filename)
    pl.savefig(filename)
    pl.show()


def fit_p3_edot5(data):
    """
        using Geoff's 1 - 1/x dependence, with Xiaoxi's correlation
    """
    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])

    xi_xs = []
    for i, p in enumerate(periods):
        val = p ** (-1.78) * (pdots[i]/1e-15)
        xi_xs.append(val)
    xi_xs = np.array(xi_xs)


    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))

    # alised values
    #p3aliased = lambda p3, n: p3 / (1  - n * p3)
    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            elif p3 == 1:
                return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                elif p3_ == 0:
                    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)


    # dependence in question
    #p3fun = lambda x, w, th: np.fabs(1 / (1 - w / (x - th)))
    def p3fun(xs, w, th):
        if type(xs) == np.float64:
            if xs != th:
                return np.fabs(1 / (1 - w / (xs - th)))
                #return np.fabs(1 / (1 - 1/ xs))
            else:
                return 1.
        else:
            res = []
            for x in xs:
                if x != th:
                    res.append(np.fabs(1 / (1 - w / (x - th))))
                    #res.append(np.fabs(1 / (1 - 1 / x)))
                else:
                    res.append(1)
            return np.array(res)

    # observed values
    def p3_obs(v, edot10):
        w = v[0]
        th = v[1]
        ym = p3fun(edot10, w, th)
        y = ym
        #print(edot10, ym)
        n = 1
        if ym < 1e-3:
            return 1e-3
        while y < 2:
            #print(edot10, ym, p3aliased(ym, n), n)
            y = np.abs(p3aliased(ym, n))
            #print(y, n)
            n += 1
        return y

    w = 0.25
    th = -1.1
    vobs = [w, th]
    xt = np.linspace(-3, 4, num=1000)
    yt = p3fun(xt, vobs[0], vobs[1])

    yo = np.empty(xt.size)
    for i, x in enumerate(xt):
        yo[i] = p3_obs(vobs, x)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    fig = pl.figure(figsize=(7.086614, 1.3*4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    pl.subplot(2,1,1)
    pl.minorticks_on()
    sc = pl.scatter(xi_xs, p3s, c="tab:blue", s=5, zorder=1)
    #pl.plot(10**xh, 10**yh, c="tab:blue")
    #pl.plot(10**xh, 10**yintr, c="tab:red")
    #pl.plot(10**xp, 10**yintr, c="tab:orange", alpha=0.3, lw=3)
    lt, = pl.plot(10**xt, yt, c="black", lw=1, ls="-")
    lo, = pl.plot(10**xt, yo, c="red", lw=2, ls="-", alpha=0.4)
    #pl.legend()
    pl.loglog()
    #pl.axhline(y=1)
    yl = pl.ylim()
    #pl.ylim([0.01, yl[1]])
    pl.ylim([0.7, 300])
    pl.xlabel(r"$(P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15} {\rm s/s })$")
    pl.ylabel(r"$P_3$ in $P$")
    xlims = pl.xlim()

    pl.figtext(0.5, 0.4, r"$f(x) = \frac{1}{1-\frac{w}{x-t}}$", size=15)

    # Sliders
    axcolor = 'lightgoldenrodyellow'
    wax = pl.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    thax = pl.axes([0.25, 0.3, 0.65, 0.03], facecolor=axcolor)

    sw = Slider(wax, 'w', 0.0, 3, valinit=w, valstep=0.01)
    sth = Slider(thax, 't', -7, 4, valinit=th, valstep=0.01)
    def update(val):
        w = sw.val
        th = sth.val
        yt = p3fun(xt, w, th)
        yo = np.empty(xt.size)
        for i, x in enumerate(xt):
            yo[i] = p3_obs([w, th], x)
        lt.set_ydata(yt)
        lo.set_ydata(yo)
        fig.canvas.draw_idle()
        #print(w, th)

    sw.on_changed(update)
    sth.on_changed(update)

    filename = "output/fit_p3_edot5.pdf"
    print(filename)
    pl.savefig(filename)
    pl.show()



def fit_p3_edot6(data):
    """
        using Geoff's 1 - 1/x dependence, with Xiaoxi's correlation fixing dep.
    """
    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])

    xi_xs = []
    for i, p in enumerate(periods):
        val = p ** (-1.78) * (pdots[i]/1e-15)
        xi_xs.append(val)
    xi_xs = np.array(xi_xs)

    """
    # YES this is the JiuJitsu way :D
    f = open("data/p3_ppdotedot.txt", "w")
    for i in range(len(p3s)):
        f.write("{} {} {} {} {}\n".format(p3s[i], ep3s[i], periods[i], pdots[i], edots[i]))
    f.close()
    return
    #"""

    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))

    # alised values
    #p3aliased = lambda p3, n: p3 / (1  - n * p3)
    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            elif p3 == 1:
                return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                elif p3_ == 0:
                    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)


    # dependence in question
    #p3fun = lambda x, w, th: np.fabs(1 / (1 - w / (x - th)))
    def p3fun(xs, w, th, r):
        if type(xs) == np.float64:
            if xs != th:
                return np.fabs(1 / (1 - (w / (xs - th))**r))
                #return np.fabs(1 / (1 - 1/ xs))
            else:
                return 1.
        else:
            res = []
            for x in xs:
                if x != th:
                    res.append(np.fabs(1 / (1 - (w / (x - th))**r)))
                    #res.append(np.fabs(1 / (1 - 1 / x)))
                else:
                    res.append(1)
            return np.array(res)

    # observed values
    def p3_obs(v, edot10):
        w = v[0]
        th = v[1]
        r = v[2]
        ym = p3fun(edot10, w, th, r)
        y = ym
        #print(edot10, ym)
        n = 1
        if ym < 1e-3:
            return 1e-3
        while y < 2:
            #print(edot10, ym, p3aliased(ym, n), n)
            y = np.abs(p3aliased(ym, n))
            #print(y, n)
            n += 1
        return y

    w = 1 #1 # 0.5 # 6.32 # 9.65
    th = -1.05 #-1 # -0.42 # -3.16 # -3.16
    r = 1
    vobs = [w, th, r]
    xt = np.logspace(np.log10(1e-3), np.log10(1e4), num=100)
    yt = p3fun(xt, vobs[0], vobs[1], vobs[2])

    yo = np.empty(xt.size)
    for i, x in enumerate(xt):
        yo[i] = p3_obs(vobs, x)

    xor = np.logspace(np.log10(1), np.log10(1e4), num=10000)
    yor = p3fun(xor, 1, 0, 1)
    yor2 = np.empty(xor.size)
    for i, x in enumerate(xor):
        yor2[i] = p3_obs([1,0, 2], x)

    # fitting line
    f, v = da.fit_line()

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    fig = pl.figure(figsize=(7.086614, 1.3*4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    pl.subplot(2,1,1)
    pl.minorticks_on()
    sc = pl.scatter(xi_xs, p3s, c="tab:blue", s=5, zorder=1)
    #pl.plot(10**xh, 10**yh, c="tab:blue")
    #pl.plot(10**xh, 10**yintr, c="tab:red")
    #pl.plot(10**xp, 10**yintr, c="tab:orange", alpha=0.3, lw=3)
    #lt, = pl.plot(10**xt, yt, c="black", lw=1, ls="-")
    #lo, = pl.plot(10**xt, yo, c="red", lw=2, ls="-", alpha=0.4)
    lt, = pl.plot(xt, yt, c="black", lw=1, ls="-")
    lo, = pl.plot(xt, yo, c="red", lw=2, ls="-", alpha=0.4)
    pl.plot(xor, yor, c="grey", lw=1, ls="--")
    pl.plot(xor, yor2, c="red", lw=1, ls="--")
    #pl.legend()
    pl.loglog()
    #pl.axhline(y=1)
    yl = pl.ylim()
    #pl.ylim([0.01, yl[1]])
    pl.ylim([0.2, 300])
    pl.xlabel(r"$(P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15} {\rm s/s })$")
    pl.ylabel(r"$P_3$ in $P$")
    xlims = pl.xlim()

    pl.figtext(0.5, 0.4, r"$f(x) = \frac{1}{1-\frac{w}{(x-t)^r}}$", size=15)

    # Sliders
    axcolor = 'lightgoldenrodyellow'
    wax = pl.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    thax = pl.axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)
    rax = pl.axes([0.25, 0.3, 0.65, 0.03], facecolor=axcolor)

    sw = Slider(wax, 'w', 0, 10, valinit=w, valstep=0.01)
    sth = Slider(thax, 't', -10, 20, valinit=th, valstep=0.01)
    sr = Slider(rax, 'r', 0.1, 4, valinit=r, valstep=0.01)
    def update(val):
        th = sth.val
        r = sr.val
        # TODO automatic w based on fit
        w = sw.val
        #w = f(v, th)
        #sw.set_val(w)
        #print(w, th)
        yt = p3fun(xt, w, th, r)
        yo = np.empty(xt.size)
        for i, x in enumerate(xt):
            yo[i] = p3_obs([w, th, r], x)
        lt.set_ydata(yt)
        lo.set_ydata(yo)
        fig.canvas.draw_idle()
        #print(w, th)

    sw.on_changed(update)
    sth.on_changed(update)
    sr.on_changed(update)

    filename = "output/fit_p3_edot6.pdf"
    print(filename)
    pl.savefig(filename)
    pl.show()


def fit_p3_edot7(data, a=-0.9, y0=0.7):
    """
        NO GOOD result! check fit_p3_edot6
        using parabolic dependence derived from linear fit in low-edot region and intrinsic p3 values in high-edot region
    """
    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])

    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))

    #p3 = lambda p3obs, n: p3obs / (n* p3obs + 1)
    # intrinsic P3 based on observed P3
    def p3_intr(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            elif p3 == 1:
                return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                elif p3_ == 0:
                    res.append(100)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    b = y0 - a * np.log10(1e31)
    p3lin = lambda x: a * x + b
    ylin = p3lin(np.log10(edots))

    # finding parabolic dependence...
    # break point
    edotbr = (np.log10(2) - b) / a
    # high-edots / aliased
    xl = []
    yl = []
    xh = []
    yh = []
    for i,edot in enumerate(edots10):
        if edot >= edotbr:
            xh.append(edot)
            yh.append(p3s10[i])
        else:
            xl.append(edot)
            yl.append(p3s10[i])
    xh = np.array(xh)
    yh = np.array(yh)
    xl = np.array(xl)
    yl = np.array(yl)
    yl_lin = p3lin(xl)
    # linear dependence for high-edot
    linear = lambda v, x: v[0] * x + v[1]
    v0 = [1, 1]
    xh2, yh2, vp = least_sq(xh, yh, linear, v0, xmax=None)
    # intrinsic p3 for high edot - based on linear dependence
    n = 1
    yintr = np.empty(xh2.size)
    for i,x in enumerate(xh2):
        yintr[i] = np.log10(p3_intr(10**yh2[i], n))
    # parabolic fit for the whole sample

    #parabolic = lambda v, x: v[0] * x ** 2 + v[1] * x + v[2]   # +  v[3] * x **3 + v[4] * x**4 + v[5] * x**5 + v[6] * x**6 + v[7] * x**7 + v[8] * x**8 + v[9] * x**9 + v[10] * x**10 + v[11] * x**11
    parabolic = lambda v, x: np.exp(- v[0] * (x ** v[2] - v[1])) # best?
    #parabolic = lambda v, x: v[0] ** (x ** v[2] - v[1])
    #parabolic = lambda v, x: v[3] * v[0] ** (v[2]* x - v[1])
    #parabolic = lambda v, x: np.exp(-x ** v[0] - v[1]) # best?

    #v0 = [1,33,1, 1,1, 1, 1, 1, 1, 1, 1, 1]
    #v0 = [1,33,1, 1,1]
    #v0 = [0.2, 33, 1]
    v0 = [1, 33, 1,1, 1]

    xp, yp, vp = least_sq(np.array(list(xl) + list(xh2)), np.array(list(yl_lin) + list(yintr)), parabolic, v0, xmax=None)
    print(vp)

    pl.scatter(xl, yl)
    pl.plot(xl, yl_lin)
    pl.scatter(xh, yh)
    pl.plot(xh2, yh2)
    pl.plot(xh2, yintr)
    pl.plot(xp, yp, lw=3)
    #pl.plot(xp, parabolic(vp, xp), lw=3)
    pl.show()


    return

    # adding nsp dependence
    nsp_fun = lambda v, x: v[0] * x + v[1] + x ** 2 * v[2]
    v0 = [1, 1, 1]
    xpoints = np.array([29, 31, 33, 35, 36.4])
    ypoints = np.array([15, 10, 5, 3, 2])
    xnsp, ynsp, vnsp = least_sq(xpoints, ypoints, nsp_fun, v0, xmax=None)
    p3minimum = 1. / ynsp

    #p3ob_al10 = lambda p3true10, n: np.log10(np.abs(10 ** p3true10 / (1 - n * 10 ** p3true10)))

    def genp3(edots, true, para, vp, nsp_fun, vnsp, sigma, edot_min):
        p3s_obs = []
        edots_obs = []
        p3s_intr = []
        edots_intr = []
        for edot in edots:
            x = np.log10(edot)
            p3i = true(x)
            p3p = para(vp, x)
            if p3i < p3p:
                p3 = p3p
            else:
                p3 = p3i
            p3m = p3

            sig = sigma
            p3 = np.random.normal(p3m, sig)
            while 10 ** p3 < 1 / nsp_fun(vnsp, x):
                p3 = np.random.normal(p3m, sig)


            if x < 32.5:
                sig = 1.0 * 10 ** p3
            else:
                sig = 0.1 * 10 ** p3
            p3_lin = np.abs(np.random.normal(10 ** p3, sig))
            while p3_lin == 0:
                p3_lin = np.abs(np.random.normal(10 ** p3, sig))
            p3 = np.log10(p3_lin)

            while 10 ** p3 < 1 / nsp_fun(vnsp, x):
                p3_lin = np.abs(np.random.normal(10 ** p3, sig))
                while p3_lin == 0:
                    p3_lin = np.abs(np.random.normal(10 ** p3, sig))
                p3 =  np.log10(p3_lin)
            alias = False
            n = 1
            p3obs = p3
            while (10 ** p3obs < 2):
                p3obs = np.log10(np.abs(10 ** p3 / (1 - n * 10 ** p3)))
                print("\t", 10 ** p3obs, n)
                if np.isnan(p3obs):
                    #p3 = np.random.normal(p3, sigma)
                    #n = 0
                    #p3obs = p3
                    print(10 ** p3)
                    print(10 ** p3obs, n)
                #print("p3", 10 ** p3)
                #print((1 - n * 10 ** p3))
                n += 1
                alias = True
                #print(n)
            if alias is True:
                p3s_intr.append(10 ** p3)
                edots_intr.append(edot)

            p3s_obs.append(10 ** p3obs)
            edots_obs.append(edot)
        return np.array(p3s_obs), np.array(edots_obs), np.array(p3s_intr), np.array(edots_intr)

    sigma = np.std(np.log10(p3s))
    print("SIGMA: ", sigma)
    p3s_model, edots_model, p3s_true, edots_true = genp3(edots, true, para, vp, nsp_fun, vnsp, sigma, edot_min)

    # randomize observations
    #"""
    for zz in range(len(p3s10)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s10_rand[zz] = np.log10(ran)
        #p3s10_rand[zz] = p3s10[zz] # not random
    #"""

    pvals, xpvals, pvals_bin, xpvals_bin = da.calculate_divergance(p3s10_rand, np.log10(p3s_model), edots10, ep3s, sample=30)
    #pvals, xpvals, pvals_bin, xpvals_bin = calculate_divergance(p3s10, p3s_model, edots10, ep3s) # Note ep3s not ep3s10 (is it fine..?)

    divergence = len([pv for pv in pvals if pv < 0.05])
    print("Divergence: ", divergence)


    print("LOW: ", vl)
    print("HIGH: ", vh)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    fig = pl.figure(figsize=(7.086614, 1.3*4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.15, right=0.99, top=0.99)
    pl.subplot(2,1,1)
    pl.minorticks_on()
    sc = pl.scatter(edots_low, p3s_low, c="tab:green", s=5, zorder=1)
    pl.plot(10**xl, 10**yl, c="tab:green")
    sc = pl.scatter(edots_high, p3s_high, c="tab:blue", s=5, zorder=1)
    pl.plot(10**xh, 10**yh, c="tab:blue")
    pl.plot(10**xf, 10**yh2, c="tab:blue", lw=2, alpha=0.4)
    pl.plot(10**xh, 10**yintr, c="tab:red")
    pl.plot(10**xp, 10**yp, c="tab:orange", alpha=0.3, lw=3)
    pl.plot(10**xnsp, p3minimum, c="black", lw=1, ls="--")
    sc = pl.plot(edots, 10**ytrue, c="tab:pink")
    sc = pl.scatter(edots_model, p3s_model, c="tab:orange", s=5)
    sc = pl.scatter(edots_true, p3s_true, c="tab:grey", s=5)
    #pl.legend()
    pl.loglog()
    #pl.axhline(y=1)
    yl = pl.ylim()
    #pl.ylim([0.01, yl[1]])
    #pl.ylim([0.7, 50])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    xlims = pl.xlim()

    pl.subplot(2,1,2)
    pl.minorticks_on()
    #pl.scatter(xpvals, pvals)
    try:
        pl.plot(xpvals_bin, pvals_bin, lw=2, c="C1")
    except:
        pass
    pl.axhline(y=0.05, c="C2", ls="--")
    pl.xlim(np.log10(xlims[0]), np.log10(xlims[1]))

    filename = "output/fit_p3_edot7.pdf"
    print(filename)
    pl.savefig(filename)
    pl.show()


def check_p3edot_geoff(data):

    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])

    xi_xs = []
    for i, p in enumerate(periods):
        val = p ** (-1.78) * (pdots[i]/1e-15)
        xi_xs.append(val)
    xi_xs = np.array(xi_xs)

    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))
    p3s_rand = np.empty(len(p3s))

    edot_break = 7e30
    p3s_low = []
    ep3s_low = []
    edots_low = []

    p3s_high = []
    ep3s_high = []
    edots_high = []

    for i, edot in enumerate(edots):
        if edot <= edot_break:
            p3s_low.append(p3s[i])
            ep3s_low.append(ep3s[i])
            edots_low.append(edots[i])
        else:
            p3s_high.append(p3s[i])
            ep3s_high.append(ep3s[i])
            edots_high.append(edots[i])

    p3s_low = np.array(p3s_low)
    ep3s_low = np.array(ep3s_low)
    edots_low = np.array(edots_low)
    p3s_high = np.array(p3s_high)
    ep3s_high = np.array(ep3s_high)
    edots_high = np.array(edots_high)

    # alised values
    #p3aliased = lambda p3, n: p3 / (1  - n * p3)
    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            elif p3 == 1:
                return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                elif p3_ == 0:
                    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    # dependence in question
    #p3fun = lambda x, w, th: np.fabs(1 / (1 - w / (x - th)))
    def p3fun(xs, w, th):
        if type(xs) == np.float64:
            if xs != th:
                return np.fabs(1 / (1 - w / (xs - th)))
            else:
                return 1.
        else:
            res = []
            for x in xs:
                if x != th:
                    res.append(np.fabs(1 / (1 - w / (x - th))))
                else:
                    res.append(1)
            return np.array(res)



    w = 1.201 # 0.767 #1 #9.02 #1.271 #1 #1 # 0.5 # 6.32 # 9.65
    th = -1.272# -0.791 #-1.057 #-3.02 #-1.373 #-1.05 #-1 # -0.42 # -3.16 # -3.16

    # get minimum/maximum value
    xt = np.logspace(np.log10(1e-3), np.log10(1e4), num=100)
    yt = p3fun(xt, w, th)
    p3m_max = np.max(yt)
    p3m_min = np.min(yt)

    # model data
    p3s_notobs = []
    xs_notobs = []
    p3s_model = np.empty(len(p3s))
    xs_model = np.empty(len(p3s))

    n = 1
    sigma = np.std(p3s10)
    print("SIGMA: ", sigma)
    #sigma = np.std(np.log10(p3s_low))
    #print("SIGMA low: ", sigma)

    # sigma fun
    sigma_fun = lambda v, x: v[0] + v[1] * x
    v0 = [1, 1]
    xpoints = np.array([-3, 4])
    ypoints = np.array([1, 0.1])
    xsig, ysig, vsig = least_sq(xpoints, ypoints, sigma_fun, v0, xmax=None)


    for k in range(len(xi_xs)):
        #print(edots10[k], "nsp=", nsp)
        p3m = p3fun(xi_xs[k], w, th)
        sig = sigma_fun(vsig, np.log10(xi_xs[k])) * sigma
        #"""
        if p3m > 2:
            sig = sigma
        else:
            sig = sigma_fun(vsig, np.log10(xi_xs[k])) * sigma
        #"""
        #print(np.log10(xi_xs[k]), sigma_fun(vsig, np.log10(xi_xs[k])))
        p3 = 10 ** np.random.normal(np.log10(p3m), sig)
        """
        if p3m > 1.3:
            p3 = 10 ** np.random.normal(np.log10(p3m), sigma)
        else:
            p3 = 10 ** np.random.normal(np.log10(p3m), 0.5 * sigma)
        """
        #p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(0.9*p3m)))
        if p3 > 2:
            p3s_model[k] = p3
            xs_model[k] = xi_xs[k]
        else:
            # check aliasing first
            for nn in range(1, 4):
                p3obs = p3aliased(p3, nn)
                if p3obs > 2:
                    aliased = True
                    break
            # generate new p3
            while p3obs < 1.3:
                p3 = 10 ** np.random.normal(np.log10(p3m), sig)
                #p3 = 10 ** np.random.normal(np.log10(p3m), 0.5 * sigma)
                #p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(0.9*p3m)))
                if p3 > 2:
                    p3obs = p3
                    aliased = False
                else:
                    for nn in range(1, 4):
                        p3obs = p3aliased(p3, nn)
                        if p3obs > 2:
                            aliased = True
                            break
            """
            kk = 1
            aliased = True
            while p3obs < 2:
                p3 = 10 ** np.random.normal(np.log10(p3m), 0.4* sigma)
                #p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(0.9*p3m)))
                if p3 < 2:
                    p3obs = p3aliased(p3, n)
                else:
                    p3obs = p3
                    aliased = False
            """

            """
            p3obs = p3aliased(p3, n)
            kk = 1
            aliased = True
            while p3obs < 2:
                p3 = 10 ** np.random.normal(np.log10(p3m), 0.4* sigma)
                #p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(0.9*p3m)))
                if p3 < 2:
                    p3obs = p3aliased(p3, n)
                else:
                    p3obs = p3
                    aliased = False
            """

            """
            while p3obs < 2:
                for nn in range(1,4):
                    p3 = 10 ** np.random.normal(np.log10(p3m), sigma)
                    #p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(0.9*p3m)))
                    if p3 < 2:
                        p3obs = p3aliased(p3, nn)
                    else:
                        p3obs = p3
                        aliased = False
                    if p3obs > 2:
                        break
            """
            """
                kk +=1
                if kk == 100:
                    break
                    p3obs = 100
            #"""
            if aliased is True:
                p3s_notobs.append(p3)
                xs_notobs.append(xi_xs[k])
            #print(p3obs)
            p3s_model[k] = p3obs
            xs_model[k] = xi_xs[k]

        """
        p3 = 10 ** np.random.normal(np.log10(p3m), sigma)
        #p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(0.9*p3m)))
        if p3 > 2:
            p3s_model[k] = p3
            xs_model[k] = xi_xs[k]
        else:
            p3obs = p3aliased(p3, n)
            while p3obs < 2:
                #p3 = 10 ** np.random.normal(np.log10(p3m), sigma)
                p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(0.9*p3m)))
                p3obs = p3aliased(p3, n)
            p3s_notobs.append(p3)
            xs_notobs.append(xi_xs[k])
            #print(p3obs)
            p3s_model[k] = p3obs
            xs_model[k] = xi_xs[k]
        #"""

    # randomize observations
    #"""
    for zz in range(len(p3s10)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s10_rand[zz] = np.log10(ran)
        p3s_rand[zz] = ran
        #print("p3 ", p3s[zz], "new p3", ran, "error", ep3s[zz])
        #p3s10_rand[zz] = p3s10[zz] # not random
    #"""

    """
    pl.figure()
    pl.minorticks_on()
    pl.scatter(xi_xs, p3s_rand, alpha=0.7, ec="None")
    pl.scatter(xs_model, p3s_model, alpha=0.7, ec="None")
    pl.scatter(xs_notobs, p3s_notobs, c="tab:grey", s=3, alpha=0.7)

    pl.loglog()
    pl.show()
    """

    pvals, xpvals, pvals_bin, xpvals_bin = da.calculate_divergance(p3s_rand, p3s_model, xi_xs, ep3s, sample=30)
    #pvals, xpvals, pvals_bin, xpvals_bin = calculate_divergance(p3s10, p3s_model, edots10, ep3s) # Note ep3s not used
    divergence = len([pv for pv in pvals if pv < 0.01])
    print("Divergence: ", divergence)

    pl.figure(figsize=(13.149606, 13.946563744))
    pl.subplot(2,1,1)
    pl.minorticks_on()
    pl.scatter(xi_xs, p3s_rand, alpha=0.5, ec="None")
    pl.scatter(xs_model, p3s_model, alpha=0.5, ec="None")
    pl.scatter(xs_notobs, p3s_notobs)
    pl.plot(xt, yt, c="black")
    #pl.plot(xline, yline)
    pl.axhline(y=2, ls=":", c="black")
    xlims = pl.xlim()
    pl.loglog()
    pl.subplot(2,1,2)
    pl.minorticks_on()
    #pl.scatter(xpvals, pvals)
    try:
        pl.plot(xpvals_bin, pvals_bin, lw=2, c="C1")
    except:
        pass
    pl.axhline(y=0.01, c="C2", ls="--")
    pl.semilogx()
    pl.ylabel("p-value")
    #pl.xlim(xlims[0], xlims[1])
    #pl.scatter(edots10, sigs, c="green")
    filename = "output/check_p3edot.pdf"
    print(filename)
    pl.savefig(filename)

    #pl.show()


""" Patrick's idea to fit Geoff's dependence """
def p3edot_distributions(data, size=1e5):

    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])

    xi_xs = []
    for i, p in enumerate(periods):
        val = p ** (-1.78) * (pdots[i]/1e-15)
        xi_xs.append(val)
    xi_xs = np.array(xi_xs)

    """
    # sort according to xi_xs
    idx = np.argsort(xi_xs)
    xi_xs = xi_xs[idx]
    p3s = p3s[idx]
    ep3s = ep3s[idx]
    edots = edots[idx]
    periods = periods[idx]
    pdots = pdots[idx]
    #"""

    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))
    p3s_rand = np.empty(len(p3s))

    edot_break = 7e30
    p3s_low = []
    ep3s_low = []
    edots_low = []

    p3s_high = []
    ep3s_high = []
    edots_high = []

    for i, edot in enumerate(edots):
        if edot <= edot_break:
            p3s_low.append(p3s[i])
            ep3s_low.append(ep3s[i])
            edots_low.append(edots[i])
        else:
            p3s_high.append(p3s[i])
            ep3s_high.append(ep3s[i])
            edots_high.append(edots[i])

    p3s_low = np.array(p3s_low)
    ep3s_low = np.array(ep3s_low)
    edots_low = np.array(edots_low)
    p3s_high = np.array(p3s_high)
    ep3s_high = np.array(ep3s_high)
    edots_high = np.array(edots_high)

    # alised values
    #p3aliased = lambda p3, n: p3 / (1  - n * p3)
    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            elif p3 == 1:
                return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                elif p3_ == 0:
                    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    # dependence in question
    #p3fun = lambda x, w, th: np.fabs(1 / (1 - w / (x - th)))
    def p3fun(xs, w, th):
        if type(xs) == np.float64:
            if xs != th:
                return np.fabs(1 / (1 - w / (xs - th)))
            else:
                return 1.
        else:
            res = []
            for x in xs:
                if x != th:
                    res.append(np.fabs(1 / (1 - w / (x - th))))
                else:
                    res.append(1)
            return np.array(res)

    w = 1.0#1 # 0.767 #1 #9.02 #1.271 #1 #1 # 0.5 # 6.32 # 9.65
    th = -1.05#-1.05 # -0.791 #-1.057 #-3.02 #-1.373 #-1.05 #-1 # -0.42 # -3.16 # -3.16

    # get minimum / maximum value
    xt = np.logspace(np.log10(1e-3), np.log10(1e4), num=100)
    yt = p3fun(xt, w, th)
    p3m_max = np.max(yt)
    p3m_min = np.min(yt)

    # model data
    p3s_notobs = []
    xs_notobs = []
    p3s_model = np.empty(len(p3s))
    xs_model = np.empty(len(p3s))
    p3s_dist = np.empty(len(p3s))

    n = 1
    sigma = np.std(p3s10)
    print("SIGMA: ", sigma)
    sigma = np.std(np.log10(p3s_low))
    #sigma = np.std(p3s_low)
    print("SIGMA low: ", sigma)

    # sigma fun
    sigma_fun = lambda v, x: v[0] + v[1] * x
    v0 = [1, 1]
    xpoints = np.array([-3, 4])
    ypoints = np.array([1, 1])
    xsig, ysig, vsig = least_sq(xpoints, ypoints, sigma_fun, v0, xmax=None)

    p3max = np.max(p3s)
    #print(p3max)
    probabilities = []

    for i in range(len(xi_xs)):
        #print(xi_xs[i], edots[i])
        p3m = p3fun(xi_xs[i], w, th)
        sig = sigma_fun(vsig, np.log10(xi_xs[i])) * sigma
        p3s_dist = 10 ** np.random.normal(np.log10(p3m), sig, size=int(size))
        # get the aliased values
        for k, p3 in enumerate(p3s_dist):
            if p3 < 2:
                for nn in range(1, 4):
                    p3obs = p3aliased(p3, nn)
                    if p3obs > 2 and p3obs < p3max:
                        break
                    else:
                        p3obs = p3
                while p3obs < 2 or p3obs > p3max:
                    p3 = 10 ** np.random.normal(np.log10(p3m), sig)
                    if p3 > 2 and p3obs < p3max:
                        p3obs = p3
                        break
                    else:
                        for nn in range(1, 4):
                            p3obs = p3aliased(p3, nn)
                            if p3obs > 2 and p3obs < p3max:
                                break
                p3s_dist[k] = p3obs
        #pl.hist(p3s_dist)
        #pl.show()

        # fit log-normal PDF
        # with outliers
        p3stofit = p3s_dist
        """
        # no outliers in fitting...
        p3ma2 = 10 ** (np.log10(p3m) + 3 * sig)
        p3stofit = []
        for p3 in p3s_dist:
            if p3 < p3ma2:
                p3stofit.append(p3)
        #"""
        p3stofit = np.array(p3stofit)
        hi = np.histogram(p3stofit, bins=100)
        yhi = hi[0]
        xhi = np.empty(len(yhi))
        for j in range(len(hi[1]) - 1):
            xhi[j] = hi[1][j] + 0.5 * (hi[1][j+1] - hi[1][j])
        #log_norm = lambda v, x: 1 / (v[0] * np.sqrt(2)) * np.exp(- (np.log(x) - v[1])**2 / (2 * v[0] ** 2)) # this one works
        #log_norm = lambda v, x: 1 / (x * v[0] * np.sqrt(2)) * np.exp(- (np.log(x) - v[1])**2 / (2 * v[0] ** 2)) # this one does not # problem with v[0]...
        #log_norm = lambda v, x: v[0] * np.exp(- (np.log(x) - v[1])**2) # also works
        #log_norm = lambda v, x: v[0] * np.exp(- (x - v[1]) ** 2) # normal distribution does not work
        log_norm = lambda v, x: v[0] / x * np.exp(- (np.log(x) - v[1])**2) # best (log-normal...)
        v0 = [1, 1]
        yhi = yhi /  np.max(yhi) # do the fitting with ymax = 1
        x, y, v = least_sq(xhi, yhi, log_norm, v0, xmax=None)

        # generate more complete function
        x_ = np.linspace(2, 500, num=1000)
        y_ = log_norm(v, x_)

        # calculate integral to normalize the distribution
        integral = np.trapz(y_, x_)
        v[0] = v[0] / integral
        # y = log_norm(v, x) # also works
        y /= integral
        y_ /= integral
        yhi /= integral

        # calculate probablility (integration)
        x_min = p3s[i] - 3 * ep3s[i]
        if x_min < 2:
            x_min = 2
        x_max = p3s[i] + 3*ep3s[i]
        xx = np.linspace(x_min, x_max, num=1000)
        yy = log_norm(v, xx)
        prob = np.trapz(yy, xx)
        #print("Probablility: ", prob)
        probabilities.append(prob)

        p3s_model[i] = p3s_dist[0]
        xs_model[i] = xi_xs[i]
        #"""
        if xi_xs[i] > -20:
            #print(v)
            #print(len(yhi))
            #print(len(xhi))
            #pl.hist(p3stofit, bins=100)
            print(probabilities[-1])
            pl.step(xhi, yhi, where='mid')
            pl.plot(xhi, yhi, c="tab:orange")
            pl.plot(x, y, c="tab:green")
            #pl.plot(x_, y_, c="tab:red")
            pl.axvline(x=p3s[i], lw=2, c="tab:red")
            pl.axvline(x=p3s[i] - 3*ep3s[i], lw=1, c="tab:red")
            pl.axvline(x=p3s[i] + 3*ep3s[i], lw=1, c="tab:red")
            pl.show()
        #"""
        #return
    # randomize observations
    #"""
    for zz in range(len(p3s10)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s10_rand[zz] = np.log10(ran)
        p3s_rand[zz] = ran
        #print("p3 ", p3s[zz], "new p3", ran, "error", ep3s[zz])
        #p3s10_rand[zz] = p3s10[zz] # not random
    #"""
    print("Mean probability: ", np.mean(probabilities))
    print("Median probability: ", np.median(probabilities))
    #return

    #print(xs_model)

    pl.figure(figsize=(13.149606, 13.946563744))
    pl.subplot(2,1,1)
    pl.minorticks_on()
    pl.scatter(xi_xs, p3s_rand, c=probabilities, alpha=0.7, ec="None")
    pl.errorbar(xi_xs, p3s, fmt='none', yerr=ep3s, color="tab:blue", zorder=2, alpha=0.3)
    pl.colorbar()
    pl.plot(xt, yt, c="black")
    #pl.plot(xline, yline)
    pl.axhline(y=2, ls=":", c="black")
    xlims = pl.xlim()
    pl.loglog()
    pl.subplot(2,1,2)
    pl.minorticks_on()
    #pl.scatter(xpvals, pvals)

    pl.scatter(xi_xs, p3s, c="tab:blue", alpha=0.7, ec="None", label="observed")
    pl.scatter(xs_model, p3s_model, c="tab:orange", alpha=0.7, ec="None", label="modeled")
    pl.loglog()
    #pl.xlim(xlims[0], xlims[1])
    #pl.scatter(edots10, sigs, c="green")
    filename = "output/p3edot_distributions.pdf"
    print(filename)
    pl.savefig(filename)
    #pl.show()


def p3edot_distributions_lin(data, size=1e3):

    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])

    xi_xs = []
    for i, p in enumerate(periods):
        val = p ** (-1.78) * (pdots[i]/1e-15)
        xi_xs.append(val)
    xi_xs = np.array(xi_xs)

    """
    # sort according to xi_xs
    idx = np.argsort(xi_xs)
    xi_xs = xi_xs[idx]
    p3s = p3s[idx]
    ep3s = ep3s[idx]
    edots = edots[idx]
    periods = periods[idx]
    pdots = pdots[idx]
    #"""

    p3s10 = np.log10(p3s)
    ep3s10 = np.log10(ep3s)
    edots10 = np.log10(edots)

    p3s10_rand = np.empty(len(p3s))
    p3s_rand = np.empty(len(p3s))

    edot_break = 7e30
    p3s_low = []
    ep3s_low = []
    edots_low = []

    p3s_high = []
    ep3s_high = []
    edots_high = []

    for i, edot in enumerate(edots):
        if edot <= edot_break:
            p3s_low.append(p3s[i])
            ep3s_low.append(ep3s[i])
            edots_low.append(edots[i])
        else:
            p3s_high.append(p3s[i])
            ep3s_high.append(ep3s[i])
            edots_high.append(edots[i])

    p3s_low = np.array(p3s_low)
    ep3s_low = np.array(ep3s_low)
    edots_low = np.array(edots_low)
    p3s_high = np.array(p3s_high)
    ep3s_high = np.array(ep3s_high)
    edots_high = np.array(edots_high)

    # alised values
    #p3aliased = lambda p3, n: p3 / (1  - n * p3)
    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            elif p3 == 1:
                return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                elif p3_ == 0:
                    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    # dependence in question
    a = -0.4
    y0 = 0.5
    x0 = 1
    b = y0 - a * np.log10(x0)
    p3fun = lambda x: a * x + b
    #p3fun = lambda x, w, th: np.fabs(1 / (1 - w / (x - th)))
    """
    def p3fun(xs, w, th):
        if type(xs) == np.float64:
            if xs != th:
                return np.fabs(1 / (1 - w / (xs - th)))
            else:
                return 1.
        else:
            res = []
            for x in xs:
                if x != th:
                    res.append(np.fabs(1 / (1 - w / (x - th))))
                else:
                    res.append(1)
            return np.array(res)
    """


    w = 1.0#1 # 0.767 #1 #9.02 #1.271 #1 #1 # 0.5 # 6.32 # 9.65
    th = -1.05#-1.05 # -0.791 #-1.057 #-3.02 #-1.373 #-1.05 #-1 # -0.42 # -3.16 # -3.16

    # get minimum / maximum value
    xt = np.logspace(np.log10(1e-3), np.log10(1e4), num=100)
    yt = 10**p3fun(np.log10(xt))
    #yt = p3fun(xt, w, th)
    p3m_max = np.max(yt)
    p3m_min = np.min(yt)

    #pl.plot(xt, yt)
    #pl.loglog()
    #pl.show()

    # model data
    p3s_notobs = []
    xs_notobs = []
    p3s_model = np.empty(len(p3s))
    xs_model = np.empty(len(p3s))
    p3s_dist = np.empty(len(p3s))

    n = 1
    sigma = np.std(p3s10)
    print("SIGMA: ", sigma)
    sigma = np.std(np.log10(p3s_low))
    #sigma = np.std(p3s_low)
    print("SIGMA low: ", sigma)

    # sigma fun
    sigma_fun = lambda v, x: v[0] + v[1] * x
    v0 = [1, 1]
    xpoints = np.array([-3, 4])
    ypoints = np.array([1, 1])
    xsig, ysig, vsig = least_sq(xpoints, ypoints, sigma_fun, v0, xmax=None)

    p3max = np.max(p3s)
    #print(p3max)


    probabilities = []

    for i in range(len(xi_xs)):
        #print(xi_xs[i], edots[i])
        p3m = 10 ** p3fun(np.log10(xi_xs[i]))
        #p3m = p3fun(xi_xs[i], w, th)
        sig = sigma_fun(vsig, np.log10(xi_xs[i])) * sigma
        p3s_dist = 10 ** np.random.normal(np.log10(p3m), sig, size=int(size))
        # get the aliased values
        for k, p3 in enumerate(p3s_dist):
            if p3 < 2:
                for nn in range(1, 4):
                    p3obs = p3aliased(p3, nn)
                    if p3obs > 2 and p3obs < p3max:
                        break
                    else:
                        p3obs = p3
                while p3obs < 2 or p3obs > p3max:
                    p3 = 10 ** np.random.normal(np.log10(p3m), sig)
                    if p3 > 2 and p3obs < p3max:
                        p3obs = p3
                        break
                    else:
                        for nn in range(1, 4):
                            p3obs = p3aliased(p3, nn)
                            if p3obs > 2 and p3obs < p3max:
                                break
                p3s_dist[k] = p3obs

        # fit log-normal PDF
        # with outliers
        p3stofit = p3s_dist
        """
        # no outliers in fitting...
        p3ma2 = 10 ** (np.log10(p3m) + 3 * sig)
        p3stofit = []
        for p3 in p3s_dist:
            if p3 < p3ma2:
                p3stofit.append(p3)
        #"""
        p3stofit = np.array(p3stofit)
        hi = np.histogram(p3stofit, bins=100)
        yhi = hi[0]
        xhi = np.empty(len(yhi))
        for j in range(len(hi[1]) - 1):
            xhi[j] = hi[1][j] + 0.5 * (hi[1][j+1] - hi[1][j])
        #log_norm = lambda v, x: 1 / (v[0] * np.sqrt(2)) * np.exp(- (np.log(x) - v[1])**2 / (2 * v[0] ** 2)) # this one works
        #log_norm = lambda v, x: 1 / (x * v[0] * np.sqrt(2)) * np.exp(- (np.log(x) - v[1])**2 / (2 * v[0] ** 2)) # this one does not # problem with v[0]...
        #log_norm = lambda v, x: v[0] * np.exp(- (np.log(x) - v[1])**2) # also works
        #log_norm = lambda v, x: v[0] * np.exp(- (x - v[1]) ** 2) # normal distribution does not work
        log_norm = lambda v, x: v[0] / x * np.exp(- (np.log(x) - v[1])**2) # best (log-normal...)
        v0 = [1, 1]
        yhi = yhi /  np.max(yhi) # do the fitting with ymax = 1
        x, y, v = least_sq(xhi, yhi, log_norm, v0, xmax=None)

        # generate more complete function
        x_ = np.linspace(2, 500, num=1000)
        y_ = log_norm(v, x_)

        # calculate integral to normalize the distribution
        integral = np.trapz(y_, x_)
        v[0] = v[0] / integral
        # y = log_norm(v, x) # also works
        y /= integral
        y_ /= integral
        yhi /= integral

        # calculate probablility (integration)
        x_min = p3s[i] - 3 * ep3s[i]
        if x_min < 2:
            x_min = 2
        x_max = p3s[i] + 3*ep3s[i]
        xx = np.linspace(x_min, x_max, num=1000)
        yy = log_norm(v, xx)
        prob = np.trapz(yy, xx)
        #print("Probablility: ", prob)
        probabilities.append(prob)

        p3s_model[i] = p3s_dist[0]
        xs_model[i] = xi_xs[i]
        """
        if xi_xs[i] > -20:
            #print(v)
            #print(len(yhi))
            #print(len(xhi))
            #pl.hist(p3stofit, bins=100)
            print(probabilities[-1])
            pl.step(xhi, yhi, where='mid')
            pl.plot(xhi, yhi, c="tab:orange")
            pl.plot(x, y, c="tab:green")
            #pl.plot(x_, y_, c="tab:red")
            pl.axvline(x=p3s[i], lw=2, c="tab:red")
            pl.axvline(x=p3s[i] - 3*ep3s[i], lw=1, c="tab:red")
            pl.axvline(x=p3s[i] + 3*ep3s[i], lw=1, c="tab:red")
            pl.show()
        #"""
        #return
    # randomize observations
    #"""
    for zz in range(len(p3s10)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s10_rand[zz] = np.log10(ran)
        p3s_rand[zz] = ran
        #print("p3 ", p3s[zz], "new p3", ran, "error", ep3s[zz])
        #p3s10_rand[zz] = p3s10[zz] # not random
    #"""
    print("Mean probability: ", np.mean(probabilities))
    print("Median probability: ", np.median(probabilities))
    #return

    #print(xs_model)

    pl.figure(figsize=(13.149606, 13.946563744))
    pl.subplot(2,1,1)
    pl.minorticks_on()
    pl.scatter(xi_xs, p3s_rand, c=probabilities, alpha=0.7, ec="None")
    pl.errorbar(xi_xs, p3s, fmt='none', yerr=ep3s, color="tab:blue", zorder=2, alpha=0.3)
    pl.colorbar()
    pl.plot(xt, yt, c="black")
    #pl.plot(xline, yline)
    pl.axhline(y=2, ls=":", c="black")
    xlims = pl.xlim()
    pl.loglog()
    pl.subplot(2,1,2)
    pl.minorticks_on()
    #pl.scatter(xpvals, pvals)

    pl.scatter(xi_xs, p3s, c="tab:blue", alpha=0.7, ec="None", label="observed")
    pl.scatter(xs_model, p3s_model, c="tab:orange", alpha=0.7, ec="None", label="modeled")
    pl.loglog()
    #pl.xlim(xlims[0], xlims[1])
    #pl.scatter(edots10, sigs, c="green")
    filename = "output/p3edot_distributions_lin.pdf"
    print(filename)
    pl.savefig(filename)

    #pl.show()

def p3edot_inflection_low(data, size=1e3):

    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])

    xi_xs = []
    for i, p in enumerate(periods):
        val = p ** (-1.78) * (pdots[i]/1e-15)
        xi_xs.append(val)
    xi_xs = np.array(xi_xs)

    #"""
    # sort according to xi_xs
    idx = np.argsort(xi_xs)
    xi_xs = xi_xs[idx]
    p3s = p3s[idx]
    ep3s = ep3s[idx]
    edots = edots[idx]
    periods = periods[idx]
    pdots = pdots[idx]
    #"""

    # poarabolic fit to whole data
    para = lambda v, x: v[0] * x ** 2 + v[1] * x + v[2]
    xp, yp, vp = least_sq(np.log10(xi_xs), np.log10(p3s), para, [1,1,1], xmax=None)

    indxs = np.argsort(10**yp)
    in_ = 10**xp[indxs[0]]
    #in_ = 3

    print("Inflection point: ", in_)

    xi_low = []
    xi_high = []
    p3s_low = []
    p3s_high = []
    ep3s_low = []
    ep3s_high = []

    for i,xi in enumerate(xi_xs):
        if xi <= in_:
            xi_low.append(xi)
            p3s_low.append(p3s[i])
            ep3s_low.append(ep3s[i])
        else:
            xi_high.append(xi)
            p3s_high.append(p3s[i])
            ep3s_high.append(ep3s[i])

    # linear fit
    lin = lambda v, x:  v[0] * x + v[1]
    xh, yh, vh = least_sq(np.log10(xi_high), np.log10(p3s_high), lin, [1,1], xmax=None)
    #xh2, yh2, vh2, err = least_sq1D(np.log10(xi_high), np.log10(p3s_high), lin, np.log10(ep3s_high), [1,1])
    xh3, yh3, vh3, err3 = odr(np.log10(xi_high), np.log10(p3s_high), np.array([1 for i in range(len(ep3s_high))]), np.log10(ep3s_high),lin,  [1,1])

    xl, yl, vl = least_sq(np.log10(xi_low), np.log10(p3s_low), lin, [1,1], xmax=None)
    #xl2, yl2, vl2, err2 = least_sq1D(np.log10(xi_low), np.log10(p3s_low), lin, np.log10(ep3s_low), [1,1])
    xl3, yl3, vl3, err3 = odr(np.log10(xi_low), np.log10(p3s_low), np.array([1 for i in range(len(ep3s_low))]), np.log10(ep3s_low),lin,  [1,1])

    ############################################################################
    # simulating points
    # alised values
    #p3aliased = lambda p3, n: p3 / (1  - n * p3)
    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            #elif p3 == 1:
            #    return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                #elif p3_ == 0:
                #    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    # dependence in question
    a = vl[0]
    b = vl[1]
    #a = -0.7
    #y0 = 0.2
    #x0 = 1
    #b = y0 - a * np.log10(x0)
    #a = 0.9
    #y0 = 0.1
    #x0 = 1
    #b = y0 - a * np.log10(x0)
    p3fun = lambda x: a * x + b

    # get minimum / maximum value
    xt = np.logspace(np.log10(1e-3), np.log10(1e4), num=100)
    yt = 10**p3fun(np.log10(xt))

    p3m_max = np.max(yt)
    p3m_min = np.min(yt)

    #pl.plot(xt, yt)
    #pl.loglog()
    #pl.show()
    #exit()

    # model data
    p3s_notobs = []
    xs_notobs = []
    p3s_model = [] #np.empty(len(p3s))
    xs_model = [] #np.empty(len(p3s))
    p3s_model_hl = []
    xs_model_hl = []

    n = 1
    sigma = np.std(np.log10(p3s_low))
    print("SIGMA: ", sigma)

    # sigma fun
    sigma_fun = lambda v, x: v[0] + v[1] * x
    v0 = [1, 1]
    xpoints = np.array([-3, 0, 4])
    ypoints = np.array([1, 0.3, 0.1])
    #ypoints = np.array([1, 1, 1])
    xsig, ysig, vsig = least_sq(xpoints, ypoints, sigma_fun, v0, xmax=None)
    #pl.plot(xsig, ysig)
    #pl.scatter(xpoints, ypoints)
    #pl.show()
    def sigma_fun(v, x):
        if x < 1:
            return 1
        else:
            return 0.7


    p3max = np.max(p3s)
    ns = []
    for i in range(len(xi_xs)):
        skip = False
        #print(xi_xs[i], edots[i])
        #p3m = 10 ** p3fun(np.log10(xi_xs[i]))
        p3m = 10 ** p3fun(np.log10(xi_xs[i]))
        #p3m = p3fun(xi_xs[i], w, th)
        sig = sigma_fun(vsig, np.log10(xi_xs[i])) * sigma
        p3 = 10 ** np.random.normal(np.log10(p3m), sig)
        # get the aliased values

        if p3 < 2:
            for nn in range(1, 100):
                p3obs = p3aliased(p3, nn)
                if p3obs > 2 and p3obs < p3max:
                    ns.append(nn)
                    break
                else:
                    p3obs = p3
            r = 0
            while p3obs < 2 or p3obs > p3max:
                p3 = 10 ** np.random.normal(np.log10(p3m), sig)
                if p3 > 2 and p3obs < p3max:
                    p3obs = p3
                    break
                else:
                    for nn in range(1, 100):
                        p3obs = p3aliased(p3, nn)
                        #print(nn)
                        if p3obs > 2 and p3obs < p3max:
                            #print(nn)
                            ns.append(nn)
                            break
                r += 1
                if r == 10:
                    skip = True
                    break
        else:
            p3obs = p3
        if skip is False:
            #p3s_model[i] = p3obs
            #xs_model[i] = xi_xs[i]
            p3s_model.append(p3obs)
            xs_model.append(xi_xs[i])
            if xi_xs[i] > in_:
                p3s_model_hl.append(p3obs)
                xs_model_hl.append(xi_xs[i])
        else:
            p3s_notobs.append(p3obs)
            xs_notobs.append(xi_xs[i])
    #print(ns)
    #print(len(p3s_model))
    #print(len(p3s))
    #print(len(p3s_notobs))

    p3s_model_hl = np.array(p3s_model_hl)
    xs_model_hl = np.array(xs_model_hl)
    xm, ym, vm = least_sq(np.log10(xs_model_hl), np.log10(p3s_model_hl), lin, [1,1], xmax=None)

    pl.rc("font", size=8)
    pl.rc("axes", linewidth=0.5)

    pl.figure()
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    #pl.subplot(2,1,1)
    pl.minorticks_on()
    pl.scatter(xi_xs, p3s, alpha=0.5, ec="None")
    pl.errorbar(xi_xs, p3s, fmt='none', yerr=ep3s, color="tab:blue", zorder=2, alpha=0.3)
    pl.scatter(xs_model, p3s_model, alpha=0.5, ec="None", color="tab:orange")
    pl.scatter(xs_model_hl, p3s_model_hl, alpha=0.5, ec="None", color="tab:red")
    pl.plot(xt, yt, c="black")
    pl.plot(10**xm, 10**ym, c="tab:red")
    #pl.plot(10**xp, 10**yp)
    #pl.plot(10**xl, 10**yl)
    #pl.plot(10**xl2, 10**yl2)
    #pl.plot(10**xl3, 10**yl3, c="tab:red")
    pl.plot(10**xh, 10**yh, c="tab:blue")
    #pl.plot(10**xh2, 10**yh2, c="pink")
    #pl.plot(10**xh3, 10**yh3, c="tab:green")
    #pl.plot(xline, yline)
    pl.axhline(y=2, ls=":", c="black")
    pl.axvline(x=in_, ls="--", c="tab:red")
    xlims = pl.xlim()
    pl.loglog()
    pl.xlabel(r"$(P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15})$")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3edot_inflection_low.pdf"
    print(filename)
    pl.savefig(filename)
    pl.show()


def p3edot_inflection_high(data, size=1e3):

    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])

    xi_xs = []
    for i, p in enumerate(periods):
        val = p ** (-1.78) * (pdots[i]/1e-15)
        xi_xs.append(val)
    xi_xs = np.array(xi_xs)

    #"""
    # sort according to xi_xs
    idx = np.argsort(xi_xs)
    xi_xs = xi_xs[idx]
    p3s = p3s[idx]
    ep3s = ep3s[idx]
    edots = edots[idx]
    periods = periods[idx]
    pdots = pdots[idx]
    #"""

    # poarabolic fit to whole data
    para = lambda v, x: v[0] * x ** 2 + v[1] * x + v[2]
    xp, yp, vp = least_sq(np.log10(xi_xs), np.log10(p3s), para, [1,1,1], xmax=None)

    indxs = np.argsort(10**yp)
    in_ = 10**xp[indxs[0]]
    #in_ = 3

    print("Inflection point: ", in_)

    xi_low = []
    xi_high = []
    p3s_low = []
    p3s_high = []
    ep3s_low = []
    ep3s_high = []

    for i,xi in enumerate(xi_xs):
        if xi <= in_:
            xi_low.append(xi)
            p3s_low.append(p3s[i])
            ep3s_low.append(ep3s[i])
        else:
            xi_high.append(xi)
            p3s_high.append(p3s[i])
            ep3s_high.append(ep3s[i])

    # linear fit
    lin = lambda v, x:  v[0] * x + v[1]
    xh, yh, vh = least_sq(np.log10(xi_high), np.log10(p3s_high), lin, [1,1], xmax=None)
    #xh2, yh2, vh2, err = least_sq1D(np.log10(xi_high), np.log10(p3s_high), lin, np.log10(ep3s_high), [1,1])
    xh3, yh3, vh3, err3 = odr(np.log10(xi_high), np.log10(p3s_high), np.array([1 for i in range(len(ep3s_high))]), np.log10(ep3s_high),lin,  [1,1])

    xl, yl, vl = least_sq(np.log10(xi_low), np.log10(p3s_low), lin, [1,1], xmax=None)
    #xl2, yl2, vl2, err2 = least_sq1D(np.log10(xi_low), np.log10(p3s_low), lin, np.log10(ep3s_low), [1,1])
    xl3, yl3, vl3, err3 = odr(np.log10(xi_low), np.log10(p3s_low), np.array([1 for i in range(len(ep3s_low))]), np.log10(ep3s_low),lin,  [1,1])

    ############################################################################
    # simulating points
    # alised values
    #p3aliased = lambda p3, n: p3 / (1  - n * p3)
    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            #elif p3 == 1:
            #    return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                #elif p3_ == 0:
                #    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    # dependence in question
    a = vh[0]
    b = vh[1]
    #a = -0.7
    #y0 = 0.2
    #x0 = 1
    #b = y0 - a * np.log10(x0)
    #a = 0.9
    #y0 = 0.1
    #x0 = 1
    #b = y0 - a * np.log10(x0)
    p3fun = lambda x: a * x + b

    # get minimum / maximum value
    xt = np.logspace(np.log10(1e-3), np.log10(1e4), num=100)
    yt = 10**p3fun(np.log10(xt))

    p3m_max = np.max(yt)
    p3m_min = np.min(yt)

    #pl.plot(xt, yt)
    #pl.loglog()
    #pl.show()
    #exit()

    # model data
    p3s_notobs = []
    xs_notobs = []
    p3s_model = [] #np.empty(len(p3s))
    xs_model = [] #np.empty(len(p3s))
    p3s_model_hl = []
    xs_model_hl = []

    n = 1
    sigma = 1.0 * np.std(np.log10(p3s_high)) # why 1.3?
    print("SIGMA: ", sigma)

    # sigma fun
    sigma_fun = lambda v, x: v[0] + v[1] * x
    v0 = [1, 1]
    xpoints = np.array([-3, 0, 4])
    #ypoints = np.array([0.1, 0.5, 1])
    ypoints = np.array([1, 1, 1])
    xsig, ysig, vsig = least_sq(xpoints, ypoints, sigma_fun, v0, xmax=None)
    #pl.plot(xsig, ysig)
    #pl.scatter(xpoints, ypoints)
    #pl.show()

    p3max = np.max(p3s)

    for i in range(len(xi_xs)):
        skip = False
        #print(xi_xs[i], edots[i])
        #p3m = 10 ** p3fun(np.log10(xi_xs[i]))
        p3m = 10 ** p3fun(np.log10(xi_xs[i]))
        #p3m = p3fun(xi_xs[i], w, th)
        sig = sigma_fun(vsig, np.log10(xi_xs[i])) * sigma
        p3 = 10 ** np.random.normal(np.log10(p3m), sig)
        # get the aliased values

        if p3 < 2:
            for nn in range(1, 100):
                p3obs = p3aliased(p3, nn)
                if p3obs > 2 and p3obs < p3max:
                    break
                else:
                    p3obs = p3
            r = 0
            while p3obs < 2 or p3obs > p3max:
                p3 = 10 ** np.random.normal(np.log10(p3m), sig)
                if p3 > 2 and p3obs < p3max:
                    p3obs = p3
                    break
                else:
                    for nn in range(1, 100):
                        p3obs = p3aliased(p3, nn)
                        if p3obs > 2 and p3obs < p3max:
                            break
                r += 1
                if r == 10:
                    skip = True
                    break

        else:
            p3obs = p3
        if skip is False:
            #p3s_model[i] = p3obs
            #xs_model[i] = xi_xs[i]
            p3s_model.append(p3obs)
            xs_model.append(xi_xs[i])
            if xi_xs[i] < in_:
                p3s_model_hl.append(p3obs)
                xs_model_hl.append(xi_xs[i])
        else:
            p3s_notobs.append(p3obs)
            xs_notobs.append(xi_xs[i])


    p3s_model_hl = np.array(p3s_model_hl)
    xs_model_hl = np.array(xs_model_hl)
    xm, ym, vm = least_sq(np.log10(xs_model_hl), np.log10(p3s_model_hl), lin, [1,1], xmax=None)


    pl.rc("font", size=8)
    pl.rc("axes", linewidth=0.5)
    #pl.rc("lines", linewidth=0.5)

    #pl.figure(figsize=(8/2.41, 8/1.618/2.41))\
    pl.figure()
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    #pl.subplot(2,1,1)
    pl.minorticks_on()
    pl.scatter(xi_xs, p3s, alpha=0.5, ec="None")
    pl.errorbar(xi_xs, p3s, fmt='none', yerr=ep3s, color="tab:blue", zorder=2, alpha=0.3)
    pl.scatter(xs_model, p3s_model, alpha=0.5, ec="None", color="tab:orange")
    pl.scatter(xs_model_hl, p3s_model_hl, alpha=0.5, ec="None", color="tab:red")
    pl.plot(xt, yt, c="black")
    pl.plot(10**xm, 10**ym, c="tab:red")
    #pl.plot(10**xp, 10**yp)
    pl.plot(10**xl, 10**yl, c="tab:blue")
    #pl.plot(10**xl2, 10**yl2)
    #pl.plot(10**xl3, 10**yl3, c="tab:red")
    #pl.plot(10**xh, 10**yh, c="tab:blue")
    #pl.plot(10**xh2, 10**yh2, c="pink")
    #pl.plot(10**xh3, 10**yh3, c="tab:green")
    #pl.plot(xline, yline)
    pl.axhline(y=2, ls=":", c="black")
    pl.axvline(x=in_, ls="--", c="tab:red")
    xlims = pl.xlim()
    pl.loglog()
    pl.xlabel(r"$(P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15})$")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3edot_inflection_high.pdf"
    print(filename)
    pl.savefig(filename)
    pl.show()


def p3edot_inflection_geoff(data, size=1e3):

    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])

    xi_xs = []
    for i, p in enumerate(periods):
        val = p ** (-1.78) * (pdots[i]/1e-15)
        xi_xs.append(val)
    xi_xs = np.array(xi_xs)

    #"""
    # sort according to xi_xs
    idx = np.argsort(xi_xs)
    xi_xs = xi_xs[idx]
    p3s = p3s[idx]
    ep3s = ep3s[idx]
    edots = edots[idx]
    periods = periods[idx]
    pdots = pdots[idx]
    #"""

    # poarabolic fit to whole data
    para = lambda v, x: v[0] * x ** 2 + v[1] * x + v[2]
    xp, yp, vp = least_sq(np.log10(xi_xs), np.log10(p3s), para, [1,1,1], xmax=None)

    indxs = np.argsort(10**yp)
    in_ = 10 ** xp[indxs[0]]
    #in_ = 3

    print("Inflection point: ", in_)

    xi_low = []
    xi_high = []
    p3s_low = []
    p3s_high = []
    ep3s_low = []
    ep3s_high = []

    for i,xi in enumerate(xi_xs):
        if xi <= in_:
            xi_low.append(xi)
            p3s_low.append(p3s[i])
            ep3s_low.append(ep3s[i])
        else:
            xi_high.append(xi)
            p3s_high.append(p3s[i])
            ep3s_high.append(ep3s[i])

    # linear fit
    lin = lambda v, x:  v[0] * x + v[1]
    xh, yh, vh = least_sq(np.log10(xi_high), np.log10(p3s_high), lin, [1,1], xmax=None)
    #xh2, yh2, vh2, err = least_sq1D(np.log10(xi_high), np.log10(p3s_high), lin, np.log10(ep3s_high), [1,1])
    xh3, yh3, vh3, err3 = odr(np.log10(xi_high), np.log10(p3s_high), np.array([1 for i in range(len(ep3s_high))]), np.log10(ep3s_high),lin,  [1,1])

    xl, yl, vl = least_sq(np.log10(xi_low), np.log10(p3s_low), lin, [1,1], xmax=None)
    #xl2, yl2, vl2, err2 = least_sq1D(np.log10(xi_low), np.log10(p3s_low), lin, np.log10(ep3s_low), [1,1])
    xl3, yl3, vl3, err3 = odr(np.log10(xi_low), np.log10(p3s_low), np.array([1 for i in range(len(ep3s_low))]), np.log10(ep3s_low),lin,  [1,1])

    ############################################################################
    # simulating points
    # alised values
    #p3aliased = lambda p3, n: p3 / (1  - n * p3)
    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            #elif p3 == 1:
            #    return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                #elif p3_ == 0:
                #    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    #"""
    def p3fun(xs, w, th):
        if type(xs) == np.float64:
            if xs != th:
                return np.fabs(1 / (1 - w / (xs - th)))
            else:
                return 1.
        else:
            res = []
            for x in xs:
                if x != th:
                    res.append(np.fabs(1 / (1 - w / (x - th))))
                else:
                    res.append(1)
            return np.array(res)
    #"""

    w = 1.0 # 1.0 # 1.5 #1 # 0.767 #1 #9.02 #1.271 #1 #1 # 0.5 # 6.32 # 9.65
    th = -1.03 # -1.03 # -1.55 #-1.05#-1.05 # -0.791 #-1.057 #-3.02 #-1.373 #-1.05 #-1 # -0.42 # -3.16 # -3.16

    # intrinsic dependence
    xt = np.logspace(np.log10(1e-3), np.log10(1e4), num=1000)
    yt = p3fun(xt, w, th)

    # get minimum / maximum value
    p3m_max = np.max(yt)
    p3m_min = np.min(yt)

    # intrinsic dependence in high xi region
    ythi = []
    xthi = []
    for i, y in enumerate(yt):
        if y < 2:
            xthi.append(xt[i])
            ythi.append(p3aliased(y, 1))
    xthi = np.array(xthi)
    ythi = np.array(ythi)

    #pl.plot(xt, yt)
    #pl.loglog()
    #pl.show()
    #exit()

    # model data
    p3s_notobs = []
    xs_notobs = []
    p3s_model = [] #np.empty(len(p3s))
    xs_model = [] #np.empty(len(p3s))
    p3s_model_hl = []
    xs_model_hl = []

    n = 1
    sigma = np.std(np.log10(p3s_low))
    print("SIGMA: ", sigma)
    # sigma fun
    sigma_fun = lambda v, x: v[0] + v[1] * x
    v0 = [1, 1]
    xpoints = np.array([-3, 0, 4])
    #ypoints = np.array([1., 0.2, 0.1])
    ypoints = np.array([1.01, 0.3, 0.05]) # best!
    #ypoints = np.array([0.1, 0.1, 0.1])
    #ypoints = np.array([1., 1, 1])
    xsig, ysig, vsig = least_sq(xpoints, ypoints, sigma_fun, v0, xmax=None)
    #pl.plot(xsig, ysig)
    #pl.scatter(xpoints, ypoints)
    #pl.show()
    """" for some test?
    sig_fun = lambda v, x: v[0] + v[1] * x
    v0 = [1, 1]
    xpoints = np.array([0, 4])
    ypoints = np.array([0.7, 0.001])
    #ypoints = np.array([0.9, 0.05])
    xsig, ysig, vsig = least_sq(xpoints, ypoints, sig_fun, v0, xmax=None)
    def sigma_fun(v, x):
        if x < -0.05:
            return 1
        else:
            return sig_fun(v, x)
    #"""
    # intrinsic dependence variation # TODO improve it!
    ytl = np.empty(len(yt))
    yth = np.empty(len(yt))

    # generate model + sigma_fun
    for i, y in enumerate(yt):
        ytl[i] = 10 ** (np.log10(y) - sigma_fun(vsig, np.log10(xt[i])) * sigma)
        yth[i] = 10 ** (np.log10(y) + sigma_fun(vsig, np.log10(xt[i])) * sigma)

    # intrinsic dependence new variation (circle stuff)
    ytl2 = np.zeros(len(yt))
    yth2 = np.zeros(len(yt))

    # circle data
    yci = [[] for i in range(len(xt))]
    for i,x in enumerate(xt):
        a = np.log10(x)
        b = np.log10(yt[i])
        r = np.fabs(sigma_fun(vsig, np.log10(x)) * sigma) # is it ok?
        #print(r)
        #print("a b r", a, b, r)
        # loop over the x
        ys = []
        xs = []
        for j, xc in enumerate(xt):
            if xc > 10 ** (a + r):
                break
            if xc >= 10 ** (a - r):
                y_1 = b - np.sqrt(-a ** 2 + 2 * a * np.log10(xc) + r ** 2 - np.log10(xc) ** 2)
                y_2 = np.sqrt(-a ** 2 + 2 * a * np.log10(xc) + r ** 2 - np.log10(xc) ** 2) + b
                ys.append(y_1)
                ys.append(y_2)
                xs.append(np.log10(xc))
                xs.append(np.log10(xc))
                yci[j].append(y_1)
                yci[j].append(y_2)
    for i in range(len(xt)):
        ytl2[i] = np.min(10 ** np.array(yci[i]))
        yth2[i] = np.max(10 ** np.array(yci[i]))
        #print(i, xt[i],  yt[i], ytl2[i], yth2[i])

    def get_sigma(xi):
        """ WOW this is weird - you have too much time?"""
        mi = 1e50
        sig = 10
        for i, x in enumerate(xt):
            diff = np.fabs(xi - x)
            if diff < mi:
                sig = (np.log10(yth2[i]) - np.log10(ytl2[i])) / 2
                mi = diff
        return sig


    p3max = np.max(p3s)

    for i in range(len(xi_xs)):
        aliased = False
        p3m = p3fun(xi_xs[i], w, th)
        #sig = sigma_fun(vsig, np.log10(xi_xs[i])) * sigma
        sig = get_sigma(xi_xs[i]) # new weird approach
        #print(sig)
        p3  = 10 ** np.random.normal(np.log10(p3m), sig)
        if p3 > 2:
            p3s_model.append(p3)
            xs_model.append(xi_xs[i])
        else:
            # check aliasing first
            for nn in range(1, 10):
                p3obs = p3aliased(p3, nn)
                if p3obs > 2:
                    aliased = True
                    break
            # generate new p3
            while p3obs < 2 or p3obs > 2 * p3max:
                p3 = 10 ** np.random.normal(np.log10(p3m), sig)
                #p3 = 10 ** np.random.normal(np.log10(p3m), 0.5 * sigma)
                #p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(0.9*p3m)))
                if p3 > 2:
                    p3obs = p3
                    aliased = False
                else:
                    for nn in range(1, 4):
                        p3obs = p3aliased(p3, nn)
                        if p3obs > 2:
                            aliased = True
                            break
            if aliased is True:
                p3s_notobs.append(p3)
                xs_notobs.append(xi_xs[i])
            #print(p3obs)
            p3s_model.append(p3obs)
            xs_model.append(xi_xs[i])

        if xi_xs[i] > in_:
            p3s_model_hl.append(p3s_model[-1])
            xs_model_hl.append(xi_xs[i])

    p3s_model_hl = np.array(p3s_model_hl)
    xs_model_hl = np.array(xs_model_hl)
    xm, ym, vm = least_sq(np.log10(xs_model_hl), np.log10(p3s_model_hl), lin, [1,1], xmax=None)


    """
    # randomize observations
    p3s_rand = np.zeros(len(p3s))
    for zz in range(len(p3s)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s_rand[zz] = ran
        #p3s10_rand[zz] = p3s10[zz] # not random
    #"""

    pl.rc("font", size=8)
    pl.rc("axes", linewidth=0.5)

    pl.figure()
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    #pl.subplot(2,1,1)
    pl.minorticks_on()
    pl.scatter(xi_xs, p3s, alpha=0.7, ec="None")
    pl.errorbar(xi_xs, p3s, fmt='none', yerr=ep3s, color="tab:blue", zorder=2, alpha=0.3)
    pl.scatter(xs_model, p3s_model, alpha=0.7, ec="None", color="tab:orange")
    pl.scatter(xs_model_hl, p3s_model_hl, alpha=0.7, ec="None", color="tab:orange")
    pl.plot(xt, yt, c="black")
    pl.plot(xthi, ythi, c="black", ls="--")
    #pl.plot(xt, ytl, c="pink", ls="--")
    #pl.plot(xt, yth, c="blue", ls="--")
    pl.fill_between(xt, ytl2, yth2, color="grey", alpha=0.3)
    #pl.plot(xt, ytl, c="pink", ls="--")
    #pl.plot(xt, yth, c="magenta", ls="--")
    pl.plot(10**xm, 10**ym, c="tab:red", lw=2)
    #pl.plot(10**xp, 10**yp)
    #pl.plot(10**xl, 10**yl)
    #pl.plot(10**xl2, 10**yl2)
    #pl.plot(10**xl3, 10**yl3, c="tab:red")
    pl.plot(10**xh, 10**yh, c="tab:blue", lw=2)
    #pl.plot(10**xh2, 10**yh2, c="pink")
    #pl.plot(10**xh3, 10**yh3, c="tab:green")
    #pl.plot(xline, yline)
    pl.axhline(y=2, ls=":", c="black")
    xlims = pl.xlim()
    pl.loglog()
    pl.xlabel(r"$(P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15})$")
    pl.ylabel(r"$P_3$ in $P$")
    pl.ylim(0.8, 400)
    #pl.axis("equal")
    filename = "output/p3edot_inflection_geoff.pdf"
    print(filename)
    pl.savefig(filename)

    pl.show()


def p3edot_inflection_geoff2(data, size=1e3):
    """ version with p3 fraction as error - the best one?! think so """

    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])

    p2s_theo = np.empty(shape=len(p3s))
    p2s = np.array(data[4])
    # check the sorting according to xi_xs below!

    xi_xs = []
    for i, p in enumerate(periods):
        val = p ** (-1.78) * (pdots[i]/1e-15)
        #val = (p ** (-1.78) * (pdots[i]/1e-15)) ** 0.5
        xi_xs.append(val)
        # p2
        # log(P2) = 0.14*( log(P^{-5.5}*Pdot) ) + 1.31
        p2 = 10 ** (0.14 * np.log10(p ** (-5.5) * (pdots[i]/1e-15)) + 1.31)
        p2s_theo[i] = p2
        p2s[i] = np.fabs(p2s[i])
    xi_xs = np.array(xi_xs)

    #"""
    # sort according to xi_xs
    idx = np.argsort(xi_xs)
    xi_xs = xi_xs[idx]
    p3s = p3s[idx]
    ep3s = ep3s[idx]
    edots = edots[idx]
    periods = periods[idx]
    pdots = pdots[idx]
    p2s = p2s[idx]
    p2s_theo = p2s_theo[idx]
    #"""

    # poarabolic fit to whole data
    para = lambda v, x: v[0] * x ** 2 + v[1] * x + v[2]
    xp, yp, vp = least_sq(np.log10(xi_xs), np.log10(p3s), para, [1,1,1], xmax=None)

    indxs = np.argsort(10**yp)
    in_ = 10 ** xp[indxs[0]]
    #in_ = 3

    print("Inflection point: ", in_)

    xi_low = []
    xi_high = []
    p3s_low = []
    p3s_high = []
    ep3s_low = []
    ep3s_high = []

    for i,xi in enumerate(xi_xs):
        if xi <= in_:
            xi_low.append(xi)
            p3s_low.append(p3s[i])
            ep3s_low.append(ep3s[i])
        else:
            xi_high.append(xi)
            p3s_high.append(p3s[i])
            ep3s_high.append(ep3s[i])

    # linear fit
    lin = lambda v, x:  v[0] * x + v[1]
    xh, yh, vh = least_sq(np.log10(xi_high), np.log10(p3s_high), lin, [1,1], xmax=None)
    #xh2, yh2, vh2, err = least_sq1D(np.log10(xi_high), np.log10(p3s_high), lin, np.log10(ep3s_high), [1,1])
    xh3, yh3, vh3, err3 = odr(np.log10(xi_high), np.log10(p3s_high), np.array([1 for i in range(len(ep3s_high))]), np.log10(ep3s_high),lin,  [1,1])

    xl, yl, vl = least_sq(np.log10(xi_low), np.log10(p3s_low), lin, [1,1], xmax=None)
    #xl2, yl2, vl2, err2 = least_sq1D(np.log10(xi_low), np.log10(p3s_low), lin, np.log10(ep3s_low), [1,1])
    xl3, yl3, vl3, err3 = odr(np.log10(xi_low), np.log10(p3s_low), np.array([1 for i in range(len(ep3s_low))]), np.log10(ep3s_low),lin,  [1,1])

    ############################################################################
    # simulating points
    # alised values
    #p3aliased = lambda p3, n: p3 / (1  - n * p3)
    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            #elif p3 == 1:
            #    return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                #elif p3_ == 0:
                #    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    #"""
    def p3fun(xs, w, th):
        if type(xs) == np.float64:
            if xs != th:
                return np.fabs(1 / (1 - w / (xs - th)))
            else:
                return 1.
        else:
            res = []
            for x in xs:
                if x != th:
                    res.append(np.fabs(1 / (1 - w / (x - th))))
                else:
                    res.append(1)
            return np.array(res)
    #"""

    w = 1.0 # 1.0 # 1.5 #1 # 0.767 #1 #9.02 #1.271 #1 #1 # 0.5 # 6.32 # 9.65
    #w = 1.2
    th = -1.05 # -1.03 # -1.55 #-1.05#-1.05 # -0.791 #-1.057 #-3.02 #-1.373 #-1.05 #-1 # -0.42 # -3.16 # -3.16

    # intrinsic dependence
    xt = np.logspace(np.log10(1e-3), np.log10(1e4), num=1000)
    yt = p3fun(xt, w, th)

    # get minimum / maximum value
    p3m_max = np.max(yt)
    p3m_min = np.min(yt)

    # intrinsic dependence in high xi region
    ythi = []
    xthi = []
    for i, y in enumerate(yt):
        n_ = 1
        ym = y
        while y < 2:
            y = np.fabs(p3aliased(ym, n_))
            n_ += 1
            if n_ == 100:
                break
            #print(y, n_)
        else:
            xthi.append(xt[i])
            ythi.append(y)
    xthi = np.array(xthi)
    ythi = np.array(ythi)

    #pl.plot(xt, yt)
    #pl.loglog()
    #pl.show()
    #exit()

    fraction_fun = lambda v, x: v[0] + v[1] * x
    v0 = [1, 1]
    xpoints = np.array([-3, 4])
    ypoints = np.array([0.5, 0.01])
    xsig, ysig, vsig = least_sq(xpoints, ypoints, fraction_fun, v0, xmax=None)
    #pl.plot(xsig, ysig)
    #pl.scatter(xpoints, ypoints)
    #pl.show()
    #exit()
    # intrinsic dependence variation # TODO improve it!
    ytl = np.empty(len(yt))
    yth = np.empty(len(yt))

    # generate model + sigma_fun
    for i, y in enumerate(yt):
        ytl[i] = 10 ** (np.log10(y) - np.abs(np.log10(y) - np.log10((1 - fraction_fun(vsig, np.log10(xt[i]))) * y)))
        yth[i] = 10 ** (np.log10(y) + np.abs(np.log10(y) - np.log10((1 - fraction_fun(vsig, np.log10(xt[i]))) * y)))


    # model data
    p3s_notobs = []
    xs_notobs = []
    p3s_model = [] #np.empty(len(p3s))
    xs_model = [] #np.empty(len(p3s))
    p3s_model_hl = []
    xs_model_hl = []

    p3max = np.max(p3s)

    for i in range(len(xi_xs)):
        aliased = False
        p3m = p3fun(xi_xs[i], w, th)
        p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(p3m) - np.log10((1 - fraction_fun(vsig, np.log10(xi_xs[i])))*p3m)))
        #if xi_xs[i] < 1:
        #    p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(p3m) - np.log10(0.5*p3m)))
        #else:
        #    p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(p3m) - np.log10(0.85*p3m)))

        #print(p3, p3m)
        #p3 = p3m
        if p3 > 2:
            p3s_model.append(p3)
            xs_model.append(xi_xs[i])
        else:
            # check aliasing first
            for nn in range(1, 10):
                p3obs = p3aliased(p3, nn)
                if p3obs > 2:
                    aliased = True
                    break
            # generate new p3
            while p3obs < 2 or p3obs > 2 * p3max:
                p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(p3m) - np.log10((1 - fraction_fun(vsig, np.log10(xi_xs[i])))*p3m)))
                if p3 > 2:
                    p3obs = p3
                    aliased = False
                else:
                    for nn in range(1, 4):
                        p3obs = p3aliased(p3, nn)
                        if p3obs > 2:
                            aliased = True
                            break
            if aliased is True:
                p3s_notobs.append(p3)
                xs_notobs.append(xi_xs[i])
            #print(p3obs)
            p3s_model.append(p3obs)
            xs_model.append(xi_xs[i])

        if xi_xs[i] > in_:
            p3s_model_hl.append(p3s_model[-1])
            xs_model_hl.append(xi_xs[i])

    p3s_model_hl = np.array(p3s_model_hl)
    xs_model_hl = np.array(xs_model_hl)
    xm, ym, vm = least_sq(np.log10(xs_model_hl), np.log10(p3s_model_hl), lin, [1,1], xmax=None)

    # P2 fits
    xp2m, yp2m, vp2m = least_sq(np.log10(xi_xs),  np.log10(p2s_theo), lin, [1, 1], xmax=None)
    xp2s, yp2s, vp2 = least_sq(np.log10(xi_xs),  np.log10(p2s), lin,  [1, 1], xmax=None)

    """
    # randomize observations
    p3s_rand = np.zeros(len(p3s))
    for zz in range(len(p3s)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s_rand[zz] = ran
        #p3s10_rand[zz] = p3s10[zz] # not random
    #"""

    pl.rc("font", size=8)
    pl.rc("axes", linewidth=0.5)

    pl.figure()
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    #pl.subplot(2,1,1)
    pl.minorticks_on()
    pl.scatter(xi_xs, p3s, alpha=0.7, ec="None")
    pl.errorbar(xi_xs, p3s, fmt='none', yerr=ep3s, color="tab:blue", zorder=2, alpha=0.3)
    pl.scatter(xs_model, p3s_model, alpha=0.7, ec="None", color="tab:orange")
    pl.scatter(xs_model_hl, p3s_model_hl, alpha=0.7, ec="None", color="tab:orange")
    pl.plot(xt, yt, c="black")
    pl.plot(xthi, ythi, c="black", ls="--")
    pl.plot(xt, ytl, c="pink", ls="--")
    pl.plot(xt, yth, c="magenta", ls="--")
    pl.fill_between(xt, ytl, yth, color="grey", alpha=0.3)
    pl.plot(10**xm, 10**ym, c="tab:red", lw=2)
    #pl.plot(10**xp, 10**yp)
    #pl.plot(10**xl, 10**yl)
    #pl.plot(10**xl2, 10**yl2)
    #pl.plot(10**xl3, 10**yl3, c="tab:red")
    pl.plot(10**xh, 10**yh, c="tab:blue", lw=2)
    #pl.plot(10**xh2, 10**yh2, c="pink")
    #pl.plot(10**xh3, 10**yh3, c="tab:green")
    #pl.plot(xline, yline)
    pl.plot(xi_xs, p2s_theo, c="tab:green", lw=0.5)
    pl.plot(10 **xp2m, 10**yp2m, c="tab:green", lw=1)
    pl.plot(10**xp2s, 10**yp2s, c="tab:brown", lw=1)
    pl.scatter(xi_xs, p2s, c="tab:green", ec="None", s=7)
    pl.axhline(y=2, ls=":", c="black")
    xlims = pl.xlim()
    pl.loglog()
    pl.xlabel(r"$(P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15})$")
    #pl.xlabel(r"$((P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15}))^{1/2}$")
    pl.ylabel(r"$P_3$ in $P$")
    pl.ylim(0.8, 400)
    #pl.axis("equal")
    filename = "output/p3edot_inflection_geoff2.pdf"
    print(filename)
    pl.savefig(filename)

    pl.show()


def aliasing():

    x1 = np.linspace(0, 1)
    f1 = x1
    f2 = 1 - x1

    x2 = np.linspace(1, 10)
    p1 = 1 / (1/x2)
    p2 = 1 - 1 / (1/x2)

    #pl.plot(x1, f1, c="C1")
    #pl.plot(x1, f2, c="C2")
    #pl.loglog()
    #pl.show()


    pl.plot(x2, p1, c="C1")
    pl.plot(x2, p2, c="C2")
    #pl.loglog()
    pl.show()



def p3p2edot(data, size=1e3):
    """  """

    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])
    p2s = np.fabs(np.array(data[4]))

    p2s_theo = np.empty(shape=len(p3s))
    rpc = np.empty(shape=len(p3s))

    xi_p2 = []
    xi_p3 = []
    for i, p in enumerate(periods):
        val = p ** (-5.5) * (pdots[i]/1e-15) # P2
        val2 = p ** (-1.78) * (pdots[i]/1e-15) # P3
        xi_p2.append(val)
        xi_p3.append(val2)
        p2 = 10 ** (0.14 * np.log10(p ** (-5.5) * (pdots[i]/1e-15)) + 1.31)
        p2s_theo[i] = p2
        rpc[i] = 150 * p ** (-0.5)
    xi_p2 = np.array(xi_p2)
    xi_p3 = np.array(xi_p3)

    #"""
    # sort according to xi_xs
    idx = np.argsort(xi_p3)
    xi_p2 = xi_p2[idx]
    xi_p3 = xi_p3[idx]
    p3s = p3s[idx]
    ep3s = ep3s[idx]
    edots = edots[idx]
    periods = periods[idx]
    pdots = pdots[idx]
    p2s = p2s[idx]
    p2s_theo = p2s_theo[idx]
    rpc = rpc[idx]
    #"""

    # poarabolic fit
    para = lambda v, x: v[0] * x ** 2 + v[1] * x + v[2]
    xpc, ypc, vpc = least_sq_samex(np.log10(xi_p3), np.log10(rpc), para, [1,1, 1])


    # linear fit
    lin = lambda v, x:  v[0] * x + v[1]

    # log(P2) = 0.14*( log(P^{-5.5}*Pdot) ) + 1.31
    xp, yp, vp = least_sq(np.log10(xi_p2), np.log10(p2s), lin, [1,1], xmax=None) # xi_p2 vs.
    xp2, yp2, vp2 = least_sq(np.log10(xi_p3), np.log10(p2s), lin, [1,1], xmax=None) # xi_p2 vs. real data
    xp3, yp3, vp3 = least_sq_samex(np.log10(xi_p3), np.log10(p2s_theo), lin, [1,1]) # xi_p3 vs. Xiaoxi's fit
    #print(vp)
    #print(vp2)
    print(vp3)

    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            #elif p3 == 1:
            #    return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                #elif p3_ == 0:
                #    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    def p3fun(xs, w, th):
        if type(xs) == np.float64:
            if xs != th:
                return np.fabs(1 / (1 - w / (xs - th)))
            else:
                return 1.
        else:
            res = []
            for x in xs:
                if x != th:
                    res.append(np.fabs(1 / (1 - w / (x - th))))
                else:
                    res.append(1)
            return np.array(res)
    w = 1.0 # 1.0 # 1.5 #1 # 0.767 #1 #9.02 #1.271 #1 #1 # 0.5 # 6.32 # 9.65
    th = -1.05 # -1.03 # -1.55 #-1.05#-1.05 # -0.791 #-1.057 #-3.02 #-1.373 #-1.05 #-1 # -0.42 # -3.16 # -3.16

    # intrinsic dependence
    xt = xi_p3 #np.logspace(np.log10(1e-3), np.log10(1e4), num=1000)
    yt = p3fun(xt, w, th)

    # intrinsic dependence in high xi region
    ythi = []
    xthi = []
    for i, y in enumerate(yt):
        n_ = 1
        ym = y
        while y < 2:
            y = np.fabs(p3aliased(ym, n_))
            n_ += 1
            if n_ == 100:
                break
            #print(y, n_)
        else:
            xthi.append(xt[i])
            ythi.append(y)
    xthi = np.array(xthi)
    ythi = np.array(ythi)

    # drift rate
    #dr = p2s_theo / yt  # P2 / P3
    dr = 10 ** yp3 / yt  # P2 / P3

    p2s_true = np.empty(shape=len(p3s))
    p2s_nsp = np.empty(shape=len(p3s))
    p2s_nsp_measured = np.empty(shape=len(p3s))
    dr2 = np.empty(shape=len(p3s))

    # true P2 (not measured) based on fit and drift rate..
    for i in range(len(p2s)):
        p2s_true[i] = 10**yp3[i] * (1 - dr[i] / 360)

    nsp_fun = lambda v, x: v[0] * x + v[1]
    xpoints = [-3, np.log10(7000)]
    ypoints = [50,  2]
    xnsp, ynsp, vnsp = least_sq(xpoints, ypoints, nsp_fun, [1, 1])

    for i in range(len(p2s)):
        p2s_nsp[i] = 360 / nsp_fun(vnsp, np.log10(xi_p3[i]))
        dr2[i] = p2s_nsp[i] / yt[i]
        p2s_nsp_measured[i] = p2s_nsp[i] / (1 - dr2[i] / 360)

    pl.rc("font", size=8)
    pl.rc("axes", linewidth=0.5)

    pl.figure()
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    #pl.subplot(2,1,1)
    pl.minorticks_on()
    #pl.scatter(xi_p2, p2s, alpha=0.7, ec="None")
    pl.scatter(xi_p3, p3s, alpha=0.7, ec="None", c="tab:blue", s=10)
    #pl.plot(10**xp, 10**yp, c="tab:red") #  xi_p2 vs.
    #pl.plot(10**xp2, 10**yp2, c="tab:blue")

    # 0
    pl.plot(10**xp3, 10**yp3, c="tab:green", label=r"$P_2$ (fit)", zorder=99)
    pl.scatter(xi_p3, p2s_theo, alpha=0.7, ec="None", s=5, c="tab:green")
    pl.plot(xi_p3, p2s_true, c="green", label = r"$P_{2}^t =  P_2^m (1 - D / 360)$", ls="--")
    #pl.plot(xi_p3, dr, c="tab:red", label = "drift rate $(P_2 / P_3)$", lw=3, alpha=0.7)

    pl.plot(xt, yt, c="black", label=r"$P_3$")
    pl.plot(xthi, ythi, c="black", ls="--", lw=0.7)
    # drift rate

    # 2
    pl.scatter(xi_p3, rpc, c="C1", label = r"$R_{\rm pc}$", s=5)
    pl.plot(10 ** xpc, 10 ** ypc, c="C1")

    # 1
    pl.plot(xi_p3, p2s_nsp_measured, c="C9", label = r"$P_{2,n_{sp}}^m = P_{2,n_{sp}}^t / (1 - D / 360) $")
    pl.plot(xi_p3, dr2, c="C5", label = "drift rate", ls="-")
    pl.plot(xi_p3, p2s_nsp, c="C9", label = r"$P_2^t =  360 / n_{sp}$", ls="--")


    xlims = pl.xlim()
    pl.loglog()
    #pl.xlabel(r"$(P / 1 {\rm s})^{-5.5} (\dot{P} / 10^{-15})$")
    pl.xlabel(r"$(P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15})$")
    pl.ylabel(r"$P_3$, $P_2$")
    pl.legend(loc= "upper center")
    #pl.ylim(0.8, 400)
    #pl.axis("equal")
    filename = "output/p3p2edot.pdf"
    print(filename)
    pl.savefig(filename)

    pl.show()



def p2fun(data, size=1e3):
    """ playing with P2 """

    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])
    p2s = np.fabs(np.array(data[4]))

    p2s_theo = np.empty(shape=len(p3s))
    rpc = np.empty(shape=len(p3s))

    xi_p2 = []
    xi_p3 = []
    for i, p in enumerate(periods):
        val = p ** (-5.5) * (pdots[i]/1e-15) # P2
        val2 = p ** (-1.78) * (pdots[i]/1e-15) # P3
        xi_p2.append(val)
        xi_p3.append(val2)
        p2 = 10 ** (0.14 * np.log10(p ** (-5.5) * (pdots[i]/1e-15)) + 1.31)
        p2s_theo[i] = p2
        rpc[i] = 150 * p ** (-0.5)
    xi_p2 = np.array(xi_p2)
    xi_p3 = np.array(xi_p3)

    #"""
    # sort according to xi_xs
    idx = np.argsort(xi_p3)
    xi_p2 = xi_p2[idx]
    xi_p3 = xi_p3[idx]
    p3s = p3s[idx]
    ep3s = ep3s[idx]
    edots = edots[idx]
    periods = periods[idx]
    pdots = pdots[idx]
    p2s = p2s[idx]
    p2s_theo = p2s_theo[idx]
    rpc = rpc[idx]
    #"""

    # poarabolic fit
    para = lambda v, x: v[0] * x ** 2 + v[1] * x + v[2]
    xpc, ypc, vpc = least_sq_samex(np.log10(xi_p3), np.log10(rpc), para, [1,1, 1])

    # linear fit
    lin = lambda v, x:  v[0] * x + v[1]

    xp, yp, vp = least_sq(np.log10(xi_p2), np.log10(p2s), lin, [1,1], xmax=None) # xi_p2 vs.
    xp2, yp2, vp2 = least_sq(np.log10(xi_p3), np.log10(p2s), lin, [1,1], xmax=None) # xi_p3 vs.

    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            #elif p3 == 1:
            #    return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                #elif p3_ == 0:
                #    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    def p3fun(xs, w, th):
        if type(xs) == np.float64:
            if xs != th:
                return np.fabs(1 / (1 - w / (xs - th)))
            else:
                return 1.
        else:
            res = []
            for x in xs:
                if x != th:
                    res.append(np.fabs(1 / (1 - w / (x - th))))
                else:
                    res.append(1)
            return np.array(res)
    w = 1.0 # 1.0 # 1.5 #1 # 0.767 #1 #9.02 #1.271 #1 #1 # 0.5 # 6.32 # 9.65
    th = -1.05 # -1.03 # -1.55 #-1.05#-1.05 # -0.791 #-1.057 #-3.02 #-1.373 #-1.05 #-1 # -0.42 # -3.16 # -3.16

    # intrinsic dependence
    xt = xi_p3 #np.logspace(np.log10(1e-3), np.log10(1e4), num=1000)
    yt = p3fun(xt, w, th)

    # intrinsic dependence in high xi region
    ythi = []
    xthi = []
    for i, y in enumerate(yt):
        n_ = 1
        ym = y
        while y < 2:
            y = np.fabs(p3aliased(ym, n_))
            n_ += 1
            if n_ == 100:
                break
            #print(y, n_)
        else:
            xthi.append(xt[i])
            ythi.append(y)
    xthi = np.array(xthi)
    ythi = np.array(ythi)


    def gamma_fun(alpha, beta, phi):
        return np.arccos(np.cos(alpha) * np.cos(alpha + beta) + np.sin(alpha) * np.sin(alpha+beta) * np.cos(phi))

    def sigma_fun(alpha, beta, phi):
        gamma = gamma_fun(alpha, beta, phi)
        arg = np.sin(alpha+beta) * np.sin(phi) / np.sin(gamma)
        return np.arcsin(arg)


    ############################################################################
    ############################################################################
    ############################################################################
    # geometry influance on P2
    #"""
    alphas = np.deg2rad([0, 2.5, 0.8, 28, 3, 45]) # inclination angle
    betas = np.deg2rad([12, 1, 2.8, 4.7, -7, 7]) # impact parameter

    phis = np.linspace(0, np.deg2rad(120), num=500) # phase
    sigmas = np.zeros([len(alphas), len(phis)])

    for i in range(len(alphas)):
        alpha = alphas[i]
        beta = betas[i]
        for j,phi in enumerate(phis):
            #print(i, j)
            si = sigma_fun(alpha, beta, phi)
            # nesty hack to avoid breaks does not work for for all?
            #"""
            if j > 1:
                dsi = si - sigmas[i, j-1]
                dsi_pr = sigmas[i, j-1] - sigmas[i, j-2]
                if np.sign(dsi) == np.sign(dsi_pr):
                    sigmas[i, j] = si
                else:
                    nsi = np.pi - si
                    if np.fabs(sigmas[i, j-1] - nsi) < np.deg2rad(5):
                        sigmas[i, j] = nsi
                    else:
                        sigmas[i, j] = si
            else:
                sigmas[i, j] = si
            #"""

    fig = pl.figure()
    ax = fig.add_subplot(111)
    pl.grid()
    pl.minorticks_on()
    pl.plot(np.rad2deg(phis), np.rad2deg(sigmas[0]), label=r"$\alpha=0$", ls="--", c="black")
    pl.plot(np.rad2deg(phis), np.rad2deg(sigmas[1]), label=r"$\alpha=2.5$", c="tab:red")
    pl.plot(np.rad2deg(phis), np.rad2deg(sigmas[2]), label=r"$\alpha=0.8$", c="tab:green")
    pl.plot(np.rad2deg(phis), np.rad2deg(sigmas[3]), label=r"$\alpha=28$", c="tab:blue")
    pl.plot(np.rad2deg(phis), np.rad2deg(np.fabs(sigmas[4])), label=r"$\alpha=3$", c="C1")
    #pl.plot(np.rad2deg(phis), np.rad2deg(np.fabs(sigmas[5])), label=r"$\alpha=45$", c="C4")
    ax.set_aspect('equal')
    pl.legend()
    pl.xlabel("pulse phase ($\phi$)")
    pl.ylabel("azimuthal angle ($\sigma$)")
    pl.show()
    #"""
    ############################################################################
    ############################################################################
    ############################################################################



    ############################################################################
    ############################################################################
    ############################################################################
    # calculations for individual pulsar
    alpha = np.deg2rad(30)
    beta = np.deg2rad(7)
    p = 1 # pulsar period
    d = 50 # emission height [in stellar radius]

    theta_max = fun.theta_max(d, p) # coordinate of the last open field line
    rho = fun.rho_sy(theta_max) # opening angle
    pulse_width = fun.pulse_width(alpha, beta, rho)
    print("opening angle: ", np.rad2deg(rho))
    print("pulse width: ", np.rad2deg(pulse_width))

    phis = np.linspace(0, pulse_width/2, num=100)
    sigmas = np.zeros(len(phis))

    for i in range(len(phis)):
        si = sigma_fun(alpha, beta, phis[i])
        # FIX that!
        if i > 1:
            dsi = si - sigmas[i-1]
            dsi_pr = sigmas[i-1] - sigmas[i-2]
            if np.sign(dsi) == np.sign(dsi_pr):
                sigmas[i] = si
            else:
                nsi = np.pi - si
                if np.fabs(sigmas[i-1] - nsi) < np.deg2rad(5):
                    sigmas[i] = nsi
                else:
                    sigmas[i] = si
        else:
            sigmas[i] = si


    fig = pl.figure()
    ax = fig.add_subplot(111)
    pl.grid()
    pl.minorticks_on()
    pl.plot(np.rad2deg(phis), np.rad2deg(sigmas), c="red")
    #ax.set_aspect('equal')
    pl.legend()
    pl.xlabel("pulse phase ($\phi$)")
    pl.ylabel("azimuthal angle ($\sigma$)")
    pl.show()

    dphis = np.zeros(len(phis)-1)
    dsigmas = np.zeros(len(phis)-1)
    for i in range(len(phis)-1):
        dphis[i] = phis[i+1] - phis[i]
        dsigmas[i] = sigmas[i+1] - sigmas[i]


    fig = pl.figure()
    ax = fig.add_subplot(111)
    pl.grid()
    pl.minorticks_on()
    pl.plot(np.rad2deg(phis[:-1]), dsigmas/dphis, c="red")
    pl.axhline(y = np.mean(dsigmas/dphis), c="black", ls="--")
    pl.text(np.mean(np.rad2deg(phis)), np.mean(dsigmas/dphis), r"$\overline{\frac{d \sigma}{d \phi}}=$"+"{}".format(np.mean(dsigmas/dphis)))

    #ax.set_aspect('equal')
    pl.xlabel("pulse phase ($\phi$)")
    pl.ylabel("($d \sigma / d \phi$)")
    pl.show()


    #return

    pl.rc("font", size=8)
    pl.rc("axes", linewidth=0.5)

    pl.figure()
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    #pl.subplot(2,1,1)
    pl.minorticks_on()

    pl.scatter(xi_p3, p3s, alpha=0.7, ec="None", c="tab:blue", s=10, label=r"$P_3$")
    pl.scatter(xi_p3, p2s, alpha=0.7, ec="None", s=7, c="tab:green", label=r"$P_2$")

    pl.plot(xt, yt, c="black", label=r"$P_3$")
    pl.plot(xthi, ythi, c="black", ls="--", lw=0.7)
    pl.plot(10**xp2, 10**yp2, c="tab:green") #  xi_p2 vs.


    xlims = pl.xlim()
    pl.loglog()
    #pl.xlabel(r"$(P / 1 {\rm s})^{-5.5} (\dot{P} / 10^{-15})$")
    pl.xlabel(r"$(P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15})$")
    pl.ylabel(r"$P_3$, $P_2$")
    pl.legend(loc= "upper center")
    #pl.ylim(0.8, 400)
    #pl.axis("equal")
    filename = "output/p2fun.pdf"
    print(filename)
    pl.savefig(filename)

    pl.show()


def p3edot_final(data, size=1e3):
    """ version with p3 fraction as error - for the paper?"""

    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])
    p3s_rs = np.zeros(len(p3s))
    p3s_rs2 = np.zeros(len(p3s))

    # check the sorting according to xi_xs below!


    nsp_fun = lambda v, x: v[0] * x + v[1]
    xpoints = [-3, np.log10(7000)]
    ypoints = [50,  50]
    xnsp, ynsp, vnsp = least_sq(xpoints, ypoints, nsp_fun, [1, 1])


    xi_xs = []
    for i, p in enumerate(periods):
        val = p ** (-1.78) * (pdots[i]/1e-15)
        #val = (p ** (-1.78) * (pdots[i]/1e-15)) ** 0.5
        xi_xs.append(val)
        nsp =  nsp_fun(vnsp, np.log10(val))
        p3s_rs[i] = 20 / nsp * (edots[i] / 5e32) ** 0.5
        #p3s_rs2[i] = 5.6 / nsp * (edots[i] / 4e31) ** 0.5 # Basu
    xi_xs = np.array(xi_xs)

    #"""
    # sort according to xi_xs
    idx = np.argsort(xi_xs)
    xi_xs = xi_xs[idx]
    p3s = p3s[idx]
    ep3s = ep3s[idx]
    edots = edots[idx]
    periods = periods[idx]
    pdots = pdots[idx]
    p3s_rs = p3s_rs[idx]
    p3s_rs2 = p3s_rs2[idx]
    #"""

    # poarabolic fit to whole data
    para = lambda v, x: v[0] * x ** 2 + v[1] * x + v[2]
    xp, yp, vp = least_sq(np.log10(xi_xs), np.log10(p3s), para, [1,1,1], xmax=None)

    indxs = np.argsort(10**yp)
    in_ = 10 ** xp[indxs[0]]
    #in_ = 3

    print("Inflection point: ", in_)

    xi_low = []
    xi_high = []
    p3s_low = []
    p3s_high = []
    ep3s_low = []
    ep3s_high = []

    for i,xi in enumerate(xi_xs):
        if xi <= in_:
            xi_low.append(xi)
            p3s_low.append(p3s[i])
            ep3s_low.append(ep3s[i])
        else:
            xi_high.append(xi)
            p3s_high.append(p3s[i])
            ep3s_high.append(ep3s[i])

    # linear fit
    lin = lambda v, x:  v[0] * x + v[1]

    # RS model
    xrs, yrs, vrs = least_sq(np.log10(xi_xs), np.log10(p3s_rs), lin, [1,1], xmax=None)

    xh, yh, vh = least_sq(np.log10(xi_high), np.log10(p3s_high), lin, [1,1], xmax=None)
    #xh2, yh2, vh2, err = least_sq1D(np.log10(xi_high), np.log10(p3s_high), lin, np.log10(ep3s_high), [1,1])
    xh3, yh3, vh3, err3 = odr(np.log10(xi_high), np.log10(p3s_high), np.array([1 for i in range(len(ep3s_high))]), np.log10(ep3s_high),lin,  [1,1])

    xl, yl, vl = least_sq(np.log10(xi_low), np.log10(p3s_low), lin, [1,1], xmax=None)
    #xl2, yl2, vl2, err2 = least_sq1D(np.log10(xi_low), np.log10(p3s_low), lin, np.log10(ep3s_low), [1,1])
    xl3, yl3, vl3, err3 = odr(np.log10(xi_low), np.log10(p3s_low), np.array([1 for i in range(len(ep3s_low))]), np.log10(ep3s_low),lin,  [1,1])

    ############################################################################
    # simulating points
    # alised values
    #p3aliased = lambda p3, n: p3 / (1  - n * p3)
    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            #elif p3 == 1:
            #    return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                #elif p3_ == 0:
                #    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    #"""
    def p3fun(xs, w, th):
        if type(xs) == np.float64:
            if xs != th:
                return np.fabs(1 / (1 - w / (xs - th)))
            else:
                return 1.
        else:
            res = []
            for x in xs:
                if x != th:
                    res.append(np.fabs(1 / (1 - w / (x - th))))
                else:
                    res.append(1)
            return np.array(res)
    #"""

    w = 1.0 # 1.0 # 1.5 #1 # 0.767 #1 #9.02 #1.271 #1 #1 # 0.5 # 6.32 # 9.65
    #w = 1.2
    th = -1.03 # -1.03 # -1.55 #-1.05#-1.05 # -0.791 #-1.057 #-3.02 #-1.373 #-1.05 #-1 # -0.42 # -3.16 # -3.16

    # intrinsic dependence
    xt = np.logspace(np.log10(1e-3), np.log10(1e4), num=1000)
    yt = p3fun(xt, w, th)

    # get minimum / maximum value
    p3m_max = np.max(yt)
    p3m_min = np.min(yt)

    # intrinsic dependence in high xi region
    ythi = []
    xthi = []
    for i, y in enumerate(yt):
        n_ = 1
        ym = y
        while y < 2:
            y = np.fabs(p3aliased(ym, n_))
            n_ += 1
            if n_ == 100:
                break
            #print(y, n_)
        else:
            xthi.append(xt[i])
            ythi.append(y)
    xthi = np.array(xthi)
    ythi = np.array(ythi)

    #pl.plot(xt, yt)
    #pl.loglog()
    #pl.show()
    #exit()

    fraction_fun = lambda v, x: v[0] + v[1] * x
    v0 = [1, 1]
    xpoints = np.array([-3, 4])
    ypoints = np.array([0.55, 0.01])
    xsig, ysig, vsig = least_sq(xpoints, ypoints, fraction_fun, v0, xmax=None)
    #pl.plot(xsig, ysig)
    #pl.scatter(xpoints, ypoints)
    #pl.show()
    #exit()
    # intrinsic dependence variation # TODO improve it!?
    ytl = np.empty(len(yt))
    yth = np.empty(len(yt))

    # generate model + sigma_fun
    for i, y in enumerate(yt):
        ytl[i] = 10 ** (np.log10(y) - np.abs(np.log10(y) - np.log10((1 - fraction_fun(vsig, np.log10(xt[i]))) * y)))
        yth[i] = 10 ** (np.log10(y) + np.abs(np.log10(y) - np.log10((1 - fraction_fun(vsig, np.log10(xt[i]))) * y)))

    # intrinsic dependence new variation (circle stuff)
    ytl2 = np.zeros(len(yt))
    yth2 = np.zeros(len(yt))

    # circle data
    yci = [[] for i in range(len(xt))]
    for i,x in enumerate(xt):
        a = np.log10(x)
        b = np.log10(yt[i])
        r = np.fabs(np.log10(yt[i]) - np.log10((1 - fraction_fun(vsig, np.log10(x))) * yt[i])) # is it ok?
        #print(r)
        #print("a b r", a, b, r)
        # loop over the x
        ys = []
        xs = []
        for j, xc in enumerate(xt):
            if xc > 10 ** (a + r):
                break
            if xc >= 10 ** (a - r):
                y_1 = b - np.sqrt(-a ** 2 + 2 * a * np.log10(xc) + r ** 2 - np.log10(xc) ** 2)
                y_2 = np.sqrt(-a ** 2 + 2 * a * np.log10(xc) + r ** 2 - np.log10(xc) ** 2) + b
                ys.append(y_1)
                ys.append(y_2)
                xs.append(np.log10(xc))
                xs.append(np.log10(xc))
                yci[j].append(y_1)
                yci[j].append(y_2)
    for i in range(len(xt)):
        ytl2[i] = np.min(10 ** np.array(yci[i]))
        yth2[i] = np.max(10 ** np.array(yci[i]))
        #print(i, xt[i],  yt[i], ytl2[i], yth2[i])

    def get_sigma(xi):
        """ WOW this is weird - you have too much time?"""
        mi = 1e50
        sig = 10
        for i, x in enumerate(xt):
            diff = np.fabs(xi - x)
            if diff < mi:
                sig = (np.log10(yth2[i]) - np.log10(ytl2[i])) / 2
                mi = diff
        return sig

    # model data
    p3s_notobs = []
    xs_notobs = []
    p3s_model = [] #np.empty(len(p3s))
    xs_model = [] #np.empty(len(p3s))
    p3s_model_hl = []
    xs_model_hl = []

    p3max = np.max(p3s)

    for i in range(len(xi_xs)):
        aliased = False
        p3m = p3fun(xi_xs[i], w, th)
        p3 = 10 ** np.random.normal(np.log10(p3m), get_sigma(xi_xs[i]))
        #p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(p3m) - np.log10((1 - fraction_fun(vsig, np.log10(xi_xs[i])))*p3m)))
        #if xi_xs[i] < 1:
        #    p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(p3m) - np.log10(0.5*p3m)))
        #else:
        #    p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(p3m) - np.log10(0.85*p3m)))

        #print(p3, p3m)
        #p3 = p3m
        if p3 > 2:
            p3s_model.append(p3)
            xs_model.append(xi_xs[i])
        else:
            # check aliasing first
            for nn in range(1, 10):
                p3obs = p3aliased(p3, nn)
                if p3obs > 2:
                    aliased = True
                    break
            # generate new p3
            while p3obs < 2 or p3obs > 2 * p3max:
                p3 = 10 ** np.random.normal(np.log10(p3m), get_sigma(xi_xs[i]))
                #p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(p3m) - np.log10((1 - fraction_fun(vsig, np.log10(xi_xs[i])))*p3m)))
                if p3 > 2:
                    p3obs = p3
                    aliased = False
                else:
                    for nn in range(1, 4):
                        p3obs = p3aliased(p3, nn)
                        if p3obs > 2:
                            aliased = True
                            break
            if aliased is True:
                p3s_notobs.append(p3)
                xs_notobs.append(xi_xs[i])
            #print(p3obs)
            p3s_model.append(p3obs)
            xs_model.append(xi_xs[i])

        if xi_xs[i] > in_:
            p3s_model_hl.append(p3s_model[-1])
            xs_model_hl.append(xi_xs[i])

    p3s_model_hl = np.array(p3s_model_hl)
    xs_model_hl = np.array(xs_model_hl)
    xm, ym, vm = least_sq(np.log10(xs_model_hl), np.log10(p3s_model_hl), lin, [1,1], xmax=None)

    """
    # randomize observations
    p3s_rand = np.zeros(len(p3s))
    for zz in range(len(p3s)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s_rand[zz] = ran
        #p3s10_rand[zz] = p3s10[zz] # not random
    #"""

    pl.rc("font", size=8)
    pl.rc("axes", linewidth=0.5)

    pl.figure()
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    #pl.subplot(2,1,1)
    pl.minorticks_on()
    pl.scatter(xi_xs, p3s_rs, s=7, alpha=1.0, ec="None", c="tab:green", zorder=999)
    pl.plot(10**xrs, 10**yrs, c="tab:green", zorder=999)
    #pl.scatter(xi_xs, p3s_rs2, s=7, alpha=1.0, ec="None", c="pink", zorder=999)
    pl.scatter(xi_xs, p3s, alpha=0.7, ec="None")
    pl.errorbar(xi_xs, p3s, fmt='none', yerr=ep3s, color="tab:blue", zorder=2, alpha=0.3)
    pl.scatter(xs_model, p3s_model, alpha=0.7, ec="None", color="tab:orange")
    pl.scatter(xs_model_hl, p3s_model_hl, alpha=0.7, ec="None", color="tab:orange")
    pl.plot(xt, yt, c="black")
    pl.plot(xthi, ythi, c="black", ls="--")
    pl.fill_between(xt, ytl2, yth2, color="grey", alpha=0.3)
    pl.plot(10**xm, 10**ym, c="tab:red", lw=2)
    #pl.plot(10**xp, 10**yp)
    #pl.plot(10**xl, 10**yl)
    #pl.plot(10**xl2, 10**yl2)
    #pl.plot(10**xl3, 10**yl3, c="tab:red")
    pl.plot(10**xh, 10**yh, c="tab:blue", lw=2)
    #pl.plot(10**xh2, 10**yh2, c="pink")
    #pl.plot(10**xh3, 10**yh3, c="tab:green")
    #pl.plot(xline, yline)
    pl.axhline(y=2, ls=":", c="black")
    xlims = pl.xlim()
    pl.loglog()
    pl.xlabel(r"$(P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15})$")
    #pl.xlabel(r"$((P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15}))^{1/2}$")
    pl.ylabel(r"$P_3$ in $P$")
    pl.ylim(0.8, 400)
    #pl.axis("equal")
    filename = "output/p3edot_final.pdf"
    print(filename)
    pl.savefig(filename)

    pl.show()


def p3edot_final2(data, size=1e3):
    """ version with p3 fraction as error - for the paper?"""

    p3s = np.array(data[0])
    ep3s = np.array(data[1])
    edots = np.array(data[9])
    periods = np.array(data[17])
    pdots = np.array(data[18])
    p3s_rs = np.zeros(len(p3s))
    p3s_rs2 = np.zeros(len(p3s))

    # check the sorting according to xi_xs below!


    nsp_fun = lambda v, x: v[0] * x + v[1]
    xpoints = [-3, np.log10(7000)]
    ypoints = [50,  50]
    xnsp, ynsp, vnsp = least_sq(xpoints, ypoints, nsp_fun, [1, 1])


    xi_xs = []
    for i, p in enumerate(periods):
        val = p ** (-1.78) * (pdots[i]/1e-15)
        #val = (p ** (-1.78) * (pdots[i]/1e-15)) ** 0.5
        xi_xs.append(val)
        nsp =  nsp_fun(vnsp, np.log10(val))
        p3s_rs[i] = 20 / nsp * (edots[i] / 5e32) ** 0.5
        #p3s_rs2[i] = 5.6 / nsp * (edots[i] / 4e31) ** 0.5 # Basu
    xi_xs = np.array(xi_xs)

    #"""
    # sort according to xi_xs
    idx = np.argsort(xi_xs)
    xi_xs = xi_xs[idx]
    p3s = p3s[idx]
    ep3s = ep3s[idx]
    edots = edots[idx]
    periods = periods[idx]
    pdots = pdots[idx]
    p3s_rs = p3s_rs[idx]
    p3s_rs2 = p3s_rs2[idx]
    #"""

    # poarabolic fit to whole data
    para = lambda v, x: v[0] * x ** 2 + v[1] * x + v[2]
    xp, yp, vp = least_sq(np.log10(xi_xs), np.log10(p3s), para, [1,1,1], xmax=None)

    indxs = np.argsort(10**yp)
    in_ = 10 ** xp[indxs[0]]
    #in_ = 3

    print("Inflection point: ", in_)

    xi_low = []
    xi_high = []
    p3s_low = []
    p3s_high = []
    ep3s_low = []
    ep3s_high = []

    for i,xi in enumerate(xi_xs):
        if xi <= in_:
            xi_low.append(xi)
            p3s_low.append(p3s[i])
            ep3s_low.append(ep3s[i])
        else:
            xi_high.append(xi)
            p3s_high.append(p3s[i])
            ep3s_high.append(ep3s[i])

    # linear fit
    lin = lambda v, x:  v[0] * x + v[1]

    # RS model
    xrs, yrs, vrs = least_sq(np.log10(xi_xs), np.log10(p3s_rs), lin, [1,1], xmax=None)

    xh, yh, vh = least_sq(np.log10(xi_high), np.log10(p3s_high), lin, [1,1], xmax=None)
    #xh2, yh2, vh2, err = least_sq1D(np.log10(xi_high), np.log10(p3s_high), lin, np.log10(ep3s_high), [1,1])
    xh3, yh3, vh3, err3 = odr(np.log10(xi_high), np.log10(p3s_high), np.array([1 for i in range(len(ep3s_high))]), np.log10(ep3s_high),lin,  [1,1])

    xl, yl, vl = least_sq(np.log10(xi_low), np.log10(p3s_low), lin, [1,1], xmax=None)
    #xl2, yl2, vl2, err2 = least_sq1D(np.log10(xi_low), np.log10(p3s_low), lin, np.log10(ep3s_low), [1,1])
    xl3, yl3, vl3, err3 = odr(np.log10(xi_low), np.log10(p3s_low), np.array([1 for i in range(len(ep3s_low))]), np.log10(ep3s_low),lin,  [1,1])

    ############################################################################
    # simulating points
    # alised values
    #p3aliased = lambda p3, n: p3 / (1  - n * p3)
    def p3aliased(p3, n):
        if type(p3) == np.float64 or type(p3) == float:
            if p3 > 1:
                return p3 / (n * p3 - 1)
            #elif p3 == 1:
            #    return 100.
            else:
                return p3 / (1 - n * p3)
        else:
            res = []
            for p3_ in p3:
                if p3_ > 1:
                    res.append(p3_ / (n * p3_ - 1))
                #elif p3_ == 0:
                #    res.append(200)
                elif p3_ < 1:
                    res.append(p3 / (1 - n * p3))
            return np.array(res)

    #"""
    def p3fun(xs, w, th):
        if type(xs) == np.float64:
            if xs != th:
                return np.fabs(1 / (1 - w / (-xs - th)))
            else:
                return 1.
        else:
            res = []
            for x in xs:
                if x != th:
                    res.append(np.fabs(1 / (1 - w / (-x - th))))
                else:
                    res.append(1)
            return np.array(res)
    #"""

    w = 1.0 # 1.0 # 1.5 #1 # 0.767 #1 #9.02 #1.271 #1 #1 # 0.5 # 6.32 # 9.65
    #w = 1.2
    th = -1.03 # -1.03 # -1.55 #-1.05#-1.05 # -0.791 #-1.057 #-3.02 #-1.373 #-1.05 #-1 # -0.42 # -3.16 # -3.16

    # intrinsic dependence
    xt = np.logspace(np.log10(1e-3), np.log10(1e4), num=1000)
    yt = p3fun(xt, w, th)

    # get minimum / maximum value
    p3m_max = np.max(yt)
    p3m_min = np.min(yt)

    # intrinsic dependence in high xi region
    ythi = []
    xthi = []
    for i, y in enumerate(yt):
        n_ = 1
        ym = y
        while y < 2:
            y = np.fabs(p3aliased(ym, n_))
            n_ += 1
            if n_ == 100:
                break
            #print(y, n_)
        else:
            xthi.append(xt[i])
            ythi.append(y)
    xthi = np.array(xthi)
    ythi = np.array(ythi)

    #pl.plot(xt, yt)
    #pl.loglog()
    #pl.show()
    #exit()

    fraction_fun = lambda v, x: v[0] + v[1] * x
    v0 = [1, 1]
    xpoints = np.array([-3, 4])
    ypoints = np.array([0.55, 0.01])
    xsig, ysig, vsig = least_sq(xpoints, ypoints, fraction_fun, v0, xmax=None)
    #pl.plot(xsig, ysig)
    #pl.scatter(xpoints, ypoints)
    #pl.show()
    #exit()
    # intrinsic dependence variation # TODO improve it!?
    ytl = np.empty(len(yt))
    yth = np.empty(len(yt))

    # generate model + sigma_fun
    for i, y in enumerate(yt):
        ytl[i] = 10 ** (np.log10(y) - np.abs(np.log10(y) - np.log10((1 - fraction_fun(vsig, np.log10(xt[i]))) * y)))
        yth[i] = 10 ** (np.log10(y) + np.abs(np.log10(y) - np.log10((1 - fraction_fun(vsig, np.log10(xt[i]))) * y)))

    # intrinsic dependence new variation (circle stuff)
    ytl2 = np.zeros(len(yt))
    yth2 = np.zeros(len(yt))

    # circle data
    yci = [[] for i in range(len(xt))]
    for i,x in enumerate(xt):
        a = np.log10(x)
        b = np.log10(yt[i])
        r = np.fabs(np.log10(yt[i]) - np.log10((1 - fraction_fun(vsig, np.log10(x))) * yt[i])) # is it ok?
        #print(r)
        #print("a b r", a, b, r)
        # loop over the x
        ys = []
        xs = []
        for j, xc in enumerate(xt):
            if xc > 10 ** (a + r):
                break
            if xc >= 10 ** (a - r):
                y_1 = b - np.sqrt(-a ** 2 + 2 * a * np.log10(xc) + r ** 2 - np.log10(xc) ** 2)
                y_2 = np.sqrt(-a ** 2 + 2 * a * np.log10(xc) + r ** 2 - np.log10(xc) ** 2) + b
                ys.append(y_1)
                ys.append(y_2)
                xs.append(np.log10(xc))
                xs.append(np.log10(xc))
                yci[j].append(y_1)
                yci[j].append(y_2)
    for i in range(len(xt)):
        ytl2[i] = np.min(10 ** np.array(yci[i]))
        yth2[i] = np.max(10 ** np.array(yci[i]))
        #print(i, xt[i],  yt[i], ytl2[i], yth2[i])

    def get_sigma(xi):
        """ WOW this is weird - you have too much time?"""
        mi = 1e50
        sig = 10
        for i, x in enumerate(xt):
            diff = np.fabs(xi - x)
            if diff < mi:
                sig = (np.log10(yth2[i]) - np.log10(ytl2[i])) / 2
                mi = diff
        return sig

    # model data
    p3s_notobs = []
    xs_notobs = []
    p3s_model = [] #np.empty(len(p3s))
    xs_model = [] #np.empty(len(p3s))
    p3s_model_hl = []
    xs_model_hl = []

    p3max = np.max(p3s)

    for i in range(len(xi_xs)):
        aliased = False
        p3m = p3fun(xi_xs[i], w, th)
        p3 = 10 ** np.random.normal(np.log10(p3m), get_sigma(xi_xs[i]))
        #p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(p3m) - np.log10((1 - fraction_fun(vsig, np.log10(xi_xs[i])))*p3m)))
        #if xi_xs[i] < 1:
        #    p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(p3m) - np.log10(0.5*p3m)))
        #else:
        #    p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(p3m) - np.log10(0.85*p3m)))

        #print(p3, p3m)
        #p3 = p3m
        if p3 > 2:
            p3s_model.append(p3)
            xs_model.append(xi_xs[i])
        else:
            # check aliasing first
            for nn in range(1, 10):
                p3obs = p3aliased(p3, nn)
                if p3obs > 2:
                    aliased = True
                    break
            # generate new p3
            while p3obs < 2 or p3obs > 2 * p3max:
                p3 = 10 ** np.random.normal(np.log10(p3m), get_sigma(xi_xs[i]))
                #p3 = 10 ** np.random.normal(np.log10(p3m), np.abs(np.log10(p3m) - np.log10((1 - fraction_fun(vsig, np.log10(xi_xs[i])))*p3m)))
                if p3 > 2:
                    p3obs = p3
                    aliased = False
                else:
                    for nn in range(1, 4):
                        p3obs = p3aliased(p3, nn)
                        if p3obs > 2:
                            aliased = True
                            break
            if aliased is True:
                p3s_notobs.append(p3)
                xs_notobs.append(xi_xs[i])
            #print(p3obs)
            p3s_model.append(p3obs)
            xs_model.append(xi_xs[i])

        if xi_xs[i] > in_:
            p3s_model_hl.append(p3s_model[-1])
            xs_model_hl.append(xi_xs[i])

    p3s_model_hl = np.array(p3s_model_hl)
    xs_model_hl = np.array(xs_model_hl)
    xm, ym, vm = least_sq(np.log10(xs_model_hl), np.log10(p3s_model_hl), lin, [1,1], xmax=None)

    """
    # randomize observations
    p3s_rand = np.zeros(len(p3s))
    for zz in range(len(p3s)):
        ran = np.random.normal(p3s[zz], ep3s[zz])
        while ran < 2:
            ran = np.random.normal(p3s[zz], ep3s[zz])
        p3s_rand[zz] = ran
        #p3s10_rand[zz] = p3s10[zz] # not random
    #"""

    pl.rc("font", size=8)
    pl.rc("axes", linewidth=0.5)

    pl.figure()
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    #pl.subplot(2,1,1)
    pl.minorticks_on()
    pl.scatter(xi_xs, p3s_rs, s=7, alpha=1.0, ec="None", c="tab:green", zorder=999)
    pl.plot(10**xrs, 10**yrs, c="tab:green", zorder=999)
    #pl.scatter(xi_xs, p3s_rs2, s=7, alpha=1.0, ec="None", c="pink", zorder=999)
    pl.scatter(xi_xs, p3s, alpha=0.7, ec="None")
    pl.errorbar(xi_xs, p3s, fmt='none', yerr=ep3s, color="tab:blue", zorder=2, alpha=0.3)
    pl.scatter(xs_model, p3s_model, alpha=0.7, ec="None", color="tab:orange")
    pl.scatter(xs_model_hl, p3s_model_hl, alpha=0.7, ec="None", color="tab:orange")
    pl.plot(xt, yt, c="black")
    pl.plot(xthi, ythi, c="black", ls="--")
    pl.fill_between(xt, ytl2, yth2, color="grey", alpha=0.3)
    pl.plot(10**xm, 10**ym, c="tab:red", lw=2)
    #pl.plot(10**xp, 10**yp)
    #pl.plot(10**xl, 10**yl)
    #pl.plot(10**xl2, 10**yl2)
    #pl.plot(10**xl3, 10**yl3, c="tab:red")
    pl.plot(10**xh, 10**yh, c="tab:blue", lw=2)
    #pl.plot(10**xh2, 10**yh2, c="pink")
    #pl.plot(10**xh3, 10**yh3, c="tab:green")
    #pl.plot(xline, yline)
    pl.axhline(y=2, ls=":", c="black")
    xlims = pl.xlim()
    pl.loglog()
    pl.xlabel(r"$(P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15})$")
    #pl.xlabel(r"$((P / 1 {\rm s})^{-1.78} (\dot{P} / 10^{-15}))^{1/2}$")
    pl.ylabel(r"$P_3$ in $P$")
    pl.ylim(0.8, 400)
    #pl.axis("equal")
    filename = "output/p3edot_final2.pdf"
    print(filename)
    pl.savefig(filename)

    pl.show()
