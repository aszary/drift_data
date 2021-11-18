import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import ttest_1samp, ks_2samp

from modules.functions import least_sq
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
