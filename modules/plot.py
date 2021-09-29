import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable


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
        sc = pl.scatter(edots[i], p3s[i], c=dps[i], s=5, cmap=cmaps[i], zorder=1)
        pl.errorbar(edots[i], p3s[i], fmt='none', yerr=ep3s[i], color=colors[i], zorder=2, label=labels[i])
        #co = pl.colorbar(sc, shrink=0.5)
        #pl.clim([0,70])
    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([0.7, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    filename = "output/p3_edot3.pdf"
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
