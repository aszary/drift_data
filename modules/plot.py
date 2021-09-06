import numpy as np
import matplotlib.pyplot as pl


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
    pl.savefig("output/p3_edot.pdf")
    #pl.show()


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


def p3_edot2(datas, labels):
    #print(da.info)

    p3s = []
    ep3s = []
    edots = []

    for da in datas:
        p3_, ep3_, edots_ = get_p3s_edots(da)
        p3s.append(p3_)
        ep3s.append(ep3_)
        edots.append(edots_)

    sets = len(datas)

    pl.rc("font", size=12)
    pl.rc("axes", linewidth=0.5)
    pl.rc("lines", linewidth=0.5)

    pl.figure(figsize=(7.086614, 4.38189))  # 18 cm x 11.13 cm # golden ratio
    pl.subplots_adjust(left=0.11, bottom=0.13, right=0.99, top=0.99)
    for i in range(sets):
        # first component and feature only
        pl.scatter(edots[i][0][0], p3s[i][0][0], color="C{}".format(i+1), label=labels[i], s=5, zorder=1)
        pl.errorbar(edots[i][0][0], p3s[i][0][0], fmt='none', yerr=np.array(ep3s[i][0][0]), color="C{}".format(i+1), zorder=2)
    pl.legend()
    pl.loglog()
    yl = pl.ylim()
    pl.ylim([1, yl[1]])
    pl.xlabel("$\dot{E}$ (ergs/s)")
    pl.ylabel(r"$P_3$ in $P$")
    pl.savefig("output/p3_edot2.pdf")
    #pl.show()
