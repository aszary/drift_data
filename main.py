#! /usr/bin/env python
import modules.data as da
import modules.plot as pl


def test():
    d = da.read_p0_edot("data/pulsars.csv")
    pl.test_plot(d)


def latex_test():
    d = da.read_highedot("data/stats.csv")


def latex_test2():
    d = da.high_edot2("data/stats.csv")
    da.latexify(d)
    d.show_in_browser(jsviewer=True)


def plot_edot():
    d1, d2 = da.drifting_p3only()
    pl.p3_edot(d1, d2)


def plot_edot2():
    d1, d2, d3 = da.positive_negative_mixed()
    pl.p3_edot2([d1, d2, d3], ["positive", "negative", "mixed"])


def plot_edot3():
    d1, d2, d3 = da.positive_negative_mixed3()
    pl.p3_edot3([d1, d2, d3], ["positive", "negative", "mixed"])


def plot_edot_sec():
    d1, d2, d3 = da.positive_negative_mixed3()
    pl.p3_edot_sec([d1, d2, d3], ["positive", "negative", "mixed"])


def plot_edot_rahul():
    d1, d2, d3 = da.positive_negative_mixed()
    pl.p3_edot_rahul([d1, d2, d3], ["positive", "negative", "mixed"])


def plot_edot_rahul2():
    d1, d2, d3 = da.positive_negative_mixed3()
    pl.p3_edot_rahul2([d1, d2, d3], ["positive", "negative", "mixed"])


def plot_p():
    d1, d2 = da.drifting_p3only()
    pl.p3_p([d1, d2], ["drifting", "P3only"])

def plot_pnew():
    d1, d2 = da.drifting_p3only2()
    pl.p3_pnew([d1, d2], ["drifting", "P3only"])
    pl.p3_psec([d1, d2], ["drifting", "P3only"])


def plot_p3mean_p():
    d1, d2 = da.drifting_p3only2()
    d3 = da.allfeaturesp3_xiaoxi()
    pl.p3_p_xcheck([d1, d2], d3)
    pl.p3mean_p([d1, d2], d3)


def plot_edotp():
    d1, d2, d3 = da.positive_negative_mixed3()
    pl.p3_edotp([d1, d2, d3], ["positive", "negative", "mixed"])


def generate_edot(edot_min=5e30, edot_max=2e31, fmod=""):
    """
        OBSOLETE. Please use check find_bestp3edot.
    """
    d1, d2, d3 = da.positive_negative_mixed3()
    pl.generate_p3_edot([d1, d2, d3], ["positive", "negative", "mixed"], edot_min=edot_min, edot_max=edot_max, fmod=fmod)


def generate_edot2(edot_min=5e30, edot_max=2e31, fmod=""):
    """
        OBSOLETE. Please use check find_bestp3edot.
    """
    # for dominant P3 features only
    d1 = da.p3dominant_driftonly()
    pl.generate_p3_edot2(d1, edot_min=edot_min, edot_max=edot_max, fmod=fmod)

def find_bestp3edot():
    d1 = da.p3dominant_driftonly()
    da.find_bestp3edot(d1)


def main():
    #test()
    #latex_test() # obsolete
    #latex_test2()
    #plot_edot()
    #plot_edot2()
    #plot_edot3()
    #plot_edot_rahul() # obsolete
    #plot_edot_rahul2()
    #plot_p() # obsolete
    #plot_pnew()
    #plot_edot_sec()
    #plot_p3mean_p()
    #plot_edotp()
    #generate_edot(edot_min=5e30, edot_max=2e31, fmod="_1")
    #generate_edot2(edot_min=5e30, edot_max=2e31, fmod="_1")
    #generate_edot2(edot_min=3e30, edot_max=2e31, fmod="_2")
    #da.table_1() # not used
    #da.table_2()
    find_bestp3edot()
    print("Bye")


if __name__ == "__main__":
    main()
