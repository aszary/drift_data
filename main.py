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


def plot_p():
    d1, d2 = da.drifting_p3only()
    pl.p3_p([d1, d2], ["drifting", "P3only"])


def main():
    #test()
    #latex_test()
    #latex_test2()
    #plot_edot()
    plot_edot2()
    #plot_p()
    print("Bye")


if __name__ == "__main__":
    main()
