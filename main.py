import modules.data as da
import modules.plot as pl


def test():
    d = da.read_p0_edot("data/pulsars.csv")
    pl.test_plot(d)

def latex_test():
    d = da.read_highedot("data/stats.csv")

def latex_test2():
    d = da.read_highedot2("data/stats.csv")
    da.latexify(d)

def main():
    #test()
    #latex_test()
    latex_test2()
    print("Bye")

if __name__ == "__main__":
    main()
