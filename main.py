import modules.data as da
import modules.plot as pl


def test():
    d = da.read_p0_edot("data/pulsars.csv")
    pl.test_plot(d)


def main():
    test()
    print("Bye")

if __name__ == "__main__":
    main()
