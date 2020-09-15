module DriftData
    include("modules/data.jl")
    include("modules/plot.jl")


    function main()
        df = Data.data()
        Plot.p3_edot(df; mod="raw")

        println("Bye")
    end

end  # module DriftData

dd = DriftData.main()
