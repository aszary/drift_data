module DriftData
    include("modules/data.jl")
    include("modules/plot.jl")


    function main()
        df = Data.data()
        Plot.p3_edot(df; mod="raw")
        Plot.p3_edot(df; p3_key="P3_LBC", mod="rahul")
        println("Bye")
    end

end  # module DriftData

dd = DriftData.main()
