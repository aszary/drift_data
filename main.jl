module DriftData
    include("modules/data.jl")
    include("modules/plot.jl")


    function first_table()
        df = Data.data()
        #Plot.p3_edot(df; p3_key="P3_LBC", mod="rahul") # old generic function
        Plot.p3_edot_raw(df)
        Plot.p3_edot_rahul(df)
        Plot.p3_edot_andrzej(df)
        Data.latex_table!(df)
    end

    function second_table()
        df = Data.data2()
    end

    function main()
        #first_table()
        second_table()
        println("Bye")
    end

end  # module DriftData

dd = DriftData.main()
