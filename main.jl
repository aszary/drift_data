module DriftData
    using DataFrames, Query

    include("modules/data.jl")
    include("modules/plot.jl")


    function first_table()
        df = Data.data()
        #Plot.p3_edot(df; p3_key="P3_LBC", mod="rahul") # old generic function
        #Plot.p3_edot_raw(df)
        Plot.p3_edot_rahul(df)
        #Plot.p3_edot_andrzej(df)
        Data.latex_table!(df)
    end

    function second_table()
        #(a, b) = (14.7, -0.49) # andrzej2
        (a, b) = (-16.0, 0.49) # andrzej5
        df = Data.data(;a=a, b=b)
        df2 = Data.data2(;a=a, b=b)

        #println(df)
        #println(df2)

        #Plot.p3_edot_simple(df2; p3_key="P3", mod="others_raw") # old generic function
        Plot.p3_edot_andrzej1(df, df2)
        #Plot.p3_edot_andrzej2(df, df2, a, b)
        #Plot.p3_edot_andrzej3(df, df2)
        #Plot.p3_edot_andrzej4(df, df2)
        #Plot.p3_edot_andrzej5(df, df2, a, b)
        Data.latex_table2!(df2)
    end

    function main()
        #first_table()
        second_table()
        println("Bye")
    end

end  # module DriftData

dd = DriftData.main()
