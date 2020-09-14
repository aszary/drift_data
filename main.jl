module DriftData
    include("modules/data.jl")


    function main()
        Data.save_data()

        println("Bye")
    end

end  # module DriftData

dd = DriftData.main()
