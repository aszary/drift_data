module Plot
    using PyPlot
    using DataFrames, Query


    function p3_edot(df; mod="1")
        #println(df)
        #println(names(df))

        q1 = @from i in df begin
            @where i.drift_dir == "ND"
            @select {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end


        q2 = @from i in df begin
            @where i.drift_dir == "PD"
            @select {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end


        q3 = @from i in df begin
            @where (i.drift_dir != "PD") && (i.drift_dir != "ND")
            @select {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end


        println(q1)
        println(q2)


        figure()
        scatter(q1[:EDOT], q1[:P3], label="ND")
        scatter(q2[:EDOT], q2[:P3], label="PD")
        scatter(q3[:EDOT], q3[:P3], label="?")
        loglog()
        legend()
        xlim(9e29, 1.1e33)
        ylim(0.9, 110)
        savefig("outdir/p3_edot_$mod.pdf")
        close()

    end


end  # module Plot
