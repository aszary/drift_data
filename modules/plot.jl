module Plot
    using PyPlot
    using DataFrames, Query


    function p3_edot(df; p3_key="P3", mod="1")
        #println(df)
        #println(names(df))

        q1 = @from i in df begin
            @where i.drift_dir == "ND"
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end


        q2 = @from i in df begin
            @where i.drift_dir == "PD"
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end


        q3 = @from i in df begin
            @where (i.drift_dir != "PD") && (i.drift_dir != "ND")
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end

        #println(q1)
        #println(q2)


        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 1.9465685427418518), frameon=true)  # 8cm x  golden ratio
        subplots_adjust(left=0.17, bottom=0.21, right=0.99, top=0.99, wspace=0., hspace=0.)
        scatter(q1[:EDOT], q1[p3_key], label="ND")
        scatter(q2[:EDOT], q2[p3_key], label="PD")
        scatter(q3[:EDOT], q3[p3_key], label="PD/ND")
        for i in 1:length(q1[:EDOT])
            text(q1[:EDOT][i], q1[p3_key][i], q1[:PSRJ][i], size=5)
        end
        loglog()
        legend()
        #xlim(9e29, 1.1e33)
        #ylim(0.9, 110)
        xlabel("\$\\dot{E}\$ (ergs/s)")
        ylabel("\$P_3\$ (in units \$P\$)")
        savefig("outdir/p3_edot_$mod.pdf")
        close()

    end


end  # module Plot
