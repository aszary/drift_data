module Plot
    using PyPlot
    using DataFrames, Query, GLM

    include("functions.jl")


    function p3_edot(df; p3_key="P3", mod="1")
        #println(df)
        #println(names(df))

        q1 = @from i in df begin
            @where i.drift_dir == "PD"
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end

        q2 = @from i in df begin
            @where (i.drift_dir == "ND")
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end

        q3 = @from i in df begin
            @where (i.drift_dir != "PD") && (i.drift_dir != "ND")
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end

        qtest = @from i in df begin
            @where i.PSRJ != "J0421-0345"
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end


        #println(q1)
        #println(q2)

        x = collect(df[:EDOT])
        y = collect(df[p3_key])
        (co, (edots, p3s)) = Functions.fit_line_log10(x, y; show_=true)
        (co1, (edots1, p3s1)) = Functions.fit_line_log10(collect(q1[:EDOT]), collect(q1[p3_key]))
        (co2, (edots2, p3s2)) = Functions.fit_line_log10(collect(q2[:EDOT]), collect(q2[p3_key]); show_=true)
        #=
        data = DataFrame(X=x, Y=y)
        ols = lm(@formula(Y ~ X), data)
        println(ols)
        println(coef(ols))

        edots = 10 .^ (range(30; stop=33, length=100))
        f(x) = 10^coef(ols)[1] * x^coef(ols)[2]
        p3s = f.(edots)
        =#
        ms = 10

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 1.9465685427418518), frameon=true)  # 8cm x  golden ratio
        subplots_adjust(left=0.15, bottom=0.22, right=0.99, top=0.99, wspace=0., hspace=0.)
        scatter(q1[:EDOT], q1[p3_key], label="PD", marker="s", s=ms, c="C0")
        scatter(q2[:EDOT], q2[p3_key], label="ND", s=ms, c="C1")
        #scatter(q3[:EDOT], q3[p3_key], label="PD/ND", s=ms, c="C2", marker="^")
        scatter(q2[:EDOT], q2[p3_key], label="aliased", s=ms, marker="+", c="grey")
        plot(edots, p3s, c="black")
        plot(edots1, p3s1, c="C0", ls="--")
        plot(edots2, p3s2, c="C1", ls=":")
        #=
        for i in 1:length(q1[:EDOT])
            text(q1[:EDOT][i], q1[p3_key][i], q1[:PSRJ][i], size=5)
        end
        for i in 1:length(q2[:EDOT])
            text(q2[:EDOT][i], q2[p3_key][i], q2[:PSRJ][i], size=5)
        end
        =#
        loglog()
        legend(fontsize=7)
        #xlim(9e29, 1.1e33)
        #ylim(0.9, 110)
        yticks([1, 10], [1, 10])
        ylim(0.9, 30)
        xlabel("\$\\dot{E}\$ (ergs/s)")
        ylabel("\$P_3\$ (in units \$P\$)")
        savefig("outdir/p3_edot_$mod.pdf")
        close()

    end


    function p3_edot_raw(df;)

        q1 = @from i in df begin
            @where i.drift_dir == "PD"
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end

        q2 = @from i in df begin
            @where (i.drift_dir == "ND")
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end

        x = collect(df[:EDOT])
        y = collect(df[:P3])
        (co, (edots, p3s)) = Functions.fit_line_log10(x, y; show_=true)
        ms = 10

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 1.9465685427418518), frameon=true)  # 8cm x  golden ratio
        subplots_adjust(left=0.15, bottom=0.22, right=0.99, top=0.99, wspace=0., hspace=0.)
        scatter(q1[:EDOT], q1[:P3], label="PD", marker="s", s=ms, c="C0")
        scatter(q2[:EDOT], q2[:P3], label="ND", s=ms, c="C1")
        plot(edots, p3s, c="black")
        loglog()
        legend(fontsize=7)
        #xlim(9e29, 1.1e33)
        #ylim(0.9, 110)
        yticks([1, 10], [1, 10])
        ylim(0.9, 30)
        xlabel("\$\\dot{E}\$ (ergs/s)")
        ylabel("\$P_3\$ (in units \$P\$)")
        savefig("outdir/p3_edot_raw.pdf")
        close()

    end


    function p3_edot_rahul(df;)

        q1 = @from i in df begin
            @where i.drift_dir == "PD"
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end

        q2 = @from i in df begin
            @where (i.drift_dir == "ND") && i.PSRJ != "J0421-0345 ignore"
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end

        q3 = @from i in df begin
            @where (i.drift_dir != "PD") && (i.drift_dir != "ND")
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end

        x = collect(df[:EDOT])
        y = collect(df[:P3_LBC])
        (co, (edots, p3s)) = Functions.fit_line_log10(x, y; show_=true)
        (co1, (edots1, p3s1)) = Functions.fit_line_log10(collect(q1[:EDOT]), collect(q1[:P3_LBC]))
        (co2, (edots2, p3s2)) = Functions.fit_line_log10(collect(q2[:EDOT]), collect(q2[:P3_LBC]); show_=true)

        ms = 10

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 1.9465685427418518), frameon=true)  # 8cm x  golden ratio
        subplots_adjust(left=0.15, bottom=0.22, right=0.99, top=0.99, wspace=0., hspace=0.)
        scatter(q1[:EDOT], q1[:P3_LBC], label="PD", marker="s", s=ms, c="C0")
        scatter(q2[:EDOT], q2[:P3_LBC], label="ND", s=ms, c="C1")
        scatter(q2[:EDOT], q2[:P3_LBC], label="aliased", s=ms, marker="+", c="grey")
        plot(edots, p3s, c="black")
        plot(edots1, p3s1, c="C0", ls="--")
        plot(edots2, p3s2, c="C1", ls=":")
        loglog()
        legend(fontsize=7)
        #xlim(9e29, 1.1e33)
        #ylim(0.9, 110)
        yticks([1, 10], [1, 10])
        ylim(0.5, 30)
        xlabel("\$\\dot{E}\$ (ergs/s)")
        ylabel("\$P_3\$ (in units \$P\$)")
        savefig("outdir/p3_edot_rahul.pdf")
        close()
    end



    function p3_edot_andrzej(df; )


        q1 = @from i in df begin
            @where i.drift_dir == "PD"
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end

        q2 = @from i in df begin
            @where (i.drift_dir == "ND") && i.PSRJ != "J0421-0345 ignore"
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end

        q3 = @from i in df begin
            @where (i.n != 0)
            @select i # {i.P3, i.EDOT, i.PSRJ, i.drift_dir}
            @collect DataFrame
        end

        println(df)
        #println(q1)
        #println(q2)

        x = collect(df[:EDOT])
        y = collect(df[:P3_MC])
        (co, (edots, p3s)) = Functions.fit_line_log10(x, y; show_=true)
        (co1, (edots1, p3s1)) = Functions.fit_line_log10(collect(q1[:EDOT]), collect(q1[:P3_MC]))
        (co2, (edots2, p3s2)) = Functions.fit_line_log10(collect(q2[:EDOT]), collect(q2[:P3_MC]); show_=true)
        ms = 10

        rc("font", size=8.)
        rc("axes", linewidth=0.5)
        rc("lines", linewidth=0.5)

        figure(figsize=(3.14961, 1.9465685427418518), frameon=true)  # 8cm x  golden ratio
        subplots_adjust(left=0.15, bottom=0.22, right=0.99, top=0.99, wspace=0., hspace=0.)
        scatter(q1[:EDOT], q1[:P3_MC], label="PD", marker="s", s=ms, c="C0")
        scatter(q2[:EDOT], q2[:P3_MC], label="ND", s=ms, c="C1")
        scatter(q3[:EDOT], q3[:P3_MC], label="aliased", s=ms, marker="+", c="grey")
        plot(edots, p3s, c="black")
        plot(edots1, p3s1, c="C0", ls="--")
        plot(edots2, p3s2, c="C1", ls=":")
        loglog()
        legend(fontsize=7)
        yticks([1, 10], [1, 10])
        ylim(0.5, 30)
        xlabel("\$\\dot{E}\$ (ergs/s)")
        ylabel("\$P_3\$ (in units \$P\$)")
        savefig("outdir/p3_edot_andrzej.pdf")
        close()
    end



end  # module Plot
