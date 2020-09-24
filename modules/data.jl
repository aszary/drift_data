module Data
    using CSV
    using DataFrames
    using Query
    using Latexify
    using Printf

    include("functions.jl")


    """ get pulsar parameters using psrcat """
    function get_pars(PSRJs, pars; psrcatdir="/home/szary/software/psrcat_tar/")
        cmd = `$(psrcatdir)/psrcat -c $pars $PSRJs -o short_csv`
        run(pipeline(cmd, stdout="data/pulsars.csv"))
        df = CSV.read("data/pulsars.csv" ; comment="#")
        #df = CSV.read("data/pulsars.csv")
        insertcols!(df, 2, "PSRJ"=>PSRJs)
        # delete column 5!
        select!(df, Not([:Column5]))
        # rename EDOT
        rename!(df, Dict("(ergs/s)" => :EDOT))
        rename!(df, Dict("(s)" => :P0))
        rename!(df, Dict(" _1" => :PSRB))
        return df
    end

    function add_p3!(df)
        p3_rahul = [13, 2.04, 14.4, 3.1, 2.36, 11.1, 4.7, 3.9, 2.15, 2.49, 3.01, 17.5, 2.452, 9.7, 4.1, 19.6, 3.05, 6.6, 6.1, 13.8, 2.75, 2.05, 2.1]
        p3_rahul_err = [1, 0.08, 0.8, 0.1, 0.01, 0.1, 0.1, 0.2, 0.01, 0.03, 0.05, 3.6, 0.006, 1.6, 0.2, 1.6, 0.09, 0.6, 0.3, 0.7, 0.04, 0.05, 0.1]
        rahul_dir = ["PD", "ND", "PD", "ND", "PD", "PD", "PD", "PD", "ND", "ND", "ND", "PD", "PD", "PD", "PD", "PD", "PD", "PD", "ND", "PD", "ND", "PD/ND", "ND" ]
        insertcols!(df, 6, :P3=>p3_rahul)
        insertcols!(df, 7, :P3_err=>p3_rahul_err)
        insertcols!(df, 8, :drift_dir=>rahul_dir)
    end

    function add_p3_modes!(df)
        psrs = ["J0034-0721", "J0034-0721", "J1555-3134", "J1727-2739", "J1822-2256", "J1921+1948", "J1921+1948", "J1946+1805", "J2305+3100"] # "J1822-2256",
        p3_rahul = [6.5, 4.0, 10.2, 5.2, 10.7, 3.8, 2.45, 6.1, 3] # 14.3,
        p3_rahul_err = [0.5, 0.5, 1.0, 0.9, 1.8, 1.1, 0.1, 0.04, 1.7, 0]
        for (ii, psr) in enumerate(psrs)
            rec = @from i in df begin
                @where i.PSRJ == psr
                @select i
                @collect DataFrame
            end
            rec = rec[1, :] # one record only, but DataFrameRow!
            rec[:P3] = p3_rahul[ii]
            rec[:P3_err] = p3_rahul_err[ii]
            push!(df, collect(rec))
        end
    end

    """ simple LBC interpretation of aliasing """
    function lbc_p3!(df)

        df[:P3_LBC] = -7.0 # magic number
        for r in eachrow(df)
            if r.drift_dir == "ND"
                r.P3_LBC = abs(Functions.p3(r.P3, -1))
            elseif r.drift_dir == "PD"
                r.P3_LBC = r.P3
            else
                r.P3_LBC = r.P3
                #=
                # not the best idea below, but used for plotting!
                r.P3_LBC = r.P3
                r.drift_dir = "PD"
                # adding it twice (strange I know, but it is in Rahul's paper)
                # TODO this is bad!
                ro = collect(r)
                ro[end] = abs(Functions.p3(r.P3, -1))
                ro[8] = "ND"
                push!(df, ro)
                =#
            end
        end
        #println(df)
        # some tests here
        #=
        rec = @from i in df begin
            @where occursin("PD", i.drift_dir)
            #@where i.drift_dir == "PD"
            @select {i.P3, i.P0}
            #@select i
            @collect DataFrame
        end
        println(rec[:P3])
        P0s = rec[:P0]
        for (i, p3_obs) in enumerate(rec[:P3])
            println("Edwards")
            for j in [-3, -2, -1, 0, 1, 2, 3]
                p3 = Functions.p3_edwards(p3_obs, P0s[i], j)
                p3b = Functions.p3(p3_obs, j)
                println("n:$j P3:$p3")
                println("n:$j P3:$p3b")
            end
            println("Gupta")
            for j in 0:6
                p3 = Functions.p3_gupta(p3_obs, P0s[i], j)
                println("n:$j P3:$p3")
            end
        end
        =#
    end


    """ the MC model interpretation of aliasing """
    function mc_p3!(df)
        df[:P3_MC] = -7.0 # magic number
        df[:n] = -7 # magic number
        for r in eachrow(df)
            (p3, n) = Functions.p3_edot(r.P3, r.EDOT)
            r.P3_MC = p3
            r.n = n
            #=
            if r.PSRJ == "J0108+6608" # used to find best p3_prediction function
                r.n = 1
                r.P3_MC = abs(Functions.p3(r.P3, r.n))
            end
            =#
        end
    end


    function data()
        #df = DataFrame(PSRJ=String[], P0=Float64[], P1=Float64[], P3=Float64[], EDOT=Float64[])
        psrs_coherent = ["J0034-0721", "J0108+6608", "J0151-0635", "J0421-0345", "J0459-0210", "J0814+7429", "J0820-1350", "J0934-5249", "J0946+0951", "J1418-3921", "J1543-0620", "J1555-3134", "J1720-2933", "J1727-2739", "J1816-2650", "J1822-2256", "J1901-0906", "J1919+0134", "J1921+1948", "J1946+1805", "J2046-0421", "J2305+3100", "J2313+4253"]
        df = get_pars(psrs_coherent, "PSRB P0 EDOT")
        add_p3!(df)
        add_p3_modes!(df)
        lbc_p3!(df)
        mc_p3!(df)
        sort!(df, [order(:PSRJ)]) # sorting
        select!(df, Not(" "))  # removing column " "
        #println(names(df))
        #println(df)

        return df
    end

    function data2()
        psrs_switching = ["J0323+3944", "J0815+0939", "J0820-4114", "J1034-3224", "J1842-0359", "J1921+2153", "J2321+6024"]
        psrs_diffuse  = ["J0152-1637", "J0304+1932", "J0525+1115", "J0630-2834", "J0823+0159", "J0944-1354", "J0959-4809", "J1041-1942", "J1703-1846", "J1720-0212", "J1741-0840", "J1840-0840", "J2018+2839", "J2046+1540", "J2317+2149"]
        psrs_lowmixed= ["J0624-0424", "J0837+0610", "J0846-3533", "J1239+2453", "J1328-4921", "J1625-4048", "J1650-1654", "J1700-3312", "J1703-3241", "J1733-2228", "J1740+1311", "J1801-2920", "J1900-2600", "J1912+2104", "J2006-0807", "J2048-1616"]

        psrs = vcat(psrs_switching, psrs_diffuse, psrs_lowmixed)

        df = get_pars(psrs, "PSRB P0 EDOT")
        #push!(df, ["J0034-0721", -7, -7, -7, -7])
        println(df)

    end

    function latex_table!(df)
        # make proper strings
        df[:P0] = [@sprintf("%.3f", p0) for p0 in df[:P0]]
        df[:EDOT] = [@sprintf("%.1f", ed/1e30) for ed in df[:EDOT]]
        df[!, :Number] =  ["" for i in axes(df, 1)]
        df[!, :P3_STR] =  ["" for i in axes(df, 1)]

        # P3 with errors
        for r in eachrow(df)
            r.P3_STR = "\$ $(r.P3) \\pm $(r.P3_err) \$" #@sprintf("\$%.3f \\pm %.3f\$", r.P3, r.P3_err)
            r.P3_LBC = parse(Float64, @sprintf("%.7f", r.P3_LBC)[1:length("$(r.P3)")]) # you need a break Who writes like that?
            r.P3_MC = parse(Float64, @sprintf("%.7f", r.P3_MC)[1:length("$(r.P3)")]) # you need a break Who writes like that?
        end

        #reorder
        permutecols!(df, [:Number, :PSRJ, :PSRB, :P0, :EDOT, :P3_STR, :drift_dir, :P3_LBC, :n, :P3_MC])

        # cleaning
        num = 1
        for (i,r) in enumerate(eachrow(df))
            if i == 1
                r.Number = "$num"
                num += 1
            else
                #println(df[:PSRJ][i-1])
                #println(r.PSRJ)
                if (df[:PSRJ][i-1] == r.PSRJ)
                    r.PSRJ = ""
                    r.PSRB = ""
                    r.P0 = ""
                    r.EDOT = ""
                    r.drift_dir = ""
                elseif ((df[:PSRJ][i-1] == "") && (df[:PSRJ][i-2] == r.PSRJ))  # three modes tops
                    r.PSRJ = ""
                    r.PSRB = ""
                    r.P0 = ""
                    r.EDOT = ""
                    r.drift_dir = ""
                else
                    r.Number = "$num"
                    num += 1
                end
            end
        end
        # change column names
        rename!(df, :P0=>"\$P\$")
        rename!(df, :P3_STR=>"\$P_3\$")
        rename!(df, :P3_LBC=>"\$P_3^{\\rm LBC }\$")
        rename!(df, :P3_MC=>"\$P_3^{\\rm MC }\$")
        rename!(df, :EDOT=>"\$\\dot{E}\$")
        rename!(df, :drift_dir=>"Type")

        println(df)
        la = latexify(df, env=:table, latex=false)
        println(la)
    end


end  # module Data
