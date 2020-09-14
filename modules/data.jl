module Data
    using CSV
    using DataFrames


    """ get pulsar parameters using psrcat """
    function get_pars(PSRJs, pars; psrcatdir="/home/szary/software/psrcat_tar/")
        cmd = `$(psrcatdir)/psrcat -c $pars $PSRJs -o short_csv`
        run(pipeline(cmd, stdout="data/pulsars.csv"))
        df = CSV.read("data/pulsars.csv" ; comment="#")
        #df = CSV.read("data/pulsars.csv")
        insertcols!(df, 2, :PSRJ=>PSRJs)
        # delete column 5!
        select!(df, Not([:Column5]))
        return df
    end

    function add_p3!(df)
        p3_rahul = [13, 2.04, 14.4, 3.1, 2.36, 11.1, 4.7, 3.9, 2.15, 2.49, 3.01, 17.5, 2.452, 9.7, 4.1, 19.6, 3.05, 6.6, 6.1, 13.8, 2.75, 2.05, 2.1]
        p3_rahul_err = [1, 0.08, 0.8, 0.1, 0.01, 0.1, 0.1, 0.2, 0.01, 0.03, 0.05, 3.6, 0.006, 1.6, 0.2, 1.6, 0.09, 0.6, 0.3, 0.7, 0.04, 0.05, 0.1]
        rahul_dir = ["ND", "PD", "ND", "PD", "ND", "ND", "ND", "ND", "PD", "PD", "PD", "ND", "ND", "ND", "ND", "ND", "ND", "ND", "PD", "ND", "PD", "PD/ND", "PD" ]

        println(length(p3_rahul))
        println(length(p3_rahul_err))
        println(length(rahul_dir))
        insertcols!(df, 6, :P_3=>p3_rahul)
        insertcols!(df, 7, :P3_err=>p3_rahul_err)
        insertcols!(df, 8, :Drift_dir=>rahul_dir)
        println(df)
    end


    function save_data()
        #df = DataFrame(PSRJ=String[], P0=Float64[], P1=Float64[], P3=Float64[], EDOT=Float64[])
        # magic number -7
        psrs_coherent = ["J0034-0721", "J0108+6608", "J0151-0635", "J0421-0345", "J0459-0210", "J0814+7429", "J0820-1350", "J0934-5249", "J0946+0951", "J1418-3921", "J1543-0620", "J1555-3134", "J1720-2933", "J1727-2739", "J1816-2650", "J1822-2256", "J1901-0906", "J1919+0134", "J1921+1948", "J1946+1805", "J2046-0421", "J2305+3100", "J2313+4253"]
        df = get_pars(psrs_coherent, "PSRB P0 EDOT")
        add_p3!(df)


        #=
        psrs_switching = ["J0323+3944", "J0815+0939", "J0820-4114", "J1034-3224", "J1842-0359", "J1921+2153", "J2321+6024"]
        get_pars(psrs_switching, "PSRB P0 EDOT")

        psrs_diffuse  = ["J0152-1637", "J0304+1932", "J0525+1115", "J0630-2834", "J0823+0159", "J0944-1354", "J0959-4809", "J1041-1942", "J1703-1846", "J1720-0212", "J1741-0840", "J1840-0840", "J2018+2839", "J2046+1540", "J2317+2149"]
        get_pars(psrs_diffuse, "PSRB P0 EDOT")

        psrs_lowmixed= ["J0624-0424", "J0837+0610", "J0846-3533", "J1239+2453", "J1328-4921", "J1625-4048", "J1650-1654", "J1700-3312", "J1703-3241", "J1733-2228", "J1740+1311", "J1801-2920", "J1900-2600", "J1912+2104", "J2006-0807", "J2048-1616"]
        get_pars(psrs_lowmixed, "PSRB P0 EDOT")
        =#
        #push!(df, ["J0034-0721", -7, -7, -7, -7])
        #println(df)

    end


end  # module Data
