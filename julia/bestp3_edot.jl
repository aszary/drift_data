module Simulation
    using LsqFit
    using PyPlot
    using Statistics
    using Random, Distributions
    using HypothesisTests
    using PyCall
    @pyimport scipy.stats as stats


    function read_data(filename; l10=true, edot_min=nothing, edot_max=nothing)
        p3s = []
        ep3s = []
        edots = []
        f = open(filename)
        for line in readlines(f)
            r = split(line)
            edot = parse(Float64, r[3])
            if (edot_min != nothing) && (edot_max != nothing)
                if (edot >= edot_min) && (edot <= edot_max)
                    push!(p3s, parse(Float64, r[1]))
                    push!(ep3s, parse(Float64, r[2]))
                    push!(edots, parse(Float64, r[3]))
                end
            elseif edot_max != nothing
                if edot <= edot_max
                    push!(p3s, parse(Float64, r[1]))
                    push!(ep3s, parse(Float64, r[2]))
                    push!(edots, parse(Float64, r[3]))
                end
            elseif edot_min != nothing
                if edot >= edot_min
                    push!(p3s, parse(Float64, r[1]))
                    push!(ep3s, parse(Float64, r[2]))
                    push!(edots, parse(Float64, r[3]))
                end
            else
                push!(p3s, parse(Float64, r[1]))
                push!(ep3s, parse(Float64, r[2]))
                push!(edots, parse(Float64, r[3]))
            end
        end
        close(f)
        if l10 == true
            p3s = log10.(p3s)
            ep3s = log10.(ep3s)
            edots = log10.(edots)
        end
        return p3s, ep3s, edots
    end


    function read_data2(filename; edot_min=nothing, edot_max=nothing)
        p3s = []
        ep3s = []
        periods = []
        pdots = []
        edots = []
        f = open(filename)
        for line in readlines(f)
            r = split(line)
            edot = parse(Float64, r[5])
            if (edot_min != nothing) && (edot_max != nothing)
                if (edot >= edot_min) && (edot <= edot_max)
                    push!(p3s, parse(Float64, r[1]))
                    push!(ep3s, parse(Float64, r[2]))
                    push!(periods, parse(Float64, r[3]))
                    push!(pdots, parse(Float64, r[4]))
                    push!(edots, parse(Float64, r[5]))
                end
            elseif edot_max != nothing
                if edot <= edot_max

                    push!(p3s, parse(Float64, r[1]))
                    push!(ep3s, parse(Float64, r[2]))
                    push!(periods, parse(Float64, r[3]))
                    push!(pdots, parse(Float64, r[4]))
                    push!(edots, parse(Float64, r[5]))
                end
            elseif edot_min != nothing
                if edot >= edot_min
                    push!(p3s, parse(Float64, r[1]))
                    push!(ep3s, parse(Float64, r[2]))
                    push!(periods, parse(Float64, r[3]))
                    push!(pdots, parse(Float64, r[4]))
                    push!(edots, parse(Float64, r[5]))
                end
            else
                push!(p3s, parse(Float64, r[1]))
                push!(ep3s, parse(Float64, r[2]))
                push!(periods, parse(Float64, r[3]))
                push!(pdots, parse(Float64, r[4]))
                push!(edots, parse(Float64, r[5]))
            end
        end
        close(f)
        return p3s, ep3s, periods, pdots, edots
    end



    function calculate_pvals(data1, data2, edots, sample=30)

        i = sortperm(edots)
        d1 = data1[i]
        d2 = data2[i]
        xd = edots[i]

        rngs = collect(range(1, step=sample, stop=length(xd)))
        if rngs[end] < length(xd)
            rngs[end] = length(xd)
        end
        pvals = Array{Float64}(undef, length(rngs)-1)
        for i in 1:length(rngs)-1
            dd1 = d1[rngs[i]:rngs[i+1]]
            dd2 = d2[rngs[i]:rngs[i+1]]
            p1 = pvalue(EqualVarianceTTest(dd1, dd2))
            #p1 = pvalue(EqualVarianceTTest(10 .^ (dd1), 10 .^ (dd2))) # less defined
            #p2 = pvalue(UnequalVarianceTTest(dd1, dd2))
            #println(p1, " ",  p2)
            pvals[i] = p1
        end
        xds = (rngs .+ sample/2)[1:end-1]
        #=

    xds = (rngs + int(sample/2))[:-1]
    bins = len(pvals) * 2
    pvals_bin = np.empty(bins)
    xpvals_bin = np.empty(bins)
    for i in range(0, bins, 2):
        pvals_bin[i] = pvals[i//2]
        pvals_bin[i+1] = pvals[i//2]

    xpvals_bin[0] = xd[0]
    xpvals_bin[-1] = xd[-1]
    j = 1
    for i in range(1, bins - 1, 2):
        xpvals_bin[i] = xd[xds[j]]
        xpvals_bin[i+1] = xd[xds[j]]
        j += 1
        =#
        return [pvals, xds]
    end


    function best_p3edot(;y=(-1.5, 1.5), a=(-1.0, 1.0), size_=100, repeat=100, filename="../data/p3_edot.txt" ) #, filename="../data/p3_edot-synthetic.txt")
        p3s, ep3s, edots = read_data(filename) #; edot_min=4e31) #; edot_min=5e32) #  edot_max!!! # HERE HERE HERE
        ys = range(y[1], y[2], length=size_)
        as = range(a[1], a[2], length=size_)
        x_th = log10(1e32)
        sigma = std(p3s)
        two10 = log10(2)

        obsnum = length(p3s)

        p3s_obs = Array{Float64}(undef, obsnum)
        p3s_model = Array{Float64}(undef, obsnum)
        edots_model = Array{Float64}(undef, obsnum)
        t_tests = Array{Float64}(undef, size_, size_, repeat)
        t_test = Array{Float64}(undef, size_, size_)
        #t_tests = zeros(Float64, (size_, size_, repeat))
        #t_test = zeros(Float64, (size_, size_))

        for i in 1:size_
            for j in 1:size_
                for k in 1:repeat
                    y = ys[i]
                    a = as[j]
                    b = y - a * x_th
                    p3fun(x) = a * x + b
                    xline = range(28, 37, length=100)
                    yline = p3fun.(xline)
                    p3s_notobs = []
                    edots_notobs = []
                    skip = false
                    for ii in 1:length(edots) # generate p3
                        p3 = rand(Normal(p3fun(edots[ii]), sigma))
                        if p3 >= two10
                            p3s_model[ii] = p3
                            edots_model[ii] = edots[ii]
                        else
                            push!(p3s_notobs, p3)
                            push!(edots_notobs, edots[ii])
                            p3obs = nothing
                            for n in 1:1000
                                p3nolog = 10 ^ p3
                                p3obss = abs(p3nolog / (1 - n * p3nolog))
                                if p3obss > 2
                                    p3obs = p3obss
                                    break
                                end
                            end
                            if p3obs != nothing
                                p3s_model[ii] = log10(p3obs) # it is ok
                                edots_model[ii] = edots[ii]
                            else
                                skip = true
                            end
                        end
                    end # generate p3

                    # randomize observations
                    for ii in 1:obsnum
                        error = 10 ^ ep3s[ii]
                        p3 = 10 ^ p3s[ii]
                        new_p3 = rand(Normal(p3, error))
                        while (new_p3 < 2)
                            new_p3 = rand(Normal(p3, error))
                        end
                        #println(p3s_obs[ii])
                        #println(p3, " ", new_p3, " ", error)
                        p3s_obs[ii] = log10(new_p3)
                        #p3s_obs[ii] = p3s[ii] # not random
                    end


                    if skip != true
                        #pvals, xpvals = calculate_pvals(p3s, p3s_model, edots)
                        pvals, xpvals = calculate_pvals(p3s_obs, p3s_model, edots)
                        #=
                        figure(figsize=(13.149606, 13.946563744))
                        subplot(2,1,1)
                        minorticks_on()
                        scatter(edots, p3s)
                        scatter(edots_model, p3s_model)
                        scatter(edots_notobs, p3s_notobs)
                        plot(xline, yline)
                        xlims = xlim()
                        subplot(2,1,2)
                        minorticks_on()
                        axhline(y=0.05, c="C2", ls="--")
                        #xlim(xlims[0], xlims[1])
                        scatter(xpvals, pvals, c="green")
                        show()
                        q = readline(stdin; keep=false)
                        if q == "q"
                            return
                        end
                        close()
                        =#
                        t_tests[i, j, k] = length([pv for pv in pvals if pv < 0.05])
                    else
                        t_tests[i, j, k] = 11
                    end
                end # k
                t_test[i, j] = trunc(mean(t_tests[i, j, :]))
                #t_test[i, j] = minimum(t_tests[i, j, :])
                #println(i, j)
                #println(t_test[i, j])
            end # j
        end # i

        print(size(t_test))

        figure(figsize=(13.149606, 13.946563744))
        minorticks_on()
        #scatter([-0.3], [-0.5], c="red", s=10)
        #imshow(transpose(t_test), origin="lower", extent=[ys[1], ys[end], as[1], as[end]])
        imshow(t_test, origin="lower", extent=[as[1], as[end], ys[1], ys[end]])
        xlabel("a")
        ylabel("y")
        ticks = range(0, stop=10, length=11)
        colorbar(ticks=ticks)
        savefig("best_p3edot.pdf")
        show()
        q = readline(stdin; keep=false)
        PyPlot.close()

    end

    """
    same as bestp3edot but with limit on minimum p3 based on number of sparks i.e. p3 > 1/n_s (P4 >P)
    """
    function best_p3edot_limited(;y=(-1.5, 1.5), a=(-1.0, 1.0), size_=10, repeat=10, filename="../data/p3_edot.txt" ) #, filename="../data/p3_edot-synthetic.txt")
        p3s, ep3s, edots = read_data(filename) #; edot_max=4e31) #; edot_min=4e31) #; edot_min=5e32) #  edot_max!!! # HERE HERE HERE
        ys = range(y[1], y[2], length=size_)
        as = range(a[1], a[2], length=size_)
        x_th = log10(1e32)
        sigma = std(p3s)
        two10 = log10(2)

        obsnum = length(p3s)

        p3s_obs = Array{Float64}(undef, obsnum)
        p3s_model = Array{Float64}(undef, obsnum)
        edots_model = Array{Float64}(undef, obsnum)
        t_tests = Array{Float64}(undef, size_, size_, repeat)
        t_test = Array{Float64}(undef, size_, size_)
        #t_tests = zeros(Float64, (size_, size_, repeat))
        #t_test = zeros(Float64, (size_, size_))

        nsp = 20

        for i in 1:size_
            for j in 1:size_
                for k in 1:repeat
                    y = ys[i]
                    a = as[j]
                    b = y - a * x_th
                    p3fun(x) = a * x + b
                    xline = range(28, 37, length=100)
                    yline = p3fun.(xline)
                    p3s_notobs = []
                    edots_notobs = []
                    skip = false
                    for ii in 1:length(edots) # generate p3
                        p3 = rand(Normal(p3fun(edots[ii]), sigma))
                        if p3 >= two10
                            p3s_model[ii] = p3
                            edots_model[ii] = edots[ii]
                        elseif (p3 < 1 / nsp)
                            p3 = rand(Normal(1 / nsp, sigma))
                            while (p3 < 1 / nsp)
                                p3 = rand(Normal(1 / nsp, sigma))
                            end
                            if p3 >= two10
                                p3s_model[ii] = p3
                                edots_model[ii] = edots[ii]
                            else
                                push!(p3s_notobs, p3)
                                push!(edots_notobs, edots[ii])
                                p3obs = nothing
                                for n in 1:1000
                                    p3nolog = 10 ^ p3
                                    p3obss = abs(p3nolog / (1 - n * p3nolog))
                                    if p3obss > 2
                                        p3obs = p3obss
                                        break
                                    end
                                end
                                if p3obs != nothing
                                    p3s_model[ii] = log10(p3obs) # it is ok
                                    edots_model[ii] = edots[ii]
                                else
                                    skip = true
                                end
                            end
                        else
                            push!(p3s_notobs, p3)
                            push!(edots_notobs, edots[ii])
                            p3obs = nothing
                            for n in 1:1000
                                p3nolog = 10 ^ p3
                                p3obss = abs(p3nolog / (1 - n * p3nolog))
                                if p3obss > 2
                                    p3obs = p3obss
                                    break
                                end
                            end
                            if p3obs != nothing
                                p3s_model[ii] = log10(p3obs) # it is ok
                                edots_model[ii] = edots[ii]
                            else
                                skip = true
                            end
                        end
                    end # generate p3

                    # randomize observations
                    for ii in 1:obsnum
                        error = 10 ^ ep3s[ii]
                        p3 = 10 ^ p3s[ii]
                        new_p3 = rand(Normal(p3, error))
                        while (new_p3 < 2)
                            new_p3 = rand(Normal(p3, error))
                        end
                        #println(p3s_obs[ii])
                        #println(p3, " ", new_p3, " ", error)
                        p3s_obs[ii] = log10(new_p3)
                        #p3s_obs[ii] = p3s[ii] # not random
                    end


                    if skip != true
                        pvals, xpvals = calculate_pvals(p3s_obs, p3s_model, edots)
                        t_tests[i, j, k] = length([pv for pv in pvals if pv < 0.05])
                        #=
                        if t_tests[i, j, k] < 2
                        figure(figsize=(13.149606, 13.946563744))
                        subplot(2,1,1)
                        minorticks_on()
                        scatter(edots, p3s)
                        scatter(edots_model, p3s_model)
                        scatter(edots_notobs, p3s_notobs)
                        plot(xline, yline)
                        xlims = xlim()
                        subplot(2,1,2)
                        minorticks_on()
                        axhline(y=0.05, c="C2", ls="--")
                        #xlim(xlims[0], xlims[1])
                        scatter(xpvals, pvals, c="green")
                        show()
                        q = readline(stdin; keep=false)
                        if q == "q"
                            return
                        end
                        close()
                        end
                        =#

                    else
                        t_tests[i, j, k] = 11
                    end
                end # k
                t_test[i, j] = trunc(mean(t_tests[i, j, :]))
                #t_test[i, j] = minimum(t_tests[i, j, :])
            end # j
        end # i

        print(size(t_test))

        figure(figsize=(13.149606, 13.946563744))
        minorticks_on()
        #scatter([-0.3], [-0.5], c="red", s=10)
        #imshow(transpose(t_test), origin="lower", extent=[ys[1], ys[end], as[1], as[end]])
        imshow(t_test, origin="lower", extent=[as[1], as[end], ys[1], ys[end]])
        xlabel("a")
        ylabel("y")
        ticks = range(0, stop=10, length=11)
        colorbar(ticks=ticks)
        savefig("best_p3edot.pdf")
        show()
        q = readline(stdin; keep=false)
        PyPlot.close()

    end

    """
    Using parabolic fit to the observed data and assuming aliassing to reproduce the fit
    """
    function best_p3edot_parabolic(;y=(-1.5, 1.5), a=(-1.0, 1.0), size_=100, repeat=100, filename="../data/p3_edot.txt" ) #, filename="../data/p3_edot-synthetic.txt")
        p3s, ep3s, edots = read_data(filename) #; edot_max=4e31) #; edot_min=4e31) #; edot_min=5e32) #  edot_max!!! # HERE HERE HERE
        ys = range(y[1], y[2], length=size_)
        as = range(a[1], a[2], length=size_)
        x_th = log10(1e31)
        sigma = std(p3s)
        two10 = log10(2)

        obsnum = length(p3s)

        p3s_obs = Array{Float64}(undef, obsnum)
        p3s_model = Array{Float64}(undef, obsnum)
        edots_model = Array{Float64}(undef, obsnum)
        t_tests = Array{Float64}(undef, size_, size_, repeat)
        t_test = Array{Float64}(undef, size_, size_)
        #t_tests = zeros(Float64, (size_, size_, repeat))
        #t_test = zeros(Float64, (size_, size_))

        emin = minimum(edots)
        emax = maximum(edots)

        # poarabolic fit to whole data
        @. para(x, p) = p[1] * x^2 + p[2] * x + p[3]
        fit = curve_fit(para, edots, p3s, [1., 1., 1.])
        p1 = coef(fit)
        err = stderror(fit)
        xp = range(emin, emax, length=100)
        yp = para(xp, p1)

        # unalised values assuming all observations are aliased
        @. p3(p3obs, n) = p3obs / (n * p3obs + 1)
        yintr = log10.(p3(10 .^ yp, 1)) # nice

        # polynominal fit to intrinsic p3
        @. poly(x, p) = p[1] * x^3 + p[2] * x^2 + p[3] * x + p[4]
        fit = curve_fit(poly, xp, yintr, [1., 1., 1., 1.])
        p2 = coef(fit)
        yp2 = poly(xp, p2)

        p3intr(x, p) = poly(x, p) # intrinsic (unaliased) p3

        # nsp dependence
        @. nsp_fun(x, p) = p[1] * x^2 + p[2] * x + p[3]
        xpoints = [29, 31, 33, 35, 36.4]
        ypoints = [15, 10, 5, 3, 2]
        #ypoints = [2, 2, 2, 2, 2]
        fit = curve_fit(nsp_fun, xpoints, ypoints, [1., 1., 1.])
        pnsp = coef(fit)
        nsps = nsp_fun(xp, pnsp)
        p3mins = log10.(1. ./ nsps)

        #=
        scatter(edots, p3s)
        plot(xp, yp)
        plot(xp, yintr)
        plot(xp, yp2)
        plot(xp, p3mins, ls="--", c="black")
        show()
        q = readline(stdin; keep=false)
        PyPlot.close()
        =#


        for i in 1:size_
            for j in 1:size_

                y0 = ys[i]
                a = as[j]
                b = y0 - a * x_th

                # dependence in question
                p3fun(x) = a * x + b
                xline = range(emin, emax, length=100)
                yline = p3fun.(xline)

                #=
                        figure(figsize=(13.149606, 13.946563744))
                        subplot(2,1,1)
                        minorticks_on()
                        scatter(edots, p3s)
                        #scatter(edots_model, p3s_model)
                        #scatter(edots_notobs, p3s_notobs)
                        plot(edots, yt)
                        xlims = xlim()
                        subplot(2,1,2)
                        minorticks_on()
                        axhline(y=0.05, c="C2", ls="--")
                        #xlim(xlims[0], xlims[1])
                        #scatter(xpvals, pvals, c="green")
                        show()
                        q = readline(stdin; keep=false)
                        if q == "q"
                            return
                        end
                        close()
                =#

                for k in 1:repeat
                    p3s_notobs = []
                    edots_notobs = []
                    skip = false
                    # generate p3
                    for ii in 1:length(edots)
                        p3f = p3fun(edots[ii])
                        p3i = p3intr(edots[ii], p2) # p2 fitting parameters
                        if  p3f >= p3i
                            p3 = p3f
                        else
                            p3 = p3i
                        end
                        p3m = p3
                        p3 = rand(Normal(p3m, sigma))
                        while 10 ^ p3 < 1. / nsp_fun(edots[ii], pnsp)
                            p3 = rand(Normal(p3m, sigma))
                            #println("p3 ", 10^p3)
                            #println("p3 limit ", 1. / nsp_fun(edots[ii], pnsp))
                        end
                        if 10 ^ p3 >= 2
                            p3s_model[ii] = p3
                            edots_model[ii] = edots[ii]
                        else
                            while 10 ^ p3 == 1
                                p3 = rand(Normal(p3m, sigma))
                            end
                            push!(p3s_notobs, p3)
                            push!(edots_notobs, edots[ii])
                            alias = false
                            n = 1
                            p3obs = p3
                            while (10 ^ p3obs < 2)
                                p3obs = log10(abs(10 ^ p3 / (1 - n * 10 ^ p3)))
                                n += 1
                                alias = true
                            end
                            p3s_model[ii] = p3obs # it is ok
                            edots_model[ii] = edots[ii]
                        end
                    end # generate p3

                    # randomize observations
                    for ii in 1:obsnum
                        error = 10 ^ ep3s[ii]
                        p3 = 10 ^ p3s[ii]
                        new_p3 = rand(Normal(p3, error))
                        while (new_p3 < 2)
                            new_p3 = rand(Normal(p3, error))
                        end
                        #println(p3s_obs[ii])
                        #println(p3, " ", new_p3, " ", error)
                        p3s_obs[ii] = log10(new_p3)
                        #p3s_obs[ii] = p3s[ii] # not random
                    end


                    if skip != true
                        pvals, xpvals = calculate_pvals(p3s_obs, p3s_model, edots)
                        t_tests[i, j, k] = length([pv for pv in pvals if pv < 0.05])
                    else
                        t_tests[i, j, k] = 11
                    end
                end # k
                t_test[i, j] = trunc(mean(t_tests[i, j, :]))
                #t_test[i, j] = minimum(t_tests[i, j, :])
            end # j
        end # i
        #return

        print(size(t_test))
        figure(figsize=(13.149606, 13.946563744))
        minorticks_on()
        #scatter([-0.3], [-0.5], c="red", s=10)
        #imshow(transpose(t_test), origin="lower", extent=[ys[1], ys[end], as[1], as[end]])
        imshow(t_test, origin="lower", extent=[as[1], as[end], ys[1], ys[end]])
        xlabel("a")
        ylabel("y")
        ticks = range(0, stop=10, length=11)
        colorbar(ticks=ticks)
        savefig("best_p3edot_parabolic.pdf")
        show()
        q = readline(stdin; keep=false)
        PyPlot.close()

    end


    """
    Using Geoff's fit
    """
    function best_p3edot_geoff(;w=(0.5, 1.5), t=(-1.5, -0.5), size_=100, repeat=100, filename="../data/p3_ppdotedot.txt")
        # values not in log scale
        p3s, ep3s, periods, pdots, edots = read_data2(filename) #; edot_max=4e31) #; edot_min=4e31) #; edot_min=5e32) #  edot_max!!! # HERE HERE HERE

        obsnum = length(p3s)
        xi_xs = Array{Float64}(undef, obsnum)
        for (i, p) in enumerate(periods)
            xi_xs[i] = p ^ (-1.78) * (pdots[i] / 1e-15)
        end

        ws = range(w[1], w[2], length=size_)
        ts = range(t[1], t[2], length=size_)
        x_th = 7e30 # low - high edot


        p3s_obs = Array{Float64}(undef, obsnum)
        p3s_model = Array{Float64}(undef, obsnum)
        edots_model = Array{Float64}(undef, obsnum)
        xs_model = Array{Float64}(undef, obsnum)
        t_tests = Array{Float64}(undef, size_, size_, repeat)
        t_test = Array{Float64}(undef, size_, size_)

        emin = minimum(edots)
        emax = maximum(edots)

        """ aliased values """
        function p3aliased(p3, n)
            if p3 > 1
                return p3 / (n * p3 - 1)
            elseif p3 == 1
                return 100.
            else
                return p3 / (1 - n * p3)
            end
        end

        # dependence in question
        function p3fun(x, w, t)
            return abs(1 / (1 - w / (x - t)))
        end

        n = 1
        sigma = std(log10.(p3s))
        for i in 1:size_
            for j in 1:size_
                w_ = ws[i]
                t_ = ts[j]
                #w_ = 1.0
                #t_ = -1.05
                xline = 10 .^ range(-3, 4, length=100)
                yline = p3fun.(xline, w_, t_)
                p3s_notobs = []
                xs_notobs = []
                for k in 1:repeat
                    p3s_notobs = []
                    xs_notobs = []
                    skip = false
                    # generate p3
                    for ii in 1:length(xi_xs)
                        p3m = p3fun(xi_xs[ii], w_, t_)
                        p3 = 10 ^ rand(Normal(log10(p3m), sigma))
                        if p3 > 2
                            p3s_model[ii] = p3
                            xs_model[ii] = xi_xs[ii]
                        else
                            p3obs = p3aliased(p3, n)
                            kk = 1
                            while p3obs < 2
                                p3 = 10 ^ rand(Normal(log10(p3m), sigma))
                                p3obs = p3aliased(p3, n)
                                kk += 1
                                if kk == 1000
                                    p3obs = 100
                                    break
                                end

                            end
                            p3s_model[ii] = p3obs
                            xs_model[ii] = xi_xs[ii]
                            push!(p3s_notobs, p3)
                            push!(xs_notobs, xi_xs[ii])
                        end
                    end # generate p3

                    # randomize observations
                    for ii in 1:obsnum
                        error = 10 ^ ep3s[ii]
                        p3 = 10 ^ p3s[ii]
                        new_p3 = rand(Normal(p3s[ii], ep3s[ii]))
                        while (new_p3 < 2)
                            new_p3 = rand(Normal(p3s[ii], ep3s[ii]))
                        end
                        p3s_obs[ii] = new_p3
                    end

                    if skip != true
                        pvals, xpvals = calculate_pvals(p3s_obs, p3s_model, edots)
                        t_tests[i, j, k] = length([pv for pv in pvals if pv < 0.05])
                    else
                        t_tests[i, j, k] = 11
                    end
                end # k

                println("w = $w_ t = $t_ ", trunc(mean(t_tests[i, j, :])))
                #=
                    figure(figsize=(13.149606, 13.946563744))
                    subplot(2,1,1)
                    minorticks_on()
                    #scatter(xi_xs, p3s)
                    scatter(xi_xs, p3s_obs, c="tab:blue")
                    scatter(xi_xs, p3s_model, c="tab:orange")
                    scatter(xs_notobs, p3s_notobs, c="tab:grey")
                    #scatter(edots_model, p3s_model)
                    #scatter(edots_notobs, p3s_notobs)
                    plot(xline, yline, c="tab:red")
                    xlims = xlim()
                    loglog()
                    subplot(2,1,2)
                    minorticks_on()
                    axhline(y=0.05, c="C2", ls="--")
                    #xlim(xlims[0], xlims[1])
                    #scatter(xpvals, pvals, c="green")
                    #show()
                    savefig("test.pdf")
                    PyPlot.close()
                    q = readline(stdin; keep=false)
                    if q == "q"
                        return
                    end
                =#
                t_test[i, j] = trunc(mean(t_tests[i, j, :]))
                #t_test[i, j] = minimum(t_tests[i, j, :])
            end # j
        end # i
        #return

        print(size(t_test))
        figure(figsize=(13.149606, 13.946563744))
        minorticks_on()
        #scatter([-0.3], [-0.5], c="red", s=10)
        #imshow(transpose(t_test), origin="lower", extent=[ys[1], ys[end], as[1], as[end]])
        imshow(t_test, origin="lower", extent=[ts[1], ts[end], ws[1], ws[end]])
        xlabel("t?")
        ylabel("w?")

        ticks = range(0, stop=10, length=11)
        colorbar(ticks=ticks)
        savefig("best_p3edot_geoff.pdf")
        show()
        q = readline(stdin; keep=false)
        PyPlot.close()

    end




    function test_python()
        # seems to work...
        #s1 = [1,43,3,43,3,2,4,53, 3,2,5,3,9,9,3,1,32,5,33,2,1,5,2,5]
        #s2 = [1,43,3,11,3,2,4,3, 3,2,5,23,9,19,3,1,3,5,3,2,1,5,2,5]
        s1 = rand(30)
        s2 = rand(30)

        pv = pvalue(EqualVarianceTTest(s1, s2))
        println(pv)
        println(stats.ttest_ind(s1, s2))
    end


    function main()
        #best_p3edot()
        #best_p3edot_limited()
        #best_p3edot_parabolic()
        best_p3edot_geoff()
        println("Bye")
    end


end

#Simulation.test_python()
Simulation.main()
