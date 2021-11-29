module Simulation
    using PyPlot
    using Statistics
    using Random, Distributions
    using HypothesisTests


    function read_data(filename; l10=true, edot_max=nothing)
        p3s = []
        ep3s = []
        edots = []
        f = open(filename)
        for line in readlines(f)
            r = split(line)
            edot = parse(Float64, r[3])
            if edot_max != nothing
                if edot < edot_max
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


    function best_p3edot(;y=(0.05, 30), a=(-1.5, 1.5), size_=100, repeat=100, filename="../data/p3_edot.txt")
        p3s, ep3s, edots = read_data(filename, edot_max=5e32)
        ys = range(log10(y[1]), log10(y[2]), length=size_)
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
                                p3obs_p = abs(p3nolog / (1 - n * p3nolog))
                                p3obs_n = abs(p3nolog / (1 - (-n) * p3nolog))
                                if p3obs_p > 2
                                    p3obs = p3obs_p
                                    break
                                end
                                if p3obs_n > 2
                                    p3obs = p3obs_n
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
                        p3s_obs[ii] = rand(Normal(p3s[ii], abs(ep3s[ii])))
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
                t_test[i, j] = mean(t_tests[i, j, :])
                #t_test[i, j] = minimum(t_tests[i, j, :])
                #println(i, j)
                #println(t_test[i, j])
            end # j
        end # i

    print(size(t_test))

    figure(figsize=(13.149606, 13.946563744))
    minorticks_on()
    imshow(transpose(t_test), origin="lower", extent=[ys[1], ys[end], as[1], as[end]])
    xlabel("y")
    ylabel("a")
    ticks = range(0, stop=2, length=10)
    colorbar(ticks=ticks)
    show()
    q = readline(stdin; keep=false)
    PyPlot.close()



    end


    function main()
        best_p3edot()

        println("Bye")
    end

end

Simulation.main()
