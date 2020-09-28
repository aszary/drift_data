module Functions

    using DataFrames, GLM


    function p3(p3_obs, n) # P3_obs in Periods
        p3 = 1 / (n + 1 / (p3_obs))
        return p3
    end

    function p3_edwards(p3_obs, p, n) # P3_obs in Periods
        p3_sec = p / (n + p / (p3_obs*p))
        p3 = p3_sec / p
        return p3
    end


    function p3_gupta(p3_obs, p, n) # P3_obs in Periods
        fobs = 1 / (p3_obs * p)
        fN = 0.5 / p
        k = trunc(Int64, (n+1) / 2)
        l = mod(n, 2)
        f = 2 * k * fN + (-1)^l * fobs
        # in seconds
        p3_sec = 1 / f
        p3 = p3_sec / p
        return p3
    end

    """ determines P3 based on P3~Edot anticorrelation """
    function p3_edot(p3_obs, edot)
        p3_pre = p3_prediction(edot)
        n_range = -20:1:20#range(-10, 10, step=1)
        min_ = 1e50
        n_best = n_range[1]
        p3_best = 1e50
        for n in n_range
            p3_ = abs(p3(p3_obs, n))
            dist = abs(log10(p3_pre) - log10(p3_))
            if dist < min_
                min_ = dist
                n_best = n
                p3_best = p3_
            end
            #println("$n $p3_")
        end
        #println("$n_best $p3_best")
        return p3_best, n_best
    end


    function p3_prediction(edot)
        return 10 ^ 19.1789 * edot ^ -0.591211 # update those numbers
        #return 10 ^ -6.0 * edot ^ 0.18 # update those numbers
    end


    """ determines P3 based on P3~Edot anticorrelation """
    function p3_edot2(p3_obs, edot)
        p3_pre = p3_prediction2(edot)
        n_range = -70:1:70#range(-10, 10, step=1)
        min_ = 1e50
        n_best = n_range[1]
        p3_best = 1e50
        for n in n_range
            p3_ = abs(p3(p3_obs, n))
            dist = abs(log10(p3_pre) - log10(p3_))
            if dist < min_
                min_ = dist
                n_best = n
                p3_best = p3_
            end
        end
        return p3_best, n_best
    end


    function p3_prediction2(edot)
        #return 10 ^ -6.0 * edot ^ 0.18
        #return 10 ^ -13.0 * edot ^ 0.4
        return 10 ^ -16.0 * edot ^ 0.5
        #return 10 ^ -22.0 * edot ^ 0.7
        #return 10 ^ -23.0 * edot ^ 0.7
    end


    function fit_line_log10(x_, y_; show_=false)
        x = log10.(x_)
        y = log10.(y_)
        data = DataFrame(X=x, Y=y)
        ols = lm(@formula(Y ~ X), data)
        co = coef(ols)
        if show_ == true
            println(ols)
        end
        (x_min, x_max) = extrema(x)
        edots = 10 .^ (range(x_min; stop=x_max, length=100))
        f(x) = 10^co[1] * x^co[2]
        p3s = f.(edots)
        return co, [edots, p3s]
    end


end  # module Functions
