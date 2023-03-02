module Functions

    using Statistics
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
            #println("\t", n, " ", p3_, " ", dist, " p3pre=", p3_pre)
            if dist < min_
                min_ = dist
                n_best = n
                p3_best = p3(p3_obs, n) # not abs
            end
            #println("$n $p3_")
        end
        #println("$n_best $p3_best")
        return p3_best, n_best
    end


    function p3_prediction(edot)
        return 10 ^ 19.0414 * edot ^ -0.58647 # update those numbers
        #return 10 ^ 19.1789 * edot ^ -0.591211 # update those numbers
        #return 10 ^ 12.6468 * edot ^ -0.385241 # the whole sample
    end


    """ determines P3 based on P3~Edot correlation """
    function p3_edot(p3_obs, edot, a, b)
        p3_pre = p3_prediction(edot, a, b)
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
        end
        return p3_best, n_best
    end

    function p3_prediction(edot, a, b)
        return 10 ^ a * edot ^ b
    end


    function p3_prediction_rahul(edot)
        return 10 ^ 16.1955 * edot ^ -0.496491
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
        x = log10.(abs.(x_))
        y = log10.(abs.(y_))
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


    function squared_error(ydata, yfit)
        return sum((yfit .- ydata) .* (yfit .- ydata))
    end


    """ https://pythonprogramming.net/how-to-program-r-squared-machine-learning-tutorial/ """
    function rsquared(ydata, yfit)
        mean_line = [mean(ydata) for y in ydata]
        squared_error_fit = squared_error(ydata, yfit)
        squared_error_mean = squared_error(ydata, mean_line)
        return 1 - (squared_error_fit / squared_error_mean)
    end


end  # module Functions
