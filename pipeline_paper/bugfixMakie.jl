import Makie.get_minor_tickvalues

function get_minor_tickvalues(i::IntervalsBetween, scale::Union{typeof(log), typeof(log2), typeof(log10)}, tickvalues, vmin, vmax)

    vals = Float64[]
    length(tickvalues) < 2 && return vals
    n = i.n

    invscale = Makie.inverse_transform(scale)

    if i.mirror
        firstinterval_scaled = scale(tickvalues[2]) - scale(tickvalues[1])
        stepsize = firstinterval_scaled / n
        prevtick = invscale(scale(tickvalues[1]) - firstinterval_scaled)
        stepsize = (tickvalues[1] - prevtick) / n
        v = tickvalues[1] - stepsize
        prepend!(vals, v:-stepsize:vmin)
    end

    for (lo, hi) in zip(@view(tickvalues[1:end-1]), @view(tickvalues[2:end]))
        interval = hi - lo
        stepsize = interval / n
        v = lo
        for i in 1:n-1
            v += stepsize
            push!(vals, v)
        end
    end

    if i.mirror
        lastinterval_scaled = scale(tickvalues[end]) - scale(tickvalues[end-1])
        nexttick = invscale(scale(Float64(tickvalues[end])) + Float64(lastinterval_scaled))
        stepsize = (nexttick - tickvalues[end]) / n
        v = tickvalues[end] + stepsize
        append!(vals, v:stepsize:vmax)
    end

    vals
end