using Distributions

thestates = [(10.0*r,θ*30.0, crs*30.0, spd, qual) for r in 1:30, 
        θ in 0:11, crs in 0:11, spd in 1:2, qual in [0, 1, 10]];

states_r = [10.0*r for r in 1:30]
states_θ = [30.0*θ for θ in 0:12]
states_crs = [30*c for c in 0:12]
states_qual = [q for q in [0, 1, 10]]
#states_hdg = [30*h for h in 1:12]
states_spd = [1,2];

grid = RectangleGrid(states_r, states_θ, states_crs, states_spd, states_qual)  	# rectangular grid

wrapmap = Dict()
wrapmap2 = Dict()

for r in states_r
    for crs in states_crs
        for spd in states_spd
            for qual in states_qual
                key = interpolants(grid, [r, 360, crs, spd, qual])[1][1]
                val = interpolants(grid, [r, 0, crs, spd, qual])[1][1]
                wrapmap[key] = val
            end
        end
    end
end

for r in states_r
    for theta in states_θ
        for spd in states_spd
            for qual in states_qual
                key = interpolants(grid, [r, theta, 360, spd, qual])[1][1]
                val = interpolants(grid, [r, theta, 0, spd, qual])[1][1]
                wrapmap2[key] = val
            end
        end
    end
end

function polar_grid(vec, grid=grid)
    polants = interpolants(grid, vec)
    for (i, p) in enumerate(polants[1])
        if haskey(wrapmap, p)
            polants[1][i] = wrapmap[p]
        end
        if haskey(wrapmap2, p)
            polants[1][i] = wrapmap[p]
        end
    end
    return polants
end

function weighted_grid_2(b::ParticleCollection)
    beta = zeros(length(grid));
    for row in particles(b)
        for (i, x) in enumerate(polar_grid(row)[1])
            beta[x] += polar_grid(row)[2][i]
        end
    end
    vv = var(beta)
    beta = zeros(length(grid));
    for row in particles(b)
        row = [float(r) for r in row]
        push!(row, vv)
        for (i, x) in enumerate(polar_grid(row)[1])
            beta[x] += polar_grid(row)[2][i]
        end
    end
    return beta
end

#function weighted_grid_2(b::ParticleCollection)
#    beta = spzeros(length(grid))
#    mapout = pmap(polar_grid, particles(b))
#    for row in mapout
#        for j in 1:length(row[1])
#            beta[row[1][j]] += row[2][j]
#        end
#    end
#    vv = var(beta)
#    beta = zeros(length(grid));
#    for row in particles(b)
#        row = [float(r) for r in row]
#        push!(row, vv)
#        for (i, x) in enumerate(polar_grid(row)[1])
#            beta[x] += polar_grid(row)[2][i]
#        end
#    end
#    return beta
#end

# returns random argument of tied maxima
# else (no ties) mimics argmax

function argmax2(thing, rng)
    options = []
    max = maximum([i for i in thing if !isnan(i)])
    for (i, x) in enumerate(thing)
        if x == max
            push!(options, i)
        end
    end
    #print(options)
    return rand(rng, options)
end

function max2(thing, rng)
    options = []
    max = maximum([i for i in thing if !isnan(i)])
    for (i, x) in enumerate(thing)
        if x == max
            push!(options, i)
        end
    end   
    return thing[rand(rng, options)]
end


function next_action(thing, epsilon, rng)
    if rand(rng) > epsilon
        return (trunc(Int, argmax2(thing, rng)),0, trunc(Int, argmax2(thing, rng)))
    end
    return (trunc(Int, rand(rng, 1:length(thing))),1, trunc(Int, argmax2(thing, rng)))
end

function softmax_action(thing, k, rng)
    thing = exp.(k*thing)/sum(exp.(k*thing))
    next = trunc(Int, sample(rng, 1:length(thing), StatsBase.weights(thing)))
    best = trunc(Int, argmax2(thing, rng))
    return(next, (next == best)*1, best)
end                                