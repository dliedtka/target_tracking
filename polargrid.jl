thestates = [(10.0*r,θ*30.0, crs*30.0, hdg*30.0, spd) for r in 1:30, 
        θ in 0:11, crs in 1:12, hdg in 1:12, spd in 1:2];

states_r = [10.0*r for r in 1:30]
states_θ = [30.0*θ for θ in 0:12]
states_crs = [30*c for c in 1:12]
states_hdg = [30*h for h in 1:12]
states_spd = [1,2];

somestates_x, somestates_y = [0., 1., 2.], [0, 90, 180, 270, 360]
grid = RectangleGrid(states_r, states_θ, states_crs, states_hdg, states_spd)  	# rectangular grid

gridData = [8., 1., 6., 3., 5., 7., 4., 9., 2.]   	# vector of value data at each cut
x = [1.5, 325]
#@show interpolants(grid, x)

wrapmap = Dict()

for r in states_r
    for crs in states_crs
        for hdg in states_hdg
            for spd in states_spd
                key = interpolants(grid, [r, 360, crs, hdg, spd])[1][1]
                val = interpolants(grid, [r, 0, crs, hdg, spd])[1][1]
                wrapmap[key] = val
            end
        end
    end
end

polants =  [12.5, 358, 270, 60, 2]

function polar_grid(vec, grid=grid)
    polants = interpolants(grid, vec)
    for (i, p) in enumerate(polants[1])
        if haskey(wrapmap, p)
            polants[1][i] = wrapmap[p]
        end    
    end
    return polants
end

function weighted_grid(b::ParticleCollection)
    beta = Dict()
    for row in particles(b)
        for (i, x) in enumerate(polar_grid(row)[1])
            if !haskey(beta, x)
                beta[x] = 0
            end
            beta[x] += polar_grid(row)[2][i]
        end
    end
    return beta
end

function weighted_grid_2(b::ParticleCollection)
    beta = zeros(length(grid));
    for row in particles(b)
        for (i, x) in enumerate(polar_grid(row)[1])
            beta[x] += polar_grid(row)[2][i]
        end
    end
    return beta
end

# returns random argument of tied maxima
# else (no ties) mimics argmax
function argmax2(thing, rng)
    options = []
    max = maximum(thing)
    for (i, x) in enumerate(thing)
        if x == max
            push!(options, i)
        end
    end
    return rand(rng, options)
end


