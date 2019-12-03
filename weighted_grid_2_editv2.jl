using Plots
using ParticleFilters
using Distributions
using StaticArrays
using LinearAlgebra
using Random
using StatsBase
#using Reel
using SparseArrays
using GridInterpolations
using DataStructures
using DataFrames
using CSV
using Distributed

include("atan2.jl")
include("obs_rel.jl")
include("polargrid_rel_qual.jl")
rng = MersenneTwister(2);

# returns vector rather than Tuple, for particle filter
function f2(x, u, rng)
    temp = [i for i in f(x, u, rng)]
    return temp
end

##########

#=

N = 500
model = ParticleFilterModel{Vector{Float64}}(f2, g)
pfilter = SIRParticleFilter(model, N);

x = [rand(rng, 25:150), rand(rng,0:359), rand(rng,0:11)*30, 1, 1];   
b = ParticleCollection([x[1:4] for i in 1:N]);

# In loop:

xp = f2(x, actions_[u], rng)
y = h(xp, rng)
b = update(pfilter, b, actions_[u], y)

=#

##########

thestates = [(10.0*r,θ*30.0, crs*30.0, spd, qual) for r in 1:30, 
        θ in 0:11, crs in 1:12, spd in 1:2, qual in [0, 1, 10]];

states_r = [10.0*r for r in 1:30]
states_θ = [30.0*θ for θ in 0:12]
states_crs = [30*c for c in 1:12]
states_qual = [q for q in [0, 1, 10]]
#states_hdg = [30*h for h in 1:12]
states_spd = [1,2];

grid = RectangleGrid(states_r, states_θ, states_crs, states_spd, states_qual)  	# rectangular grid

wrapmap = Dict()

### Essentially used to map down the 360 degree parameters.

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

function polar_grid(vec, grid=grid)
    polants = interpolants(grid, vec)
    for (i, p) in enumerate(polants[1])
        if haskey(wrapmap, p)
            polants[1][i] = wrapmap[p]
        end    
    end
    return polants
end

function weighted_grid_2_edit(b::ParticleCollection)
    
    #=
    beta = zeros(length(grid));
    for row in particles(b)
        for (i, x) in enumerate(polar_grid(row)[1])
            beta[x] += polar_grid(row)[2][i]
        end
    end
    =#
    
    # Creating the degree-radian conversion factor.
    
    deg_rad_conv = pi/180;
    
    # Creating the initial data structure to filter out the "b_particle = particles(b)" 
    
    b_particle = particles(b);
    row_num = length(b_particle);
    col_num = length(b_particle[1]);
    b_particle_filtered = zeros(row_num, col_num);

    # Filtering out the b_particle "Array{Array}" into an "Array{Int64}" matrix.
    
    for row_index in 1:row_num
        
        for col_index in 1:col_num
            
            b_particle_filtered[row_index, col_index] = b_particle[row_index][col_index];
            
        end
        
    end
    
    # Finding the mean of the 'r', 'θ', & crs values.
    
    r_mean = 0;

    cos_count_θ = 0;
    sin_count_θ = 0;
    
    cos_count_crs = 0;
    sin_count_crs = 0;
    
    for mean_index in 1:row_num
        
        r = b_particle_filtered[mean_index, 1];
        θ = b_particle_filtered[mean_index, 2];
        crs = b_particle_filtered[mean_index, 3];
        
        r_mean += (r / row_num);
        
        cos_count_θ += (1/row_num) * cos(θ * deg_rad_conv);
        sin_count_θ += (1/row_num) * sin(θ * deg_rad_conv);
        
        cos_count_crs += (1/row_num) * cos(crs * deg_rad_conv);
        sin_count_crs += (1/row_num) * sin(crs * deg_rad_conv);
        
    end
    
    r_mean = round(r_mean);
    θ_mean = round(atan2(cos_count_θ,sin_count_θ) / deg_rad_conv);
    crs_mean = round(atan2(cos_count_crs,sin_count_crs) / deg_rad_conv);
    
    # If "θ_mean" or "crs_mean" are less than 0 degrees, the following code adds 360 degrees.
    
    if(θ_mean < 0.0)
        
        θ_mean += 360
        
    end
    
    if(crs_mean < 0.0)
        
        crs_mean += 360
        
    end

    # Source for mean radius operation: http://ballistipedia.com/index.php?title=Mean_Radius
    # Source for mean angle operation: https://rosettacode.org/wiki/Averages/Mean_angle
    # Some other information on the normal distributions of r & θ: http://ballistipedia.com/index.php?title=Closed_Form_Precision#Mean_Radius_.28MR.29
    
    # Finding the standard deviation of the crs values.
    
    sum_squares_r = 0;
    sum_squares_θ = 0;
    sum_squares_crs = 0;
    
    for std_index in 1:row_num
        
        r_i = b_particle_filtered[std_index,1];
        θ_i = b_particle_filtered[std_index,2];
        crs_i = b_particle_filtered[std_index,3];
        
        sum_squares_r += (r_i - r_mean)^2;
        sum_squares_θ += (θ_i - θ_mean)^2;
        sum_squares_crs += (crs_i - crs_mean)^2;
        
    end
    
    std_r = sqrt(sum_squares_r / (row_num - 1) );
    std_θ = sqrt(sum_squares_θ / (row_num - 1) );
    std_crs = sqrt(sum_squares_crs / (row_num - 1) );
    
    ##########
    
    #=

    # Converting all of the 'r' & 'θ' values into 'x' and 'y' coordinates.
    # The 'r' & 'θ' values are subsequently replaced with the 'x' and 'y' coordinates.
    
    for row_index_convert in 1:row_num
        
        r = b_particle_filtered[row_index_convert, 1];
        θ = b_particle_filtered[row_index_convert, 2];
        
        x = r * cos(θ * deg_rad_conv);
        y = r * sin(θ * deg_rad_conv);
        
        b_particle_filtered[row_index_convert, 1] = x;
        b_particle_filtered[row_index_convert, 2] = y;
        
    end
    
    x_center = mean(b_particle_filtered[:,1]);
    y_center = mean(b_particle_filtered[:,2]);
    
    x_var = std(b_particle_filtered[:,1]);
    y_var = std(b_particle_filtered[:,2]);
    
    r_center = sqrt(x_center^2 + y_center^2);
    
    r_var = sqrt(x_var^2 + y_var^2);
    # This isn't right because you first need to find the covariance matrix, ...
    # ... then use the singular value decomposition (SVD), and more...
    # Let's try to avoid this.
    
    # Covariance: https://en.wikipedia.org/wiki/Covariance
    # SVD: https://stackoverflow.com/questions/16585980/variance-matrix-from-polar-to-cartesian-coordinates
    
    =# 
    
    ##########
    
    kk = 25
    qual_ = max(11-exp((std_r+std_θ+std_crs)/kk),0)
    
    b_individual = [r_mean, θ_mean, crs_mean, b_particle[1][4], qual_]; # std_r, std_θ,
    
    b_individual_print = [r_mean, std_r, θ_mean, std_θ, crs_mean, std_crs, b_particle[1][4]];
    println(b_individual)
    
    beta = zeros(length(grid));
    for (i, x) in enumerate(polar_grid(b_individual)[1])
        beta[x] += polar_grid(b_individual)[2][i]
    end
            
    return beta
    
    ##########
    
    #=
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
    =#
    
    ##########
    
end

particles(b)[1]

N = 500;
model = ParticleFilterModel{Vector{Float64}}(f2, g);
pfilter = SIRParticleFilter(model, N);

x = [rand(rng, 25:150), rand(rng,0:359), rand(rng,0:11)*30, 1, 1];   
b = ParticleCollection([x[1:4] for i in 1:N]);

b = ParticleCollection([[rand(rng, 19:35), rand(rng,1:2), rand(rng,7:11)*30, 1, 1] for i in 1:N])

test = weighted_grid_2_edit(b);

[transpose(particles(b)[1])];
particles(b)

count_test = 0;
test_val = 0.00;

for counter_test in (1:length(test))
    
    if (test[counter_test]) > test_val
        
        print(grid[counter_test]);
        print(", ");
        println(test[counter_test]);
        
        count_test += 1;
        
    end
    
end

println(length(test))
println(length(grid))


