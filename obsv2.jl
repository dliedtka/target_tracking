SENSOR_RANGE = 150
prob = 0.3

function obs(state)
    #initialize relative bearing and range from state
    rel_brg = state[2]-state[4]
    range = state[1]
    
    #convert negative bearing into positive between 0-360
    if rel_brg < 0 rel_brg += 360 end
    
    #observation zone 1
    if ((60 < rel_brg < 90) || (270 < rel_brg < 300)) && (range < SENSOR_RANGE)
        return [1-prob, prob, 0, 0]::Array{Float64,1} 
    end
    
    #observation zone 2
    if ((90 <= rel_brg < 120) || (240 < rel_brg <= 270)) && (range < SENSOR_RANGE)
        return [1-prob, 0, prob, 0]::Array{Float64,1} 
    end
    if (120 <= rel_brg <= 240) && (range < SENSOR_RANGE/2)
        return [1-prob, 0, prob, 0]::Array{Float64,1} 
    end
    return [1.0, 0, 0, 0]::Array{Float64, 1}
end        

function g(x, a, xp, o)
    return obs(xp)[o+1]
end

function h(x, rng)
    #weights = [obs0(x), obs1(x), obs2(x), obs3(x)]
    weights = obs(x)
    obsers::Array{Int64,1} = [0, 1, 2, 3]
    return sample(obsers, Weights(weights))::Int64
end