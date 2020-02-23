# rules-based method for choosing next action based on estimate of target location
# action_space = ((-30,1), (-30, 2), (0, 1), (0, 2), (30, 1), (30, 2))
# corresponds to    1         2        3       4         5       6
# negative is left, positive is right


function rules_action(state, risk_parameter=0.5)
    rel_brg = state[2]-state[4]
    range = state[1]
    r1 = 15 + risk_parameter*150
    r2 = 200
    
    #convert negative bearing into positive between 0-360
    if rel_brg < 0 rel_brg += 360 end
        
    #behavior if target range is "distant"
    if range >= r2
        if (rel_brg > 340) || (rel_brg < 20)
            return 4
        elseif (rel_brg > 180)
            return 6
        return 2
        end
    end
    
    #behavior if target range is "long" WIP WIP
    if range >= r1
        
end