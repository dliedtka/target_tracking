# rules-based method for choosing next action based on estimate of target location
# action_space = ((-30,1), (-30, 2), (0, 1), (0, 2), (30, 1), (30, 2))
# corresponds to    1         2        3       4         5       6
# negative is right, positive is left

debugging = false

function rules_action(state, risk_parameter=0.5)
    rel_brg = state[2]-state[4]
    range = state[1]
    r1 = 15 + risk_parameter*100
    r2 = 215 - risk_parameter*100
    
    if (range < 10) && (debugging == true)
        println("CLOSE -- Range: ", range, "   Bearing: ", rel_brg)
    elseif range > 250 && debugging == true
        println("FAR -- Range: ", range, "   Bearing: ", rel_brg)
    elseif debugging == true 
        println("FINE -- Range: ", range, "   Bearing: ", rel_brg) 
    end
    
    #convert negative bearing into positive between 0-360
    if rel_brg < 0 rel_brg += 360 end
        
    #behavior if target range is "distant" -- turn toward target and move quickly
    if range >= r2
        if debugging print("distant ") end
        if (rel_brg > 330) || (rel_brg < 30)
            return 4
        elseif (rel_brg > 180)
            return 6
        end
        return 2
    
    #behavior if target range is "long" 
    elseif range >= r1
        if debugging print("long ") end
        #if target is ahead, place target off closest beam
        if 30 >= rel_brg >=0
            return 6
        elseif rel_brg >= 330
            return 2
        #if target is forward of beam, maintain course and draw target down beam
        elseif (90 >= rel_brg > 30) || (330 > rel_brg >= 270)
            return 4
        #if target is abaft of beam, maintain course and match speed with target
        elseif (120 >= rel_brg > 90)
            return 1
        elseif (270 > rel_brg >= 240)
            return 5
        #if target is aft, place target off closest beam
        elseif 240 > rel_brg > 120
            return 2
        return 6
        end
    end
    if debugging print("close ") end
    
    #behavior if target range is "close"
       
    #if target is ahead, place target off closest beam
    if 90 >= rel_brg >=0
        return 5
    elseif rel_brg >= 270
        return 1
    #if target is abaft of beam, maintain course and speed up
    elseif 270 > rel_brg > 90
        return 4       
    end
    println("Range: ", range, "   Bearing: ", rel_brg)
    return 3    
end