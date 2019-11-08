function atan2(x,y)
    if x > 0
        return atan(y/x)
    elseif (x < 0) & (y >= 0)
        return atan(y/x) + π
    elseif (x < 0) & (y < 0)
        return atan(y/x) - π
    elseif (x == 0) & (y > 0)
        return π/2
    elseif (x == 0) & (y < 0)
        return -π/2
    elseif (x == 0) & (y == 0)
        return false
    end
end