module BioprocessControllers

############## CLASSICAL CONTROL ####################
function PI_control(spt, val, kp, ki, int_last, t_S)
    err = spt - val
    int = clamp(int_last + ki * t_S * err, 0, 1)
    u = clamp(kp * err + int, 0, 1) # clamping u between 0 and 1
    return u, int
end

function PI_control_b(spt, val, kp, ki, int_last, t_S)
    err = spt - val
    int = clamp(int_last + ki * t_S * err, -1, 1)
    u = clamp(kp * err + int, -1, 1) # clamping u between 0 and 1
    return u, int
end

############### CLOSED LOOP CONTROL ##########################

# pO2 Controls 
function pO2_controls(pO2_Spt, pO2, kp, ki, int_last, NSt_min, NSt_max, t_S)
    (u, int) = PI_control(pO2_Spt, pO2, kp, ki, int_last, t_S)
    NSt = NSt_min + u * (NSt_max - NSt_min) # translating u into stirrer speed
    return NSt, int
end

function pump_controls(F_Rw, F_R, slope, offset, kp, ki, int_last, t_S)
    # slope in % / (L/h)
    (u, int) = PI_control_b(F_Rw, F_R, kp, ki, int_last, t_S)
    PUMP = F_Rw * slope + offset + u * 100 # feedforward feedback
    return PUMP, int
end

############### FEEDFORWARD ##################
qS(t, m_X0, q_Sw, c_SR, Y_XS) = m_X0*q_Sw/c_SR * exp(Y_XS*q_Sw*t)
qS(q_Sw, m_X, c_SR) = m_X * q_Sw / c_SR

function impulse(status, period, integral, timespan, t_S, tick, started)
    out = 0
    if status > 0
        if tick < (timespan / t_S)
            started = 1
        end
        if tick > (period / t_S)
            started = 1
            tick = 0
        end
        if started > 0
            if tick < (timespan / t_S)
                out = integral / timespan
            else
                started = 0
            end
        end
        tick = tick + 1
    else
        started = 0
        tick = 0
    end
    return out, tick, started
end

end