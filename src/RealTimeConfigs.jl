include("./Lucapi.jl")
include("./observers.jl")
include("./BioprocessControllers.jl")

import .Lucrest
import .Observers
import .BioprocessControllers as control

fermenter_name = "Fermenter 10"
process_name = "DWA_Feed_tank_2_20230422"

function kalman_flow_rate_startup(rtl::Lucrest.RealTimeLoop) # hier vielleicht ein Makro draus machen?
    reactor_id = Lucrest.get_reactor_ids()[fermenter_name] # diese am besten in die start_loop function auslagern
    signal_df = Lucrest.get_process_signal_info(process_name) # diese am besten in die start_loop function auslagern
    rtl.state = [[-1111,0], [10 0; 0 10]]
    println("Kalman states xk and Pk initialized")
    return (reactor_id=reactor_id, signal_df=signal_df)
end

function kalman_flow_rate_timer(rtl::Lucrest.RealTimeLoop)
    port_names = ["OPC_WeightFeed1"]
    data = Lucrest.get_current_values(rtl.startup_data.reactor_id, rtl.startup_data.signal_df, port_names)
    (volume, flow, rtl.state[1], rtl.state[2]) = Observers.kalman_flow_rate(rtl.tS_hours,
    data["OPC_WeightFeed1"],1000,0,rtl.state[1],rtl.state[2])
    println("WeightFeed1: ", data["OPC_WeightFeed1"], "\nVolume : $volume\nFlow : $flow\nx_k : ", rtl.state[1], "\nP_K : ", rtl.state[2])
    return (port_names = port_names, volume=data, flow=flow)
end

function K2S1_startup(rtl::Lucrest.RealTimeLoop) # hier vielleicht ein Makro draus machen?
    reactor_id = Lucrest.get_reactor_ids()[fermenter_name] # diese am besten in die start_loop function auslagern
    signal_df = Lucrest.get_process_signal_info(process_name) # diese am besten in die start_loop function auslagern
    rtl.state = [[0.4,22.5,0,0]]
    println("K2S1 state initialized")
    return (reactor_id=reactor_id, signal_df=signal_df)
end

function K2S1_timer(rtl::Lucrest.RealTimeLoop)
    port_names = ["OPC_WeightFeed1"]
    data = Lucrest.get_current_values(rtl.startup_data.reactor_id, rtl.startup_data.signal_df, port_names)
    (volume, flow, rtl.state[1], rtl.state[2]) = Observers.K2S1_dynamical()
    println("WeightFeed1: ", data["OPC_WeightFeed1"], "\nVolume : $volume\nFlow : $flow\nx_k : ", rtl.state[1], "\nP_K : ", rtl.state[2])
    return (port_names = port_names, volume=data, flow=flow)
end

# Volume: WeightReactor -> Biomare_estimate_VL
# Feed: OPC_WeightFeed2 -> 
#
# 

Feed1_KalmanLoop = ()->Lucrest.start_loop(kalman_flow_rate_startup,kalman_flow_rate_timer;interval=5)
t = Feed1_KalmanLoop()
close(t)