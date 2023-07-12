module Lucrest

using HTTP
using JSON
using DataFrames
using Interpolations
using Dates

rest_url = "http://128.131.132.179:8080/lpims/rest/v1/"
user = "student"
pass = "bioVT"
auth = HTTP.base64encode(user * ":" * pass)

url_tails = Dict(
    :running_proc => "processes?running=true",
    :running_reac => "reactors?running=true",
    :reactors => "reactors",
    :alarms => "alarms?running=true&active=true",
)

signal_port_subset(signal_df, port_names) = signal_df[map(x -> x in port_names, signal_df.portName), :]
signal_port_ids(signal_df, port_names) = Dict(signal_port_subset(signal_df, port_names)[!,"portName"] .=> signal_port_subset(signal_df, port_names)[!,"portID"])
signal_signal_ids(signal_df, port_names) = Dict(signal_port_subset(signal_df, port_names)[!,"portName"] .=> signal_port_subset(signal_df, port_names)[!,"signalID"])

function format_response_output(response)
    return JSON.parse(String(response.body))["data"]
end

function get_running_processes()
    calling_str = rest_url * url_tails[:running_proc]
    response = HTTP.request("GET", calling_str, headers = ["Authorization" => "Basic $(auth)"])
    return format_response_output(response)
end

function get_running_reactors()
    calling_str = rest_url * url_tails[:running_reac]
    response = HTTP.request("GET", calling_str, headers = ["Authorization" => "Basic $(auth)"])
    return format_response_output(response)
end

function get_reactor_ids()
    calling_str = rest_url * "reactors/"
    response = HTTP.request("GET", calling_str, headers = ["Authorization" => "Basic $(auth)"])
    jsondata = JSON.parse(String(response.body))["data"]
    name = [reactor["name"] for reactor in jsondata]
    id = [reactor["id"] for reactor in jsondata]
    return Dict(zip(name,id))
end

function get_process_id(name)
    calling_str = rest_url * "processes?name=" * name
    response = HTTP.request("GET", calling_str, headers = ["Authorization" => "Basic $(auth)"])
    return format_response_output(response)[1]["id"]
end

function get_port_id(name)
    calling_str = rest_url * "ports?name=" * name
    response = HTTP.request("GET", calling_str, headers = ["Authorization" => "Basic $(auth)"])
    return JSON.parse(String(response.body))["data"][1]["id"]
end

function get_reactor_id(name)
    calling_str = rest_url * "reactors=" * name
    response = HTTP.request("GET", calling_str, headers = ["Authorization" => "Basic $(auth)"])
    return format_response_output(response)[1]["id"]
end

function get_process_signal_info(process_name)
    calling_str = rest_url * "signals?processId=" * string(get_process_id(process_name))
    response = HTTP.request("GET", calling_str, headers = ["Authorization" => "Basic $(auth)"])
    data = format_response_output(response)
    df = DataFrame()
    df.subDevice = [signal["subDevice"]["name"] for signal in data]
    df.device = [signal["device"]["name"] for signal in data]
    df.portName = [signal["port"]["name"] for signal in data]
    df.signalID = [signal["id"] for signal in data]
    df.portID = [signal["port"]["id"] for signal in data]
    return df
end

function get_signals(process_name, port_names::Vector{String}; interval=0, interpolate=false)
    signal_df = get_process_signal_info(process_name)
    signal_ids = signal_signal_ids(signal_df, port_names)
    signal_str = join(values(signal_ids), ",")
    calling_str = rest_url * "signals?ids=" * signal_str * "&interval=" * string(interval) * "&"
    response = HTTP.request("GET", calling_str, headers = ["Authorization" => "Basic $(auth)"])
    data = JSON.parse(String(response.body))["data"]
    signals = []
    for signal in data
        df = DataFrame()
        df[!,"time"] = [datapoint[1] for datapoint in signal["values"]]
        df[!,signal["port"]["name"]] = [datapoint[2] for datapoint in signal["values"]]
        push!(signals,df)
    end
    if interpolate
        return interpolate_signals(signals)
    end
    return signals
end

function get_signals(process_name, signal_refs::Vector{Tuple{String, Int64}}; interval=0, interpolate=false)
    #signal_df = get_process_signal_info(process_name)
    #signal_ids = signal_signal_ids(signal_df, port_names)
    signal_ids = [signal_ref[2] for signal_ref in signal_refs]
    signal_str = join(signal_ids, ",")
    calling_str = rest_url * "signals?ids=" * signal_str * "&interval=" * string(interval) * "&"
    response = HTTP.request("GET", calling_str, headers = ["Authorization" => "Basic $(auth)"])
    data = JSON.parse(String(response.body))["data"]
    signals = []
    for signal in data
        df = DataFrame()
        df[!,"time"] = [datapoint[1] for datapoint in signal["values"]]
        df[!,signal["port"]["name"]] = [datapoint[2] for datapoint in signal["values"]]
        push!(signals,df)
    end
    if interpolate
        return interpolate_signals(signals)
    end
    return signals
end

function interpolate_signals(signals)
    _,signal_idx = findmax([nrow(signal) for signal in signals])
    interp_linears = [LinearInterpolation(signal[:,1], signal[:,2],extrapolation_bc=Line()) for signal in signals]
    unified_time_vector = signals[signal_idx][:,1]
    unified_signals = [interp_linear.(unified_time_vector) for interp_linear in interp_linears]
    df = DataFrame([names(signal)[2] for signal in signals].=>unified_signals)
    df.time = unified_time_vector
    return df
end

function get_signals(process_name, port_references::Dict; interval=0, tuples=false)
    port_names = [value for value in values(port_references)]
    signal_df = get_signals(process_name, port_names; interval=interval, interpolate=true)
    if ~tuples
        signal_df = [signal_df.time signal_df[!, port_names]]
        rename!(signal_df,vcat(:time, Symbol.(keys(port_references)))) 
    else
        port_names_2 = [value[1] for value in port_names]
        signal_df = [signal_df.time signal_df[!, port_names_2]]
        rename!(signal_df,vcat(:time, Symbol.(keys(port_references)))) 
    end
    return signal_df
end



function get_port_ids(port_names) return [get_port_id(port_name) for port_name in port_names] end

function get_current_values(reactor_name::String, process_name::String, port_names)
    reactor_id = get_reactor_ids()[reactor_name]
    signal_df = get_process_signal_info(process_name)
    port_ids = signal_port_ids(signal_df, port_names)
    port_str = join(values(port_ids), ",")
    calling_str = rest_url * "reactors/" * string(reactor_id) * "?currentValues=" * port_str
    response = HTTP.request("GET", calling_str, headers = ["Authorization" => "Basic $(auth)"])
    data = (JSON.parse(String(response.body))["data"])["process"]["currentValues"]
    names = [port["name"] for port in data]
    currentVals = [port["value"] for port in data]
    return Dict(zip(names,currentVals))
end

# Calling the get_cuurent_values with signal_df and reactor_id is faster than with reactor_name and process name,
# because only one call needs to be made.
function get_current_values(reactor_id::Int, signal_df::DataFrame, port_names)
    port_ids = signal_port_ids(signal_df, port_names)
    port_str = join(values(port_ids), ",")
    calling_str = rest_url * "reactors/" * string(reactor_id) * "?currentValues=" * port_str
    response = HTTP.request("GET", calling_str, headers = ["Authorization" => "Basic $(auth)"])
    data = (JSON.parse(String(response.body))["data"])["process"]["currentValues"]
    names = [port["name"] for port in data]
    currentVals = [port["value"] for port in data]
    return Dict(zip(names,currentVals))
end

function set_current_values(reactor_id::Int, signal_df::DataFrame, write_dict)
    signal_ids = signal_signal_ids(signal_df, collect(keys(write_dict)))
    #port_str = join(values(port_ids), ",")
    for (key, value) in signal_ids
        calling_str = rest_url * "signals/" * string(value)
        #{"Content-Type":"application/json"}
        #response = HTTP.request("PUT", calling_str, "{\"ReactorTemperature_SptExt\" : 20.01}", headers = ["Content-Type" => "application/json","Authorization" => "Basic $(auth)"])
        response = HTTP.put(calling_str, ["Content-Type" => "application/json","Authorization" => "Basic $(auth)"], JSON.json(Dict("currentValue"=>write_dict[key])))
        #data = (JSON.parse(String(response.body))["data"])["process"]["currentValues"]
        #names = [port["name"] for port in data]
        #currentVals = [port["value"] for port in data]
    end
    return nothing
end


############################ REAL TIME CAPABILITIES #########################################################
mutable struct RealTimeLoop
    NOW::DateTime 
    tS_seconds::Float64 
    tS_hours::Float64 
    interval::Float64
    startup_fun::Function
    startup_data::NamedTuple
    timer_fun::Function
    timer_data::NamedTuple
    state::Vector{Any}
end

function timing(rtl::RealTimeLoop)
    rtl.timer_data = rtl.timer_fun(rtl)
    timer_fun = rtl.timer_fun
    tS = (now()-rtl.NOW).value/1000.0
    NOW = rtl.NOW
    println("$NOW: $timer_fun ran with sampling time = $tS seconds\n")
    rtl.tS_seconds = tS
    rtl.tS_hours = tS/3600
    rtl.NOW = now()
    return rtl
end

function start_loop(startup_fun, timer_fun; interval=5)
    # timer = start_loop(TEST_startup,TEST_timerfun;interval=5)
    # close(t3)
    rtl = RealTimeLoop(now(), 0, 0, interval, startup_fun, NamedTuple(), timer_fun, NamedTuple(),[])
    rtl.startup_data = rtl.startup_fun(rtl)
    println("$startup_fun executed")
    start = now()
    sleep(interval)
    rtl.tS_seconds = (now()-start).value/1000.0
    rtl.tS_hours = rtl.tS_seconds/3600
    return Timer(_->timing(rtl), 0, interval=interval)
end

function align_time(signals_online, signals_offline)
    process_start = signals_offline.time[1]
    signals_offline.time = signals_offline.time.-process_start # shift time to start at 0
    signals_online.time = signals_online.time.-process_start
    signals_online = signals_online[signals_online[:,:time].>=0,:]
    signals_online = signals_online[signals_online[:,:time].<=signals_offline.time[end],:]
    return signals_online, signals_offline
end

########### TESTFUNCTIONS FOR CONTROL LOOP ###########################################################
function TEST_startup(rtl::RealTimeLoop)
    reactor_id = Lucrest.get_reactor_ids()["Fermenter 10"]
    signal_df = Lucrest.get_process_signal_info("DWA_Feed_tank_2_20230422")
    rtl.state[1] = 100
    println("PID state initialized at $rtl.state\n")
    return (reactor_id=reactor_id, signal_df=signal_df)
end

function TEST_timerfun(rtl::RealTimeLoop)
    port_names = ["AgitatorSpeed","AIR N2 IN","OxygenIn","DO2Redox"]
    timer_data = Lucrest.get_current_values(rtl.startup_data.reactor_id, rtl.startup_data.signal_df, port_names)
    println(timer_data)
    rtl.state[1] = rtl.state[1] + timer_data["AgitatorSpeed"]/10
    println(rtl.state)
    return (port_names = port_names, data_dict=timer_data)
end



end



