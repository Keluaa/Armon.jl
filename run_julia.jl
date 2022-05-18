
using Printf


scheme = :GAD_minmod
riemann = :acoustic
iterations = 4
cfl = 0.6
Dt = 0.
maxtime = 0.0
maxcycle = 500
euler_projection = false
cst_dt = false
ieee_bits = 64
silent = 2
write_output = false
use_ccall = false
use_threading = true
use_simd = true
interleaving = false
use_gpu = false

tests = []
cells_list = []

base_file_name = ""
gnuplot_script = ""


i = 1
while i <= length(ARGS)
    arg = ARGS[i]
    if arg == "-s"
        global scheme=Symbol(replace(ARGS[i+1], '-' => '_'))
        global i += 1
    elseif arg == "--ieee"
        global ieee_bits = parse(Int, ARGS[i+1])
        global i += 1
    elseif arg == "--cycle"
        global maxcycle = parse(Int, ARGS[i+1])
        global i += 1
    elseif arg == "--riemann"
        global riemann = Symbol(replace(ARGS[i+1], '-' => '_'))
        global i += 1
    elseif arg == "--iterations"
        global iterations = parse(Int, ARGS[i+1])
        global i += 1
    elseif arg == "--verbose"
        global silent = parse(Int, ARGS[i+1])
        global i += 1
    elseif arg == "--time"
        global maxtime = parse(Float64, ARGS[i+1])
        global i += 1
    elseif arg == "--cfl"
        global cfl = parse(Float64, ARGS[i+1])
        global i += 1
    elseif arg == "--dt"
        global Dt = parse(Float64, ARGS[i+1])
        global i += 1
    elseif arg == "--write-output"
        global write_output = parse(Bool, ARGS[i+1])
        global i += 1
    elseif arg == "--use-ccall"
        global use_ccall = parse(Bool, ARGS[i+1])
        global i += 1
    elseif arg == "--use-threading"
        global use_threading = parse(Bool, ARGS[i+1])
        global i += 1
    elseif arg == "--use-simd"
        global use_simd = parse(Bool, ARGS[i+1])
        global i += 1
    elseif arg == "--interleaving"
        global interleaving = parse(Bool, ARGS[i+1])
        global i += 1
    elseif arg == "--use-gpu"
        global use_gpu = parse(Bool, ARGS[i+1])
        global i += 1
    elseif arg == "--euler"
        global euler_projection = parse(Bool, ARGS[i+1])
        global i += 1
    elseif arg == "--cst-dt"
        global cst_dt = parse(Bool, ARGS[i+1])
        global i += 1

    # Additionnal params
    elseif arg == "--gpu"
        gpu = ARGS[i+1]
        global i += 1
        if gpu == "ROCM"
            ENV["USE_ROCM_GPU"] = "true"
        elseif gpu == "CUDA"
            ENV["USE_ROCM_GPU"] = "false"
        else
            println("Unknown gpu: ", gpu)
            exit(1)
        end
        global use_gpu = true
    elseif arg == "--block-size"
        block_size = parse(Int, ARGS[i+1])
        global i += 1
        ENV["GPU_BLOCK_SIZE"] = block_size
    elseif arg == "--tests"
        global tests = Symbol.(split(ARGS[i+1], ','))
        global i += 1
    elseif arg == "--cells-list"
        list = split(ARGS[i+1], ',')
        global cells_list = convert.(Int, parse.(Float64, list))
        global i += 1
    elseif arg == "--data-file"
        global base_file_name = ARGS[i+1]
        global i += 1
    elseif arg == "--gnuplot-script"
        global gnuplot_script = ARGS[i+1]
        global i += 1
    elseif arg == "--use-std-threads"
        use_std_threads = parse(Bool, ARGS[i+1])
        if use_std_threads
            ENV["USE_STD_LIB_THREADS"] = "true"
        else
            ENV["USE_STD_LIB_THREADS"] = "false"
        end
        global i += 1
    else
        println("Wrong option: ", arg)
        exit(1)
    end
    global i += 1
end


println("Loading...")
include("armon_module_gpu.jl")
using .Armon


println("Compiling...")
for test in tests
    armon(ArmonParameters(; ieee_bits, test=test, riemann, scheme, iterations, nbcell=10000, cfl, Dt, euler_projection, cst_dt, maxtime, maxcycle=1, silent=5, write_output=false, use_ccall, use_threading, use_simd, interleaving, use_gpu))
end


for test in tests
    data_file_name = base_file_name * string(test) * ".csv"

    for cells in cells_list
        @printf(" - %s, %-10g cells: ", test, cells)

        cells_per_sec = armon(ArmonParameters(; ieee_bits, test, riemann, scheme, iterations, nbcell=cells, cfl, Dt, euler_projection, cst_dt, maxtime, maxcycle, silent, write_output, use_ccall, use_threading, use_simd, interleaving, use_gpu))
        
        @printf("%.2g Giga cells/sec\n", cells_per_sec)

        # Append the result to the data file
        open(data_file_name, "a") do data_file
            println(data_file, cells, ", ", cells_per_sec)
        end

        if !isempty(gnuplot_script)
            # We redirect the output of gnuplot to null so that there is no warning messages displayed
            run(pipeline(`gnuplot $(gnuplot_script)`, stdout=devnull, stderr=devnull))
        end
    end
end
