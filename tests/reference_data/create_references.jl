

include("../../armon_2D_MPI_async.jl")
using .Armon


include("reference_functions.jl")


function create_reference_data()
    for type in (Float32, Float64), test in (:Sod, :Sod_y, :Sod_circ, :Bizarrium)
        ref_params = get_reference_params(test, type)
        dt, cycles, data = run_armon_reference(ref_params)
        ref_file_name = get_reference_data_file_name(test, type)
        open(ref_file_name, "w") do ref_file
            write_reference_data(ref_params, ref_file, data, dt, cycles)
        end
    end
end


create_reference_data()
