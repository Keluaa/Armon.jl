
struct CSVSolverIO{IO_t} <: AbstractSolverIO
    file          :: IO_t
    vars          :: Tuple{Vararg{Symbol}}
    precision     :: Int
    all_ghosts    :: Bool
    global_ghosts :: Bool
    for_3D        :: Bool
end

supports_mpi(::ObjOrType{CSVSolverIO}) = false
supports_threads(::ObjOrType{CSVSolverIO}) = false
supports_temporal_data(::ObjOrType{CSVSolverIO}) = false
file_extension(::ObjOrType{CSVSolverIO}) = ".csv"
format_from_name(::Val{:csv}) = CSVSolverIO


function CSVSolverIO(
    file::IO, ::Type{T};
    vars=saved_vars(), precision=nothing, all_ghosts=false, global_ghosts=false, for_3D=true
) where {T}
    if isnothing(precision)
        # Exact decimal output by default
        if T === Float16
            precision = 4
        elseif T === Float32
            precision = 9
        elseif T === Float64
            precision = 17
        else
            error("no default precision for $T")
        end
    end
    return CSVSolverIO{typeof(file)}(file, vars, precision, all_ghosts, global_ghosts, for_3D)
end


function domain_writer(::Type{CSVSolverIO}, filename::AbstractString, params::ArmonParameters{T}, grid::BlockGrid; kwargs...) where {T}
    file_path = build_file_path(CSVSolverIO, filename, params, grid.global_dt.cycle)
    file = open(file_path, "w")
    return CSVSolverIO(file, T; kwargs...)
end


function domain_writer(::Type{CSVSolverIO}, file::IO, params::ArmonParameters{T}, grid::BlockGrid; kwargs...) where {T}
    return CSVSolverIO(file, T; kwargs...)
end


function domain_reader(::Type{CSVSolverIO}, filename::AbstractString, params::ArmonParameters{T}, cycle=nothing; kwargs...) where {T}
    file_path = build_file_path(CSVSolverIO, filename, params, cycle)
    file = open(file_path, "r")
    return CSVSolverIO(file, T; kwargs...)
end


function domain_reader(::Type{CSVSolverIO}, file::IO, params::ArmonParameters{T}, cycle=nothing; kwargs...) where {T}
    return CSVSolverIO(file, T; kwargs...)
end


Base.close(csv::CSVSolverIO) = close(csv.file)


function write_domain_to_file(csv::CSVSolverIO, ::ArmonParameters, grid::BlockGrid; row_iter_params=())
    p = csv.precision
    (; global_ghosts, all_ghosts) = csv

    format = Printf.Format(join(repeat(["%#$(p+7).$(p)e"], length(csv.vars)), ", ") * "\n")

    # Write cells in the correct ascending (X, Y, Z) order, combining the cells of all blocks
    prev_row_idx = nothing
    for (blk, row_idx, row_range) in BlockRowIterator(grid, row_iter_params...; global_ghosts, all_ghosts)
        row_idx = row_idx[2:end]
        if csv.for_3D && prev_row_idx != row_idx && !isnothing(prev_row_idx)
            println(csv.file)  # Separate rows to use pm3d plotting with gnuplot
        end

        blk_vars = var_arrays(blk, csv.vars; on_device=false)
        for idx in row_range
            Printf.format(csv.file, format, getindex.(blk_vars, idx)...)
        end

        prev_row_idx = row_idx
    end
end


function read_domain_from_file(csv::CSVSolverIO, ::ArmonParameters{T}, grid::BlockGrid; row_iter_params=()) where {T}
    (; global_ghosts, all_ghosts) = csv
    for (blk, _, row_range) in BlockRowIterator(grid, row_iter_params...; global_ghosts, all_ghosts)
        blk_vars = var_arrays(blk, csv.vars; on_device=false)
        for idx in row_range
            for var in blk_vars[1:end-1]
                var[idx] = parse(T, readuntil(csv.file, ','))
            end
            blk_vars[end][idx] = parse(T, readuntil(csv.file, '\n'))
        end
    end
end
