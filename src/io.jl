
"""
    AbstractSolverIO

Base type for file input/output formats.
"""
abstract type AbstractSolverIO end


"""
    supports_mpi(::AbstractSolverIO)
    supports_mpi(::Type{AbstractSolverIO})

If the file format supports MPI, then all sub-domains will write to the same file at the same time.
Otherwise, one file is written per sub-domain, with `"_P=\$(join(params.cart_coords))"` appended
to the end of the file name.
"""
function supports_mpi end


"""
    supports_threads(::AbstractSolverIO)
    supports_threads(::Type{AbstractSolverIO})

If the file format supports multiple threads concurrently reading or writing to the same file.
If not, then blocks are written one by one to the file.
"""
function supports_threads end


"""
    supports_temporal_data(::AbstractSolverIO)
    supports_temporal_data(::Type{AbstractSolverIO})

If the file format supports writing temporal data, i.e. storing data of different cycle to the same file.
Otherwise, one file is written per cycle, with `"_c=\$cycle"` appended to the end of the file name
(after the rank's coordinates, if needed).
"""
function supports_temporal_data end


"""
    file_extension(::AbstractSolverIO)
    file_extension(::Type{AbstractSolverIO})

The extension used by this IO format.
"""
function file_extension end


"""
    format_from_name(name::Symbol)
    format_from_name(::Val{name})

The [`AbstractSolverIO`](@ref) type associated with that name.
"""
format_from_name(name::Symbol) = format_from_name(Val(name))


"""
    domain_writer(
        format::Union{Symbol, Type{<:AbstractSolverIO}}, file,
        params::ArmonParameters, grid::BlockGrid;
        kwargs...
    )

Create a new file for the given `format`, under the prefix `file` (if it is a `String`), to write
the domain represented by `params` and `grid`.

`file` can also be an already opened file object of type `Base.IO` (or of the file type used by the
format).

`vars` are the cell variables to write, it is a `Tuple` of `Symbol`s.

`kwargs` are specific to the `format`.
"""
domain_writer(format::Symbol, file, params, grid; kwargs...) =
    domain_writer(format_from_name(format), file, params, grid; kwargs...)


"""
    domain_reader(
        format::Union{Symbol, Type{<:AbstractSolverIO}}, file,
        params::ArmonParameters, cycle::Union{Int, Nothing}=nothing;
        kwargs...
    )

Open a file of the given `format`, under the prefix `file` (if it is a `String`), matching the domain
represented by `params` at `cycle`.

`file` can also be an already opened file object of type `Base.IO` (or of the file type used by the
format).

`cycle` is either the solver cycle to read from, or `nothing` if there is none.

`vars` are the cell variables to read, it is a `Tuple` of `Symbol`s.

`kwargs` are specific to the `format`.
"""
domain_reader(format::Symbol, file, params, cycle=nothing; kwargs...) =
    domain_reader(format_from_name(format), file, params, cycle; kwargs...)


"""
    close(io::AbstractSolverIO)

Close the file(s) associated with `io`, which then becomes invalid.
"""
Base.close(::AbstractSolverIO) = nothing


"""
    write_domain_to_file(io::AbstractSolverIO, params::ArmonParameters, grid::BlockGrid; kwargs...)

Write the whole `grid` to `io`.
"""
function write_domain_to_file end


"""
    write_block_to_file(io::AbstractSolverIO, params::ArmonParameters, blk::LocalTaskBlock; kwargs...)

Write the `blk` to `io`. Only called when `supports_threads(io) == true`.
"""
function write_block_to_file end


"""
    read_domain_from_file(io::AbstractSolverIO, params::ArmonParameters, grid::BlockGrid; kwargs...)

Read the whole `grid` from `io`.
"""
function read_domain_from_file end


"""
    read_block_from_file(io::AbstractSolverIO, params::ArmonParameters, blk::LocalTaskBlock; kwargs...)

Read the `blk` from `io`. Only called when `supports_threads(io) == true`.
"""
function read_block_from_file end


"""
    write_sub_domain_file(params::ArmonParameters, grid::BlockGrid, file_name::String; options...)
    write_sub_domain_file(format, params::ArmonParameters, grid::BlockGrid, file_name::String; options...)

Write `grid` to `file_name` with the given `format` (defaults to `params.io_format`).
`options` are specific to the format.
"""
function write_sub_domain_file(format, params::ArmonParameters, grid::BlockGrid, file_name::AbstractString; options...)
    if isnothing(params.io_writer)
        writer = domain_writer(format, file_name, params, grid; params.io_options..., options...)
        write_domain_to_file(writer, params, grid)
        close(writer)
    else
        write_domain_to_file(params.io_writer, params, grid)
    end
    return grid
end

write_sub_domain_file(params::ArmonParameters, grid::BlockGrid, file_name::AbstractString; options...) =
    write_sub_domain_file(params.io_format, params, grid, file_name; options...)


"""
    read_sub_domain_file!(params::ArmonParameters, grid::BlockGrid, file_name::String; options...)
    read_sub_domain_file!(format, params::ArmonParameters, grid::BlockGrid, file_name::String; options...)

Read `grid` from `file_name` with the given `format` (defaults to `params.io_format`).
`options` are specific to the format.
"""
function read_sub_domain_file!(format, params::ArmonParameters, grid::BlockGrid, file_name::AbstractString; options...)
    reader = domain_reader(format, file_name, params, grid.global_dt.cycle; params.io_options..., options...)
    read_domain_from_file(reader, params, grid)
    close(reader)
    return grid
end

read_sub_domain_file!(params::ArmonParameters, grid::BlockGrid, file_name::AbstractString; options...) =
    read_sub_domain_file!(params.io_format, params, grid, file_name; options...)


function build_file_path(io_format::ObjOrType{AbstractSolverIO}, file_path::AbstractString, params::ArmonParameters, cycle)
    dir = dirname(file_path)
    if !isempty(dir) && !isdir(dir)
        mkpath(dir)
    end

    file_path, file_ext = splitext(file_path)
    if !isempty(file_ext) && file_ext != file_extension(io_format)
        error("file path '$file_path$file_ext' uses another extension than the format's: '$(file_extension(io_format))'")
    end

    if params.use_MPI && !supports_mpi(io_format)
        file_path *= "_P=" * join(params.cart_coords, 'x')
    end

    if !isnothing(cycle) && !supports_temporal_data(io_format)
        file_path *= "_c=$cycle"
    end

    return file_path * file_extension(io_format)
end
