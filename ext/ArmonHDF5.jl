module ArmonHDF5

using Armon
import Armon: ObjOrType
using MPI
import HDF5


struct HDF5BlockGridInfo{D} <: Armon.AbstractSolverIO
    file          :: HDF5.File
    box_offset    :: Int            # Number of boxes before this MPI rank 
    boxes         :: HDF5.Dataset   # AMR boxes bounds in each dataset of `cells`
    boxes_offsets :: Array{Int, D}  # Offsets of each block in the `cells` datasets
    cell_offset   :: Int            # Number of global cells before this MPI rank 
    cells         :: Dict{Symbol, HDF5.Dataset}
    # Steps are only present for temporal HDF5 files
    steps         :: Union{Nothing, @NamedTuple{
        nsteps            :: HDF5.Attribute,  # number of steps
        last_cycle        :: HDF5.Attribute,  # last cycle written to the file
        values            :: HDF5.Dataset,    # time step values
        box_offsets       :: HDF5.Dataset,    # offset of each step in `boxes`
        num_boxes         :: HDF5.Dataset,    # number of boxes for each step in `boxes`
        cell_vars_offsets :: Dict{Symbol, HDF5.Dataset}  # offsets of `cells` for each step
    }}
end

Armon.supports_mpi(::ObjOrType{HDF5BlockGridInfo}) = true  # Suppose `HDF5.has_parallel()` returns `true`
Armon.supports_threads(::ObjOrType{HDF5BlockGridInfo}) = HDF5.API.h5_is_library_threadsafe()
Armon.supports_temporal_data(::ObjOrType{HDF5BlockGridInfo}) = true
Armon.file_extension(::ObjOrType{HDF5BlockGridInfo}) = ".vtkhdf"
Armon.format_from_name(::Val{:hdf5}) = HDF5BlockGridInfo


function Armon.domain_writer(
    ::Type{HDF5BlockGridInfo}, filename::String,
    params::ArmonParameters, grid::BlockGrid;
    kwargs...
)
    file_path = Armon.build_file_path(HDF5BlockGridInfo, filename, params, grid.global_dt.cycle)
    file = create_file(file_path, params)
    return write_domain_header(file, params, grid; kwargs...)
end


function Armon.domain_writer(
    ::Type{HDF5BlockGridInfo}, file::HDF5.File,
    params::ArmonParameters, grid::BlockGrid;
    kwargs...
)
    return write_domain_header(file, params, grid; kwargs...)
end


function Armon.domain_reader(
    ::Type{HDF5BlockGridInfo}, filename::AbstractString,
    params::ArmonParameters{T}, cycle=nothing;
    kwargs...
)
    file_path = Armon.build_file_path(HDF5BlockGridInfo, filename, params, cycle)
    error("NYI")  # TODO
    return
end


function Armon.domain_reader(
    ::Type{HDF5BlockGridInfo}, file::IO,
    params::ArmonParameters{T}, cycle=nothing;
    kwargs...
) where {T}
    error("NYI")  # TODO
    return
end


Base.close(info::HDF5BlockGridInfo) = close(info.file)


function create_file(filename::AbstractString, params::ArmonParameters)
    # TODO: mark metadata access as collective for higher performance
    #  => H5Pset_coll_metadata_write
    #  => H5Pset_all_coll_metadata_ops
    # TODO: mark `AMRBox` as `dxpl_mpio=:collective`, since all processes write only once to that dataset
    # TODO: maybe mark all variable datasets as `dxpl_mpio=:collective` (at least when we know that all processes have the same number of blocks)
    return if params.use_MPI
        h5open(filename; driver=HDF5.Drivers.MPIO(params.cart_comm), dxpl_mpio=:collective)
    else
        h5open(filename)
    end
end


function Armon.write_domain_to_file(info::HDF5BlockGridInfo, params::ArmonParameters, grid::BlockGrid)
    if !isnothing(info.steps)
        # Temporal file: if we are writing a new step, we need to update the offsets datasets and
        # expand the cells datasets beforehand.
        prev_cycle = read(info.steps.last_cycle)
        write(info.steps.last_cycle, grid.global_dt.cycle)
        new_step = prev_cycle != -1 && prev_cycle != grid.global_dt.cycle

        nsteps = read(info.steps.nsteps)

        if new_step
            nsteps += 1
            write(info.steps.nsteps, nsteps)

            # TODO: we are iterating `Dict`s, which have a random order in all ranks, is this a problem?

            # Extents are metadata: they must be changed by all processes
            HDF5.set_extent_dims(info.steps.values, (nsteps,))
            HDF5.set_extent_dims(info.steps.box_offsets, (nsteps,))
            HDF5.set_extent_dims(info.steps.num_boxes, (nsteps,))
            for var_dataset in values(info.steps.cell_vars_offsets)
                HDF5.set_extent_dims(var_dataset, (nsteps,))
            end

            # The offset of the step's cell data is the number of cells of all the previous steps
            dims, _ = values(info.cells) |> first |> HDF5.get_extent_dims
            step_offset = last(dims)

            if params.is_root
                # The new values only need to be written once by the root process
                info.steps.values[nsteps] = grid.global_dt.current_dt
                info.steps.box_offsets[nsteps] = 0  # the number of boxes doesn't change (for now)
                info.steps.num_boxes[nsteps] = 0    # idem

                for var_dataset in values(info.steps.cell_vars_offsets)
                    var_dataset[nsteps] = step_offset
                end
            end

            # Since the AMR boxes (our blocks) didn't change, no need to update `boxes` and `boxes_offsets`

            # Extend the cells' datasets by the number of cells
            total_cells = prod(params.global_grid)
            for var_dataset in values(info.cells)
                dims, _ = HDF5.get_extent_dims(var_dataset)
                HDF5.set_extent_dims(var_dataset, (dims[1], step_offset + total_cells))
            end
        else
            # Reuse the offset of the previous step, or 0 if it is the first one.
            # This would allow to overwrite the previous step.
            total_cells = prod(params.global_grid)
            dims, _ = values(info.cells) |> first |> HDF5.get_extent_dims
            step_offset = last(dims) - total_cells
        end
    else
        step_offset = 0
        prev_cycle = -1
    end

    if prev_cycle == -1
        # The first time we write to the file, we also need to write the positions of all blocks
        write_AMR_boxes!(info.boxes, params, grid)
    end
    write_domain!(info.cells, grid; offset=step_offset)
end


function Armon.read_domain_from_file(info::HDF5BlockGridInfo, params::ArmonParameters, grid::BlockGrid)
    error("NYI")  # TODO
end


Armon.write_block_to_file(info::HDF5BlockGridInfo, ::ArmonParameters, blk::Armon.LocalTaskBlock) = write_block!(info, blk)
Armon.read_block_from_file(info::HDF5BlockGridInfo, ::ArmonParameters, blk::Armon.LocalTaskBlock) = read_block!(blk, info)


function compute_global_offsets(params::ArmonParameters, grid::BlockGrid)
    block_offset = 0
    cells_offset = 0
    total_blocks = prod(grid.grid_size)
    if params.use_MPI
        block_offset, cells_offset = MPI.Exscan([total_blocks, prod(params.N)], +, params.cart_comm)
        total_blocks = MPI.Allreduce(total_blocks, +, params.cart_comm)
        if params.rank == 0
            # Exscan results are undefined on rank 0
            cells_offset = 0
            block_offset = 0
        else
            cells_offset -= 1
            block_offset -= 1
        end
    end
    return block_offset, cells_offset, total_blocks
end


function compute_boxes_offsets(grid::BlockGrid)
    boxes_offsets = zeros(Int, grid.grid_size)
    offset = 0
    for block_pos in CartesianIndices(grid.grid_size)
        blk::Armon.LocalTaskBlock = Armon.block_at(grid, block_pos)
        boxes_offsets[block_pos] = offset
        offset += prod(Armon.real_block_size(blk))
    end
    return boxes_offsets
end


function read_domain_header(file::HDF5.File, params::ArmonParameters, grid::BlockGrid)
    l0 = file["VTKHDF"]["Level0"]
    boxes_bb = l0["AMRBox"]
    cells = l0["CellData"]
    cell_vars = Dict(
        Symbol(var_name) => var_dataset
        for (var_name, var_dataset) in pairs(cells)
    )

    if haskey(file["VTKHDF"], "Steps")
        # Temporal file
        steps_group = file["VTKHDF"]["Steps"]
        steps_attrs = HDF5.attributes(steps_group)
        steps_l0 = steps_group["Level0"]
        steps = (;
            nsteps = steps_attrs["NSteps"],
            last_cycle = steps_attrs["ArmonLastCycle"],
            values = steps_group["Values"],
            box_offsets = steps_l0["AMRBoxOffsets"],
            num_boxes = steps_l0["NumberOfAMRBoxes"],
            cell_vars_offsets = Dict(
                Symbol(var_name) => var_dataset
                for (var_name, var_dataset) in pairs(steps_l0["CellDataOffsets"])
            ),
        )
    else
        steps = nothing
    end

    block_offset, cells_offset, _ = compute_global_offsets(params, grid)
    boxes_offsets = compute_boxes_offsets(grid)

    return HDF5BlockGridInfo{D}(file, block_offset, boxes_bb, boxes_offsets, cells_offset, cell_vars, steps)
end


function write_domain_header(
    file::HDF5.File, params::ArmonParameters, grid::BlockGrid{T, D}; vars=Armon.saved_vars(), temporal=true
) where {T, D}
    vtk_root = HDF5.create_group(file, "VTKHDF")

    dim = length(grid.grid_size)
    dim > 3 && error("VTKHDF does not support $(dim)D data")

    top_attrs = HDF5.attributes(vtk_root)
    top_attrs["Type"] = "OverlappingAMR"
    top_attrs["Version"] = [2, 2]

    origin = zeros(Float64, 3)  # [Ox, Oy, Oz]
    origin[1:dim] .= params.origin
    top_attrs["Origin"] = origin

    # Datasets and attributes which aren't part of the VTKHDF format, used mainly to read
    # back the file afterward.
    custom_data = HDF5.create_group(file, "Armon")
    custom_attributes = HDF5.attributes(custom_data)
    # Number of ranks in each direction: since it determines the order in which boxes and
    # cells are stored, it is mandatory when reading back the file.
    custom_attributes["MPI_grid"] = params.proc_dims
    custom_attributes["Dimension"] = D  # dimension of the domain

    # OverlappingAMR is built by box levels. All of our blocks have the same spacing
    # between cells therefore we only use a single level.
    l0 = HDF5.create_group(vtk_root, "Level0")

    l0_attrs = HDF5.attributes(l0)
    spacing = ones(Float64, 3)  # [Δx, Δy, Δz]
    spacing[1:dim] .= params.domain_size ./ params.global_grid
    l0_attrs["Spacing"] = spacing

    block_offset, cells_offset, total_blocks = compute_global_offsets(params, grid)

    # TODO: add chunking to this dataset?
    # `AMRBox` stores a 3D bounding box of indices for each block
    boxes_bb = HDF5.create_dataset(l0, "AMRBox", Int64, (6, total_blocks))
    boxes_offsets = compute_boxes_offsets(grid)

    # TODO: idea to solve the uneven domain issue, which prevent smart chunking:
    #   => compute the maximum amount of cells per domain
    #   => `domain_file_size = chunk_size * cld(max_cells, chunk_size)`
    #   => then each domain writes to its portion of the dataset: `(1:my_cell_count) .+ (domain_file_size * (rank-1) - 1)`
    #   => some cells will remain undefined in the file, but the `AMRBox`es will not reference them
    #   => then we know that two ranks will never access the same chunk, reducing contention? improving performance?
    #   => the extra cost in file size should be small, as long as the chunk size isn't too big
    #   => this could also remove the need for some `MPI_Exscan` operations
    cells = HDF5.create_group(l0, "CellData")
    vars = setdiff(vars, (:x, :y, :z))
    var_size = (var in Armon.dim_vars() ? dim : 1 for var in vars)
    total_cells = prod(params.global_grid)
    cell_vars = Dict(
        # TODO: chunking would be very beneficial as we are writing N-D data in a 1D array
        var => HDF5.create_dataset(cells, string(var), T, (size, total_cells))
        for (var, size) in zip(vars, var_size)
    )

    if temporal
        # VTKHDF temporal format.
        # The data layout stays the same, each step are stored contiguously and the "Steps" group
        # gives the offsets of each time step in those datasets.
        # See https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#temporal-data
        # TODO: we may need to have the cells datasets with `UNLIMITED` dimensions in order to append data to them
        steps_group = HDF5.create_group(vtk_root, "Steps")
        steps_attrs = HDF5.attributes(steps_group)
        steps_attrs["NSteps"] = 1  # each new step will increment this counter
        steps_attrs["ArmonLastCycle"] = -1  # the last cycle written to the file. -1 for none.

        offsets_space = HDF5.dataspace((1,), (-1,))  # single element, but can be appended to by an unlimited amount of times

        values = HDF5.create_dataset(steps_group, "Values", T, offsets_space)  # time value of each step

        l0_offsets     = HDF5.create_group(steps_group, "Level0")
        l0_box_offsets = HDF5.create_dataset(l0_offsets, "AMRBoxOffsets", offsets_space)
        l0_num_boxes   = HDF5.create_dataset(l0_offsets, "NumberOfAMRBoxes", offsets_space)

        cell_offsets = HDF5.create_group(l0_offsets, "CellDataOffsets")
        cell_vars_offsets = Dict(
            var => HDF5.create_dataset(cell_offsets, string(var), T, offsets_space)
            for (var, size) in zip(vars, var_size)
        )

        if params.is_root
            # Only the root has to init the offsets. All offsets are 0 for the first time step.
            # Next steps can have the same AMRBoxOffset if the boxes do not change.
            values[1] = grid.global_dt.current_dt
            l0_box_offsets[1] = 0
            l0_num_boxes[1] = total_blocks
            for var_offsets in values(cell_vars_offsets)
                var_offsets[1] = 0
            end
        end

        steps = (;
            nsteps = steps_attrs["NSteps"],
            last_cycle = steps_attrs["ArmonLastCycle"],
            values,
            box_offsets = l0_box_offsets,
            num_boxes = l0_num_boxes,
            cell_vars_offsets,
        )
    else
        steps = nothing
    end

    return HDF5BlockGridInfo{D}(file, block_offset, boxes_bb, boxes_offsets, cells_offset, cell_vars, steps)
end


function write_AMR_boxes!(info::HDF5BlockGridInfo, params::ArmonParameters, grid::BlockGrid)
    box_indices = Matrix{Int64}(undef, (6, prod(grid.grid_size)))
    for (box_idx, block_pos) in enumerate(CartesianIndices(grid.grid_size))
        global_pos = Armon.block_origin(grid, block_pos) .+ params.N_origin .- 1
        block_size = Armon.block_size_at(grid, block_pos)
        x_min = get(global_pos, 1, 1)
        y_min = get(global_pos, 2, 1)
        z_min = get(global_pos, 3, 1)
        x_max = get(block_size, 1, 1) + x_min - 1
        y_max = get(block_size, 2, 1) + y_min - 1
        z_max = get(block_size, 3, 1) + z_min - 1
        box_indices[:, box_idx] .= (x_min, x_max, y_min, y_max, z_min, z_max) .- 1
    end
    info.boxes[:, axes(box_indices, 2) .+ info.box_offset] = box_indices
end


function write_domain!(info::HDF5BlockGridInfo, grid::BlockGrid; offset=0)
    for block_pos in CartesianIndices(grid.grid_size)
        blk = Armon.block_at(grid, block_pos)
        write_block!(info, blk; offset)
    end
end


function read_AMR_boxes(info::HDF5BlockGridInfo)
    # TODO
    # => produce a block size and grid size
    # => compute the offsets of the blocks
end


function read_domain!(grid::BlockGrid, info::HDF5BlockGridInfo)
    # TODO
end


read_block!(blk::Armon.LocalTaskBlock, info::HDF5BlockGridInfo; offset=0) =
    write_block!(blk::Armon.LocalTaskBlock, info::HDF5BlockGridInfo; read=true, offset=0) 

function write_block!(info::HDF5BlockGridInfo, blk::Armon.LocalTaskBlock; read=false, offset=0)
    block_size = Armon.block_size(blk)
    real_size  = Armon.real_block_size(blk)
    real_cells_count = prod(real_size)

    # `cells_range` is the range of cells of our block/AMRBox in the dataset
    cells_range = HDF5.BlockRange(Base.OneTo(real_cells_count) .+ (offset + info.cell_offset + info.boxes_offsets[blk.pos]))

    # `real_data_slice` are the N-D indices of the real cells in our data
    one_idx = one(CartesianIndex{ndims(blk)})
    real_data_slice = CartesianIndices(real_size) .+ (one_idx * Armon.ghosts(blk))
    subview_range = ntuple(ndims(blk)) do d
        # This is equivalent to: `view(reshape(var_array, block_size), real_data_slice)`
        # as a tuple of `BlockRange` which are then used to create an hyperslab of the data.
        HDF5.BlockRange(real_data_slice.indices[d])
    end

    blk_data = Armon.block_data(blk; on_device=false)
    for (var, dataset) in info.cells
        var_arrays = Armon.var_arrays(blk_data, var)
        for (var_i, var_array) in enumerate(var_arrays)
            file_indices = (var_i, cells_range)
            var_matrix = reshape(var_array, block_size)
            # Equivalent to `dataset[var_i, cell_range] = view(var_matrix, real_data_slice)`
            write_array_to_1D!(dataset, var_matrix, subview_range, file_indices; read)
        end
    end
end


function write_array_to_1D!(dset::HDF5.Dataset, X::Array, src_indices::Tuple, dest_indices::Tuple; read=false)
    # Adapted from HDF5's `Base.setindex!(dset::Dataset, X::Array{T}, I::IndexType...)`
    # This avoids an unneccesary copy of `X` when it is a subview, and supports `BlockRange`
    # as a dataspace index.
    HDF5.get_jl_type(dset) !== eltype(X) && error("mismatched data types")

    filetype = HDF5.datatype(dset)
    memtype = HDF5._memtype(filetype, eltype(X))
    close(filetype)

    dspace = HDF5.dataspace(dset)
    stype = HDF5.API.h5s_get_simple_extent_type(dspace)
    stype == HDF5.API.H5S_NULL && error("attempting to write to null dataspace")

    dspace = HDF5.hyperslab(dspace, dest_indices...)

    memspace = HDF5.dataspace(X)
    memspace = HDF5.hyperslab(memspace, src_indices...)

    if HDF5.API.h5s_get_select_npoints(dspace) != HDF5.API.h5s_get_select_npoints(memspace)
        error("number of elements in src and dest arrays must be equal")
    end

    try
        if read
            HDF5.API.h5d_read(dset, memtype, memspace, dspace, dset.xfer, X)
        else
            HDF5.API.h5d_write(dset, memtype, memspace, dspace, dset.xfer, X)
        end
    finally
        close(memtype)
        close(memspace)
        close(dspace)
    end

    return X
end

end
