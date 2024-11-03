module ArmonHDF5

using Armon
using MPI
import HDF5


struct HDF5BlockGridInfo{D} <: Armon.AbstractSolverIO
    file          :: HDF5.File
    box_offset    :: Int  # Number of boxes before this MPI rank 
    boxes         :: HDF5.Dataset
    boxes_offsets :: Array{Int, D}  # Offsets of each block in the `cells` datasets
    cell_offset   :: Int  # Number of global cells before this MPI rank 
    cells         :: Dict{Symbol, HDF5.Dataset}
end

Armon.supports_mpi(::HDF5BlockGridInfo) = true  # Suppose `HDF5.has_parallel()` returns `true`
Armon.supports_threads(::HDF5BlockGridInfo) = HDF5.API.h5_is_library_threadsafe()
Armon.supports_temporal_data(::HDF5BlockGridInfo) = true  # TODO


function Armon.domain_writer(::Val{:hdf5}, filename::String, params::ArmonParameters, grid::BlockGrid; kwargs...)
    file = create_file(filename, params)
    return write_domain_header(file, params, grid; kwargs...)
end

function Armon.domain_writer(::Val{:hdf5}, file::HDF5.File, params::ArmonParameters, grid::BlockGrid; kwargs...)
    return write_domain_header(file, params, grid; kwargs...)
end


function create_file(filename::AbstractString, params::ArmonParameters)
    return if params.use_MPI
        h5open(filename; driver=HDF5.Drivers.MPIO(params.cart_comm), dxpl_mpio=:collective)
    else
        h5open(filename)
    end
end


function write_domain(file::HDF5.File, params::ArmonParameters, grid::BlockGrid)
    writer_info = if haskey(file, "VTKHDF")
        get_writer_info(file, params, grid) 
    else
        write_domain_header(file, params, grid)
    end

    write_AMR_boxes!(writer_info.boxes, params, grid)
    write_domain!(writer_info.cells, grid)
end


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


function get_writer_info(file::HDF5.File, params::ArmonParameters, grid::BlockGrid)
    l0 = file["VTKHDF"]["Level0"]
    boxes_bb = l0["AMRBox"]
    cells = l0["CellData"]
    cell_vars = Dict(
        Symbol(var_name) => var_dataset
        for (var_name, var_dataset) in pairs(cells)
    )

    block_offset, cells_offset, _ = compute_global_offsets(params, grid)
    boxes_offsets = compute_boxes_offsets(grid)

    return HDF5BlockGridInfo{D}(file, block_offset, boxes_bb, boxes_offsets, cells_offset, cell_vars)
end


function write_domain_header(
    file::HDF5.File, params::ArmonParameters, grid::BlockGrid{T, D}; vars=Armon.saved_vars()
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
    boxes_offsets = boxes_offsets(grid)

    cells = HDF5.create_group(l0, "CellData")
    vars = setdiff(vars, (:x, :y, :z))
    var_size = (var in Armon.dim_vars() ? dim : 1 for var in vars)
    total_cells = prod(params.global_grid)
    cell_vars = Dict(
        # TODO: chunking would be very beneficial as we are writing N-D data in a 1D array
        var => HDF5.create_dataset(cells, string(var), T, (size, total_cells))
        for (var, size) in zip(vars, var_size)
    )

    return HDF5BlockGridInfo{D}(file, block_offset, boxes_bb, boxes_offsets, cells_offset, cell_vars)
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


function write_domain!(info::HDF5BlockGridInfo, grid::BlockGrid)
    for block_pos in CartesianIndices(grid.grid_size)
        blk = Armon.block_at(grid, block_pos)
        write_block!(info, blk)
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


read_block!(blk::Armon.LocalTaskBlock, info::HDF5BlockGridInfo) =
    write_block!(blk::Armon.LocalTaskBlock, info::HDF5BlockGridInfo; read=true) 

function write_block!(info::HDF5BlockGridInfo, blk::Armon.LocalTaskBlock; read=false)
    block_size = Armon.block_size(blk)
    real_size  = Armon.real_block_size(blk)
    real_cells_count = prod(real_size)

    # `cells_range` is the range of cells of our block/AMRBox in the dataset
    cells_range = HDF5.BlockRange(Base.OneTo(real_cells_count) .+ (info.cell_offset + info.boxes_offsets[blk.pos]))

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
