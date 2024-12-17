
using CSV
import HDF5

@info "HDF5 library supports MPI: $(HDF5.has_parallel())"
@info "HDF5 library is thread-safe: $(HDF5.h5_is_library_threadsafe())"


function check_vtkhdf_format(filename)
    @test HDF5.ishdf5(filename)
    HDF5.h5open(filename, "r") do file
        @test haskey(file, "VTKHDF")
        vtkhdf = file["VTKHDF"]

        vtk_attrs = HDF5.attributes(vtkhdf)
        @test haskey(vtk_attrs, "Version")
        @test vtk_attrs["Type"] == "OverlappingAMR"
        @test length(vtk_attrs["Origin"]) == 3

        is_temporal = haskey(file, "Steps")
        steps = is_temporal ? file["Steps"] : nothing
        # TODO: steps: values + NSteps

        levels = count(startswith("Level"), keys(file))
        @test levels > 0

        for i_level in 0:levels-1
            level = file["Level$i_level"]

            level_attrs = HDF5.attributes(level)
            @test length(level_attrs["Spacing"]) == 3
            @test count(!iszero, level_attrs["Spacing"]) == 3

            # TODO: AMRBox dataset
            # TODO: {Point|Cell|Field}Data datasets

            if is_temporal
                level_offsets = steps["Level$i_level"]
                # TODO: AMRBoxOffsets, NumberOfAMRBoxes, {Point|Cell|Field}DataOffsets datasets
            end
        end
    end
end


@testset "IO" begin
    @testset "CSV" begin
        
    end

    @testset "HDF5" begin

    end
end
