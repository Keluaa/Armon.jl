
import Armon: Axis, Side, Neighbours

@testset "Dimension agnostism" begin

@testset "Axis" begin
    @test Integer(Axis.T(0x1)) == UInt8(Axis.T(1)) == 0x1
    @test typeof(Axis.T(0x1)) == Axis.T
    @test typeof(Integer(Axis.T(0x1))) == UInt8

    @test sizeof(Axis.T) == sizeof(UInt8)
    @test Int(typemin(Axis.T)) == 1
    @test Int(typemax(Axis.T)) == 127
    @test typeof(Base.cconvert(UInt8, Axis.X)) == UInt8

    @test_throws "Axis index" Axis.T(128)
    @test_throws "Axis index" Axis.T(0)
    @test_throws "should not be iterated" instances(Axis.T)
    @test_throws "should not be iterated" Base.Enums.namemap(Axis.T)

    @test length(Armon.axes_of(8)) == 8

    @test hash(Axis.X) != hash(1) != hash(Axis.Y)
    @test length(unique(hash.(Armon.axes_of(127)))) == 127

    @test UInt8(Axis.X) == UInt8(Axis.T(1)) == 0x1
    @test UInt8(Axis.Y) == UInt8(Axis.T(2)) == 0x2
    @test UInt8(Axis.Z) == UInt8(Axis.T(3)) == 0x3

    @test Symbol(Axis.X) === :X
    @test Symbol(Axis.Y) === :Y
    @test Symbol(Axis.Z) === :Z
    @test Symbol(Axis.T(4)) === :Axis_4
    @test Symbol(Axis.T(127)) === :Axis_127

    @test count('\n', sprint(Base.show, MIME"text/plain"(), Axis.T)) == 8
    @test sprint(Base.show, MIME"text/plain"(), Axis.X) == "X::Axis.T = 1"
end


@testset "Side" begin
    @test Integer(Side.T(0x1)) == UInt8(Side.T(1)) == 0x1
    @test typeof(Side.T(0x1)) == Side.T
    @test typeof(Integer(Side.T(0x1))) == UInt8

    @test sizeof(Side.T) == sizeof(UInt8)
    @test Int(typemin(Side.T)) == 1
    @test Int(typemax(Side.T)) == 254
    @test typeof(Base.cconvert(UInt8, Side.Left)) == UInt8

    @test_throws "Side index" Side.T(255)
    @test_throws "Side index" Side.T(0)
    @test_throws "should not be iterated" instances(Side.T)
    @test_throws "should not be iterated" Base.Enums.namemap(Side.T)

    @test length(Armon.sides_of(8)) == 16

    @test hash(Side.Left) != hash(1) != hash(Side.Right)
    @test length(unique(hash.(Armon.sides_of(127)))) == 254

    @test UInt8(Side.Left)   == UInt8(Side.T(1)) == 0x1
    @test UInt8(Side.Right)  == UInt8(Side.T(2)) == 0x2
    @test UInt8(Side.Bottom) == UInt8(Side.T(3)) == 0x3
    @test UInt8(Side.Top)    == UInt8(Side.T(4)) == 0x4
    @test UInt8(Side.Back)   == UInt8(Side.T(5)) == 0x5
    @test UInt8(Side.Front)  == UInt8(Side.T(6)) == 0x6

    @test Symbol(Side.Left)   === :Left
    @test Symbol(Side.Right)  === :Right
    @test Symbol(Side.Bottom) === :Bottom
    @test Symbol(Side.Top)    === :Top
    @test Symbol(Side.Back)   === :Back
    @test Symbol(Side.Front)  === :Front
    @test Symbol(Side.T(7))   === :Side_7
    @test Symbol(Side.T(254)) === :Side_254

    @test count('\n', sprint(Base.show, MIME"text/plain"(), Side.T)) == 8
    @test sprint(Base.show, MIME"text/plain"(), Side.Left) == "Left::Side.T = 1"
end


@testset "Axes and Sides" begin
    @test Armon.sides_along(Axis.X) == (Side.Left, Side.Right)
    @test Armon.sides_along(Axis.Y) == (Side.Bottom, Side.Top)

    @testset "$(D)D" for D in (1, 2, 3, 5, 16, 127)
        VD = Val{D}()
        A = Axis.T(D)
        FS = Armon.first_side(A)
        LS = Armon.last_side(A)
        
        @test length(Armon.axes_of(D))  == length(Armon.axes_of(VD))  == D
        @test length(Armon.sides_of(D)) == length(Armon.sides_of(VD)) == 2D

        @test Armon.next_axis(A, D) == Armon.next_axis(A, VD) == Axis.X
        @test Armon.offset_to(A, D) == Armon.offset_to(A, VD)
        @test count(iszero, Armon.offset_to(A, D)) == D - 1

        @test Armon.sides_along(A) == (FS, LS)

        @test Armon.first_sides(D) == Armon.first_sides(VD)
        @test Armon.last_sides(D)  == Armon.last_sides(VD)
        @test length(Armon.first_sides(D)) == length(Armon.last_sides(D)) == D
        @test (Side.Left,  FS) ⊆ Armon.first_sides(D)
        @test (Side.Right, LS) ⊆ Armon.last_sides(D)

        @test  Armon.first_side(FS)
        @test !Armon.first_side(LS)
        @test  Armon.last_side(LS)
        @test !Armon.last_side(FS)

        @test all(Armon.axis_of.(Armon.sides_along(A)) .== A)

        @test Armon.opposite_of(FS) == LS
        @test Armon.opposite_of(LS) == FS

        @test Armon.offset_to(FS, D) == Armon.offset_to(FS, VD)
        @test sum(Armon.offset_to(FS, D)) == -1
        @test findall(==(-1), Armon.offset_to(FS, 127)) == [D]
        @test sum(Armon.offset_to(LS, D)) == 1
        @test findall(==(1), Armon.offset_to(LS, 127)) == [D]

        @test Armon.side_from_offset(Armon.offset_to(FS, D)) == FS
        @test Armon.side_from_offset(Armon.offset_to(LS, D)) == LS
        @test Armon.side_from_offset(Armon.offset_to(FS, 127)) == FS
        @test Armon.side_from_offset(Armon.offset_to(LS, 127)) == LS
    end
end


@testset "Neighbours" begin
    @test typeof(Neighbours(((1, 2),))) == Neighbours{Int, 1}

    @test length(Neighbours(Returns(0), 0)) == 0
    @test length(Neighbours(Returns(0), 1)) == 1
    @test length(Neighbours(Returns(0), 9)) == 9

    n2D = Neighbours((_, s) -> s, 2)
    @test size(n2D) == (2,)
    @test length(n2D) == 2
    @test eltype(eachindex(n2D)) == Int
    @test all(n2D.val .== collect(n2D) .== ((Side.Left, Side.Right), (Side.Bottom, Side.Top)))
    @test n2D[1] == n2D[begin] == n2D[Axis.X] == n2D.X == (Side.Left, Side.Right)
    @test n2D[2] == n2D[end]   == n2D[Axis.Y] == n2D.Y == (Side.Bottom, Side.Top)
    @test n2D[Side.Left] == n2D.Left == Side.Left
    @test n2D[Side.Top]  == n2D.Top  == Side.Top
    @test all(getindex.(Ref(n2D), Armon.sides_of(2)) .== Armon.sides_of(2))
    @test all(getindex.(Ref(n2D), Armon.axes_of(2))  .== n2D.val)

    @test_throws BoundsError n2D.Z
    @test_throws BoundsError n2D.Front
end

end
