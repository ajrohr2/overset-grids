using OversetGrids, CurvilinearGrids, Test, StaticArrays

@testset verbose = true "UnitTests" begin
    include("unit/full_examples.jl")
    include("unit/test_intersection2d.jl")
    include("unit/test_intersection3d.jl")
    include("unit/test_interpolation.jl")
end
