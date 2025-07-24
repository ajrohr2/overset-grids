module OversetGrids

using CurvilinearGrids
using WriteVTK
using LinearAlgebra
using StaticArrays
using NearestNeighbors

include("types.jl")
include("helper_functions.jl")

include("components2d.jl")
include("components3d.jl")
export create_components

include("save_vtk.jl")
export save_vtk_with_threshold

include("intersection2d.jl")
include("intersection3d.jl")
include("slice_and_mark3d.jl")
include("slice_and_mark2d.jl")
export mark_interpolation_cells!

include("interpolation.jl")
export interpolate_to_grid!, error_estimate

end
