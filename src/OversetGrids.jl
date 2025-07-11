module OversetGrids

using CurvilinearGrids
using WriteVTK
using LinearAlgebra
using StaticArrays

include("components2d.jl")
include("components3d.jl")
export create_components

include("save_vtk.jl")
export save_vtk_with_threshold

include("intersection2d.jl")
include("intersection3d.jl")
include("slice_and_mark3d.jl")
include("slice_and_mark2d.jl")
export slice_interior!, mark_interpolation_cells!

end
