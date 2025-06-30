module OversetGrids

using CurvilinearGrids
using WriteVTK

include("components.jl")
export create_components

include("save_vtk.jl")
export save_vtk_with_threshold

include("intersection.jl")
include("slice_and_mark.jl")
export slice_interior!, mark_interpolation_cells!

end
