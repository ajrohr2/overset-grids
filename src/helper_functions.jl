function reassociate(meshes::Dict{Int64, Vector{T}}, grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid}}) where T <: ComponentMesh
    return_dict = Dict{Int64, Vector{NamedTuple{(:grid, :component_mesh), Tuple{CurvilinearGrids.AbstractCurvilinearGrid, T}}}}()
    
    idx = 0
    for (key, value) in meshes
        return_dict[key] = Vector{NamedTuple{(:grid, :component_mesh), Tuple{CurvilinearGrids.AbstractCurvilinearGrid, T}}}()
        for mesh in value
            idx += 1
            push!(return_dict[key], (grid=grids[mesh.grid_index], component_mesh=mesh))
        end
    end
    return return_dict
end
