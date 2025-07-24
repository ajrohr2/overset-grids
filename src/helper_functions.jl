function reassociate(meshes::Dict{Int64, Vector{T}}, grids::Tuple) where {T <: ComponentMesh}
    return_dict = Dict{Int64, Vector}()
    
    idx = 0
    for (key, value) in meshes
        return_dict[key] = Vector()
        for mesh in value
            idx += 1
            push!(return_dict[key], (grid=grids[mesh.grid_index], component_mesh=mesh))
        end
    end
    return return_dict
end
