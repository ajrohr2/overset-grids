function save_vtk_with_threshold(mesh, blank, fn; centroids=true)
    @info "Writing to $fn.vti"

    xyz_n = CurvilinearGrids.coords(mesh)
    domain = mesh.iterators.cell.domain

    @views vtk_grid(fn, xyz_n) do vtk
        vtk["J", VTKCellData()] = mesh.cell_center_metrics.J[domain]

        # vtk["volume", VTKCellData()] = CurvilinearGrids.cellvolume.(Ref(mesh), domain)

        vtk["xi", VTKCellData(), component_names=["x1", "x2", "t"]] = (
            mesh.cell_center_metrics.ξ.x₁[domain],
            mesh.cell_center_metrics.ξ.x₂[domain],
            mesh.cell_center_metrics.ξ.t[domain],
        )

        vtk["eta", VTKCellData(), component_names=["x1", "x2", "t"]] = (
            mesh.cell_center_metrics.η.x₁[domain],
            mesh.cell_center_metrics.η.x₂[domain],
            mesh.cell_center_metrics.η.t[domain],
        )

        vtk["dx_di", VTKCellData(), component_names=["xi", "eta"]] = (
            mesh.cell_center_metrics.x₁.ξ[domain], mesh.cell_center_metrics.x₁.η[domain]
        )

        vtk["dy_di", VTKCellData(), component_names=["xi", "eta"]] = (
            mesh.cell_center_metrics.x₂.ξ[domain], mesh.cell_center_metrics.x₂.η[domain]
        )
              
        if centroids
            vtk["blank", VTKCellData()] = blank
        else
            vtk["blank", VTKPointData()] = blank
        end
    end
end
