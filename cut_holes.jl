function cut_holes(mesh1, mesh2; overlap_cell_num=5)
    # For now, we assume that mesh1 is completely contained in mesh2

    cut_toward_center = true

    x_m1, y_m1 = CurvilinearGrids.coords(mesh1)
    x_m2, y_m2 = CurvilinearGrids.coords(mesh2)

    domain_m1, domain_m2 = mesh1.iterators.cell.domain, mesh2.iterators.cell.domain

    overlap_points_m2 = []

    centroid_m1_x = mesh1.centroid_coordinates.x[domain_m1]
    centroid_m1_y = mesh1.centroid_coordinates.y[domain_m1]
    boundary_of_object = mesh1.centroid_coordinates.x[domain_m1[1,:]]

    blanked_cells = shift(domain_m1[1:overlap_cell_num,:], -5)

    for c_m1 in blanked_cells, c_m2 in CartesianIndices((1:size(x_m2)[1], 1:size(x_m2)[2]))
        tolerance_x = max(abs(mesh1.cell_center_metrics.x₁.ξ[domain_m1][c_m1]), abs(mesh1.cell_center_metrics.x₁.η[domain_m1][c_m1])) # * 0.6
        tolerance_y = max(abs(mesh1.cell_center_metrics.x₂.ξ[domain_m1][c_m1]), abs(mesh1.cell_center_metrics.x₂.η[domain_m1][c_m1])) # * 0.6
        if abs(centroid_m1_x[c_m1]-x_m2[c_m2]) ≤ tolerance_x && abs(centroid_m1_y[c_m1]-y_m2[c_m2]) ≤ tolerance_y
            push!(overlap_points_m2, c_m2)
        end
    end

    return unique(overlap_points_m2)
end
