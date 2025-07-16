using OversetGrids, CurvilinearGrids

inner_sphere = rthetaphi_grid((0.25, 0, 0), (0.75, π, π), (50, 50, 50), :meg6);
outer_sphere = rthetaphi_grid((0.25, 0, 0), (1, π, π), (40, 40, 40), :meg6);

meshes = create_components((inner_sphere, outer_sphere); num_interp_points=3)

# Then save to vtk files to view in paraview. An example function call is below.
# save_vtk_with_threshold(meshes[1][1].grid, meshes[1][1].blank_mask, "inner_sphere")
