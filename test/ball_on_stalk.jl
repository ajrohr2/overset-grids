using OversetGrids, CurvilinearGrids

ball = rtheta_grid((0.75, (-π/2)+0.2), (1.15, (3π/2)-0.2), (80, 170), :meg6)
stalk_right = RectilinearGrid2D((0.05, -2), (0.75, -0.75), (80, 80), :meg6)
stalk_left = RectilinearGrid2D((-0.75, -2), (-0.05, -0.75), (80, 80), :meg6)
background = RectilinearGrid2D((-1.5, -2), (1.5, 2), (100, 100), :meg6)

meshes = create_components((ball, stalk_left, stalk_right, background); num_interp_points=5)

# Then save to vtk files to view in paraview. An example function call is below.
# save_vtk_with_threshold(meshes[1][1].grid, meshes[1][1].blank_mask, "ball")
