@testset "Interpolation" begin
    square = UniformGrid2D((-1, -1), (1, 1), 0.02, :meg6)
    square_domain = square.iterators.cell.domain
    s_centroids_x, s_centroids_y = square.centroid_coordinates.x[square_domain], square.centroid_coordinates.y[square_domain]

    circle = rtheta_grid((0.25, 0), (0.75, 2π), (200, 200), :meg6)
    circle_domain = circle.iterators.cell.domain
    c_centroids_x, c_centroids_y = circle.centroid_coordinates.x[circle_domain], circle.centroid_coordinates.y[circle_domain]
    meshes = create_components((square, circle), (2, 1))

    f(x,y) = exp(sin(x) + cos(y))
    F_square = similar(s_centroids_x)
    F_square .= f.(s_centroids_x, s_centroids_y)

    F_circle = similar(c_centroids_x)
    F_circle_exact = similar(F_circle)
    F_circle_exact .= f.(c_centroids_x, c_centroids_y)

    interpolate_to_grid!(meshes[2][1], meshes[1][1], F_square, F_circle)
    error_matrix = error_estimate.(F_circle_exact, F_circle)
    error_matrix[CartesianIndices((6:200-5, 6:200-5))] .= 0.0 

    @test maximum(error_matrix) ≤ 1e-5
end
