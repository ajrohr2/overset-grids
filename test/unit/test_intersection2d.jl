@testset "Orientation" begin
    A = (0, 0)
    B = (0, 1)
    C = (1, 0)

    @test OversetGrids.orientation(A, B, C) == -1
end

@testset "Intersection" begin
    # Intersecting lines
    int_line1 = OversetGrids.Line((0, 0), (0, 1))
    int_line2 = OversetGrids.Line((1, 0), (-1, 1))

    @test OversetGrids.determine_intersection(int_line1, int_line2) == true

    # Non-intersecting lines
    nint_line1 = OversetGrids.Line((0, 0), (1, 0))
    nint_line2 = OversetGrids.Line((0, 1), (0, 2))

    @test OversetGrids.determine_intersection(nint_line1, nint_line2) == false

    # Intersecting, ABC colinear
    abc_line1 = OversetGrids.Line((0, 0), (1, 0))
    abc_line2 = OversetGrids.Line((0.5, 0), (0, 1))

    @test OversetGrids.determine_intersection(abc_line1, abc_line2) == true

    # Intersecting, ABD colinear
    abd_line1 = OversetGrids.Line((0, 0), (1, 0))
    abd_line2 = OversetGrids.Line((0, 1), (0.5, 0))

    @test OversetGrids.determine_intersection(abd_line1, abd_line2) == true

    # Intersecting, CDA colinear
    cda_line1 = OversetGrids.Line((0, 0), (1, 0))
    cda_line2 = OversetGrids.Line((0, 1), (0, -1))

    @test OversetGrids.determine_intersection(cda_line1, cda_line2) == true

    # Intersecting, CDB colinear
    cdb_line1 = OversetGrids.Line((1, 0), (0, 0))
    cdb_line2 = OversetGrids.Line((0, 1), (0, -1))

    @test OversetGrids.determine_intersection(cdb_line1, cdb_line2) == true

    # All points colinear, shouldn't Intersect 
    all_line1 = OversetGrids.Line((0, 0), (1, 0))
    all_line2 = OversetGrids.Line((-1, 0), (0.5, 0))

    @test OversetGrids.determine_intersection(all_line1, all_line2) == false
end

@testset "Boundary Creation" begin
    grid = UniformGrid2D((0, 0), (1, 1), 0.1, :meg6)
    bounds, polygon = OversetGrids.get_boundary(grid)

    @test bounds == (0, 1, 0, 1)

    correct_polygon_left = [OversetGrids.Line((0, i), (0, i+0.1)) for i in 0:0.1:0.9]
    correct_polygon_right = [OversetGrids.Line((1, i), (1, i+0.1)) for i in 0:0.1:0.9]
    correct_po1ygon_top = [OversetGrids.Line((i, 1), (i+0.1, 1)) for i in 0:0.1:0.9]
    correct_polygon_bottom = [OversetGrids.Line((i, 0), (i+0.1, 0)) for i in 0:0.1:0.9]

    correct_polygon = vcat(correct_po1ygon_top, correct_polygon_left, correct_polygon_right, correct_polygon_bottom)

    @test length(correct_polygon) == length(polygon)

    for line in polygon
        @test any(l -> isapprox(l.point_1[1], line.point_1[1]), correct_polygon) && any(l -> isapprox(l.point_1[2], line.point_1[2]), correct_polygon) && any(l -> isapprox(l.point_0[1], line.point_0[1]), correct_polygon) && any(l -> isapprox(l.point_0[2], line.point_0[2]), correct_polygon)
    end
end
