function plot_meshes(mesh1, mesh2)
    x_m1, y_m1 = CurvilinearGrids.coords(mesh1)
    x_m1_flat, y_m1_flat = reshape(x_m1, prod(size(x_m1))), reshape(y_m1, prod(size(y_m1)))
           
    x_m2, y_m2 = CurvilinearGrids.coords(mesh2)
    x_m2_flat, y_m2_flat = reshape(x_m2, prod(size(x_m2))), reshape(y_m2, prod(size(y_m2)))
    
    p = plot(aspect_ratio=:equal)
    
    scatter!(p, x_m1_flat, y_m1_flat, markersize=0.25, label="")
    scatter!(p, x_m2_flat, y_m2_flat, markersize=0.25, label="")

    for i in 1:size(x_m1)[1]
        for j in 1:size(x_m1, 2)-1  # Loop over each pair of consecutive points in a row
            plot!(p, [x_m1[i, j], x_m1[i, j+1]], [y_m1[i, j], y_m1[i, j+1]], line=:solid, label="", color=:black)
        end
    end
    
    for j in 1:size(x_m1)[2]
        for i in 1:size(x_m1, 1)-1  # Loop over each pair of consecutive points in a column
            plot!(p, [x_m1[i, j], x_m1[i+1, j]], [y_m1[i, j], y_m1[i+1, j]], line=:solid, label="", color=:black)
        end
    end
    
    for i in 1:size(x_m2)[1]
        for j in 1:size(x_m2, 2)-1  # Loop over each pair of consecutive points in a row
            plot!(p, [x_m2[i, j], x_m2[i, j+1]], [y_m2[i, j], y_m2[i, j+1]], line=:solid, label="", color=:red)
        end
    end
    
    for j in 1:size(x_m2)[2]
        for i in 1:size(x_m2, 1)-1  # Loop over each pair of consecutive points in a column
            plot!(p, [x_m2[i, j], x_m2[i+1, j]], [y_m2[i, j], y_m2[i+1, j]], line=:solid, label="", color=:red)
        end
    end
    
    return p
end
