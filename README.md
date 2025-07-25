# OversetGrids.jl

## Description:
The purpose of this package is to make overlaying grids automatic and easy. It uses a raycasting technique to determine regions of overlap. Depending on the user input, it will then generate an array of -1/0/1/2 values that represent interpolation/untouched/removed/interior cells. This allows solvers to work seamlessly on multiple grids.

## How to use:
1. Generate the grids you'd like to use using `CurvilinearGrids.jl`.
2. Optionally specify the z-order of your grids (i.e. which grids will cut which. A z-order of n will cut all grids with z-order >n).
3. Call `create_components` (see the docstring for more information on how to use this function). This will return a dictionary of (z-order, NamedTuple(:grid, :component_mesh)). The `component_mesh` value of the NamedTuple will store the generated metadata, including the `blank_mask` (this is the array described in the description above).

That's it! Your grids have now been cut!

## Use in a solver
`OversetGrids.jl` provides an interpolation function `interpolate_to_grid!` that will handle interpolating a solution from one grid to another. To use this function, you will need the dictionary that `create_components` creates as well as an array of function values for both grids. See the docstring for specifics.

Examples of both `create_components` and `interpolate_to_grid!` can be found in the test folder.

## A key of blank_mask values
The `blank_mask` array can have 4 different entries.
* -1 denotes an interpolation cell. This is a cell that overlaps the cutting grid, but is still solved on to allow the interpolation function to move the solution between grids.
* 0 denotes an untouched cell. This is just a standard cell.
* 1 denotes a "removed" cell. This is a cell that overlaps the cutting grid and is not an interpolation cell. It's not necessary to solve for this region of a grid.
* 2 denotes an "interior" cell. This is a cell that doesn't overlap the cutting grid, but falls in a region that likely doesn't need to be solved. The prototypical example of such a situation is a closed circular grid with a hole on top of a rectangular grid. The hole in the circular grid likely doesn't need to be solved for on the rectangular grid. Marking these cells is optional, specified by the `mark_interior` argument of `create_components`.
