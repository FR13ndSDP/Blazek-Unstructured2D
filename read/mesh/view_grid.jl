using WriteVTK
using DelimitedFiles

points = readdlm("coords.txt", Float64)
connection = readdlm("./connection.txt", Int64)

cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, connection[i, :]) for i in 1:size(connection)[1]]

vtk_grid("grid", points', cells) do vtk
    # add datasets...
end