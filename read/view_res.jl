using WriteVTK
using DelimitedFiles

res = readdlm("./res.txt", Float64)
connection = readdlm("./res_connection.txt", Int64)

points = @view res[:, 1:2]

cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, connection[i, :].+1) for i in 1:size(connection)[1]]

vtk_grid("result", points', cells) do vtk
    vtk["rho"] = @view res[:, 3]
    vtk["u"] = @view res[:, 4]
    vtk["v"] = @view res[:, 5]
    vtk["p"] = @view res[:, 6]
    vtk["T"] = @view res[:, 7]
    vtk["M"] = @view res[:, 8]
end