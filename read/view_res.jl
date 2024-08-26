# Usage: julia view_res.jl [path_to iso.v2d]

using WriteVTK

f = open(ARGS[1], "r")
case_name = readline(f)
none = readline(f)
none = readline(f)

# read number of variables
var_number = readline(f)
_,nvar = parse.(Int64, split(var_number, " "))

# read variable names
name = Vector{String}(undef, nvar)
for i in 1:nvar
    name[i] = readline(f)
end

# read number of nodes and triangles
none = readline(f)
nNodes,nTria = parse.(Int64, split(readline(f), " "))
none = readline(f)

# read result as a matrix, res[nvar, nNodes]
res = Matrix{Float64}(undef, nvar, nNodes)
for i in 1:nNodes
    res[:, i] = parse.(Float64, split(readline(f), " "))
end

# read triangle connections
connection = Matrix{Int64}(undef, 3, nTria)
for i in 1:nTria
    connection[:, i] = parse.(Int64, split(readline(f), " "))
end


# construct cells, onnection index must start from 1, so +1 here 
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, connection[:, i].+1) for i in 1:nTria]
points = @view res[1:2, :]

vtk_grid(case_name, points, cells) do vtk
    for i in 1:nvar
        vtk[name[i]] = @view res[i, :]
    end
end