import Base: size

# only support N=3 for now.
struct BZGrid{N,T}
   # map reducible point to irreducible point (index in irrpos)
   map2ir::Array{Int,N}
   # symmetry operation that map irreducible point to current one.
   map2ir_symop::Array{Int,N}
   # crystal coordinates of the irreducible points.
   irrpos::Vector{SVector{N,T}}
   # symmetry operations used to compute the irreducible points.
   symop::Vector{SMatrix{N,N,T}}
   # tetrahedra, if N=2 (2D case), we need to change 4 to 3.
   tetra::Vector{SVector{4,Int}}
end

# get the size of the BZ grid.
size(g::BZGrid) = size(g.map2ir)

BZGrid(grid_dim::AbstractVector{<:Integer}, symop::AbstractVector{<:SMatrix}) =
   BZGrid(NTuple{3,Int}(grid_dim), symop)

function BZGrid(grid_dim::NTuple{3,Int}, symop::Vector{<:SMatrix{3,3,T}}) where {T}
   map2ir = zeros(Int, grid_dim)
   map2ir_symop = zeros(Int, grid_dim)

   irr_idx = 0
   for pt in CartesianIndices(grid_dim)
      # skip the already labelled reducible points.
      map2ir[pt] > 0 && continue
      # the current point is a irreducible one
      irr_idx += 1
      map2ir[pt] = irr_idx
      ir_pts = pos2pts(grid_dim, pt.I)
      #find reducible points of the current irreducible one
      for (ip, opt) in Iterators.enumerate(symop)
         # applying the symmetry operation
         c_pts = opt * ir_pts
         rpt = pts2pos(grid_dim, c_pts)
         # skip if rpt is nothing, meaning it isn't on the grid.
         rpt === nothing && continue
         # skip if the reducible point is already labelled. It occurs when multiple 
         # symmetry operations connect the irreducible points to the current reducible one 
         rpt = CartesianIndex(rpt)
         map2ir[rpt] > 0 && continue
         # label the current reducible point
         map2ir[rpt] = irr_idx
         map2ir_symop[rpt] = ip
      end
   end
   # collect the crystal coordinate of all the irreducible points.
   # for irreducible points, their corresponding map2ir_symop is 0
   irrpos = Vector{SVector{3,T}}(undef, irr_idx)
   for (i, v) in Iterators.enumerate(findall(iszero, map2ir_symop))
      irrpos[i] = pos2pts(grid_dim, v.I)
   end
   # generate tetrahedra
   tetra = get_tetra(grid_dim)

   return BZGrid{3,T}(map2ir, map2ir_symop, irrpos, symop, tetra)
end

@inline function pos2pts(grid_dim::NTuple{3,Int}, pos::NTuple{3,Int})
   # shift to start from 0
   pos = pos .- 1
   # fold to the Gamma-centered FBZ, -0.5 <= kpt[i] < 0.5, in crystal coordinate
   i_fold = round.(pos ./ grid_dim, RoundNearestTiesUp)
   pos_fold = pos .- grid_dim .* i_fold

   return SVector{3}(pos_fold ./ grid_dim)
end

@inline function pts2pos(grid_dim::NTuple{3,Int}, pts::SVector{3,T}) where {T}
   xpt = pts .* grid_dim
   ipt = round.(xpt)
   # return nothing if the pts isn't on the grid.
   norm(xpt .- ipt) > eps(Float32) && return nothing

   return NTuple{3,Int}((ipt .% grid_dim .+ grid_dim) .% grid_dim .+ 1)
end

function get_tetra(grid_dim::NTuple{3,Int})
   @assert all(grid_dim .> 0) "Negative dimension in grid size: $(grid_dim)!"
   npts = prod(grid_dim)
   tetra = Matrix{SVector{4,Int}}(undef, 6, npts)

   for (i, v) in Iterators.enumerate(CartesianIndices(grid_dim))
      t = _cube_to_tetra(grid_dim, v.I)
      for j in 1:6
         @inbounds tetra[j, i] = t[j]
      end
   end
   return vec(tetra)
end

@inline function _cube_to_tetra(grid_dim::NTuple{3,Int}, start::NTuple{3,Int})
   ## Note that this is different from the fortran version of perturbo
   ## here we use the convention consistent with CartesianIndex in Julia
   ## 1) i, j, k start from 1;  2) here i is fast index, instead of k.
   # i=1, nk1; j=1, nk2; k=1, nk3;
   # num = 1 + (i-1) + (j-1)*nk1 + (k-1)*nk1*nk2; map (i,j,k) to num
   nk1, nk2, nk3 = grid_dim
   i, j, k = start .- 1
   # construct a cube: (i,j,k) -> (ip1, jp1, kp1)
   ip1 = (i + 1) % nk1
   jp1 = (j + 1) % nk2
   kp1 = (k + 1) % nk3
   # n1-n8 are the indices of k-point 1-8 forming a cube
   # num = LinearIndices(A)[i+1, j+1, k+1]
   n1 = 1 + i + j * nk1 + k * nk1 * nk2
   n2 = 1 + i + j * nk1 + kp1 * nk1 * nk2
   n3 = 1 + i + jp1 * nk1 + k * nk1 * nk2
   n4 = 1 + i + jp1 * nk1 + kp1 * nk1 * nk2
   n5 = 1 + ip1 + j * nk1 + k * nk1 * nk2
   n6 = 1 + ip1 + j * nk1 + kp1 * nk1 * nk2
   n7 = 1 + ip1 + jp1 * nk1 + k * nk1 * nk2
   n8 = 1 + ip1 + jp1 * nk1 + kp1 * nk1 * nk2

   # along 1-8, this is a better choice when i, j, k are symmetric
   return (
      #tetra 1: 2, 4, 1, 8;
      SVector{4}(n2, n4, n1, n8),
      #tetra 2: 4, 1, 3, 8;
      SVector{4}(n4, n1, n3, n8),
      #tetra 3: 2, 1, 6, 8;
      SVector{4}(n2, n1, n6, n8),
      #tetra 4: 1, 3, 8, 7;
      SVector{4}(n1, n3, n8, n7),
      #tetra 5: 1, 8, 5, 7;
      SVector{4}(n1, n8, n5, n7),
      #tetra 6: 1, 6, 8, 5;
      SVector{4}(n1, n6, n8, n5),
   )
end