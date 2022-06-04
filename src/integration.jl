"""
Brillouin Zone Integration using the tetrahedron method
"""
function compute_dos(
   kgrid::BZGrid,
   ene_bands::AbstractArray{T,2},
   egrid::AbstractVector{T},
) where {T}
   numb, nkpts = size(ene_bands)
   @assert nkpts == length(kgrid.irrpos) "num of kpts does not match !"

   dos = zeros(T, size(egrid))
   tetra, map2ir = kgrid.tetra, kgrid.map2ir

   fnk = ones(SVector{4,T})
   for tet in tetra, ib in 1:numb
      enk = SVector{4,T}(ene_bands[ib, map2ir[ik]] for ik in tet)
      tetra_int!(enk, fnk, egrid, dos)
   end
   return dos ./ length(tetra)
end