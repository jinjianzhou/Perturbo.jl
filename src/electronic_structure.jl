using StaticArrays, LinearAlgebra


function bands(el::ElecHam, kpt::SVector)
   nband = el.basic_data[:num_wann]
   Hk = zeros(ComplexF64, nband, nband)
   _band!(Hk, el, kpt)
end

function bands(el::ElecHam, kpts::AVec{<:SVector})
   nband = el.basic_data[:num_wann]
   Hk = zeros(ComplexF64, nband, nband)

   eigs = zeros(Float64, nband, length(kpts))

   for (i, k) in enumerate(kpts)
      eigs[:,i] = _band!(Hk, el, k)
   end
   return eigs .* Rydberg2eV
end

function bands(el::ElecHam, kpts::AbstractArray{T,2}) where T
   size(kpts, 1) != 3 && throw(error("1-dimension of kpts should be 3"))
   bands(el, vec(reinterpret(SVector{3,T}, kpts)))
end

function _band!(Hk::AbstractMatrix, el::ElecHam, kpt::SVector)
   hamiltonian!(Hk, el, kpt)
   eigvals( Hermitian(Hk) ) .* Rydberg2eV
end

exp_ikr(rvecs::AbstractVector{<:SVector}, kpt::SVector) = cis.(Ref(kpt') .* rvecs)

_ham_elem(hopping, rvec_idx, e_ikr) = sum( e_ikr[val] * hopping[n] for (n, val) in enumerate(rvec_idx) )

function hamiltonian!(H::AbstractMatrix, el::ElecHam, kpt::SVector)
   nband = LinearAlgebra.checksquare(H)
   nband == el.basic_data[:num_wann] || throw(DimensionMismatch("H is not square: dimensions are $(size(H))"))
   
   #compute exp(i 2\pi*K*R)
   e_ikr = exp_ikr(el.ws_rvecs, 2π*kpt)

   # compute upper triangle of ham
   k = 0
   @inbounds for j = 1:nband, i = 1:j
      k += 1
      H[i,j] = _ham_elem(el.hopping[k], el.rvec_idx[k], e_ikr)
   end
end

#function bands(ph::LatticeIFC, qpt::SVector)
#   nsize = 3*ph.basic_data[:nat]
#   dynmat = zeros(ComplexF64, nsize, nsize)
#   _dynamic_matrix!(dynmat, ph, qpt)
#end

#function _dynamic_matrix!(dynmat::AbstractMatrix, ph::LatticeIFC, qpt::SVector) where T
#   n, nsize = size(dynmat)
#   (n == nsize && nsize == 3*ph.basic_data[:nat]) || 
#      throw(DimensionMismatch("dynmat dimensions are $(size(dynmat))"))
#
#   #compute exp(i 2\pi*K*R)
#   e_ikr = exp_ikr(ph.ws_rvecs, 2π*qpt)
#
#   k = 0
#   R = CartesianIndices((-2:0, -2:0))
#   mass = ph.basic_data[:mass]
#   @inbounds for j = 1:ph.basic_data[:nat], i = 1:j
#      k += 1
#      mfactor = 1.0 / sqrt( mass[i]*mass[j] )
#      dynmat[R .+ CartesianIndex(3*i,3*j)] .= 
#         _dyn_mat(ph.forceconst[k], ph.rvec_idx[k], e_ikr) .* mfactor
#   end
#end
#
#
#function _dyn_mat(ifc, rvec_idx, e_ikr)
#   dmat = zeros(eltype(e_ikr), 3, 3)
#   for (n, val) in enumerate(rvec_idx)
#      dmat .+= e_ikr[val] .* view(ifc, :, :, n)
#   end
#   return dmat
#end
#