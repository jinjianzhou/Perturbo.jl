using StaticArrays, LinearAlgebra

function bands(el::ElecHam, kpt::SVector)
   nband = el.basic_data[:num_wann]
   Hk = zeros(ComplexF64, nband, nband)
   _band!(Hk, el, kpt)
end

bands(el::ElecHam, kpt::AbstractVector{T}) where T = bands(el, SVector{3,T}(kpt))

function bands(el::ElecHam, kpts::AVec{<:SVector})
   nband = el.basic_data[:num_wann]
   Hk = zeros(ComplexF64, nband, nband)

   eigs = zeros(Float64, nband, length(kpts))

   for (i, k) in enumerate(kpts)
      eigs[:,i] = _band!(Hk, el, k)
   end
   return eigs
end

function bands(el::ElecHam, kpts::AbstractArray{T,2}) where T
   size(kpts, 1) != 3 && throw(error("1-dimension of kpts should be 3"))
   bands(el, vec(reinterpret(SVector{3,T}, kpts)))
end

function _band!(Hk::AbstractMatrix, el::ElecHam, kpt::SVector)
   hamiltonian!(Hk, el, kpt)
   eigvals( Hermitian(Hk) )
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