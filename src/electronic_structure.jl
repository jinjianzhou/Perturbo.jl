
struct Hopping{T}
   t::Vector{Complex{T}}
   r_idx::Vector{Int}
end
Hopping(t::Vector{Complex{T}}, r_idx::Vector{Int}) where T = Hopping{T}(t, r_idx)

struct ElecHam{T}
   orbit_pos::Vector{ SVector{3,T} }
   ws_rvecs::Vector{ SVector{3,T} }
   hr::Vector{ Hopping{T} }
end

ElecHam(op::Vector{SVector{3,T}}, rvec::Vector{SVector{3,T}}, hr::Vector{Hopping{T}}) where T = ElecHam{T}(op, rvec, hr)

exp_ikr(rvecs::AbstractVector{<:SVector}, kpt::SVector) = cis.(Ref(kpt') .* rvecs)
_ham_elem(hr::Hopping, e_ikr::AbstractVector) = sum(e_ikr[val] * hr.t[n] for (n, val) in enumerate(hr.r_idx))

function _band!(Hk::AbstractMatrix, el::ElecHam, kpt::SVector)
   hamiltonian!(Hk, el, kpt)
   eigvals!( Hermitian(Hk) )
end

function hamiltonian!(H::AbstractMatrix, el::ElecHam, kpt::SVector)
   nband = LinearAlgebra.checksquare(H)
   nband == length(el.orbit_pos) || throw(DimensionMismatch("H: $(size(H)) vs $(length(el.orbit_pos))"))
   #compute exp(i 2\pi*K*R)
   e_ikr = exp_ikr(el.ws_rvecs, 2π*kpt)

   # compute upper triangle of ham
   k = 0
   @inbounds for j = 1:nband, i = 1:j
      H[i,j] = _ham_elem(el.hr[ k+=1 ], e_ikr)
   end
end

function bands(el::ElecHam, kpt::SVector)
   nband = length(el.orbit_pos)
   Hk = zeros(ComplexF64, nband, nband)
   _band!(Hk, el, kpt)
end

bands(el::ElecHam{T}, kpt::AbstractVector{T}) where T = bands(el, SVector{3,T}(kpt))

function bands(el::ElecHam{T}, kpts::AVec{<:SVector}) where T
   nband = length(el.orbit_pos)
   Hk = zeros(Complex{T}, nband, nband)

   eigs = zeros(T, nband, length(kpts))

   for (i, k) in enumerate(kpts)
      eigs[:,i] = _band!(Hk, el, k)
   end
   return eigs
end

function bands(el::ElecHam{T}, kpts::AbstractArray{T,2}) where T
   size(kpts, 1) != 3 && throw(error("1-dimension of kpts should be 3"))
   bands(el, vec(reinterpret(SVector{3,T}, kpts)))
end