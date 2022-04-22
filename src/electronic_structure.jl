
struct Hopping{T}
   t::Vector{Complex{T}}
   r_idx::Vector{Int}
end

struct ElecHam{T}
   orbit_pos::Vector{ SVector{3,T} }
   ws_rvecs::Vector{ SVector{3,T} }
   hr::Vector{ Hopping{T} }
end

ElecHam(fn::String) = ElecHam( (h5open(fn,"r") do fid 
   bdata = load_basic_data(fid)
   #
   w_center = bdata[:wannier_center_cryst]
   lattice = bdata[:at]
   rdim = bdata[:kc_dim]
   num_wann = bdata[:num_wann]
   #
   hopping = load_electron_wannier(fid, num_wann)

   wscell = WSCell(SVector{3,Int}(rdim), SMatrix{3,3}(lattice))
   orb_pos = vec(collect( reinterpret(SVector{3,eltype(w_center)}, w_center) ))
   ws_rvecs, rvec_idx = wiger_seitz_cell(wscell, orb_pos)

   hr = [Hopping(hopping[i], rvec_idx[i]) for i in eachindex(hopping, rvec_idx)]
   
   return orb_pos, ws_rvecs, hr
end)... )

exp_ikr(rvecs::AVec{<:SVector}, kpt::SVector) = cis.(Ref(kpt') .* rvecs)
_ham_elem(hr::Hopping, e_ikr::AVec) = sum(e_ikr[val] * hr.t[n] for (n, val) in enumerate(hr.r_idx))

function hamiltonian!(H::AbstractMatrix, el::ElecHam, kpt::SVector)
   nband = LinearAlgebra.checksquare(H)
   nband == length(el.orbit_pos) || throw(DimensionMismatch("H: $(size(H)) vs $(length(el.orbit_pos))"))
   #compute exp(i 2\pi*K*R)
   e_ikr = exp_ikr(el.ws_rvecs, 2π*kpt)

   # compute upper triangle of H
   k = 0
   @inbounds for j = 1:nband, i = 1:j
      H[i,j] = _ham_elem(el.hr[ k+=1 ], e_ikr)
   end
end

function _band!(Hk::AbstractMatrix, el::ElecHam, kpt::SVector)
   hamiltonian!(Hk, el, kpt)
   eigvals!( Hermitian(Hk) )
end

bands(el::ElecHam{T}, kpt::AVec{T}) where T = bands(el, SVector{3,T}(kpt))

function bands(el::ElecHam{T}, kpt::SVector{3,T}) where T
   nband = length(el.orbit_pos)
   Hk = zeros(Complex{T}, nband, nband)
   _band!(Hk, el, kpt)
end

function bands(el::ElecHam{T}, kpts::AbstractArray{SVector{3,T}}) where T
   nband = length(el.orbit_pos)
   Hk = zeros(Complex{T}, nband, nband)

   eigs = zeros(T, nband, length(kpts))

   for (i, k) in enumerate(kpts)
      eigs[:,i] = _band!(Hk, el, k)
   end
   return eigs
end

function bands(el::ElecHam{T}, kpts::AbstractArray{T,2}) where T
   size(kpts, 1) == 3 || throw(error("dim 1 of kpts should be 3"))
   bands(el, reinterpret(SVector{3,T}, kpts))
end