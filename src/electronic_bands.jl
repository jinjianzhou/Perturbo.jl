using StaticArrays
using LinearAlgebra
using HDF5

struct Hopping{S<:AbstractVector, T<:AbstractVector}
   ws_cell::S
   hop::T
end

struct ElecHam{S<:Signed, T<:AbstractMatrix, H<:AbstractVector{<:Hopping}}
   nwan::S
   real_latt::T
   recip_latt::T
   hoppings::H
end

ElecHam(fn::String) = ElecHam( (h5open(fn,"r") do fid load_electron_ham(fid) end)... )

function load_electron_ham(fid)
   gid = fid["basic_data"]
   var_names = ["num_wann", "kc_dim", "at", "bg", "wannier_center_cryst"]
   nbnd, ndim, at, bg, tau = map(x->read(gid, x), var_names)
   
   T  = eltype(at)
   Tw = Vector{ SVector{3,T} }
   Th = Vector{ Complex{T} }
   hoppings = Vector{Hopping{Tw, Th}}(undef, (nbnd*(nbnd+1))>>1 )
   
   k = 0
   @inbounds for j = 1:nbnd, i = 1:j
      k += 1
      rep = read(fid, "electron_wannier/hopping_r" * string(k) )
      imp = read(fid, "electron_wannier/hopping_i" * string(k) )

      hop = complex.(rep, imp)
      ws_cell = WignerSeitzCell.wigner_seitz_cell(ndim, at, tau[:,i], tau[:,j])
      
      hoppings[k] = Hopping(ws_cell, hop)
   end

   return Int(nbnd), SMatrix{3,3}(at), SMatrix{3,3}(bg), hoppings
end

function bands(el::ElecHam, kpts::AbstractArray{T,2}) where T
   size(kpts, 1) != 3 && throw(error("inconsistent dimension"))
   eigs = zeros(T, el.nwan, size(kpts,2))
   Hk = zeros(Complex{T}, el.nwan, el.nwan)

   for (i, k) in enumerate(Iterators.partition(kpts,3))
      hamiltonian!(Hk, el.hoppings, SVector{3}(k))
      eigs[:,i] = eigvals( Hermitian(Hk) )
   end

   return eigs
end

function bands(el::ElecHam, kpts::AbstractVector{<:SVector{3,T}}) where T
   eigs = zeros(T, el.nwan, length(kpts))
   Hk = zeros(Complex{T}, el.nwan, el.nwan)

   for (i, k) in enumerate(kpts)
      hamiltonian!(Hk, el.hoppings, k)
      eigs[:,i] = eigvals( Hermitian(Hk) )
   end

   return eigs
end

function hamiltonian!(H::AbstractMatrix, hoppings::AbstractVector{<:Hopping}, kpt::SVector)
   m = LinearAlgebra.checksquare(H)
   length(hoppings) != ( (m*(m+1))>>1 ) && throw(error("inconsistent dimension"))

   k = 0
   @inbounds for j = 1:m, i = 1:j
      H[i,j] = _ham_elem(hoppings[ k += 1 ], kpt)
   end

   return nothing
end

_ham_elem(hopping::Hopping, kpt::SVector) = _ham_elem(hopping.ws_cell, hopping.hop, kpt)

# @btime gives: 8.264 μs (0 allocations: 0 bytes)
function _ham_elem(ws_cell::AbstractVector{<:SVector{L,S}}, hop::AbstractVector{Complex{S}}, kpt::SVector{L,S}) where {L, S}
   ham = zero( eltype(hop) )
   for i in eachindex(ws_cell, hop)
      ham += cis( dot(kpt, ws_cell[i]) ) * hop[i]
   end
   return ham
end