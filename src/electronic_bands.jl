using StaticArrays
using StructArrays
using LinearAlgebra
using HDF5
#using TimerOutputs

struct Hopping{T<:Number}
   ws_cell::AbstractVector{SVector{3,T}}
   hop::AbstractVector{Complex{T}}
end

struct ElecHam{T<:Number}
   nwan::Integer
   hoppings::AbstractVector{Hopping{T}}
end

ElecHam(fn::String) = ElecHam((h5open(fn,"r") do fid load_electron_ham(fid) end)...)

function load_electron_ham(fid)
   nbnd = read(fid, "basic_data/num_wann")
   tau = read(fid, "basic_data/wannier_center_cryst")
   at = read(fid, "basic_data/at")
   ndim = Tuple( read(fid, "basic_data/kc_dim") )

   hoppings = Hopping{eltype(at)}[]
   k = 0
   @inbounds for j = 1:nbnd, i = 1:j
      k += 1
      dname = "electron_wannier/hopping_r" * string(k)
      rep = read(fid, dname)
      dname = "electron_wannier/hopping_i" * string(k)
      imp = read(fid, dname)

      #complex.(rep, imp)
      hop = StructArray{Complex{eltype(at)}}( (rep, imp) )
      ws_cell = WignerSeitzCell.wigner_seitz_cell(ndim, at, tau[:,i], tau[:,j])
      
      push!(hoppings, Hopping(ws_cell, hop))
   end
   return nbnd, hoppings
end

function bands(el::ElecHam{T}, kpts::AbstractArray{T,2}) where T
   size(kpts, 1) != 3 && throw(error("inconsistent dimension"))
   eigs = zeros( T, el.nwan, size(kpts,2) )
   Hk = zeros(Complex{T}, el.nwan, el.nwan)
   for (i, k) in enumerate(Iterators.partition(kpts,3))
      hamiltonian!(Hk, el.hoppings, SVector{3}(k))
      eigs[:,i] = eigvals( Hermitian(Hk) )
   end
   return eigs
end

function bands(el::ElecHam{T}, kpts::AbstractVector{<:SVector{3,T}}) where T
   eigs = zeros(T, el.nwan, length(kpts))
   Hk = zeros(Complex{T}, el.nwan, el.nwan)
   #reset_timer!()
   for (i, k) in enumerate(kpts)
      hamiltonian!(Hk, el.hoppings, k)
      eigs[:,i] = eigvals( Hermitian(Hk) )
   end
   #print_timer()
   return eigs
end

function hamiltonian!(H::Matrix{Complex{T}}, hoppings::AbstractVector{<:Hopping{T}}, kpt::SVector{3,T}) where T
   m = LinearAlgebra.checksquare(H)
   length(hoppings) != ( (m*(m+1))>>1 ) && throw(error("inconsistent dimension"))
   k = 0
   @inbounds for j = 1:m, i = 1:j
      H[i,j] = _ham_elem(hoppings[ k += 1 ], kpt)
   end
   return nothing
end

_ham_elem(hopping::Hopping{S}, kpt::SVector{3,S}) where S = _ham_elem(hopping.ws_cell, hopping.hop, kpt)

function _ham_elem(ws_cell::AbstractVector{<:SVector{L,S}}, hop::AbstractVector{Complex{S}}, kpt::SVector{L,S}) where {L, S<:Number}
   ham = zero(Complex{S})
   for i in eachindex(ws_cell, hop)
      ham += cis( dot(kpt, ws_cell[i]) ) * hop[i]
   end
   return ham
end
