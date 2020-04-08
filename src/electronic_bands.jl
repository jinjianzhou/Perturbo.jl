using StaticArrays
using LinearAlgebra
using HDF5
#for plot bandstructure
using RecipesBase

const Rydberg2eV = 13.605698066

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

#for plot bandstructure
@userplot Plotbands

@recipe function f(bs::Plotbands; nkpoints=30)
   tb, kpath = bs.args[1:2]
   tics, kloc, kpts = gen_kpts(kpath, tb.recip_latt, nkpoints)
   eigs = bands(tb, kpts) .* Rydberg2eV
   
   tics_label = length(bs.args) >= 3 ? bs.args[3] : string.(collect(eachindex(tics)))
   xticks --> (tics, tics_label)
   xlims --> (kloc[1], kloc[end])
   seriescolor --> :blue
   ygrid --> false
   legend --> false
   framestyle --> :box
   yguide --> "Energy (eV)"
   aspect_ratio --> 0.5

   #remove npoints from plotattributes, avoid unexpected results
   delete!(plotattributes, :nkpoints)

   kloc, eigs'
end

function gen_kpts(kpath, recip_latt, nk)
   size(kpath, 1) != 3 && throw(error("size of the first dimension is not 3"))
   kpts = [SVector{3}(kpath[:,1])]
   kloc = [0.0]
   npath = size(kpath, 2)
   ktic = zeros(npath)

   calc_len = WignerSeitzCell.length_func(recip_latt)
   @inbounds for i in 2:npath
      dk = SVector{3}(kpath[:,i]-kpath[:,i-1])
      dk_len = calc_len(dk)
      nkt = i == 2 ? nk : Int(round(kloc[i] * nk / kloc[2]))
      iter = Iterators.drop(range(0.,1., length=nkt), 1)
      ktic[i] = dk_len + ktic[i-1]
      append!(kloc, dk_len .* iter .+ ktic[i-1])
      append!(kpts, [dk] .* iter .+ [SVector{3}(kpath[:,i-1])] )
   end

   return ktic, kloc, kpts
end