const ws_search_range = 3
const ws_sr = ws_search_range #a short name for ws_search_range

struct WSCell{T}
   WS_size::SVector{3,Int}
   latt_vectors::SMatrix{3,3,T,9}
   rvecs::Vector{ Vector{SVector{3,T}} }
   cutoff::T

   function WSCell{T}(wsize::SVector{3,Int}, latvec::SMatrix{3,3,T,9}) where {T<:Real}
      get_length = vec -> norm( latvec * vec)
      cutoff = set_cutoff(get_length, wsize)
      nr1, nr2, nr3 = wsize

      all_vecs = Vector{ Vector{SVector{3,T}} }(undef, prod(wsize))
      for (i, v) in enumerate( Iterators.product(0:nr1-1, 0:nr2-1, 0:nr3-1) )
         all_vecs[i] = valid_cells(get_length, SVector{3}(v), wsize, cutoff)
      end
      return new(wsize, latvec, all_vecs, cutoff)
   end
end

function set_cutoff(get_length::Function, rdim; lcutoff = false)
   ndim = lcutoff ? rdim : (rdim .÷ 2 .+ 1)
   max_len = 0.0
   for vec in Iterators.product( (-1, 1), (-1, 1), (-1, 1) )
      dist = get_length( SVector{3}(float.(vec .* ndim)) )
      dist > max_len && (max_len = dist)
   end
   return max_len
end

function valid_cells(get_length::Function, v0::SVector, rdim::SVector, cutoff::T) where T
   [vec for i in Iterators.product(-ws_sr:ws_sr, -ws_sr:ws_sr, -ws_sr:ws_sr) 
        for vec = Ref(SVector{3,T}(i .* rdim .+ v0)) if get_length(vec) < cutoff]
end

function wiger_seitz_cell(ws::WSCell{T}, tau::Vector{SVector{3,T}}) where T
   get_length = vec -> norm( ws.latvec * vec)
   #
   delta_tau = [ (tau[j] - tau[i]) for j in 1:length(tau)  for i in 1:j ]

   ws_vecs = Vector{ Vector{SVector{3,T}} }(undef, length(ws.rvecs))
   idx_tmp = Vector{ Vector{Int} }(undef, length(delta_tau))
   r_idx = [ Vector{Int}() for _ in 1:length(delta_tau) ]

   pre_len = 0
   for (vid, vecs) in enumerate(ws.rvecs)
      labels = zeros(Int, length(vecs))
      dist = zeros(T, length(vecs))

      for (it, dt) in enumerate(delta_tau)
         dist .= get_length.( Ref(dt) .+ vecs )
         mindist = minimum(dist) + 1.0E-6  #a small threshold
         idx = findall(<(mindist), dist)
         
         labels[idx] .+= 1
         idx_tmp[it] = idx
      end

      #collect vectors
      sel_idx = findall(labels .> 0)
      #re-label
      labels[sel_idx] = 1:length(sel_idx)
      # store relavent vectors
      ws_vecs[vid] = vecs[sel_idx]

      #update index of idx_tmp
      for it in eachindex(delta_tau)
         #map the old index to the re-labeled index
         idx_tmp[it] .= pre_len .+ labels[ idx_tmp[it] ]
         append!(r_idx[it], idx_tmp[it])
      end
      #
      pre_len += length(sel_idx)
   end
   return [v for vecs in ws_vecs for v in vecs], r_idx
end

#module WignerSeitzCell
#   export wigner_seitz_cell
#   #
#   using LinearAlgebra: norm
#   using StaticArrays
#
#   
#   # for performance, please use SMatrix{3,3} for at.
#   #return function to get length of a vector 
#   length_func(at) = vec -> norm( at * vec )
#
#   function wigner_seitz_cell(rdim, lattice::AbstractMatrix, tau::AbstractMatrix)
#      get_length = length_func( SMatrix{3,3}(lattice) )
#      cutoff = set_cutoff(get_length, rdim)
#
#      ntmp, ntau = size(tau)
#      ntmp != 3 && throw(DimensionMismatch("tau dimensions are $(size(tau))"))
#
#      delta_tau = [SVector{3}(view(tau,:,j) .- view(tau,:,i)) for j in 1:ntau for i in 1:j]
#
#      ws_rvecs = Vector{ SVector{3, Float64} }()
#      rvec_idx_tmp = Vector{ Vector{Int} }(undef, length(delta_tau))
#      rvec_idx = Vector{ Vector{Int} }(undef, length(delta_tau))
#
#      (nr1, nr2, nr3) = Tuple(rdim)
#      for rv in Iterators.product(0:nr1-1, 0:nr2-1, 0:nr3-1)
#         vecs = valid_cells(get_length, rv, rdim, cutoff)
#
#         vecs_label = zeros(Int, length(vecs))
#         for (it, dt) in enumerate(delta_tau)
#            dist = [get_length( v + dt) for v in vecs]
#            mindist = minimum(dist)
#            idx = findall(x -> ((x-mindist) < 1.0E-6), dist)
#            #
#            vecs_label[idx] .+= 1
#            rvec_idx_tmp[it] = idx
#         end
#
#         #collect vectors
#         col_idx = findall(x->(x > 0), vecs_label)
#         #re-label
#         vecs_label[col_idx] = 1:length(col_idx)
#         #
#         prev_len = length(ws_rvecs)
#         append!(ws_rvecs, vecs[col_idx])
#
#         #update rvec_idx
#         for it in eachindex(delta_tau)
#            rvec_idx_tmp[it] .= prev_len .+ vecs_label[ rvec_idx_tmp[it] ]
#            try 
#               append!(rvec_idx[it], rvec_idx_tmp[it])
#            catch error
#               error isa UndefRefError ? rvec_idx[it]=rvec_idx_tmp[it] : rethrow(error)
#            end
#         end
#      end
#      return ws_rvecs, rvec_idx
#   end
#
#end #module 
#   #lattice:  'at' in QE
#   function wigner_seitz_cell(ndim, lattice, tau_a, tau_b)
#      #
#      get_length = length_func( SMatrix{3,3}(lattice) )
#      cutoff = set_cutoff(get_length, ndim)
#      
#      r_vectors = rvec_image(get_length, ndim, cutoff)
#   
#      wscells = Vector{ SVector{3, Float64} }()
#      #
#      for vecs in r_vectors
#         #
#         dist = zeros( length(vecs) )
#         for (i, v) in enumerate(vecs)
#            dist[i] = get_length( SVector{3}(float.(v .+ tau_b .- tau_a)) )
#         end
#         #
#         mindist = minimum(dist)
#         idx = findall(x -> ( (x-mindist) < 1.0E-6), dist)
#         #
#         append!(wscells, vecs[idx])
#      end
#      #
#      return wscells
#   end #wiger_seitz_cell

#   let rvecs = Vector{Vector{ SVector{3,Float64} }}(), ndim = (0, 0, 0)
#      #
#      global function rvec_image(get_length::Function, rdim, cutoff)
#         Tuple(rdim) != ndim && rvec_image!(get_length, rdim, cutoff, rvecs)
#         return copy(rvecs)
#      end
#      #
#      function rvec_image!(get_length::Function, rdim, cutoff, rvecs)
#         empty!(rvecs)
#         (nr1, nr2, nr3) = ndim = Tuple(rdim)
#      
#         for i in Iterators.product(0:nr1-1, 0:nr2-1, 0:nr3-1)
#            v0 = SVector(float.(i))
#            push!(rvecs, valid_cells(get_length, v0, rdim, cutoff))
#         end
#         return nothing
#      end
#   end #let block
