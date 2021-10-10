
module WignerSeitzCell
   export wigner_seitz_cell
   #
   using LinearAlgebra: norm
   using StaticArrays

   const ws_search_range = 3
   
   # for performance, please use SMatrix{3,3} for at.
   #return function to get length of a vector 
   length_func(at) = vec -> norm( at * vec )

   function wigner_seitz_cell(rdim, lattice::AbstractMatrix, tau::AbstractMatrix)
      get_length = length_func( SMatrix{3,3}(lattice) )
      cutoff = set_cutoff(get_length, rdim)

      ntmp, ntau = size(tau)
      ntmp != 3 && throw(DimensionMismatch("tau dimensions are $(size(tau))"))

      delta_tau = [SVector{3}(view(tau,:,j) .- view(tau,:,i)) for j in 1:ntau for i in 1:j]

      ws_rvecs = Vector{ SVector{3, Float64} }()
      rvec_idx_tmp = Vector{ Vector{Int} }(undef, length(delta_tau))
      rvec_idx = Vector{ Vector{Int} }(undef, length(delta_tau))

      (nr1, nr2, nr3) = Tuple(rdim)
      for rv in Iterators.product(0:nr1-1, 0:nr2-1, 0:nr3-1)
         vecs = valid_cells(get_length, rv, rdim, cutoff)

         vecs_label = zeros(Int, length(vecs))
         for (it, dt) in enumerate(delta_tau)
            dist = [get_length( v + dt) for v in vecs]
            mindist = minimum(dist)
            idx = findall(x -> ((x-mindist) < 1.0E-6), dist)
            #
            vecs_label[idx] .+= 1
            rvec_idx_tmp[it] = idx
         end

         #collect vectors
         col_idx = findall(x->(x > 0), vecs_label)
         #re-label
         vecs_label[col_idx] = 1:length(col_idx)
         #
         prev_len = length(ws_rvecs)
         append!(ws_rvecs, vecs[col_idx])

         #update rvec_idx
         for it in eachindex(delta_tau)
            rvec_idx_tmp[it] .= prev_len .+ vecs_label[ rvec_idx_tmp[it] ]
            try 
               append!(rvec_idx[it], rvec_idx_tmp[it])
            catch error
               error isa UndefRefError ? rvec_idx[it]=rvec_idx_tmp[it] : rethrow(error)
            end
         end
      end
      return ws_rvecs, rvec_idx
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

   function valid_cells(get_length::Function, v0, rdim, cutoff)
      vectors = Vector{SVector{3,Float64}}()
      nr = ws_search_range
      for i in Iterators.product(-nr:nr, -nr:nr, -nr:nr)
         vec = SVector{3}(float.(i .* rdim .+ v0))
         get_length(vec) < cutoff && push!(vectors, vec)
      end
      return vectors
   end
end #module 
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
