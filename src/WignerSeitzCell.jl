
module WignerSeitzCell
   export wigner_seitz_cell
   #
   using LinearAlgebra
   using StaticArrays

   const ws_search_range = 3
   
   #lattice:  'at' in QE
   function wigner_seitz_cell(ndim, lattice, tau_a, tau_b)
      #
      get_length = length_func( SMatrix{3,3}(lattice) )
      cutoff = set_cutoff(get_length, ndim)
      
      r_vectors = rvec_image(get_length, ndim, cutoff)
   
      wscells = Vector{ SVector{3, Float64} }()
      #
      for vecs in r_vectors
         #
         dist = zeros( length(vecs) )
         for (i, v) in enumerate(vecs)
            dist[i] = get_length( SVector{3}(float.(v .+ tau_b .- tau_a)) )
         end
         #
         dist .= dist .- minimum(dist)
         idx = findall(x -> (x < 1.0E-6), dist)
         #
         append!(wscells, vecs[idx])
      end
      #
      return wscells
   end #wiger_seitz_cell
   
   # for performance, please use SMatrix{3,3} for at.
   #return function to get length of a vector 
   length_func(at) = vec -> norm( at * vec )

   function set_cutoff(get_length::Function, rdim; lcutoff = false)
      ndim = lcutoff ? rdim : (rdim .÷ 2 .+ 1)
     
      max_len = 0.0
      for vec in Iterators.product( (-1, 1), (-1, 1), (-1, 1) )
         dist = get_length( SVector{3}(float.(vec .* ndim)) )
         dist > max_len && (max_len = dist)
      end
      return max_len
   end

   function valid_cells(get_length::Function, v0, ws_dim, cutoff)
      vectors = Vector{SVector{3,Float64}}()
      nr = ws_search_range
      for i in Iterators.product(-nr:nr, -nr:nr, -nr:nr)
         vec = SVector{3}(float.(i .* ws_dim .+ v0))
         get_length(vec) < cutoff && push!(vectors, vec)
      end
      return vectors
   end

   let rvecs = Vector{Vector{ SVector{3,Float64} }}(), ndim = (0, 0, 0)
      #
      global function rvec_image(get_length::Function, rdim, cutoff)
         Tuple(rdim) != ndim && rvec_image!(get_length, rdim, cutoff, rvecs)
         return copy(rvecs)
      end
      #
      function rvec_image!(get_length::Function, rdim, cutoff, rvecs)
         empty!(rvecs)
         (nr1, nr2, nr3) = ndim = Tuple(rdim)
      
         for i in Iterators.product(0:nr1-1, 0:nr2-1, 0:nr3-1)
            v0 = SVector(float.(i))
            push!(rvecs, valid_cells(get_length, v0, rdim, cutoff))
         end
         return nothing
      end
   end #let block
end #module 
