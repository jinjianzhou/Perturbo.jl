const ws_search_range = 3
const ws_sr = ws_search_range #a short name for ws_search_range

struct WSCell{S<:Integer,T<:Real}
   ws_size::NTuple{3,S}
   rvecs::Vector{Vector{SVector{3,S}}}
   latt_vectors::SMatrix{3,3,T,9}
   cutoff::T
end

function WSCell(wsize::NTuple{3,S}, latvec::SMatrix{3,3,T,9}) where {S<:Integer,T<:Real}
   get_length = vec -> norm(latvec * vec)
   cutoff = set_cutoff(get_length, wsize)
   nr1, nr2, nr3 = wsize

   all_vecs = Vector{Vector{SVector{3,Int}}}(undef, prod(wsize))
   for (i, v) in enumerate(Iterators.product(0:nr1-1, 0:nr2-1, 0:nr3-1))
      all_vecs[i] = valid_cells(get_length, v, wsize, cutoff)
   end
   return WSCell{S,T}(wsize, all_vecs, latvec, cutoff)
end

WSCell(wsize::SVector{3}, latvec) = WSCell(NTuple{3}(wsize), latvec)

function set_cutoff(get_length::Function, rdim; lcutoff = false)
   ndim = lcutoff ? rdim : (rdim .÷ 2 .+ 1)
   max_len = 0.0
   for vec in Iterators.product((-1, 1), (-1, 1), (-1, 1))
      dist = get_length(SVector{3,Float64}(vec .* ndim))
      dist > max_len && (max_len = dist)
   end
   return max_len
end

function valid_cells(
   get_length::Function,
   v0::NTuple{3,S},
   rdim::NTuple{3,S},
   cutoff::Real,
) where {S<:Integer}
   [
      vec for i in Iterators.product(-ws_sr:ws_sr, -ws_sr:ws_sr, -ws_sr:ws_sr) for
      vec = Ref(SVector{3,S}(i .* rdim .+ v0)) if get_length(vec) < cutoff
   ]
end

function wiger_seitz_cell(
   ws::WSCell{S,T},
   delta_tau::AVec{<:SVector{3}},
) where {S<:Integer,T<:Real}
   get_length = vec -> norm(ws.latt_vectors * vec)
   #
   #delta_tau = [(tau[j] - tau[i]) for j in 1:length(tau) for i in 1:j]

   ws_vecs = Vector{Vector{SVector{3,S}}}(undef, length(ws.rvecs))
   idx_tmp = Vector{Vector{S}}(undef, length(delta_tau))
   r_idx = [Vector{S}() for _ in 1:length(delta_tau)]

   pre_len = 0
   for (vid, vecs) in enumerate(ws.rvecs)
      labels = zeros(S, length(vecs))
      dist = zeros(T, length(vecs))

      for (it, dt) in enumerate(delta_tau)
         dist .= get_length.(Ref(dt) .+ vecs)
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
         idx_tmp[it] .= pre_len .+ labels[idx_tmp[it]]
         append!(r_idx[it], idx_tmp[it])
      end
      #
      pre_len += length(sel_idx)
   end
   return [v for vecs in ws_vecs for v in vecs], r_idx
end