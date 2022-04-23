
function gen_kpts(get_length::Function, kpath::AbstractVector{SVector{3,T}}, num_kpts::Integer) where T
   npath = length(kpath)
   npath > 1 || throw(error("At least two k-points should be provided !"))

   len = [get_length(kpath[i+1]-kpath[i]) for i in 1:(npath-1)]
   deltak = len[1] / num_kpts
   nk = Int[ round(val/deltak) for val in len ]
   total_nk = sum(nk) + 1

   kpts = Vector{ SVector{3,T} }(undef, total_nk)
   kpts[1] = kpath[1]
   xloc = zeros(T, total_nk)
   
   prev = 1
   for i in 1:(npath-1)
      idx_range = prev+1 : prev+nk[i]
      kpts[idx_range] .= Iterators.drop( range(kpath[i], stop=kpath[i+1], length=nk[i]+1), 1)
      xloc[idx_range] .= Iterators.drop( range(zero(T), len[i], length=nk[i]+1), 1) .+ xloc[prev]
      prev = prev+nk[i]
   end
   return xloc, kpts
end