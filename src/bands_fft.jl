
function bands_fft(
   tb::ElecHam{T},
   fft_size::NTuple{3,Int},
   kpt::SVector{3,T} = zero(SVector{3,T}),
) where {T}
   hk_fft = hamk_fft(tb, fft_size, kpt)
   nbnd = length(tb.orbit_pos)
   diag_hamk(hk_fft, Val(nbnd))
end

function hamk_fft(tb::ElecHam{T}, fft_size::NTuple{3,Int}, kpt::SVector{3,T}) where {T}
   #get minimum fft_size
   min_size = grid_size(tb.ws_rvecs)
   any(fft_size .< min_size) && error("fft_size should be larger than $(min_size)!")
   nelem = length(tb.hr)

   hr_fft = zeros(Complex{T}, (fft_size..., nelem))
   #compute exp(i 2\pi*K*R)
   e_ikr = exp_ikr(tb.ws_rvecs, 2π * kpt)
   #map the H(R) data to an array for fftw transformation.
   for i in 1:nelem, (n, val) in Iterators.enumerate(tb.hr[i].r_idx)
      rvec = tb.ws_rvecs[val]
      fft_pos = get_fft_pos(fft_size, rvec)
      hr_fft[CartesianIndex(fft_pos), i] = tb.hr[i].t[n] * e_ikr[val]
   end

   #perform FFT 
   p = plan_fft(hr_fft, 1:3, flags = FFTW.PATIENT)
   hk_fft = similar(hr_fft)
   mul!(hk_fft, p, hr_fft)
   #
   return hk_fft
end

function diag_hamk(hk_fft::AbstractArray{Complex{T},4}, ::Val{D}) where {T,D}
   fft_size = size(hk_fft)[1:3]
   eigs = zeros(T, (D, fft_size...))
   Hk = zeros(Complex{T}, D, D)

   for ik in CartesianIndices(fft_size)
      n = 0
      for j in 1:D, i in 1:j
         Hk[i, j] = hk_fft[ik, n+=1]
      end
      H = D < 5 ? SMatrix{D,D}(Hermitian(Hk)) : Hk
      eigs[:, ik] = eigvals(Hermitian(H))
   end
   return eigs
end

@inline function grid_size(rvecs::AbstractVector{SVector{3,T}}) where {T}
   rsize = zeros(T, 3)
   for i in 1:3
      rmin, rmax = extrema(getindex.(rvecs, i))
      @inbounds rsize[i] = rmax - rmin + one(T)
   end
   return rsize
end

@inline exp_ikr(rvecs::AbstractVector{<:SVector}, kpt::SVector) = cis.(Ref(kpt') .* rvecs)

@inline function get_fft_pos(fft_size::NTuple{3,S}, rvec::SVector{3,S}) where {S<:Integer}
   return NTuple{3}((rvec .% fft_size .+ fft_size) .% fft_size .+ one(S))
end

function get_fft_kpts(fft_size::NTuple{3,Int})
   kpts = Array{SVector{3,Float64},3}(undef, fft_size)
   for i in CartesianIndices(fft_size)
      kpts[i] = pos2pts(fft_size, i.I)
   end
   return vec(kpts)
end