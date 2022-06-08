
function elph_fft(
   ep::ElecPHon{T,NA,NB},
   kgrid::NTuple{3,Int},
   qgrid::NTuple{3,Int},
   kpt::SVector{3,T} = zero(SVector{3,T}),
   qpt::SVector{3,T} = zero(SVector{3,T}),
)
   #check size
   min_size = grid_size(ep.ws_rvecs_e)
   any(kgrid .< min_size) && error("kgrid should be larger than $(min_size)!")
   min_size = grid_size(ep.ws_rvecs_p)
   any(qgrid .< min_size) && error("qgrid should be larger than $(min_size)!")

   ep_fft = zeros(Complex{T}, (3, kgrid..., qgrid...))
   #Array of Array to hold all the results
   ephwan = Array{Array{Complex{T},7},3}(undef, NB, NB, NA)

   e_ikr = exp_ikr(ep.ws_rvecs_e, 2π * kpt)
   e_iqr = exp_ikr(ep.ws_rvecs_p, 2π * qpt)

   e_pos = [CartesianIndex(get_fft_pos(kgrid, rv)) for rv in ep.ws_rvecs_e]
   p_pos = [CartesianIndex(get_fft_pos(qgrid, rv)) for rv in ep.ws_rvecs_p]

   h5open(ep.path, "r") do fid
      for i in CartesianIndices((NB, NB, NA))
         #read data from HDF5 file, TODO:
         g = load_eph_data(fid, i.I)
         #map to ep_fft array to do FFT 
         for (np, pval) in Iterators.enumerate(ep.eph[i].rp_idx),
            (ne, eval) in Iterators.enumerate(ep.eph[i].re_idx),
            j in 1:3

            ep_fft[j, e_pos[eval], p_pos[pval]] = g[j, ne, np] * e_ikr[eval] * e_iqr[pval]
         end

         #perform FFT
         p = plan_fft(ep_fft, 2:7)
         ephwan[i] = similar(ep_fft)
         mul!(ephwan[i], p, ep_fft)
      end
      return ephwan
   end
end

function load_eph_data(fid, ep_id::NTuple{3,Int})
   id = join(string.(reverse(ep_id)), "_")
   eph_r = read(fid, "eph_matrix_wannier/ep_hop_r_" * id)
   eph_i = read(fid, "eph_matrix_wannier/ep_hop_i_" * id)
   complex.(eph_r, eph_i)
end