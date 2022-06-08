
function ElecPhon(
   fn::String,
   ::Val{NA},
   ::Val{NB},
   atom_pos::AbstractVector{SVector{3,T}},
   orb_pos::AbstractVector{SVector{3,T}},
   ws_rvec_e,
   ws_rvec_p,
   eph::AbstractArray{EphElement{S},3},
) where {T,S,NA,NB}
   ElecPhon{T,S,NA,NB}(
      fn,
      SVector{NA}(atom_pos),
      SVector{NB}(orb_pos),
      ws_rvec_e,
      ws_rvec_p,
      eph,
   )
end

function ElecPhon(fn::String)
   nb, na, orb_pos, atom_pos, ws_rvec_e, ws_rvec_p, eph = h5open(fn, "r") do fid
      bdata = load_basic_data(fid)

      w_center = bdata[:wannier_center_cryst]
      lattice = bdata[:at]
      recip_latt = bdata[:bg]
      kc_dim = bdata[:kc_dim]
      qc_dim = bdata[:qc_dim]
      atom_pos_cart = bdata[:tau]
      atom_pos_crys = recip_latt * atom_pos_cart

      orb_pos = reinterpret(SVector{3,Float64}, w_center) |> vec |> collect
      atom_pos = reinterpret(SVector{3,Float64}, atom_pos_crys) |> vec |> collect
      nband = length(orb_pos)
      natom = length(atom_pos)

      wscell = WSCell(SVector{3,Int}(kc_dim), SMatrix{3,3}(lattice))
      delta_tau = [j - i for j in orb_pos for i in orb_pos]
      ws_rvecs_el, rvec_idx_el = wiger_seitz_cell(wscell, delta_tau)
      rvec_idx_el = reshape(rvec_idx_el, (nband, nband))

      wscell = WSCell(SVector{3,Int}(qc_dim), SMatrix{3,3}(lattice))
      delta_tau = [j - i for j in atom_pos for i in orb_pos]
      ws_rvecs_ph, rvec_idx_ph = wiger_seitz_cell(wscell, delta_tau)
      rvec_idx_ph = reshape(rvec_idx_ph, (nband, natom))

      eph = Array{EphElement{Int},3}(undef, (nband, nband, natom))
      for i in CartesianIndices(eph)
         iw, jw, ia = i.I
         eph[i] = EphElement(rvec_idx_el[iw, jw], rvec_idx_ph[iw, ia])
      end
      return nband, natom, orb_pos, atom_pos, ws_rvecs_el, ws_rvecs_ph, eph
   end

   ElecPhon(fn, Val(na), Val(nb), atom_pos, orb_pos, ws_rvec_e, ws_rvec_p, eph)
end