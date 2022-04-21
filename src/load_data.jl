using HDF5, StaticArrays

ElecHam(fn::String) = ElecHam( (h5open(fn,"r") do fid 
   bdata = load_basic_data(fid)
   #
   w_center = bdata[:wannier_center_cryst]
   lattice = bdata[:at]
   rdim = bdata[:kc_dim]
   num_wann = bdata[:num_wann]
   #
   hopping = load_electron_wannier(fid, num_wann)

   wscell = WSCell(SVector{3}(rdim), SMatrix{3,3}(lattice))
   orb_pos = vec( reinterpret(SVector{3,eltype(w_center)}, w_center) )
   ws_rvecs, rvec_idx = wiger_seitz_cell(wscell, orb_pos)

   hr = [Hopping(hopping[i], rvec_idx[i]) for i in eachindex(hopping, rvec_idx)]

   return orb_pos, ws_
end)... )

function load_basic_data(fid)
   gid = fid["basic_data"]
   var_names = ("alat", "at", "bg", "epsil", "kc_dim", "mass", "nat", "nsym", 
      "num_wann", "polar_alpha", "qc_dim", "spinor", "symop", "system_2d", "tau", 
      "volume", "wannier_center", "wannier_center_cryst", "zstar")
   
   #read and store basic data in a NamedTuple
   (; zip(Symbol.(var_names), map(x->read(gid, x), var_names))...)
end
   
function load_electron_wannier(fid, num_wann)
   #only stored upper triangle part: for j = 1:nbnd, i = 1:j ==> (i,j)
   nelem = (num_wann * (num_wann+1) )>>1
   hopping = Vector{ Vector{ComplexF64} }(undef, nelem)
   for k in 1:nelem
      hop_r = read(fid, "electron_wannier/hopping_r" * string(k) )
      hop_i = read(fid, "electron_wannier/hopping_i" * string(k) )
      hopping[k] = complex.(hop_r, hop_i)
   end
   return hopping
end