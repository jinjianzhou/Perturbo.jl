using HDF5

ElecHam(fn::String) = ElecHam( (h5open(fn,"r") do fid 
   bdata = load_basic_data(fid)
   #
   num_wann = bdata[:num_wann]
   w_center = bdata[:wannier_center_cryst]
   lattice = bdata[:at]
   rdim = bdata[:kc_dim]
   #
   hopping = load_electron_wannier(fid, num_wann)
   ws_rvecs, rvec_idx = WignerSeitzCell.wigner_seitz_cell(rdim, lattice, w_center)
   return bdata, ws_rvecs, rvec_idx, hopping
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