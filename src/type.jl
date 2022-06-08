
struct ElecPhon{T, S, NA, NB}
   path::String
   atom_pos::SVector{NA, SVector{3,T}}
   orbit_pos::SVector{NB, SVector{3,T}}
   ws_rvecs_e::Vector{SVector{3,Int}}
   ws_rvecs_p::Vector{SVector{3,Int}}
   eph::Array{EphElement{S},3}
end

struct EphElement{S}
   re_idx::Vector{S}
   rp_idx::Vector{S}
end