
"""
Evaluate integral weight for the four coner of an tetra.

The weight are based on appendix B: since delta{e} = d step_func{e} / de
and appendix B give the formula to evalute weight for step_func{e_fermi},
replace e_fermi with e, the derivative of e is the weihgt we needed.
ps: the weights evaluate here are without V_T / V_G.

Refs: Peter Blochl, "Improved tetrahedron method or Brillouin-zone
integrations" PRB 49 16 333 (1994). 

   Check also the Fortran version: weight_dos.f90
"""
function tetra_int!(
   enk_t::AbstractVector{T},
   fnk_t::AbstractVector{S},
   e_dos::AbstractVector{T},
   t_sum::AbstractVector{S},
) where {T,S}
   @assert length(fnk_t) == length(enk_t) == 4
   "both enk_t and fnk_t should have the length of 4!"
   @assert length(e_dos) == length(t_sum)
   "e_dos and t_sum should have the same length!"

   # get sorted enk_t and their corresponding fnk_t
   idx = SVector{4}(sortperm(enk_t))
   e_sorted = enk_t[idx]
   f_sorted = fnk_t[idx]

   for ie in eachindex(e_dos, t_sum)
      # update tsum with contribution from the current tetrahedron
      @inbounds t_sum[ie] += _weight_dos_sorted(e_sorted, f_sorted, e_dos[ie])
   end
   return nothing
end

# internal function, do not export or used elsewhere.
# Assume e is in sorted order:  e[1] <= e[2] <= e[3] <= e[4]
@inline function _weight_dos_sorted(e::SVector{4,T}, f::SVector{4,S}, e_dos::T) where {T,S}
   # evaluate weight of the four corners, adopted from weight_dos.f90
   if e_dos < e[1]
      return zero(S)
   elseif e_dos < e[2] || (e_dos == e[2] == e[3] == e[4] && e[2] > e[1])
      et1 = e_dos - e[1]
      et2 = e[2] - e[1]
      et3 = e[3] - e[1]
      et4 = e[4] - e[1]
      factor = et1 * et1 / (et2 * et3 * et4)
      w2 = factor * et1 / et2
      w3 = factor * et1 / et3
      w4 = factor * et1 / et4
      w1 = 3.0 * factor - w2 - w3 - w4
   elseif e_dos < e[3]
      et1 = e_dos - e[1]
      et2 = e_dos - e[2]
      et3 = e[3] - e_dos
      et4 = e[4] - e_dos
      e31 = e[3] - e[1]
      e32 = e[3] - e[2]
      e41 = e[4] - e[1]
      e42 = e[4] - e[2]
      dc1 = 0.5 * et1 / (e41 * e31)
      c1 = 0.5 * et1 * dc1
      c2 = 0.25 / (e41 * e32 * e31)
      dc2 = c2 * (et2 * et3 + et1 * et3 - et1 * et2)
      c2 = c2 * et1 * et2 * et3
      c3 = 0.25 / (e42 * e32 * e41)
      dc3 = c3 * (2.0 * et2 * et4 - et2 * et2)
      c3 = c3 * (et2 * et2 * et4)
      w1 = dc1 + (dc1 + dc2) * et3 / e31 + (dc1 + dc2 + dc3) * et4 / e41
      w1 = w1 - (c1 + c2) / e31 - (c1 + c2 + c3) / e41
      w2 = dc1 + dc2 + dc3 + (dc2 + dc3) * et3 / e32 + dc3 * et4 / e42
      w2 = w2 - (c2 + c3) / e32 - c3 / e42
      w3 = (dc1 + dc2) * et1 / e31 + (dc2 + dc3) * et2 / e32
      w3 = w3 + (c1 + c2) / e31 + (c2 + c3) / e32
      w4 = (dc1 + dc2 + dc3) * et1 / e41 + dc3 * et2 / e42
      w4 = w4 + (c1 + c2 + c3) / e41 + c3 / e42
   elseif e_dos < e[4]
      et1 = e[4] - e[1]
      et2 = e[4] - e[2]
      et3 = e[4] - e[3]
      et4 = e[4] - e_dos
      factor = et4 * et4 / (et1 * et2 * et3)
      w1 = factor * et4 / et1
      w2 = factor * et4 / et2
      w3 = factor * et4 / et3
      w4 = 3.0 * factor - w1 - w2 - w3
   else
      return zero(S)
   end
   return convert(S, w1 * f[1] + w2 * f[2] + w3 * f[3] + w4 * f[4])
end