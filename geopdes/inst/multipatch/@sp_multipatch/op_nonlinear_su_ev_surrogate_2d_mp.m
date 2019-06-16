% OP_NONLINEAR_SU_EV_SURROGATE_2D_MP: assemble the matrix A = 
%   mat = op_nonlinear_su_ev_surrogate_2d_mp (spu, msh, Stress, DStress, u_old, M, q, [patches]);
%
% INPUT:
%    
%   spu:        object representing the space of trial and test functions (see sp_multipatch)
%   msh:        object that defines the domain partition and the quadrature rule (see msh_multipatch)
%   lambda, mu: function handles to compute the Lame' coefficients
%   M:          skip parameter
%   q:          surrogate interpolation degree (1 or 3)
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat:    assembled matrix
% 
% Copyright (C) 2015 Rafael Vazquez
% Copyright (C) 2019 Daniel Drzisga
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function A = op_nonlinear_su_ev_surrogate_2d_mp (spu, msh, Stress, DStress, u_old, M, q, patch_list)

  if (nargin < 8)
    patch_list = 1:msh.npatch;
  end

  if ((spu.npatch ~= spu.npatch) || (spu.npatch ~= msh.npatch))
    error ('op_nonlinear_su_ev_surrogate_2d_mp: the number of patches does not coincide')
  end
  
  ncounter = 0;
  for iptc = patch_list

    u_loc = u_old(spu.gnum{iptc});
    
    [rs, cs, vs] = op_nonlinear_su_ev_surrogate_2d (spu.sp_patch{iptc}, msh.msh_patch{iptc}, Stress, DStress, u_old, M, q);
    rows(ncounter+(1:numel (rs))) = spu.gnum{iptc}(rs);
    cols(ncounter+(1:numel (rs))) = spu.gnum{iptc}(cs);
    vals(ncounter+(1:numel (rs))) = vs;
    ncounter = ncounter + numel (rs);
  end

  A = sparse (rows, cols, vals, spu.ndof, spu.ndof);
  clear rows cols vals rs cs vs

end