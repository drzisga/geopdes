% OP_U_V_ELEMENTSURROGATE_2D_MP: assemble the mass matrix M = [m(i,j)], m(i,j) = (mu u_j, v_i),
% in a multipatch domain using the element surrogate method.
%
%   mat = op_u_v_elementsurrogate_2d_mp (spu, msh, coeff, M, q, patch_list);
%
% INPUT:
%
%  spu:     object representing the space of test functions (see sp_multipatch)
%  msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%  coeff:   function handle to compute the reaction coefficient
%  M:       skip parameter
%  q:       surrogate interpolation degree (1 or 3)
%  patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%  mat:    assembled surrogate mass matrix
% 
% Copyright (C) 2015, 2016 Rafael Vazquez
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

function A = op_u_v_elementsurrogate_2d_mp (spu, msh, coeff, M, q, patch_list)

  if (nargin < 6)
    patch_list = 1:msh.npatch;
  end

  if ((spu.npatch ~= spu.npatch) || (spu.npatch ~= msh.npatch))
    error ('op_u_v_elementsurrogate_2d_mp: the number of patches does not coincide')
  end
  
  ncounter = 0;
  for iptc = patch_list
    [rs, cs, vs] = op_u_v_elementsurrogate_2d (spu.sp_patch{iptc}, msh.msh_patch{iptc}, coeff, M, q);
    rows(ncounter+(1:numel (rs))) = spu.gnum{iptc}(rs);
    cols(ncounter+(1:numel (rs))) = spu.gnum{iptc}(cs);

    if (~isempty (spu.dofs_ornt))
      vs = spu.dofs_ornt{iptc}(rs)' .* vs;
    end
    if (~isempty (spu.dofs_ornt))
      vs = vs .* spu.dofs_ornt{iptc}(cs)';
    end
    
    vals(ncounter+(1:numel (rs))) = vs;
    ncounter = ncounter + numel (rs);
  end

  A = sparse (rows, cols, vals, spu.ndof, spu.ndof);
  clear rows cols vals rs cs vs

end
