% OP_F_V_MP2: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i), in a multipatch domain.
%
%   rhs = op_f_v_mp (spv, msh, coeff, [patches]);
%
% INPUT:
%     
%   spv:     object representing the function space (see sp_multipatch)
%   msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%   coeff:   function handle to compute the source function
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   rhs: assembled right-hand side
% 
% Copyright (C) 2015 Rafael Vazquez
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

function rhs = op_f_v_mp2 (space, msh, u, f, space_f, patch_list)

  if (nargin < 6)
    patch_list = 1:msh.npatch;
  end

  if (space.npatch ~= msh.npatch)
    error ('op_f_v_mp: the number of patches does not coincide')
  end
  
  rhs = zeros (space.ndof, 1);
  for iptc = patch_list

    if (isempty (space_f.dofs_ornt))
      u_ptc = u(space_f.gnum{iptc});
    else
      u_ptc = u(space_f.gnum{iptc}) .* space_f.dofs_ornt{iptc}.';
    end

    f_loc = @(sp_col, msh_col) f(u_ptc, sp_col, msh_col);
    rhs_loc = op_f_v_tp2 (space.sp_patch{iptc}, msh.msh_patch{iptc}, f_loc, space_f.sp_patch{iptc});
    
    if (~isempty (space.dofs_ornt))
      rhs_loc = space.dofs_ornt{iptc}(:) .* rhs_loc(:);
    end
    rhs(space.gnum{iptc}) = rhs(space.gnum{iptc}) + rhs_loc;
  end

end
