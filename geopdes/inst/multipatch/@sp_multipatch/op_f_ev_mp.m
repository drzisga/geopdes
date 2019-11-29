% OP_F_EV_MP: assemble the right-hand side vector r
%
%   rhs = op_f_ev_mp (spv, msh, coeff, [patches]);
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

function rhs = op_f_ev_mp (space, msh, u, f, patch_list)

  if (nargin < 5)
    patch_list = 1:msh.npatch;
  end

  if (space.npatch ~= msh.npatch)
    error ('op_f_v_mp: the number of patches does not coincide')
  end
  
  rhs_cell = cell(length(patch_list), 1);
  
  parfor iptc = patch_list
    
    if (isempty (space.dofs_ornt))
      u_ptc = u(space.gnum{iptc});
    else
      u_ptc = u(space.gnum{iptc}) .* space.dofs_ornt{iptc}.';
    end
    
    f_ = @(space,msh) f(u_ptc, space, msh);
    
    rhs_loc = op_f_ev_tp_vec (space.sp_patch{iptc}, msh.msh_patch{iptc}, f_);
    
    if (~isempty (space.dofs_ornt))
      rhs_loc = space.dofs_ornt{iptc}(:) .* rhs_loc(:);
    end
    
    rhs_cell{iptc} = rhs_loc;
  end
  
  rhs = zeros (space.ndof, 1);
  for iptc = patch_list
    rhs(space.gnum{iptc}) = rhs(space.gnum{iptc}) + rhs_cell{iptc};
  end  

end