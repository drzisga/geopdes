% OP_F_EV: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, eps(v_i)).
%
%   rhs = op_f_ev (spv, msh, coeff);
%
% INPUT:
%
%   spv:   structure representing the function space (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   coeff: source function evaluated at the quadrature points
%
% OUTPUT:
%
%   rhs: assembled right-hand side
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2017 Rafael Vazquez
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

function rhs = op_f_ev (spv, msh, FS)

 rhs   = zeros (spv.ndof, 1);
 gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
 
 ndir = size (gradv, 2);

 for iel = 1:msh.nel
   if (all (msh.jacdet(:,iel)))
     jacdet_weights = reshape (msh.jacdet(:, iel) .* msh.quad_weights(:, iel), 1, 1, msh.nqn);
     jacdet_weights_coeff_iel = jacdet_weights .* FS(:,:,:,iel);
     
     gradv_iel = gradv(:,:,:,:,iel);
     gradv_iel = reshape (gradv_iel, spv.ncomp, ndir, msh.nqn, spv.nsh_max);
     gradv_iel = reshape(gradv_iel, [spv.ncomp*ndir, msh.nqn, spv.nsh_max, 1]);
     
     jacdet_weights_coeff_iel = reshape(jacdet_weights_coeff_iel, [spv.ncomp*ndir, msh.nqn, 1, 1]);

     aux_val = jacdet_weights_coeff_iel .* gradv_iel;
     rhs_loc = sum (sum (aux_val, 1), 2);

     indices = find (spv.connectivity(:,iel));
     rhs_loc = rhs_loc(indices); conn_iel = spv.connectivity(indices,iel);
     rhs(conn_iel) = rhs(conn_iel) + rhs_loc(:); 
   else
     warning ('geopdes:jacdet_zero_at_quad_node', 'op_f_v: singular map in element number %d', iel)
   end
 end
 
end
