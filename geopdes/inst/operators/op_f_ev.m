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

function rhs = op_f_ev (spv, msh, u, lambda, mu)

 rhs   = zeros (spv.ndof, 1);
 gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
 
 ndir = size (gradv, 2);
 
 jacdet_weights = msh.jacdet .* msh.quad_weights;

 for iel = 1:msh.nel
   if (all (msh.jacdet(:, iel)))
      gradv_iel = reshape (gradv(:,:,:,:,iel), spv.ncomp, ndir, msh.nqn, spv.nsh_max);

      jacdet_weights_iel = jacdet_weights(:,iel);
      
      epsv_iel = zeros(size(gradv_iel));
      for i=1:spv.nsh_max
        for qp=1:msh.nqn
          epsv_iel(:,:,qp,i) = 0.5 * (gradv_iel(:,:,qp,i) + gradv_iel(:,:,qp,i).');
        end
      end

      elementary_values = zeros(spv.nsh_max, spv.nsh_max);
      
      indices = find (spv.connectivity(:,iel));
      u_loc = u(spv.connectivity(:,iel));
      
      grad_u = zeros(2, 2, msh.nqn);
      for qp = 1:msh.nqn
        for j = 1:spv.nsh_max
          grad_u(:,:,qp) = grad_u(:,:,qp) + u_loc(j) * gradv_iel(:,:,qp,j);
        end
      end
      
      rhs_loc = zeros(size(u_loc));
      
      for i=1:spv.nsh_max
        for qp=1:msh.nqn

          %N = grad_u(:,:,qp).' * grad_u(:,:,qp);
          E = 0.5 * (grad_u(:,:,qp) + grad_u(:,:,qp).');

          lhs = lambda * trace(E) * eye(2) + 2 * mu * E;
          lhs = jacdet_weights_iel(qp) * lhs;

          val = lhs .* epsv_iel(:,:,qp,i);
          val = sum(sum(val));

          rhs_loc(i) = rhs_loc(i) + val;
        end
      end
      rhs_loc = rhs_loc(indices); conn_iel = spv.connectivity(indices,iel);
      rhs(conn_iel) = rhs(conn_iel) + rhs_loc(:); 
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_general_su_ev: singular map in element number %d', iel)
   end
 end
 
end
