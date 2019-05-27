% OP_NONLINEAR_SU_EV:
%
%   mat = op_nonlinear_su_ev (spu, spv, msh, lambda, mu, gradientU);
%   [rows, cols, values] = op_nonlinear_su_ev (spu, spv, msh, lambda, mu, gradientU);
%
% INPUT:
%    
%   spu: structure representing the space of trial functions (see sp_vector/sp_evaluate_col)
%   spv: structure representing the space of test functions (see sp_vector/sp_evaluate_col)
%   msh: structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   lambda, mu: Lame parameters
%   gradientU: gradient of u evaluated at quadrature points
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
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

function varargout = op_nonlinear_su_ev (spu, spv, msh, lambda, mu, gradientU)

  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights;
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      gradu_iel = reshape (gradu(:,:,:,:,iel), spu.ncomp, ndir, msh.nqn, spu.nsh_max);
      gradv_iel = reshape (gradv(:,:,:,:,iel), spv.ncomp, ndir, msh.nqn, spv.nsh_max);

      jacdet_weights_iel = reshape(jacdet_weights(:,iel), 1, 1, msh.nqn, 1);
      grad_iel = gradientU(:,:,:,iel);
      grad_iel_T = permute(grad_iel, [2 1 3 4]);
      
      E = zeros(2, 2, msh.nqn, spv.nsh_max);
      
      for j = 1:spv.nsh_max
        for qp = 1:msh.nqn
          E(:,:,qp,j) = grad_iel_T(:,:,qp) * gradu_iel(:,:,qp,j);
        end
      end
      
      E = E + permute(E, [2 1 3 4]);
      
      epsv_iel = 0.5 * (gradv_iel + permute(gradv_iel, [2 1 3 4]));
      epsv_iel = reshape(epsv_iel, [spv.ncomp*ndir, msh.nqn, spv.nsh_max, 1]);
      
      trE = E(1,1,:,:) + E(2,2,:,:);
      trEId = bsxfun(@times, trE, eye(2));
      
      lhs = mu * E + (0.5 * lambda) * trEId;
      lhs = bsxfun(@times, jacdet_weights_iel, lhs);
      lhs = reshape(lhs, [spv.ncomp*ndir, msh.nqn, 1, spv.nsh_max]);
      
      aux_val = bsxfun(@times, lhs, epsv_iel);
      elementary_values = sum(sum(aux_val, 1), 2);

      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
      indices = rows_loc & cols_loc;
      rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc(indices);
      cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc(indices);
      values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = elementary_values(indices);
      ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_general_su_ev: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_general_su_ev: wrong number of output arguments')
  end

end
