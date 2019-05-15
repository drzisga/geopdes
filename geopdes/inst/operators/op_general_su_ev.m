% OP_GENERAL_SU_EV: assemble the matrix A = [a(i,j)], a(i,j) = 1/2 (C : epsilon (u_j), epsilon (v_i)).
%
%   mat = op_general_su_ev (spu, spv, msh, C);
%   [rows, cols, values] = op_general_su_ev (spu, spv, msh, C);
%
% INPUT:
%    
%   spu: structure representing the space of trial functions (see sp_vector/sp_evaluate_col)
%   spv: structure representing the space of test functions (see sp_vector/sp_evaluate_col)
%   msh: structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   C:   fourth-order stiffness tensor
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

function varargout = op_general_su_ev (spu, spv, msh, C)

  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights;
  
  jacdet_weights_C = C;
  for i=1:spv.ncomp
    for j = 1:spv.ncomp
      for k = 1:spu.ncomp
        for l = 1:spu.ncomp
          jacdet_weights_C(:,:,i,j,k,l) = jacdet_weights .* C(:,:,i,j,k,l);
        end
      end
    end
  end
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      gradu_iel = reshape (gradu(:,:,:,:,iel), spu.ncomp, ndir, msh.nqn, spu.nsh_max);
      gradv_iel = reshape (gradv(:,:,:,:,iel), spv.ncomp, ndir, msh.nqn, spv.nsh_max);

      jacdet_weights_C_iel = squeeze(jacdet_weights_C(:,iel,:,:,:,:));
      elementary_values = zeros(spv.nsh_max, spu.nsh_max);
      
      for i=1:spv.ncomp
        for j = 1:spv.ncomp
          for k = 1:spu.ncomp
            for l = 1:spu.ncomp
              for qp=1:msh.nqn
                elementary_values = elementary_values + jacdet_weights_C_iel(qp,i,j,k,l) * squeeze(gradv_iel(i,j,qp,:)) * squeeze(gradu_iel(k,l,qp,:))';
              end
            end
          end
        end
      end

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