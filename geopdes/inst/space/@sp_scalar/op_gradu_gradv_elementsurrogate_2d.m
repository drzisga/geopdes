% OP_GRADU_GRADV_ELEMENTSURROGATE_2D: assemble the stiffness matrix A = [a(i,j)],
% a(i,j) = (epsilon grad u_j, grad v_i) using the element surrogate method.
%
%   mat = op_gradu_gradv_elementsurrogate_2d (space, msh, coeff, M, q);
%   [rows, cols, values] = op_gradu_gradv_elementsurrogate_2d (space, msh, coeff, M, q);
%
% INPUT:
%
%   space: structure representing the space of trial and test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   coeff: diffusion coefficient
%   M:     skip parameter
%   q:     surrogate interpolation degree (1 to 5)
%
% OUTPUT:
%
%   mat:    assembled surrogate stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2017 Rafael Vazquez
% Copyright (C) 2019       Daniel Drzisga
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

function varargout = op_gradu_gradv_elementsurrogate_2d(space, msh, coeff, M, q)

plots = false;

if plots
    element_usage_assemble = zeros(msh.nel_dir(1) * msh.nel_dir(1),1);
end

if q == 1
  method = 'linear';
elseif q == 3
  method = 'spline';
else
  error('op_gradu_gradv_elementsurrogate_2d: q = %d is not supported', q);
end

elems_all{1} = 1:msh.nel_dir(1);
elems_inner{1} = elems_all{1}(3:(end-2));
elems_bdry{1} = elems_all{1}([1:2, (end-1):end]);
elems_sample{1} = unique(elems_inner{1}([1:M:end, end]));

h = 1/msh.nel_dir(1);

x = linspace(0.5 * h, 1-0.5*h, msh.nel_dir(1));
x_inner = x(3:(end-2));
x_sampled = unique(x_inner([1:M:end, end]));

iga_degree = min(space.degree);

num_basis_per_element_1d = (iga_degree+1)^2;
num_basis_per_element_2d = num_basis_per_element_1d^2;


local_matrices = zeros(length(x_sampled), length(x_sampled), num_basis_per_element_2d);

% Integrate active elements
j = 1;
for jel = elems_sample{1}
    
    msh_col = msh_evaluate_col(msh, jel, 'element_mask', elems_sample);
    sp_col = sp_evaluate_col(space, msh_col, 'value', false, 'gradient', true, 'element_mask', elems_sample);
    
    if plots
        element_usage_assemble(msh_col.elem_list) = 1;
    end

    if (nargin == 4)
      for idim = 1:msh.rdim
        x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
      end
      coeffs = coeff (x{:});
    else
      coeffs = ones (msh_col.nqn, msh_col.nel);
    end

    [~, ~, vals] = op_gradu_gradv(sp_col, sp_col, msh_col, coeffs);
    
    i = 1;
    for iel = elems_sample{1}
      local_matrices(i, j, :) = vals(((i-1) * num_basis_per_element_2d + 1):(i * num_basis_per_element_2d));
      i = i + 1;
    end
    
    j = j + 1;
end

[X_sampled, Y_sampled] = meshgrid(x_sampled, x_sampled);
[X, Y] = meshgrid(x_inner, x_inner);

K_surr = spalloc(space.ndof, space.ndof, 3*space.ndof);

for jel = elems_all{1}
  if ismember(jel, elems_bdry{1})
    msh_col = msh_evaluate_col(msh, jel, 'element_mask', elems_all);
    sp_col = sp_evaluate_col(space, msh_col, 'value', false, 'gradient', true, 'element_mask', elems_all);
  else
    msh_col = msh_evaluate_col(msh, jel, 'element_mask', elems_bdry);
    sp_col = sp_evaluate_col(space, msh_col, 'value', false, 'gradient', true, 'element_mask', elems_bdry);
  end

  if plots
      element_usage_assemble(msh_col.elem_list) = 1;
  end

  if (nargin == 4)
    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end
    coeffs = coeff (x{:});
  else
    coeffs = ones (msh_col.nqn, msh_col.nel);
  end

  K_surr = K_surr + op_gradu_gradv(sp_col, sp_col, msh_col, coeffs);
end

msh_col = msh_evaluate_col(msh, jel, 'element_mask', elems_inner);
sp_col = sp_evaluate_col(space, msh_col, 'value', false, 'element_mask', elems_inner);

rows = zeros (length(elems_inner{1}) * msh_col.nel * sp_col.nsh_max * sp_col.nsh_max, 1);
cols = zeros (length(elems_inner{1}) * msh_col.nel * sp_col.nsh_max * sp_col.nsh_max, 1);
values = zeros (msh_col.nel * sp_col.nsh_max * sp_col.nsh_max, length(elems_inner{1}));

for sfunc=1:num_basis_per_element_2d
  sfunc_sampled = local_matrices(:, :, sfunc);
  tmp = interp2(X_sampled, Y_sampled, sfunc_sampled, X, Y, method);
  
  for i=1:length(elems_inner{1})
    values(sfunc:num_basis_per_element_2d:end, i) = tmp(:,i);
  end
end

ncounter = 0;

for jel = elems_inner{1}
  msh_col = msh_evaluate_col_stripped(msh, jel, elems_inner);
  sp_col = sp_evaluate_col(space, msh_col, 'value', false, 'element_mask', elems_inner);

  if plots
    element_usage_assemble(msh_col.elem_list) = 2;
  end
  
  for iel = 1:msh_col.nel
      [rows_loc, cols_loc] = meshgrid(sp_col.connectivity(:,iel), sp_col.connectivity(:,iel));
      rows(ncounter+(1:sp_col.nsh(iel)*sp_col.nsh(iel))) = rows_loc;
      cols(ncounter+(1:sp_col.nsh(iel)*sp_col.nsh(iel))) = cols_loc;
      ncounter = ncounter + sp_col.nsh(iel)*sp_col.nsh(iel);
  end
end

K_surr = K_surr + sparse(rows(1:ncounter), cols(1:ncounter), values(1:ncounter), space.ndof, space.ndof);

% Enforce zero row-sum property
K_surr(logical(speye(size(K_surr)))) = diag(K_surr) - sum(K_surr, 2);

if plots
    element_usage_assemble = reshape(element_usage_assemble, msh.nel_dir(1), msh.nel_dir(1)).';
    figure;
    imagesc(element_usage_assemble);
    hold on;
    for i = 1:msh.nel_dir(1)
       plot([.5,msh.nel_dir(1)+0.5],[i-.5,i-.5],'k-');
       plot([i-.5,i-.5],[.5,msh.nel_dir(1)+0.5],'k-');
    end
    colorbar;
    title('Actually assembled elements');
    hold off;
end

if (nargout == 1)
  varargout{1} = K_surr;
elseif (nargout == 3)
  [rows, cols, vals] = find (K_surr);
  varargout{1} = rows;
  varargout{2} = cols;
  varargout{3} = vals;
end

end
