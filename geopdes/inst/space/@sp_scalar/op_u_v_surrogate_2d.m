% OP_U_V_SURROGATE_2D: assemble the mass matrix M = [m(i,j)], m(i,j) = (mu u_j, v_i), exploiting the tensor product structure
% using the surrogate method
%
%   mat = op_u_v_surrogate_2d (space, msh, coeff, M, q);
%   [rows, cols, values] = op_u_v_surrogate_2d (space, msh, coeff, M, q);
%
% INPUT:
%
%   space: structure representing the space of trial and test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   coeff: diffusion coefficient
%   M:     skip parameter
%   q:     surrogate interpolation degree (1 or 3)
%
% OUTPUT:
%
%   mat:    assembled surrogate mass matrix
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

function varargout = op_u_v_surrogate_2d(space, msh, coeff, M, q)

if msh.ndim ~= 2
  error('op_u_v_surrogate_2d: This function only supports 2D');
end

if ~all(space.ndof_dir == space.ndof_dir(1))
  error('op_u_v_surrogate_2d: Dofs in each dimension must be equal');
end

if q == 1
  method = 'linear';
elseif q == 3
  method = 'spline';
else
  error('op_u_v_surrogate_2d: q = %d is not supported', q);
end

iga_degree = min(space.degree);
num_1D_basis = space.ndof_dir(1);
num_1D_basis_inner = num_1D_basis - 4 * iga_degree;

% Obtain the rows corresponding to the cardinal basis functions
row_indices = 1:num_1D_basis;
row_indices = repmat(row_indices, num_1D_basis, 1);
for i = 2:num_1D_basis
  row_indices(i,:) = row_indices(i-1,:) + num_1D_basis;
end
ind = 2*iga_degree+1:num_1D_basis-2*iga_degree;
row_indices = row_indices(ind, ind);

% Create sample points
x = linspace(0, 1, num_1D_basis_inner);
[X, Y] = meshgrid(x);
ind = unique([1:M:num_1D_basis_inner, num_1D_basis_inner]);
X_sample  = X(ind, ind);
Y_sample  = Y(ind, ind);

% Take only the necessary rows 
row_indices_subset = row_indices(ind, ind);
row_indices_subset = row_indices_subset(:);
row_indices_subset = row_indices_subset(row_indices_subset <= (2*iga_degree + 1) * num_1D_basis);

% Create mask for active elements
element_mask{1} = sp_get_cells(space, msh, row_indices_subset.');
element_mask{1} = element_mask{1}(element_mask{1} <= (iga_degree + 1) * msh.nel_dir(1));
element_mask{1} = element_mask{1} - iga_degree * msh.nel_dir(1);
element_mask{1} = [(1:(2*iga_degree)).'; element_mask{1}; ((msh.nel_dir(1)-(2*iga_degree-1)):msh.nel_dir(1)).'];
element_mask{1} = unique(element_mask{1});

% Create mask for boundary elements
boundary_mask{1} = [1:(2*iga_degree), (msh.nel_dir(1)-(2*iga_degree-1)):msh.nel_dir(1)];
boundary_mask{1} = unique(boundary_mask{1});

% Allocate matrix memory
K_surr = spalloc(space.ndof, space.ndof, (2*iga_degree + 1)^2*space.ndof);

% Assemble only necessary elements
for iel = 1:msh.nel_dir(1)
  if ismember(iel, boundary_mask{1}) % Case corresponding to the red elements in the article
    msh_col = msh_evaluate_col(msh, iel);
    sp_col = sp_evaluate_col(space, msh_col);
  elseif ismember(iel, element_mask{1}) % Case corresponding to the blue elements in the article
    msh_col = msh_evaluate_col(msh, iel, 'element_mask', element_mask);
    sp_col = sp_evaluate_col(space, msh_col, 'element_mask', element_mask);
  else % Case corresponding to the green elements in the article
    msh_col = msh_evaluate_col(msh, iel, 'element_mask', boundary_mask);
    sp_col = sp_evaluate_col(space, msh_col, 'element_mask', boundary_mask);
  end

  for idim = 1:msh.rdim
    coords{idim} = reshape(msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
  end
  coeffs = coeff(coords{:});

  K_surr = K_surr + op_u_v(sp_col, sp_col, msh_col, coeffs);
end

sp_i = [];
sp_j = [];
sp_v = [];

% Apply local interpolation surrogate approach
for i=-iga_degree:iga_degree
  for j=-iga_degree:iga_degree
      
    shift = i + num_1D_basis * j;
        
    % Skip if not upper part or on the diagonal
    if shift < 0
      continue;
    end
    
    % Obtain stencil function values for the sample points
    stencilfunc = K_surr(sub2ind(size(K_surr), row_indices(:), row_indices(:) + shift));
    stencilfunc = reshape(stencilfunc, num_1D_basis_inner, num_1D_basis_inner);
    sf_sample = full(stencilfunc(ind, ind));
        
    % Interpolate missing values
    tmp = interp2(X_sample, Y_sample, sf_sample, X, Y, method);
        
    % Add contribution to sparse vectors
    % Add diagonal only once
    if shift == 0
      sp_i = [sp_i, row_indices(:).'];
      sp_j = [sp_j, row_indices(:).'];
      sp_v = [sp_v, tmp(:).'];
    else
      sp_i = [sp_i, row_indices(:).', row_indices(:).' + shift];
      sp_j = [sp_j, row_indices(:).' + shift, row_indices(:).'];
      sp_v = [sp_v, tmp(:).', tmp(:).'];
    end
  end
end

% Combine surrogate matrix and standard matrix
K_interp = sparse(sp_i, sp_j, sp_v, length(K_surr), length(K_surr));
idx_interp = find(K_interp);
idx_surr = find(K_surr);
idx_intersect = intersect(idx_interp,idx_surr);
K_surr(idx_intersect) = 0;
K_surr = K_surr + K_interp;

if (nargout == 1)
  varargout{1} = K_surr;
elseif (nargout == 3)
  [rows, cols, vals] = find (K_surr);
  varargout{1} = rows;
  varargout{2} = cols;
  varargout{3} = vals;
end

end
