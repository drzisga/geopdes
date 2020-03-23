% OP_NONLINEAR_SU_EV_SURROGATE_2D: 
%
%   mat = op_nonlinear_su_ev_surrogate_2d (spu, spv, msh, Stress, DStress, u_old);
%   [rows, cols, values] = op_nonlinear_su_ev_surrogate_2d (spu, spv, msh, Stress, DStress, u_old);
%
% INPUT:
%    
%   spu:     object representing the space of trial functions (see sp_vector)
%   spv:     object representing the space of test functions (see sp_vector)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_cartesian)
%   Stress:  stress functional
%   DStress: linearized stress functional
%   u:       solution of previous iteration
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
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

function varargout = op_nonlinear_su_ev_surrogate_2d (space, msh, Stress, DStress, u_old, surrogate_opts)

  for icomp = 1:space.ncomp_param
    for idim = 1:msh.ndim
      size1 = size (space.scalar_spaces{icomp}.sp_univ(idim).connectivity);
      size2 = size (space.scalar_spaces{icomp}.sp_univ(idim).connectivity);
      if (size1(2) ~= size2(2) || size1(2) ~= msh.nel_dir(idim))
        error ('One of the discrete spaces is not associated to the mesh')
      end
    end
  end
  
  if msh.ndim ~= 2
    error('op_nonlinear_su_ev_surrogate_2d: This function only supports 2D');
  end

  if ~all(space.ndof_dir == space.ndof_dir(1))
    error('op_nonlinear_su_ev_surrogate_2d: Dofs in each dimension must be equal');
  end

  iga_degree = min(space.scalar_spaces{1}.degree);
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
  ind = unique(round(linspace(1, num_1D_basis_inner, ceil((num_1D_basis_inner-1)/surrogate_opts.skip) + 1)));
  x_sample = x(ind);

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

  K_surr = spalloc (space.ndof, space.ndof, 5*space.ndof);

  for iel = 1:msh.nel_dir(1)
    if ismember(iel, boundary_mask{1}) % Case corresponding to the red elements in the article
      msh_col = msh_evaluate_col(msh, iel);
      sp_col = sp_evaluate_col(space, msh_col, 'value', false, 'gradient', true);
    elseif ismember(iel, element_mask{1}) % Case corresponding to the blue elements in the article
      msh_col = msh_evaluate_col(msh, iel, 'element_mask', element_mask);
      sp_col = sp_evaluate_col(space, msh_col, 'value', false, 'gradient', true, 'element_mask', element_mask);
    else % Case corresponding to the green elements in the article
      msh_col = msh_evaluate_col(msh, iel, 'element_mask', boundary_mask);
      sp_col = sp_evaluate_col(space, msh_col, 'value', false, 'gradient', true, 'element_mask', boundary_mask);
    end

    K_surr = K_surr + op_nonlinear_su_ev (sp_col, sp_col, msh_col, Stress, DStress, u_old);
  end
  
  deltas1D = -iga_degree:iga_degree;
  
%   if (nargout == 2 || nargout == 4)
%     cur_sf = cell(space.ncomp_param, space.ncomp_param, length(deltas1D), length(deltas1D));
%   end
  
  % Extract sub-blocks
  blocklength = space.scalar_spaces{1}.ndof;
  for I = 1:space.ncomp_param
    for J = 1:space.ncomp_param
      
      K_sub = K_surr(((I-1)*blocklength+1):(I*blocklength), ((J-1)*blocklength+1):(J*blocklength));
      sp_i = [];
      sp_j = [];
      sp_v = [];

      % Apply local interpolation surrogate approach
      for i=deltas1D
        for j=deltas1D

          shift = i + num_1D_basis * j;

          % Obtain stencil function values for the sample points
          stencilfunc = K_sub(sub2ind(size(K_sub), row_indices(:), row_indices(:) + shift));
          stencilfunc = reshape(stencilfunc, num_1D_basis_inner, num_1D_basis_inner);
          sf_sample = full(stencilfunc(ind, ind));

          % Interpolate missing values
          
%           if ~isempty(surrogate_opts.cur_sf)
%             interpolationSplines = surrogate_opts.cur_sf(I,J,i+iga_degree+1,j+iga_degree+1);
%             interpolationSplines = interpolationSplines{1};
%           else
            interpolationSplines = spapi({surrogate_opts.degree+1, surrogate_opts.degree+1}, {x_sample, x_sample}, sf_sample);
%           end
          
%           if (nargout == 2 || nargout == 4)
%             cur_sf{I,J,i+iga_degree+1,j+iga_degree+1} = interpolationSplines;
%           end
          
          tmp = spval(interpolationSplines, {x, x});
          
%           [X,Y] = meshgrid(x, x);
%           
%           surf(X,Y,tmp);
%           pause;

          % Add contribution to sparse vectors
          sp_i = [sp_i, row_indices(:).'];
          sp_j = [sp_j, row_indices(:).' + shift];
          sp_v = [sp_v, tmp(:).'];
        end
      end

      % Combine surrogate matrix and standard matrix
      K_interp = sparse(sp_i, sp_j, sp_v, length(K_sub), length(K_sub));
      idx_interp = find(K_interp);
      idx_surr = find(K_sub);
      idx_intersect = intersect(idx_interp,idx_surr);
      K_sub(idx_intersect) = 0;
      K_sub = K_sub + K_interp;
      K_surr(((I-1)*blocklength+1):(I*blocklength), ((J-1)*blocklength+1):(J*blocklength)) = K_sub;
    end
  end
  
  % Enforce zero row-sum property
  K_surr(logical(speye(size(K_surr)))) = diag(K_surr) - sum(K_surr, 2);
  
%   if (nargout == 2 || nargout == 4)
%     surrogate_opts.cur_sf = cur_sf;
%     surrogate_opts.sf{length(surrogate_opts.sf)+1} = cur_sf;
%   end

  if (nargout == 1)
    varargout{1} = K_surr;
  elseif (nargout == 2)
    varargout{1} = K_surr;
    varargout{2} = surrogate_opts;
  elseif (nargout == 3)
    [rows, cols, vals] = find (K_surr);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  elseif (nargout == 4)
    [rows, cols, vals] = find (K_surr);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
    varargout{4} = surrogate_opts;
  end

end
