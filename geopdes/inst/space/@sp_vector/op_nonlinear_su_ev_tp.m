% OP_NONLINEAR_SU_EV_TP: 
%
%   mat = op_nonlinear_su_ev_tp (spu, spv, msh, lambda, mu, u);
%   [rows, cols, values] = op_nonlinear_su_ev_tp (spu, spv, msh, lambda, mu, u);
%
% INPUT:
%    
%   spu:     object representing the space of trial functions (see sp_vector)
%   spv:     object representing the space of test functions (see sp_vector)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_cartesian)
%   lambda, mu: function handles to compute the Lame' coefficients
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

function varargout = op_nonlinear_su_ev_tp (space, msh, lambda, mu, u)

  for icomp = 1:space.ncomp_param
    for idim = 1:msh.ndim
      size1 = size (space.scalar_spaces{icomp}.sp_univ(idim).connectivity);
      size2 = size (space.scalar_spaces{icomp}.sp_univ(idim).connectivity);
      if (size1(2) ~= size2(2) || size1(2) ~= msh.nel_dir(idim))
        error ('One of the discrete spaces is not associated to the mesh')
      end
    end
  end

  A = spalloc (space.ndof, space.ndof, 5*space.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col = sp_evaluate_col (space, msh_col, 'value', false, 'gradient', true);

    gradientU = sp_eval_msh(u, sp_col, msh_col, 'gradient');

    A = A + op_nonlinear_su_ev (sp_col, sp_col, msh_col, lambda, mu, gradientU);
  end

  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end
