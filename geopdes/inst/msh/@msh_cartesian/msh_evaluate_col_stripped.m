% MSH_EVALUATE_COL_STRIPPED: stripped down evaluation of the parameterization in one column of the mesh.
%
%     msh_col = msh_evaluate_col_stripped (msh, colnum, element_mask)
%
% INPUTS:
%
%    msh:          mesh object (see msh_cartesian)
%    colnum:       number of the "column", i.e., the element in the first parametric direction.
%    element_mask: list of elements which are assembled
%
% OUTPUT:
%
%     msh_col: structure containing the quadrature rule in one column of the physical domain, which contains the following fields
%
%     FIELD_NAME    (SIZE)                    DESCRIPTION
%     ndim          (scalar)                  dimension of the parametric space
%     rdim          (scalar)                  dimension of the physical space
%     colnum        (scalar)                  number of the column
%     nel           (scalar)                  number of elements in the column
%     elem_list     (nel vector)              indices of the elements in the column
%     nel_dir       (1 x ndim vector)         number of elements in each parametric direction for the entire mesh
%     nqn           (scalar)                  number of quadrature nodes per element
%     nqn_dir       (1 x ndim vector)         number of quadrature nodes per element in each parametric direction for the entire mesh
%     quad_nodes    (ndim x nqn x nel vector) coordinates of the quadrature nodes in parametric space
%     quad_weights  (nqn x nel vector)        weights associated to the quadrature nodes
%
%  For more details, see the documentation
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
% Copyright (C) 2014 Elena Bulgarello, Carlo de Falco, Sara Frizziero
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

function msh_col = msh_evaluate_col_stripped (msh, colnum, element_mask)

  msh_col.ndim = msh.ndim;
  msh_col.rdim = msh.rdim;

  msh_col.colnum = colnum;
  msh_col.nel_dir = msh.nel_dir;
  msh_col.nel_dir(1) = 1;

  if (msh.ndim == 1)
    msh_col.elem_list = colnum;
  elseif (msh.ndim == 2)
    msh_col.elem_list = colnum + msh.nel_dir(1)*(0:msh.nel_dir(2)-1);
    msh_col.elem_list = msh_col.elem_list(element_mask{1});
    msh_col.nel_dir = [1, length(msh_col.elem_list)];
  elseif (msh.ndim == 3)
    indu = colnum * ones(msh.nel_dir(2), msh.nel_dir(3));
    indv = repmat ((1:msh.nel_dir(2))', 1, msh.nel_dir(3));
    indw = repmat ((1:msh.nel_dir(3)), msh.nel_dir(2), 1);
    elem_list = sub2ind ([msh.nel_dir(1), msh.nel_dir(2), msh.nel_dir(3)], indu, indv, indw);
    elem_list = elem_list(element_mask{1}, element_mask{2});
    msh_col.nel_dir = [1, length(element_mask{1}), length(element_mask{2})];
    msh_col.elem_list = elem_list(:);
  end
  
  msh_col.nel = numel(msh_col.elem_list);
%   msh_col.nel  = msh.nelcol;

  msh_col.nqn_dir = msh.nqn_dir;
  msh_col.nqn  = msh.nqn;

  msh_col.qn = msh.qn; 
  msh_col.qn{1} = msh.qn{1}(:,colnum);
  for idim=2:msh.ndim
    msh_col.qn{idim} = msh.qn{idim}(:,element_mask{idim-1});
  end
  
  if (~isempty (msh.qw))
    msh_col.qw = msh.qw; 
    msh_col.qw{1} = msh.qw{1}(:,colnum);
    for idim=2:msh.ndim
      msh_col.qw{idim} = msh.qw{idim}(:,element_mask{idim-1});
    end
  end

end
