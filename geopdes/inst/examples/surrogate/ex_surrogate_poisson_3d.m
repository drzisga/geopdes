% EX_SURROGATE_POISSON_3D: Example which solves the Poisson problem in 3D
% using the surrogate method and compares it to the standard assembly.
%
% This examples considers the 3D benchmark with a constant sampling length
% discussed in Section 7.3 of
% 
%   Daniel Drzisga, Brendan Keith, and Barbara Wohlmuth.
%   The surrogate matrix methodology: Low-cost assembly for isogeometric
%   analysis
%
% An explanation of this reference implementation may be found in
% 
%   Daniel Drzisga, Brendan Keith, and Barbara Wohlmuth.
%   The surrogate matrix methodology: A reference implementation for
%   low-cost assembly in isogeometric analysis
%
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

p = 2;  % IGA degree
q = 3;  % Surrogate interpolation degree (1 or 3)
M = 10; % Sampling parameter
N = 40; % Number of knots in one dimension

% Construct manufactured solution
fprintf('Initializing problem...\n');
solution      = @(x,y,z) sin(20*pi*x).*sin(20*pi*y).*sin(20*pi*z);
solution_dx   = @(x,y,z) 20*pi*cos(20*pi*x).*sin(20*pi*y).*sin(20*pi*z);
solution_dy   = @(x,y,z) 20*pi*cos(20*pi*y).*sin(20*pi*x).*sin(20*pi*z);
solution_dz   = @(x,y,z) 20*pi*cos(20*pi*z).*sin(20*pi*x).*sin(20*pi*y);
grad_solution = @(x,y,z) cat (1, ...
                  reshape(feval(solution_dx, x, y, z), [1, size(x)]), ...
                  reshape(feval(solution_dy, x, y, z), [1, size(x)]), ...
                  reshape(feval(solution_dz, x, y, z), [1, size(x)]));
coeff = @(x,y,z) 1;
f = @(x,y,z) 1200*pi^2*sin(20*pi*x).*sin(20*pi*y).*sin(20*pi*z);

% Load volume
volume = GeomBentTwistedBox(p, N);

% Set up geometry and discrete spaces
geometry = geo_load(volume);
rule     = msh_gauss_nodes(geometry.nurbs.order);
[qn,qw]  = msh_set_quad_nodes(geometry.nurbs.knots, rule);  
msh      = msh_cartesian (geometry.nurbs.knots, qn, qw, geometry); 
space    = sp_nurbs (geometry.nurbs, msh);

% Setting up the right hand side vector
rhs = op_f_v_tp (space, msh, f);

% Extract Dirichlet boundary dofs
drchlt_dofs = [];
for i = 1:6
  drchlt_dofs = union (drchlt_dofs, space.boundary(i).dofs);
end

% The remining ones are internal DOFs
int_dofs = setdiff(1:space.ndof, drchlt_dofs);

% Initialize the solution vector
u = zeros(space.ndof,1);

% Set the Dirichlet-DOFs
drchlt_sides = [1 2 3 4 5 6];
bc_solution = @(x, y, z, ind) solution(x,y,z);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, bc_solution, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

% Assembly (by numerical integration) of the stiffness matrix
fprintf('Assembling standard IGA matrix...\n');
tic;
K = op_gradu_gradv_tp (space, space, msh, coeff);
time_std = toc;
fprintf('Standard assembly time: %f s\n', time_std);
rhs_std = rhs - K*u;

% Solve systems
fprintf('Solving standard IGA problem...\n');
tic
u(int_dofs) = K(int_dofs, int_dofs)\rhs_std(int_dofs);
time = toc;
fprintf('Standard solve time: %f s\n', time);

u_surr = zeros(space.ndof,1);
u_surr(drchlt_dofs) = u_drchlt;

% Surrogate assembly
fprintf('Assembling surrogate IGA matrix...\n');
tic;
K_surr = op_gradu_gradv_surrogate_3d(space, msh, coeff, M, q);
time_surr = toc;
fprintf('Surrogate assembly time: %f s\n', time_surr);
rhs_surr = rhs - K_surr * u_surr;

% Surrogate solve
fprintf('Solving surrogate IGA problem...\n');
tic
u_surr(int_dofs) = K_surr(int_dofs, int_dofs)\rhs_surr(int_dofs);
time = toc;
fprintf('Surrogate solve time: %f s\n', time);

% Error computation
fprintf('Computing errors...\n');
[errh1, errl2]  = sp_h1_error (space, msh, u, solution, grad_solution);
[normh1, norml2] = sp_h1_error (space, msh, zeros(size(u)), solution, grad_solution);

errh1 = errh1 / normh1;
errl2 = errl2 / norml2;

fprintf('Relative error in standard IGA\n');
fprintf('  L2-norm: %d \n', errl2);
fprintf('  H1-norm: %d \n\n', errh1);

[err_surr_h1, err_surr_l2]  = sp_h1_error (space, msh, u_surr, solution, grad_solution);

err_surr_h1 = err_surr_h1 / normh1;
err_surr_l2 = err_surr_l2 / norml2;

fprintf('\n\\|A-\\tilde{A}\\|_{\\max} = %d\n\n', full(max(max(abs(K-K_surr)))));
fprintf('Relative error in surrogate IGA  \n');
fprintf('  L2-norm: %d \n',err_surr_l2);
fprintf('  H1-norm: %d \n',err_surr_h1);

fprintf('\nAssembly speed-up: %.2f%%\n', 100 * (time_std / time_surr - 1));
