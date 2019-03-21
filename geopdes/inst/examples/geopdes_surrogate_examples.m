% GEOPDES_SURROGATE_EXAMPLES: Run surrogate method examples.
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

function [] = geopdes_surrogate_examples ()

iopt = 1; 
while (iopt > 0)
  clc;
  fprintf (1, ...
           ['GeoPDEs examples menu:\n', ...
            '----------------------\n', ...
            '\n', ...
            '   (1) Surrogate method example in 2D: Poisson problem.\n \n', ...
            '   (2) Surrogate method example in 3D: Poisson problem.\n \n']);

  iopt = input ('Please choose a number from above or press <Enter> to return: ');
  clc;

  if (~isempty (iopt))
    switch (iopt)
      case 1
         fprintf('You can have a look at the source file: ex_surrogate_poisson_2d.m \n \n');
         ex_surrogate_poisson_2d;
         input ('Press <Enter> to continue: ');
       case 2
         fprintf('You can have a look at the source file: ex_surrogate_poisson_3d.m \n \n');
         ex_surrogate_poisson_3d;
         input ('Press <Enter> to continue: ');
       otherwise
         continue
    end
  end

end
end

