% GEOMQUARTERANNULUSWITHBUMPS: Generates a surface of a quarter annulus
% with bumps
%
%   surface = GeomQuarterAnnulusWithBumps (degree, num_knots);
%
% INPUT:
%
%   degree:    IGA space degree
%   num_knots: number of knots in one dimension
%
% OUTPUT:
%
%   surface: surface of a quarter annulus with bumps
% 
% Copyright (C) 2019 Brendan Keith
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

function [surface] = GeomQuarterAnnulusWithBumps(degree, num_knots)

bottom = nrbcirc(sqrt(2)*0.25, [0.75, 0.25], 5*pi/4, 7*pi/4);
left   = nrbcirc(sqrt(2)*0.25, [-0.25, 0.75], -pi/4, pi/4);
outer  = nrbcirc(1.0, [0.0, 0.0], pi/2, -2*pi);
inner  = nrbcirc(0.5, [0.0, 0.0], pi/2, -2*pi);

surface = nrbcoons(inner, outer, left, bottom);

% Degree elevation
elevate = degree-2;
if elevate > 0
  surface = nrbdegelev(surface, elevate*[1 1]);
end

% Insert knots
knots = linspace(0, 1, num_knots);
inner_knots = knots(2:end-1);
surface = nrbkntins(surface,{inner_knots, inner_knots});

end

