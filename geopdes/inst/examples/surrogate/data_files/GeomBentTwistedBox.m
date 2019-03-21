% GEOMBENTTWISTEDBOX: Generates a volume of a twisted box
%
%   vol = GeomBentTwistedBox (degree, num_knots);
%
% INPUT:
%
%   degree:    IGA space degree
%   num_knots: number of knots in one dimension
%
% OUTPUT:
%
%   vol: volume of a twisted box
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

function [vol] = GeomBentTwistedBox(degree, num_knots)

coefs(:,1,1,1) = [ 1.0  0.0  0.0  1.0];
coefs(:,2,1,1) = [ 2.0  0.0  0.0  1.0];
coefs(:,1,2,1) = [ 0.75 0.75 0.25 1.0];
coefs(:,2,2,1) = [ 1.75 1.75 0.25 1.0];
coefs(:,1,3,1) = [ 0.0  1.0  2.0  1.0];
coefs(:,2,3,1) = [ 0.0  2.0  2.0  1.0];


coefs(:,1,1,2) = [ 1.0  0.0  1.0  1.0];
coefs(:,2,1,2) = [ 2.0  0.0  1.0  1.0];
coefs(:,1,2,2) = [ 1.0  0.75 1.25 1.0];
coefs(:,2,2,2) = [ 1.75 1.75 1.25 1.0];
coefs(:,1,3,2) = [ 1.0  1.0  2.0  1.0];
coefs(:,2,3,2) = [ 1.0  2.0  2.0  1.0];

knots{1} = [0 0 1 1];
knots{2} = [0 0 0 1 1 1];
knots{3} = [0 0 1 1];

vol = nrbmak(coefs, knots);

% Degree elevation
elevate = degree-1;
if elevate > 0
  vol = nrbdegelev(vol, elevate-1*[0 1 0]);
else
 error('GeomBendedTwistedBox: needs a larger polynomial degree.');
end

% Insert knots
knots = linspace(0, 1, num_knots);
inner_knots = knots(2:end-1);
vol = nrbkntins(vol,{inner_knots inner_knots inner_knots});

end

