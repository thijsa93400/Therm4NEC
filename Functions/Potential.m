function [V] = Potential( r, J, D_e, Beta, R_e, mu)
% POTENTIAL evaluates the morse potential including centrifugal component
%   POTENTIAL( r, J, D_e, Beta, R_e, mu)
%   Input:
%       - r:        Position in morse potential (m)
%       - J:        Rotational level
%       - D_e:      Disocation energy measured from equilibrium (J)
%       - beta:     Anharmonicity (m^-1)
%       - R_e:      Equilibrium distance (m)
%   Output:
%       - V:        Potential at r (J)

%   This code is part of Therm4NEC
%   Copyright (C) 2021  T. Hazenberg
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU Lesser General Public License as 
%   published by the Free Software Foundation, either version 3 of the 
%   License, or (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser General Public License for more details.
% 
%   You should have received a copy of the GNU Lesser General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

% constants
h = 6.62607004e-34; % (m^2kg/s)

% compute
ebr_r = exp(-Beta*( r- R_e));
Vv    = D_e*(1-ebr_r).^2;
VJ    = h^2./(8*pi^2*mu.*r.^2) * J*(J+1);
V     = Vv + VJ;

end