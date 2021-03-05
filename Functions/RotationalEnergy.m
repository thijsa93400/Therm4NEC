function [Erot] = RotationalEnergy( J, v, B_e, alpha_e, D_e, Beta_e)
% ROTATIONALENERGY computes the rotational energy
%   ROTATIONALENERGY( J, v, B_e, alpha_e, D_e, Beta_e)
%   Input:
%       - J:        Rotational quantum number
%       - v:        Vibrational quantum number
%       - B_e:      Rotational constant at equilibrium (m^-1)
%       - alpha_e:  Rotational constant first term (m^-1)
%       - D_e:      Centrifugal distortion (m^-1)
%       - Beta_e:   Centrifugal distortion first term(m^-1)
%   Output:
%       - Erot:     Energy due to rotations (J)

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
c = 299792458;      % (m/s)

% Rotational constant out of equilibrium
B = B_e - alpha_e.*(v + 0.5);
D = D_e - Beta_e.*(v + 0.5);

% Rotational energy
Erot = B.*J.*(J+1) - D.*J.^2.*(J+1).^2;
Erot = Erot .* h .*c;

end