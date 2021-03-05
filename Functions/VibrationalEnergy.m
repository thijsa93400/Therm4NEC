function [Evib] = VibrationalEnergy( v, omega_e, omegaX_e, omegaY_e)
% VIBRATIONALENERGY computes the vibrational energy
%   VIBRATIONALENERGY( v, omega_e, omegaX_e, omegaY_e)
%   Input:
%      - v:        Vibrational quantum number
%      - omega_e:  First vibrational constant  (m^-1)
%      - omegaX_e: Second vibrational constant (m^-1)
%      - omegaY_e: Third vibrational constant  (m^-1)
%   Output:
%       - Evib:     Energy due to vibrations   (J)

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

Evib = omega_e.*(v+0.5) - omegaX_e.*(v+0.5).^2 + omegaY_e.*(v+0.5).^3;
Evib = Evib.*h*c;

end