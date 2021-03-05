function [D_e, Beta, R_e] = GetMorsePotential( B_e, omega_e, omegaX_e, ...
                                               mu, D_e)
% GETMORSEPOTENTIAL construct morse potential from spectroscopic constants
%   GETMORSEPOTENTIAL( B_e, omega_e, omegaX_e, mu, D_e)
%   Input:
%       - B_e:      Rotational constant at equilibrium (m^-1)
%       - omega_e:  First vibrational constant (m^-1)
%       - omegaX_e: Second vibrational constant (m^-1)
%       - mu:       Reduced mass (kg/mol)
%   Output:
%       - D_e:      Disocation energy measured from equilibrium (J)
%       - beta:     Anharmonicity (m^-1)
%       - R_e:      Equilibrium distance (m)

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
Na   = 6.022140857e23;   % molecules/mol
u    = 1/Na;             % mol/molecule
h    = 6.62607004e-34;   % (m^2kg/s)
c    = 299792458;        % (m/s)
v0   = omega_e.*c;       % m^-1 > 1/s
v0Xe = omegaX_e.*c;      % m^-1 > 1/s

% Use rotational constant to estimate equilibrium
R_e  = sqrt( h / (8*pi^2*c*mu*B_e*u) );

% Extrapolate to disocation
if nargin == 4
    D_e  = h*v0^2/(4*v0Xe);
end
    
% Use vibrational constants + equilibrium to estimate anharmonicity
Beta = v0*pi*sqrt(2*mu/D_e);

end