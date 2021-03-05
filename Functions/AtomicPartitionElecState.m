function [Z] = AtomicPartitionElecState( T, State, Debug )
% ATOMICPARITIONELECSTATE computes Parition function and 
% its moments of a single potential
%	ATOMICPARITIONELECSTATE(T, State, Debug)
%   Input:
%       - T:             Temperature range of interest
%       - State.Elec:
%           - .Te:       Electronic energy (m^-1)
%           - .S:        Spin quantum number (not used)
%           - .L         Orbital momentum quantum number (not used)
%           - .J         Total angular momentum (not used)
%       - State.Degen:
%           - .ge        Degnercy of electronic state
%       - Debug:
%           - .Pot       Debug information for atomic state
%   Output:
%       - Z(:,1):        Parition function
%       - Z(:,2):        1st moment of partition function
%       - Z(:,3):        2nd moment of partition function

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
%% constants
h  = 6.62607004e-34;  % (m^2kg/s)
c  = 299792458;       % (m/s)
kb = 1.38064852e-23;  % (m^2 kg s^-2 K^-1)
eV = 1.602176565e-19; % J/eV

% preallocate
Z     = zeros(size(T,2),3);


% Potenital
if Debug.Pot
    fprintf("State            : %s\n",State.Name);
    fprintf("Zero point energy: %f\n",(Eelec)./eV);
    fprintf("Degenercy        : %f\n", State.Degen.ge);
end

% compute
beta   = (State.Te.*h.*c)./(kb.*T');
ebeta  = exp(-beta);
Z(:,1) = State.Degen.ge .* ebeta;
Z(:,2) = State.Degen.ge .* beta .* ebeta;
Z(:,3) = State.Degen.ge .* beta.^2 .* ebeta;

end

