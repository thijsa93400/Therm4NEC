function [H,S,Cp] = Partition2Thermo(T,Z,mass)
% PARITION2THERMO computes thermodynamic properties from partition function
%   PARTITION2THERMO
%   Input:
%       - T:  Temperature        [K]
%       - Z:  Partition function and moments
%       - mass: molecule mass    [g/mol]
%   Ouput:
%       - H:  Enthalpy        H  [J/mol]
%       - S:  Entropy         S  [J/mol/K]
%       - Cp: Heat capacity   Cp [J/mol/K]

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

% read data and copy
if size(T,1) == 1 
    T = T';
end

% Constants
Ru   = 8.31446261815324;
Patm = 1.01325e5;
h    = 6.62606957e-34;
Na   = 6.02214129e23;
u    = 1/Na/1000;
kb   = Ru/Na;

          % Internal
H  = Ru * T .* ( Z(:,2)./Z(:,1) ...
   + 5/2 );
          % Internal
Cp = Ru * ( Z(:,3)./Z(:,1) - (Z(:,2)./Z(:,1)).^2 ...
   + 5/2);
          % Internal
S  = Ru * ( Z(:,2)./Z(:,1) + log(Z(:,1)) ...
   + 5/2 ...
   + log( (2 * pi * mass * u * kb .* T / h.^2 ).^1.5 .* kb .* T ./ Patm) );
end
