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
%% clean workspace
clear all;
clc;
close all;

% ensure path is loaded
addpath(strcat(pwd,'/Functions/'));

%% Run example
% Here we provide an example with a single excited state.

% The new species will use spectroscopic data from oxygen.
% All data has been obtained from: https://webbook.nist.gov/
MolFile = 'Molecules/O2.json';

% We name the species
SpName  = 'O2(1,a,a)';

% the new species contains a1D electronic state included in MolFile
Elecs = 'a^1\Delta_g^-';
% the new species contains all vibrational states
Vibs  = [0 inf];
% the new species contains all rotational states
Rots  = [0 inf];

% A refrence species is required to ensure matching formation enthalpy
% Burcats' database is not included!
% The database can be obtained from: http://garfield.chem.elte.hu/Burcat/burcat.html
% copy the THERM.DAT file to 'Molecules/THERM.DAT'
% IMPORTANT the included chemkin reader requires species without space
RefIn = 'Reference/THERM.DAT';
RefSp = 'O2';

% Compute thermo properties
% write out potential information
[Sp, Spref, Mol] = ThermoMake( MolFile, SpName, Elecs, Vibs, Rots, ...
                               RefIn, RefSp, 'None');
                           
%% write to screen
[chemkin] = WriteChemkin(Sp);
for i = 1:4
    fprintf("%s\n",chemkin(i,:));
end
                           
%% Plot the result
figure
ax = subplot(3,1,[1 2]);
hold on;

% make graph
T = 250:10:6000;
plot(T, CpT( Sp.pol, T, Sp.Tswitch));
plot(T, CpT( Spref.pol, T, Spref.Tswitch));

xlim([250 6000]);
ylabel('$C_p$ [J/mol/K]','Interpreter','latex');
grid on

ax = subplot(3,1,3);
hold on;
plot(T, dHT( Sp.pol, T, Sp.Tswitch));
plot(T, dHT( Spref.pol, T, Spref.Tswitch));

xlim([250 6000]);
xlabel('$T$ [K]','Interpreter','latex');
ylabel('$h$ [J/mol]','Interpreter','latex')
grid on
