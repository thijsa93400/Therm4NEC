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
% Here we provide an example of a single vibrational band and a range

% The new species will use spectroscopic data from N2.
% All data has been obtained from: https://webbook.nist.gov/
MolFile = 'Molecules/N2.json';

% We name the species
SpName  = 'N2(0,0,a)';

% the new species contains only the ground electronic state
Elecs = [0];
% the new species contains only the lowest vibrational band
Vibs  = [0 0];
% the new species contains all rotational states
Rots  = [0 inf];

% A refrence species is required to ensure matching formation enthalpy
% Burcats' database is not included!
% The database can be obtained from: http://garfield.chem.elte.hu/Burcat/burcat.html
% copy the THERM.DAT file to 'Molecules/THERM.DAT'
% IMPORTANT the included chemkin reader requires species without space
RefIn = 'Reference/THERM.DAT';
RefSp = 'N2';

% Compute thermo properties
% write out potential information
[Sp(1), Spref, Mol] = ThermoMake( MolFile, SpName, Elecs, Vibs, Rots, ...
                               RefIn, RefSp, 'None');

% the new species contains only the ground electronic state
Elecs = [0];
% the new species contains only the remaining vibrational states
Vibs  = [1 inf];

% We name the species
SpName  = 'N2(0,1-,a)';

% Compute thermo properties
% write out potential information                          
[Sp(2), ~, ~] = ThermoMake( MolFile, SpName, Elecs, Vibs, Rots, ...
                          RefIn, RefSp, 'None');
                           
%% write to screen
[chemkin] = WriteChemkin(Sp(1));
for i = 1:4
    fprintf("%s\n",chemkin(i,:));
end

[chemkin] = WriteChemkin(Sp(2));
for i = 1:4
    fprintf("%s\n",chemkin(i,:));
end
                           
%% Plot the result
figure
ax = subplot(3,1,[1 2]);
hold on;

% make graph
T = 250:10:6000;
plot(T, CpT( Sp(1).pol, T, Sp(1).Tswitch));
plot(T, CpT( Sp(2).pol, T, Sp(2).Tswitch));
plot(T, CpT( Spref.pol, T, Spref.Tswitch));
legend('N2(0,0,a)','N2(0,1-,a)','N2(a,a,a)','Location','East')

xlim([250 6000]);
ylabel('$C_p$ [J/mol/K]','Interpreter','latex');
grid on

ax = subplot(3,1,3);
hold on;
plot(T, dHT( Sp(1).pol, T, Sp(1).Tswitch));
plot(T, dHT( Sp(2).pol, T, Sp(2).Tswitch));
plot(T, dHT( Spref.pol, T, Spref.Tswitch));

xlim([250 6000]);
xlabel('$T$ [K]','Interpreter','latex');
ylabel('$h$ [J/mol]','Interpreter','latex')
grid on
