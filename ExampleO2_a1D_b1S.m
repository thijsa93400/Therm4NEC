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
% Here we provide an example with two combined electronic states.

% The new species will use spectroscopic data from oxygen.
% All data has been obtained from: https://webbook.nist.gov/
MolFile = 'Molecules/O2.json';

% We name the species
SpName  = 'O2(1-2,a,a)';

% the new species contains a1D and b1S electronic states
Elecs = [1 2];
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
[Sp(1), Spref, Mol] = ThermoMake( MolFile, SpName, Elecs, Vibs, Rots, ...
                               RefIn, RefSp, 'None');

% a second species with all states but a1D and b1S electronic states
Elecs = [0 3 4 5 6];

% Species name, unlike the article chemkin uses / as negate symbol
% Using / in text was
SpName  = 'O2(/1-2,a,a)';

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
legend('O2(1-2,a,a)','O2(/1-2,a,a)','O2(a,a,a)','Location','SouthEast')

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
