function [Sp, Spref, Mol] = ThermoMake( MolFile, SpecName, Elecs, Vibs, ...
                                       Rots, RefFile, Refspec, Debug)
% THERMOMAKE computes the vibrational energy
%   THERMOMAKE( v, omega_e, omegaX_e, omegaY_e)
%   Input:
%      - MolFile:  File containing spectroscopic of molecule
%      - SpecName: Name of the newly defined species
%      - Elecs:    Which electronic states to include see SpeciesPartition
%      - Vibs:     Which vibrational states to include see SpeciesPartition
%      - Rots:     Which rotational states to include see SpeciesPartition
%      - RefFile:  Chemkin format reference file
%      - Refspec:  Which species to compare to from reference file
%      - Debug:    What debug information to print
%   Output:
%       - Sp:      Newly defined species
%       - Spref:   The reference species
%       - Mol:     Molecule information see DiatomicPartionElecState

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

%% Set the temprature range
Tref    = 298;             % temperature at which to match refspec
Tlow    = 200;             % lowest range of nasa poly
Tswitch = 1000;            % switch temperature
Thigh   = 6000;            % highest temperature

T       = Tlow:10:Thigh;   % points at which to evaluate partition

%% Obtain information of reference species
Ru   = 8.31446261815324;

% load species from burcat for hf298
Spref = ReadChemkin(RefFile, Refspec);
Hf298 = dHT(Spref.pol, Tref, Spref.Tswitch);  % Refrence enthalpy

% Compute the same from out method, this includes all internal states
[Zref,Mol]  = SpeciesPartition( Tref, MolFile);
H298        = Partition2Thermo( Tref, Zref, Mol.Mass); % our enthalpy


%% Compute for the new species
[Z,Mol]   = SpeciesPartition( T, MolFile, 'Elecs', Elecs, 'Vibs', Vibs, 'Rots', Rots, 'Debug', Debug);
[H,S,Cp]  = Partition2Thermo(T, Z ,Mol.Mass);
% offset such that H(298) = 0 if alles states are included
H         = H - H298;   

% Fit nasa format polynomial
[pol,Rsqr] = FitNasa( T, Cp, S, H, Tswitch, Hf298);
Sp           = Spref;
Sp.pol       = pol';
Sp.Ts        = Tswitch;
Sp.Trange(1) = T(1);
Sp.Trange(2) = T(end);
Sp.Trange(3) = Tswitch;

%% check results
if any( (1-Rsqr(:,1)) >1e-3)
    ShowError = true;
    fprintf("Warning fit not that good\n");
    fprintf("                Cp 1-Rsqr: %e\n",  1-Rsqr(1,1));
    fprintf("                H  1-Rsqr: %e\n",  1-Rsqr(2,1));
    fprintf("                S  1-Rsqr: %e\n\n",1-Rsqr(3,1));
else
    ShowError = false;
end

if ShowError || false
    figure
    subplot(3,1,1);
    plot( T, (CpT(Sp.pol, T, Sp.Tswitch) - Cp')/Ru);

    subplot(3,1,2);
    plot( T, (dHT(Sp.pol, T, Sp.Tswitch) - (H'+Hf298))./(Ru*T) );

    subplot(3,1,3);
    plot(T, (ST(Sp.pol, T, Sp.Tswitch) - S')/Ru );
end

%% prepare output
% set name
Sp.Name     = SpecName;
% set composition
Names = fieldnames(Mol.Composition);
for i = 1:numel(Names)
    Sp.CompName{i} = Names{i};
    Sp.Comp(i)     = Mol.Composition.(Names{i});
end



end