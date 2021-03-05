function [Mol] = ReadMoleculeFile(file)
% READMOLECULEFILE Reads file containing spectroscopic constants
%   READMOLECULEFILE( file)
%   Input:
%       - file: Path to json input file
%   Output:
%       - Mol:  Struct containign molecule information

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

% read file
Mol = jsondecode( fileread(file) );

%% compute mass
Names = fieldnames(Mol.Composition);
Nele  = numel(Names);
Natm  = 0;
mass  = 0;
mu    = 0;
for i = 1:Nele
    N    = Mol.Composition.(Names{i});
    Natm = Natm + N;
    
    % compute mass
    temp = GetElementProperty(Names{i},'mass');
    mass = mass + N*temp;
    mu   = mu   + N/temp;  
end

%% set species properties
Mol.nAtom = Natm;
Mol.Mass  = mass;
Mol.mu    = 1/mu;

if Mol.nAtom == 2 && Nele == 1
    Mol.Symetric = true;
elseif Mol.nAtom == 1
    Nstates = numel(Mol.States);
    for i = 1:Nstates
        % convert cm^-1 => m^-1
        State = Mol.States(i);
        Mol.States(i).Te = State.Te.*100;
    end
    return
else
    Mol.Symetric = false;
end

%% run over states check potential and degeneracy
Nstates = numel(Mol.States);

for i = 1:Nstates
    State = Mol.States(i);
    
    % convert cm^-1 to m^1
    Mol.States(i).Te = State.Te.*100;
    
    Mol.States(i).Vibs.omega_e  = State.Vibs.omega_e.*100;
    Mol.States(i).Vibs.omegaX_e = State.Vibs.omegaX_e.*100;
    Mol.States(i).Vibs.omegaY_e = State.Vibs.omegaY_e.*100;

    Mol.States(i).Rots.B_e      = State.Rots.B_e.*100;
    Mol.States(i).Rots.alpha_e  = State.Rots.alpha_e.*100;
    Mol.States(i).Rots.D_e      = State.Rots.D_e.*100;
    Mol.States(i).Rots.beta_e   = State.Rots.beta_e.*100;
    
    % add the reduced mass, will be handy later on g/mol => kg/mol
    Mol.States(i).mu = Mol.mu*1e-3;
    
    % Compute morse potential parameters if not provided
    if ~isfield(State,'Potential') || isempty(State.Potential)
        [D_e, Beta, R_e] = GetMorsePotential( Mol.States(i).Rots.B_e, Mol.States(i).Vibs.omega_e, Mol.States(i).Vibs.omegaX_e, Mol.States(i).mu);
        Mol.States(i).Potential.D_e  = D_e;
        Mol.States(i).Potential.Beta = Beta;
        Mol.States(i).Potential.R_e  = R_e;
    end
    % Compute degeneracy from, molecular term symbols.
    % This does not account for spin-orbit coupling
    % Moreover degenercy MUST be given for atoms
    if ~isfield(State,'Degen') || isempty(State.Degen)
        [ge,gIo,gIe] = GetDegeneracy(State.S,State.lambda,Mol.Symetric);
        Mol.States(i).Degen.ge  = ge;
        Mol.States(i).Degen.gIo = gIo;
        Mol.States(i).Degen.gIe = gIe;
    end
end

end