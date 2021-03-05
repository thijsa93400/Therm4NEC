function [Z,Mol] = SpeciesPartition(T, File, varargin)
% SPECIESPARTITION compute the partition function of a series of potentials
%   SPECIESPARTITION( T, File, vararging)
%   Input:
%       - T:             Temperature range of interest
%       - File           Molecule input file
%   Output:
%       - Z(:,1):        Parition function
%       - Z(:,2):        1st moment of partition function
%       - Z(:,3):        2nd moment of partition function
%
%   List of key value pairs:
%   SPECIESPARTITION( T, File, 'Debug', Debug)
%   Input:
%       - Debug information to provide options are:
%           'Rots' for debug of rotational levels
%           'Vibs' for debug of vibrational levels
%           'Potential' for debug of potential
%           'All'  Enable all of the above debug info
%           'ReturnMol' Only read molecule file and return
%   
%   SPECIESPARTITION( T, File, 'Elecs', Elecs)
%       - Electronic states to include, options are:
%           ['a' 'b' ...] in which 'a' and 'b' are term symbols
%           [i j ...]     in which 0 and 1 are state indexes
%
%   SPECIESPARTITION( T, File, 'Vibs', Vibs)
%       - Range of vibrational states to include:
%           [i j] Lowest and highest vibrational level to include
%
%   SPECIESPARTITION( T, File, 'Rots', Rots)
%       - Range of rotational states to include:
%           [i j] Lowest and highest rotational level to include

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

%% Input handeling
p           = inputParser;

% required args
validTrange = @(x) validateattributes(x,{'numeric'},{'increasing','nonnegative','row'});
validStr    = @(x) ischar(x);
addRequired(p,'T',validTrange);
addRequired(p,'File',validStr);

% default params params
defaultElec  = 'All';
defaultVibs  = [0 inf];
defaultRots  = [0 inf];
defaultPrint = 'None';

% debug
expectedPrints = {'None','Rots','Vibs','Potential','All','ReturnMol'};
addParameter( p, 'Debug', defaultPrint, @(x) any(validatestring(x,expectedPrints)));

% Elecs to include
addParameter( p, 'Elecs', defaultElec);

% vibs to include
validVrange = @(x) validateattributes(x,{'numeric'},{'nondecreasing','numel',2});
addParameter( p, 'Vibs', defaultVibs, validVrange);

% Rots to include
addParameter( p, 'Rots', defaultRots, validVrange);

%% parse
parse( p, T, File, varargin{:});

% load
Mol = ReadMoleculeFile(File);
if strcmp(p.Results.Debug,'ReturnMol')
    Z = [];
    return
end


%% parse elecs
Elecs = p.Results.Elecs;
if ischar(Elecs)
    if strcmp(Elecs, 'All')
        % Add all
        e = 1:length(Mol.States);
    else
        % Add single
        for i = 1:numel(Mol.States)
            if strcmp(Elecs,Mol.States(i).Name)
                e = [i];
                break;
            end
        end
    end
elseif iscell(Elecs)
    e = zeros(size(Elecs));
    for j = 1:length(Elecs)
        for i = 1:numel(Mol.States)
            if strcmp(Elecs(j),Mol.States(i).Name)
                e(j) = i;
            end
        end
    end
elseif isvector(Elecs)
    e     = Elecs + 1;
end

%% parse vibs and rots
Range.Vibs = p.Results.Vibs;
Range.Rots = p.Results.Rots;

%% Parse debug
if strcmp(p.Results.Debug, expectedPrints{1})
    Debug.Rots = false;
    Debug.Vibs = false;
    Debug.Pot  = false;
elseif strcmp(p.Results.Debug, expectedPrints{2})
    Debug.Rots = true;
    Debug.Vibs = false;
    Debug.Pot  = false;
elseif strcmp(p.Results.Debug, expectedPrints{3})
    Debug.Rots = false;
    Debug.Vibs = true;
    Debug.Pot  = false;
elseif strcmp(p.Results.Debug, expectedPrints{4})
    Debug.Rots = false;
    Debug.Vibs = false;
    Debug.Pot  = true;
elseif strcmp(p.Results.Debug, expectedPrints{5})
    Debug.Rots = true;
    Debug.Vibs = true;
    Debug.Pot  = true;
end

%% Compute
Z = zeros(size(T,2),3);
for i = 1:length(e)
    if Mol.nAtom == 1
        Z = Z + AtomicPartitionElecState( T, Mol.States(e(i)), Debug);
    else
        Z = Z + DiatomicPartionElecState( T, Mol.States(e(i)), Range, Debug);
    end
end

end