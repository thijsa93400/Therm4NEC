function [Therm] = ReadChemkin(fname,spName)
%% read thermo data from senkin format
% Input:
%   - fname: name of thermo file
%   - spName: name of species to seek
% Output:
%   -Therm: thermo structure

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

% %% Setup an struct of element masses
% el(1).Name   = "H ";
% el(1).Mass   = 1.00797;
% 
% el(2).Name   = "O ";
% el(2).Mass   = 15.9994;
% 
% el(3).Name   = "N ";
% el(3).Mass   = 14.0067;
% 
% el(4).Name   = "C ";
% el(4).Mass   = 12.0107;

%% open file
fid   = fopen(fname);

Therm.Trange = zeros(3,1);
Therm.lcoefs = zeros(7,1);
Therm.hcoefs = zeros(7,1);
Therm.Mass   = 0;
coefs        = zeros(14,1);

%% seek for the species
iline = 0;
while true
    tline = fgetl(fid);
    iline = iline + 1;
    if ~ischar(tline)
        return
    end
    C     = strsplit(tline,' ');
    if strcmp(C{1},spName)
        found = true;
        break;
    end
end
Therm.Name      = C{1};

%% record the range
Therm.Trange(1) = str2double(tline(46:55));
Therm.Trange(2) = str2double(tline(56:65));
Therm.Trange(3) = str2double(tline(66:75));

Lines{1} = tline;
Lines{2} = fgetl(fid);
Lines{3} = fgetl(fid);
Lines{4} = fgetl(fid);

k = 0;
for j = 1:3
    for i = 1:5
        istart = 15*(i-1) + 1;
        iend   = istart + 15 - 1;
        k      = k + 1;
        if (k > 14)
            break;
        else
            coefs(k) = str2double( Lines{j+1}(istart:iend) );
        end
    end
end
Therm.pol(1,:) = coefs(8:14);
Therm.pol(2,:) = coefs(1:7);
Therm.iline    = iline;
Therm.Tswitch  = Therm.Trange(3);

end

