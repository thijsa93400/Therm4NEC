function [chemkin] = WriteChemkin(Sp)
%WRITECHEMKIN writes the thermodynamics in chemkin format
%   WRITECHEMKIN(Sp)
%   - Input:
%       Sp a struct containing species information
%       - Sp.Name        The name of the species
%       - Sp.CompName{i} The name of the elements that species is made of. 
%       - Sp.Comp(i)     how often element i occurs in the species.
%       - Sp.Trange(1)   The lower temperature bound
%       - Sp.Trange(2)   The higher temperature bound
%       - Sp.Trange(3)   The switch temperature
%       - Sp.pol(1,1:7)  The higher range coeficients
%       - Sp.pol(2,1:7)  The lower range coeficients
%
%   Example use:
%   fid = fopen('thermo.dat',w);
%   Sp(1).Name        = 'O2';
%   Sp(1).CompName{1} = {'O'};
%   Sp(1).Comp(1)     = 2
%   Sp(1).Trange      = [200; 6000; 1000];
%   Sp(1).pol(1,:)    = [3.7825 -0.0030 9.8473e-06 -9.6813e-09 3.2437e-12 -1.0639e+03 3.6577];
%   Sp(1).pol(2,:)    = [3.6610 6.5637e-04 -1.4115e-07 2.0580e-11 -1.2991e-15 -1.2160e+03 3.4154];
%   thermo = write_CHEMKIN(Sp(1));
%   for i = 1:4
%       fprintf(fid,"%s\n",thermo(i));
%   end

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

chemkin = char(ones(4,80) * ' ');

% set composition
comp    = char(ones(5,5) * ' ');
for i = 1:4
     comp(i,:) = char(sprintf("%s%3d",'  ',0));
end
for i = 1:length(Sp.CompName)
    % So sprintf is altough used to write string not designed to input
    % string, its buggy:D
    tmp1    = char(Sp.CompName{i});
    if length(tmp1) == 1
        tmp1(2) = ' ';
    end
    comp(i,:) = char(sprintf("%s%3d",tmp1,Sp.Comp(i)));
    clear tmp1
end

% name trickery
name = char(ones(1,18) * ' ');
iend = length(Sp.Name);
name(1:iend) = Sp.Name;

% line 1
ref          = 'EXC 21'; % This is the comment part
chemkin(1,:) = sprintf("%s%s%s%s%s%sG%10.2f%10.2f%8.2f%s 1", name, ref,...
                        comp(1,:), comp(2,:), comp(3,:), comp(4,:), ...
                        Sp.Trange(1), Sp.Trange(2), Sp.Trange(3), comp(5,:));
% coefs                    
coefs    = [Sp.pol(2,:) Sp.pol(1,:)];
chemkin(2,:) = sprintf('%#15.8e%#15.8e%#15.8e%#15.8e%#15.8e    2',coefs(1),coefs(2),coefs(3),coefs(4),coefs(5));
chemkin(3,:) = sprintf('%#15.8e%#15.8e%#15.8e%#15.8e%#15.8e    3',coefs(6),coefs(7),coefs(8),coefs(9),coefs(10));
chemkin(4,:) = sprintf('%#15.8e%#15.8e%#15.8e%#15.8e                   4',coefs(11),coefs(12),coefs(13),coefs(14));

end

