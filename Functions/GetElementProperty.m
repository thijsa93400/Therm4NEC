function [val] = GetElementProperty(name,prop)
% GETELEMENTPROPERTY returns property of element with name
%   GETELEMENTPROPERTY(name, prop)
%   Input:
%       - name: Name of the element
%       - prop: Property to be return
%   Output:
%       - val:  Requested property

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

%% If more elements are needed add to this table
% currently the nuclear spin is not used
Table{1}.Symbol = 'H';  Table{1}.Name = 'hydrogen';   Table{1}.Mass = 1.00784;  Table{1}.Spin = 1/2;
Table{2}.Symbol = 'He'; Table{2}.Name = 'helium';     Table{2}.Mass = 4.002602; Table{2}.Spin = 0;
Table{3}.Symbol = 'Li'; Table{3}.Name = 'lithium';    Table{3}.Mass = 6.94;     Table{3}.Spin = 3/2;
Table{4}.Symbol = 'Be'; Table{4}.Name = 'beryllium';  Table{4}.Mass = 9.012183; Table{4}.Spin = 3/2;
Table{5}.Symbol = 'B';  Table{5}.Name = 'boron';      Table{5}.Mass = 10.81;    Table{5}.Spin = 5/2;
Table{6}.Symbol = 'C';  Table{6}.Name = 'carbon';     Table{6}.Mass = 12.0107;  Table{6}.Spin = 0;
Table{7}.Symbol = 'N';  Table{7}.Name = 'nitrogen';   Table{7}.Mass = 14.0067;  Table{7}.Spin = 1;
Table{8}.Symbol = 'O';  Table{8}.Name = 'oxygen';     Table{8}.Mass = 15.999;   Table{8}.Spin = 0;
% ....

if strcmp(prop,'mass')
    for i = 1:numel(Table)
        if strcmpi(name,Table{i}.Symbol)
            val = Table{i}.Mass;
            return
        elseif strcmpi(name,Table{i}.Name)
            val = Table{i}.Mass;
            return
        end
    end
    error("Element not stored in the table")
elseif strcmp(prop,'name')
    for i = 1:numel(Tabel)
        if strcmpi(name,Table{i}.Symbol)
            val = Table{i}.Name;
            return
        end
    end
elseif strcmp(prop,'number')
    for i = 1:numel(Tabel)
        if strcmpi(name,Table{i}.Symbol)
            val = i;
            return
        elseif strcmpi(name,Table{i}.Name)
            val = i;
            return
        end
    end   
    error("Element not stored in the table")
elseif strcmp(prop,'spin')
    for i = 1:numel(Tabel)
        if strcmpi(name,Table{i}.Symbol)
            val = Table{i}.Spin;
            return
        elseif strcmpi(name,Table{i}.Name)
            val = Table{i}.Spin;
            return
        end
    end
    error("Element not stored in the table")
end

end