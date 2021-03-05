function [ge,gIo,gIe] = GetDegeneracy(S,lambda,symetric)
% GETDEGENERACY used to get the nuclear and rotational degenercy
%   GETDEGENERACY(S, lambda, symetric)
%   Input:
%       - S:        Spin
%       - lambda:   Orbital Angular momentum
%       - symetric: true if molecule is symetric
%   Output:
%       - ge:       Electronic degenercy
%       - gIo:      odd nuclear degenercy
%       - gIe:      even nuclear degenercy

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

% For now we ignore the nuclear degenercy
if lambda ~= 0
    ge  = (2*S + 1)*2;
    gIo = 1;
    gIe = 1;
else
    ge  = (2*S + 1);
    gIo = 1;
    gIe = 1;
end
if symetric % Homonuclear diatomics
    ge  = 0.5*ge;
end

end