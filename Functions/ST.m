function [S]  = ST( x, T, Tswitch)
% ST computes the enthalpy
%   ST(x,T) computes the heat capacity of a species at the temperatures T
%   a single set of nasa-7 coeficients must be provided
%   Input:
%       - x:    Coeficients vector of length 7
%       - T:    Temperature can be a vector
%   Output:
%       - ST:   Entropy as computed from nasa format J/mol
%
%   ST(x,T,Tswitch)  computes the heat capacity of a species at the 
%   temperatures T. Two sets of nasa-7 coeficients must be provided and 
%   the switch temperature
%   Input:
%       - x:       Matrix of 2 length 7 vectors
%       - T:       Temperature can be a vector
%       - Tswitch: Temperature to switch polynomial at
%   Output:
%       - ST:      Entropy as computed from nasa format J/mol

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
Ru = 8.31446261815324;
if exist('Tswitch','var')
    S = zeros(size(T));
    for i = 1:length(S)
        if T(i) < Tswitch
            S(i) = x(1,1)*log(T(i)) + x(1,2)*T(i) + x(1,3)/2*T(i).^2 ...
                 + x(1,4)/3*T(i).^3 + x(1,5)/4*T(i).^4 + x(1,7);
        else
            S(i) = x(2,1)*log(T(i)) + x(2,2)*T(i) + x(2,3)/2*T(i).^2 ...
                 + x(2,4)/3*T(i).^3 + x(2,5)/4*T(i).^4 + x(2,7);
        end
    end
else
    S  = x(1)*log(T) + x(2)*T + x(3)/2*T.^2 + x(4)/3*T.^3 + x(5)/4*T.^4 + x(7);
end
S  = S*Ru;
end