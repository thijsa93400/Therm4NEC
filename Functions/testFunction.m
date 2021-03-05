%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020-2021 Eindhoven University of Technology.             %
%                                                                         %
% This code is free software, you can redistribute it and/or modify it    %
% under the terms of the GNU General Public License; either version 3.0   %
% of the License, or (at your option) any later version. See LICENSE.md   %
% for details, or see <https://www.gnu.org/licenses/>                     %
%                                                                         %
% This code is distributed in the hope that it will be useful, but        %
% WITHOUT ANY WARRANTY; without even the implied warranty of              %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    %
%                                                                         %
% Author:  Thijs Hazenberg (t.hazenberg@tue.nl)                           %
%                                                                         %
% Contact: Thijs Hazenberg  (t.hazenberg@tue.nl)                          %
%          Jesper Janssen   (j.f.janssen@tue.nl)                          %
%          Jan van Dijk     (j.v.dijk@tue.nl)                             %
%          Jeroen van Oijen (j.a.v.oijen@tue.nl)                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi = testFunction(k,j,T)
% TESTFUNCTION computes the test function phi for the least squares fitting
%   TESTFUNCTION( k, j, T)
%   Input:
%       - k:   k'th test function (can be 0 to 2)
%       - j:   j'th coeficient    (can be 0 to inf)
%       - T:   Temperature        (scaled)
%   Output:
%       - phi: contribution
%
% Note: The test function has to be extend for the 9 coefficient NASA-format
% see: FitNasa.m for more information

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

if k == -1 % derivative of Cp/R
    if j == 0
        phi = 0;
    else
        phi = j*T^(j-1);
    end
elseif k == 0  % In Cp/R equation
    phi = T.^j;
elseif k == 1 % In H/RT equation
    phi = T.^j/(j+1); 
elseif k == 2
    if j == 0
        phi = log(T);
    else
        phi = T.^j/j;
    end
else
    error("Error: k can not be larger than 2");
end

end