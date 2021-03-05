function [R2,SSE] = Rsquared(y, y_est)
% RSQUARED function to compute r-squared and sum of squared errors
%   RSQUARED( y, y_est)
%   Input:
%       - y:     measurements
%       - y_est: estimations of measurements
%   Ouput:
%       - R2:    R-squared
%       - SSE:   Sum 

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

r     = y-y_est;
normr = norm(r);
SSE   = normr.^2;
SST   = norm(y-mean(y))^2;
R2    = 1 - SSE/SST;
end