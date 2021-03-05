function [pol,Rsqr] = FitNasa(T_in, Cp_in, S_in, dH_in, T_switch, Hf298)
% FITNASA fits a nasa format polynomial to thermodynamic data
%   FITNASA(T_in, Cp_in, S_in, dH_in, T_switch, Hf298)
%   Input:
%       - T_in:      vector of temperatures in K
%       - Cp_in:     vector of heat capacities in J/mol/K
%       - S_on:      vector of entropies in J/mol
%       - dH_in:     vector of enthalpies in J/mol
%       - T_switch:  Temperature to switch polynomial in K
%       - Hf298:     Enthalpy at 298K
%   Output:
%       - pol(1:7,1) Lower range polynomial
%       - pol(1:7,2) Higher range polynomial
%       - Rsqr(1,1)  R-squared of Cp
%       - Rsqr(1,2)  SEE of CP
%       - Rsqr(2,1)  R-squared of Enthalpy
%       - Rsqr(2,2)  SEE of Enthalpy
%       - Rsqr(3,1)  R-squared of Entropy
%       - Rsqr(3,2)  SEE of enthalpy
%
% Explanation of algorithm
% A least square method is used on all 3 variables simultaneously.
% The simultaneously fitting procedure is based on: 
%  F. J. Zeleznik and S. Gordon, 'Simultaneous equations least squares 
%    approximation of a function and its first integrals with application 
%    to thermodynamic data'
%
% Both the lower and higher range are fitted at the same time.
% Using Lagrangian multipliers the upper and lower ranges are matched.
% By default the algorithm sets 4 constraints:
% 1) Cp_low(T_switch) - Cp_high(T_switch) = 0
% 2) d/dT [Cp_low(T_switch)] - d/dT [Cp_high(T_switch)] = 0
% 3) H_low(T_switch) - H_high(T_switch) = 0
% 4) S_low(T_switch) - S_high(T_switch) = 0
% Note, this implies enthalpy is twice differentiable
%
%  Mx = b
%
%  With M symetric:
%      | A_1  BA_1 0    0     LA_1 |
%      | AB_1 B_1  0    0     LB_1 |
%  M = | 0    0    A_2  BA_2  LA_2 |
%      | 0    0    AB_2 B_2   LB_2 |
%      | A_1L B_1L A_2L B_2L  0    |
%
%  and x is solution vector
%      | a_1 |
%      | b_1 |
%  x = | a_2 | 
%      | b_2 | 
%      | l   |
%  _1 indicates for the lower range
%  _2 indicates for the higher range
%  L matrices are for the lagrangian multipliers

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

%% number of measurement points
m = length(T_in);

%% Initialize sub-matrices and sub-vetros
A  = zeros(5,5,2);
AB = zeros(2,5,2);  
B  = zeros(2,2,2);
cA = zeros(5,2);
cB = zeros(2,2);

% constraints
icont = 2;           % c0 and c1 continity on  Cp
icons = icont + 2;
l     = zeros(icons,1);
AL    = zeros(icons,5,2);
BL    = zeros(icons,2,2);
L     = zeros(icons,icons);

%% Scale measurements
T_scale  = 1000;
Ru       = 8.31446261815324;
T_switch = T_switch./T_scale;
T_in     = T_in./T_scale;
dH_in    = dH_in./T_scale;

%% Loop over T vector
for i = 1:m
    % Weights of equations
    c(1) = 1.0;
    c(2) = 1.0e-0;
    c(3) = 1.0e-0;
    
    % local measurements
    T    = T_in(i);
    Cp   = Cp_in(i)/Ru;
    S    = S_in(i)/Ru;
    dH   = dH_in(i)/(Ru*T);
    
    % Reset local contributions to zero
    cphi = zeros(5,3);
    phi  = zeros(3,5);
    rphi = zeros(2,5);
    r2   = zeros(2,2);
    hphi = zeros(5,1);
    rh   = zeros(2,1);
    
    % Obtain contributions
    for k = 1:3
        for j = 1:5
            phi(k,j)  = testFunction(k-1,j-1,T);
            cphi(j,k) = c(k)*phi(k,j);
        end
        if k == 1 % heat cap
            % RHS
            hphi(:)     = hphi + cphi(:,k)*Cp;
        elseif k == 2 % enthalpy
            rphi(k-1,:) = c(k)*phi(k,:)/T;
            r2(k-1,k-1) = c(k)/T.^2;
            
            % RHS
            hphi(:)     = hphi + cphi(:,k)*dH;
            rh(k-1)     = c(k)/T*dH;
        elseif k == 3 % entropy
            rphi(k-1,:) = c(k)*phi(k,:);
            r2(k-1,k-1) = c(k);
            
            
            % RHS
            hphi(:)     = hphi + cphi(:,k)*S;
            rh(k-1)     = c(k)*S;
        end
    end
    % take correct range
    if (T < T_switch)
        irange = 1;
    else 
        irange = 2;
    end
    
    % add to system sub-matrix
    A(:,:,irange)  = A(:,:,irange)  + cphi*phi;
    AB(:,:,irange) = AB(:,:,irange) + rphi;
    B(:,:,irange)  = B(:,:,irange)  + r2;
    
    % add to constants
    cA(:,irange)   = cA(:,irange) + hphi;
    cB(:,irange)   = cB(:,irange) + rh;
end

%% apply lagrangian multipliers
i = 0;
if icont >= 1
    i = i + 1;
    % c0 continuity
    % Cp_low(T_switch) - Cp_high(T_switch) = 0
    for j = 1:5
        AL(i,j,1) = c(1) *  testFunction(0,j-1,T_switch);
        AL(i,j,2) = c(1) * -testFunction(0,j-1,T_switch);
    end
end

if icont >= 2
    i = i + 1;
    % c1 continuity
    % d/dT [Cp_low(T_switch)] - d/dT [Cp_high(T_switch)] = 0
    for j = 1:5
        AL(i,j,1) = c(1) * testFunction(-1,j-1,T_switch);
        AL(i,j,2) = c(1) * -testFunction(-1,j-1,T_switch);
    end
end

if true
    % Constrain H_low(T_switch) = H_high(T_switch)
    i = i + 1;
    for j = 1:5
        AL(i,j,1) = c(2) * testFunction(1,j-1,T_switch);
    end
    for j = 1:5
        AL(i,j,2) = c(2) * -testFunction(1,j-1,T_switch);
    end
    BL(i,1,1) = c(2)/T_switch;
    BL(i,1,2) = -c(2)/T_switch;
end

if true
    % Constrain S_low(T_switch) = S_high(T_switch)
    i = i + 1;
    for j = 1:5
        AL(i,j,1) = c(3) * testFunction(2,j-1,T_switch);
    end
    for j = 1:5
        AL(i,j,2) = c(3) * -testFunction(2,j-1,T_switch);
    end
    BL(i,2,1) = c(3);
    BL(i,2,2) = -c(3);
end

%% construct system matrix
% Construct the sub matrices
BA   = permute(AB,[2 1 3]);
LA   = permute(AL,[2 1 3]);
LB   = permute(BL,[2 1 3]);
A0   = zeros(size(A(:,:,1)));
BA0  = zeros(size(BA(:,:,1)));
AB0  = BA0';
B0   = zeros(size(B(:,:,1)));
  
SyS = [A(:,:,1)  BA(:,:,1)  A0        BA0       LA(:,:,1) ;  ...
       AB(:,:,1) B(:,:,1)   AB0       B0        LB(:,:,1) ;  ...
       A0        BA0        A(:,:,2)  BA(:,:,2) LA(:,:,2) ;  ...
       AB0       B0         AB(:,:,2) B(:,:,2)  LB(:,:,2) ;  ...
       AL(:,:,1) BL(:,:,1)  AL(:,:,2) BL(:,:,2) L         ];
const = [cA(:,1)   ; ...
         cB(:,1)   ; ...
         cA(:,2)   ; ...
         cB(:,2)   ; ...
         l         ];

x  = linsolve(SyS,const);

%% scale ouput
% low range
pol(1,1) = x(1);
pol(2,1) = x(2) / T_scale;
pol(3,1) = x(3) / (T_scale^2);
pol(4,1) = x(4) / (T_scale^3);
pol(5,1) = x(5) / (T_scale^4);
pol(6,1) = x(6) * T_scale + Hf298/Ru;
pol(7,1) = x(7) - x(1) * log(T_scale);

% high range
pol(1,2) = x(8);
pol(2,2) = x(9) / T_scale;
pol(3,2) = x(10) / (T_scale^2);
pol(4,2) = x(11) / (T_scale^3);
pol(5,2) = x(12) / (T_scale^4);
pol(6,2) = x(13) * T_scale + Hf298/Ru;
pol(7,2) = x(14) - x(8) * log(T_scale);

% test fit quality
T_in     = T_in*T_scale;
T_switch = T_switch*T_scale;
dH_in    = dH_in*T_scale;
Cptest   = [CpT( pol(:,1), T_in(T_in<T_switch) ) CpT( pol(:,2), T_in(T_in>=T_switch) )];
Htest    = [dHT( pol(:,1), T_in(T_in<T_switch) ) dHT( pol(:,2), T_in(T_in>=T_switch) )];
Stest    = [ST( pol(:,1), T_in(T_in<T_switch) ) ST( pol(:,2), T_in(T_in>=T_switch) )];

[Rsqr(1,1),Rsqr(1,2)]  = Rsquared( Cp_in/Ru                 , Cptest'/Ru);
[Rsqr(2,1),Rsqr(2,2)]  = Rsquared( (dH_in+Hf298)./(Ru*T_in) , Htest'./(Ru*T_in));
[Rsqr(3,1),Rsqr(3,2)]  = Rsquared( S_in/Ru                  , Stest'/Ru);
   
end