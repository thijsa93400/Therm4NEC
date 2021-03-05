function [Z] = DiatomicPartionElecState( T, State, Range, Debug)
% DIATOMICPARTIONELECSTATE computes Parition function and its moments 
% of a single potential
%   DIATOMICPARTIONELECSTATE( T, State, Range, Debug)
%   Input:
%       - T:             Temperature range of interest
%       - State.Elec:
%           - .Te:       Electronic energy (m^-1)
%           - .S:        Spin quantum number
%           - .lambda    Orbital momentum quantum number
%           - .mu        Reduced mass (g/mol)
%       - State.Vibs:
%           - .omega_e:  First vibrational constant (m^-1)
%           - .omegaX_e: Second vibrational constant (m^-1)
%           - .omegaY_e: Third vibrational constant (m^-1)
%       - State.Rots:
%           - .B_e:      Rotational constant at equilibrium (m^-1)
%           - .alpha_e:  Rotational constant first term (m^-1)
%           - .D_e:      Centrifugal distortion (m^-1)
%           - .beta_e:   Centrifugal distortion first term(m^-1)
%       - State.Degen:
%           - .ge        Electronic state degeneracy
%           - .gIo       Odd angular momentum degeneracy
%           - .gIe       Even angular momentum degeneracy
%       - Range:
%           - .Vibs      Vibrational levels to include
%           - .Rots      Rotational levels to include
%       - Debug:
%           - .Vibs      Debug information for vibrational levels
%           - .Rots      Debug information for rotational levels
%           - .Pot       Debug information for potential
% Output:
%   - Z(:,1):        Parition function
%   - Z(:,2):        1st moment of partition function
%   - Z(:,3):        2nd moment of partition function

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


%% constants
h  = 6.62607004e-34;  % (m^2kg/s)
c  = 299792458;       % (m/s)
kb = 1.38064852e-23;  % (m^2 kg s^-2 K^-1)
eV = 1.602176565e-19; % J/eV
Ang= 1e-10;           % m/ang

%% Compute partition function
Eelec = State.Te * h * c;
% The minimum and maximum radius at which the rotational barrier should reside
Rmin  = 0.01/State.Potential.Beta + State.Potential.R_e;  % botom + offset
Rmax  = 10/State.Potential.Beta + State.Potential.R_e;    % (D_e - V) = exp(-10)
Z     = zeros(size(T,2),3);


% compute the vibrational contribution to the zero-point energy
Evib0 = VibrationalEnergy( 0                   , ...
                           State.Vibs.omega_e  , ...
                           State.Vibs.omegaX_e , ...
                           State.Vibs.omegaY_e );
                       
% compute the rotational contribution to the zero-point energy
Erot0 = RotationalEnergy( 0                   , ...
                          0                   , ...
                          State.Rots.B_e      , ...
                          State.Rots.alpha_e  , ...
                          State.Rots.D_e      , ...
                          State.Rots.beta_e   );

% Potenital
if Debug.Pot
    fprintf("State            : %s\n",State.Name);
    fprintf("Zero point energy: %f\n",(Evib0+Erot0+Eelec)./eV);
    fprintf("Degenercy        : %f\n", State.Degen.ge);
    fprintf("Morse potential:\n");
    fprintf("             De: %e\n",   State.Potential.D_e./eV );
    fprintf("           Beta: %e\n",   State.Potential.Beta.*Ang);
    fprintf("             Re: %e\n\n\n", State.Potential.R_e./Ang);
end

%% Compute the partition function for single electronic state
Evib_old = 0;
v        = 0;
% loop over v states
while true
    Evib = VibrationalEnergy( v                   , ...
                               State.Vibs.omega_e  , ...
                               State.Vibs.omegaX_e , ...
                               State.Vibs.omegaY_e );
                           
    % check if vibrational level should be excluded
    if (v < Range.Vibs(1))
        v = v + 1;
        continue;
    elseif (v > Range.Vibs(2))
        break;
    end

    % break if Evib above dissociation
    if Evib > State.Potential.D_e
        break;
    % break if taylor expension reverses
    elseif Evib < Evib_old
        if Debug.Vibs
            fprintf("Reversal of Evib\n");
        end
        break;
    end
    if Debug.Vibs
        fprintf("v: %d\n",v)
    end

    Erot_old = 0;
    J        = 0;
    
    % loop over J states
    while true
        Erot = RotationalEnergy( J                   , ...
                                 v                   , ...
                                 State.Rots.B_e      , ...
                                 State.Rots.alpha_e  , ...
                                 State.Rots.D_e      , ...
                                 State.Rots.beta_e   );
                             
        % check skip
        if (J < Range.Rots(1))
            J = J + 1;
            continue;
        elseif (J > Range.Rots(2))
            break;
        end

        if Erot + Evib > State.Potential.D_e
            
            % define minimization function
            V = @(x) -Potential( x                    , ...
                                 J                    , ...
                                 State.Potential.D_e  , ...
                                 State.Potential.Beta , ...
                                 State.Potential.R_e  , ...
                                 State.mu             );
                             
            % use bounded min search from 
            r_max = fminbnd(V,Rmin,Rmax);
            Emax  = Potential( r_max                , ...
                               J                    , ...
                               State.Potential.D_e  , ...
                               State.Potential.Beta , ...
                               State.Potential.R_e  , ...
                               State.mu             );

            % break if Erot+Evib > Emax or reversal
            if (Erot+Evib) > Emax
                break;
            elseif Erot < Erot_old
                if Debug.Rots
                    fprintf("Reversal of Erot\n");
                end
                break;
            end
        elseif Erot < Erot_old
            % In case reversal occurs before Erot + Evib > D_e
            if Debug.Rots   
                fprintf("Reversal of Erot\n");
            end
            break;
        end
        
        % even/odd degenercy is implemented but not used in the article
        % at the temperatures of interst its effect on thermodynamic 
        % properties is insignificant.
        beta  = (Erot + Evib + Eelec)./(kb.*T');
        ebeta = exp(-beta);
        gJ    = (2*J+1);
        if (mod(J,2) == 0) %even
            Z(:,1) = Z(:,1) + gJ .* State.Degen.ge .* State.Degen.gIe .* ebeta;
            Z(:,2) = Z(:,2) + gJ .* State.Degen.ge .* State.Degen.gIe .* beta .* ebeta;
            Z(:,3) = Z(:,3) + gJ .* State.Degen.ge .* State.Degen.gIe .* beta.^2 .* ebeta;
        else               %uneven
            Z(:,1) = Z(:,1) + gJ .* State.Degen.ge .* State.Degen.gIo .* ebeta;
            Z(:,2) = Z(:,2) + gJ .* State.Degen.ge .* State.Degen.gIo .* beta .* ebeta;
            Z(:,3) = Z(:,3) + gJ .* State.Degen.ge .* State.Degen.gIo .* beta.^2 .* ebeta;
        end
    
        % update quantum number and stored values
        Erot_old = Erot;
        J        = J + 1;
    end
    if Debug.Rots
        fprintf("Jmax: %d\n",J);
    end
    % update quantum number and stored values
    Evib_old = Evib;
    v        = v + 1;
end
if Debug.Vibs
    fprintf("Vmax: %d\n",v)
end

end