clc;clear;close all;

%input for ratio is (gamma #, Mach #, 'M')
%input for Mach is (gamma #, ratio #)


%T_over_Tt = isentropic_temperature_relation(gamma, 1.4, 'M');
%M = isentropic_temperature_relation(1.4, 0.7);

%P_over_Pt = isentropic_pressure_relation(gamma, 1.4, 'M');
%M = isentropic_pressure_relation(gamma, P_over_Pt);


%A_over_Astar = isentropic_area_relation(gamma, M, 'M');
%M = isentropic_area_relation(gamma, Area Ratio);


function result = isentropic_temperature_relation(varargin)
    if nargin == 2  % If only gamma and T/Tt are provided.
        gamma = varargin{1};
        T_over_Tt = varargin{2};
        % Solve for Mach number
        M = sqrt(2*((1/T_over_Tt) - 1)/(gamma - 1));
        result = M;
    elseif nargin == 3  % If gamma, Mach and a string 'M' are provided.
        gamma = varargin{1};
        M = varargin{2};
        % Calculate T/Tt
        T_over_Tt = 1 / (1 + (gamma - 1)/2 * M^2);
        result = T_over_Tt;
    else
        error('Invalid number of input arguments.');
    end
end


function result = isentropic_pressure_relation(varargin)

    if nargin == 2  % If only gamma and P/Pt are provided.
        gamma = varargin{1};
        P_over_Pt = varargin{2};
        
        % Solve for Mach number using the isentropic pressure relation
        M = sqrt(2*((P_over_Pt^(-((gamma-1)/gamma))) - 1)/(gamma - 1));
        result = M;

    elseif nargin == 3  % If gamma, Mach, and a string 'M' are provided.
        gamma = varargin{1};
        M = varargin{2};
        
        % Calculate P/Pt using the isentropic pressure relation
        P_over_Pt = (1 + (gamma - 1)/2 * M^2)^(-gamma/(gamma-1));
        result = P_over_Pt;
        
    else
        error('Invalid number of input arguments.');
    end
end


function result = isentropic_area_relation(varargin)

    if nargin == 2  % If only gamma and A/A* are provided.
        gamma = varargin{1};
        A_over_Astar = varargin{2};
        
        % Initial guess for Mach number
        M_guess = 0.5;  
        
        % Define function for the area relation
        f = @(M) (1/M) * (2/(gamma + 1) + (gamma - 1)/(gamma + 1) * M^2)^( (gamma + 1) / (2 * (gamma - 1))) - A_over_Astar;
        
        % Solve for Mach number using fsolve
        M = fsolve(f, M_guess);
        result = M;

    elseif nargin == 3  % If gamma, Mach, and a string 'M' are provided.
        gamma = varargin{1};
        M = varargin{2};
        
        % Calculate A/A* using the isentropic area relation
        A_over_Astar = (1/M) * (2/(gamma + 1) + (gamma - 1)/(gamma + 1) * M^2)^( (gamma + 1) / (2 * (gamma - 1)));
        result = A_over_Astar;
        
    else
        error('Invalid number of input arguments.');
    end
end

