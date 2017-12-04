%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Tyler Ressler
% Date Created: Nov. 16, 2017
% Course: Electromagnetic Fields & Waves (ECE 3712), Temple University - Fall 2017
% Magnetostatics Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 1: Variable Definitions
% a1            z-direction dimension of solenoid (cm)
% a2            y-direction dimension of solenoid (cm)
% N             turns
% b             length of solenoid in x-direction (cm)
% p             pitch = b/N
% u             permeability = u0 (in air)
% I             current - last 3 digits of TUID (326) in mA
% dL            differential length of coil
% i,j,k,m,n     loop index
% H             resultant H-field
% Htot          total H-field
% loc           location of point on turn
% R             distance vector from point to turn
% aR            unit vector from point on turn to test point
% Rmag          magnitude of R
%% Section 2: Variable Initialization & Flags

%clear;clc;

% Initialize Values
a1 = 23;                                    % [cm]
a2 = 1;                                     % [cm]
b = ((a1+a2)/2);                            % [cm]
birthYear = 1997;
N = birthYear - 1900;                       % [turns]
I = 326e-3;                                 % [A]
p = (b/N)*(10^-2);                          % [m]

% FLAG VARIABLES
% Variables for text prompt
Xmin = -2*b;                              % string var
Xmax = 3*b;                               % ""
Ymin = -a2;                               % ""
Ymax = 2*a2;                              % ""
Zmin = -a1;                               % ""
Zmax = 2*a1;                              % ""

% Display text prompt
disp('Enter the point P(x,y,z) with:');
text = [num2str(Xmin), ' <= x <= ', num2str(Xmax)];
disp(text);                                 % display  x limits
text = [num2str(Ymin), ' <= y <= ', num2str(Ymax)];
disp(text);                                 % display  y limits
text = [num2str(Zmin), ' <= z <= ', num2str(Zmax)];
disp(text);                                 % display z limits

prompt = 'Enter the point P(x,y,z) [cm] values in the following format "[x y z]":\n';
P = input(prompt);                          % prompt user input for P coordinates

x = P(1);
y = P(2);
z = P(3);

% Conditional Testing of Inputs
if (x <= Xmin) || (x >= Xmax)           % test x val
    t = ['ERROR: Please enter an x value between ',num2str(Xmin), ' and ', num2str(Xmax), '.'];
    disp(t);
end

if (y <= Ymin) || (y >= Ymax)           % test y val
    t = ['ERROR: Please enter an y value between ',num2str(Ymin), ' and ', num2str(Ymax), '.'];
    disp(t);
end

if (z <= Zmin) || (z >= Zmax)           % test z val
    t = ['ERROR: Please enter an z value between ',num2str(Zmin), ' and ', num2str(Zmax), '.'];
    disp(t);
end


%% Section 3: H Field Calculation

% Unit conversion cm to m
x = x*(10^-2);
y = y*(10^-2);
z = z*(10^-2);
a1 = a1*(10^-2);
a2 = a2*(10^-2);
b = b*(10^-2);

dL = 0.1e-1;                              % [m]
LdL = dL/2;                               % [m]
H = 0;                                    % H-field initialization
turn = 0;
loc = [0 LdL 0];                          % initialize location on turn
Htot = 0;                                 % initialize total H-field

for i=1:N                               % iterate through each turn
       
    for j=1:(a2/dL)            % iterate through side 1
        
        % Define vector from P to location on side 1
       
        R = [(x - loc(1)) (y-loc(2)) (z - loc(3))];
        Rx = R(1);
        Ry = R(2);
        Rz = R(3);
        
        aR = R / sqrt((Rx^2)+(Ry^2)+(Rz^2));
        Rmag = sqrt((Rx^2)+(Ry^2)+(Rz^2));  % magnitude of R
        dLaL = dL*[0 1 0];
        
        H = (I*cross(dLaL, aR)) / (4*pi*(Rmag)^2);
        Htot = Htot + H; 
        loc(2) = loc(2) + dL;
    end
    
    % initialize y & z locations to top left corner of turn
    loc(2) = a2;
    loc(3) = LdL;
    
    for k=1:(a1/dL)             % iterate through side 2
        
        R = [(x - loc(1)) (y-loc(2)) (z - loc(3))];
        Rx = R(1);
        Ry = R(2);
        Rz = R(3);
        
        aR = R / sqrt((Rx^2)+(Ry^2)+(Rz^2));
        Rmag = sqrt((Rx^2)+(Ry^2)+(Rz^2));  % magnitude of R
        dLaL = dL*[0 0 1];
        
        H = (I*cross(dLaL, aR)) / (4*pi*(Rmag)^2);
        Htot = Htot + H; 
        loc(3) = loc(3) + dL;
        
    end
    
    % initialize y & z locations to top right corner of turn
    loc(2) = a2-LdL;
    loc(3) = a1;                         

    for m=1:(a2/dL)            % iterate through side 3
        
        R = [(x - loc(1)) (y-loc(2)) (z - loc(3))];
        Rx = R(1);
        Ry = R(2);
        Rz = R(3);
        
        aR = R / sqrt((Rx^2)+(Ry^2)+(Rz^2));
        Rmag = sqrt((Rx^2)+(Ry^2)+(Rz^2));  % magnitude of R
        dLaL = dL*[0 -1 0];
        
        H = (I*cross(dLaL, aR)) / (4*pi*(Rmag)^2);
        Htot = Htot + H; 
        
        loc(2) = loc(2) - dL;   
    end
    
    % initialize y & z locations to bottom right corner of turn
    loc(2) = 0;
    loc(3) = a1;
    
    for n=1:(a1/dL)             % iterate through side 4
        
        R = [(x - loc(1)) (y-loc(2)) (z - loc(3))];
        Rx = R(1);
        Ry = R(2);
        Rz = R(3);
        
        aR = R / sqrt((Rx^2)+(Ry^2)+(Rz^2));
        Rmag = sqrt((Rx^2)+(Ry^2)+(Rz^2));  % magnitude of R
        dLaL = dL*[0 0 -1];
        
        H = (I*cross(dLaL, aR)) / (4*pi*(Rmag)^2);
        Htot = Htot + H; 
        
        loc(3) = loc(3) - dL;  
        
    end
    
    turn = turn+p;
    loc = [turn LdL 0];       % initialize location for next turn turn
end

% Display final answer
format long;
Htot 

figure(2);semilogx(delL*10,HdelL*(10^-4),'-or', 'LineWidth', 3);
xlabel('\Delta L [m]')
title('Magnetic Field Intensity vs. \Delta L')
ylabel('H-field Intensity [ A/m ]')
