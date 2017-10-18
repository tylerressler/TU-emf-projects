% MATLAB Project on Electrostatics
% ECE 3712 Electromagnetic Fields & Waves
% Tyler Ressler
%% Section 1: Definition/Initialization
% Plate 1
% [Width] X1
% [Length] Z1
% [Surface Area 1] S1 = X1*Z1;
% [Surface Charge Density S1] rho1 = 0.1 (C/cm^2)
%
% Plate 2
% [Width] X1
% [Height] Y = Y1 + Z1*tan(theta)
% [Surface Area 2] S2 = X1*(Z1*sqrt(1+(tan(theta))^2))
% [Surface Charge Density S2] rho2 = 0.4 (C/cm^2)
% -----------------------------------------------------
%% Section 2: Determine & Test Algorithm for E

clear;
clc;

% Initialize Values

X1 = 1;                                     % [cm]
Z1 = 23;                                    % [cm]
Y1 = ((X1+Z1)/2);                           % [cm]
birthYear = 1997;
theta = degtorad((birthYear - 1985) * 2);             % angle theta [rads]
rho1 = 0.1e2;                                  % [C/m^2]
rho2 = 0.4e2;                                  % [C/m^2]
e0 = 8.854e-12;                              % epsilon

% Variables for text prompt
X1min = 0.2*X1;                             % string var
X1max = 0.7*X1;                             % ""
Y1min = 0.2*Y1;                             % ""
Y1max = 0.7*Y1;                             % ""
Z1min = 0.2*Z1;                             % ""
Z1max = 0.7*Z1;                             % ""

% Display text prompt
disp('Enter the point P(x,y,z) with:');
text = [num2str(X1min), ' <= x <= ', num2str(X1max)];
disp(text);                                 % display  x limits
text = [num2str(Y1min), ' <= y <= ', num2str(Y1max)];
disp(text);                                 % display  y limits
text = [num2str(Z1min), ' <= z <= ', num2str(Z1max)];
disp(text);                                 % display z limits

prompt = 'Enter the point P(x,y,z) [cm] values in the following format "[x y z]":\n';
P = input(prompt);                          % prompt user input for P coordinates

x = P(1);
y = P(2);
z = P(3);

% Conditional Testing of Inputs
if (x <= 0.2*X1) || (x >= 0.7*X1)           % test x val
    t = ['ERROR: Please enter an x value between ',num2str(X1min), ' and ', num2str(X1max), '.'];
    disp(t);
end

if (y <= 0.2*Y1) || (y >= 0.7*Y1)           % test y val
    t = ['ERROR: Please enter an y value between ',num2str(Y1min), ' and ', num2str(Y1max), '.'];
    disp(t);
end

if (z <= 0.2*Z1) || (z >= 0.7*Z1)           % test z val
    t = ['ERROR: Please enter an z value between ',num2str(Z1min), ' and ', num2str(Z1max), '.'];
    disp(t);
end

% Calculations

% Convert cm to m
x = P(1)*(10^-2);
y = P(2)*(10^-2);
z = P(3)*(10^-2);

% Calculate E field by summation of areas for sheet 1
dx = 0.1e-2;           % [m]
dz = 0.1e-2;           % [m]

n = 0;               % initialize counter for row
m = 0;               % initialize counter for column

Ex = 0;              % initialize running E-field sum in x direction
Ey = 0;              % initialize running E-field sum in y direction
Ez = 0;              % initialize running E-field sum in z direction

zdz = dz/2;          % initialize to first point in z direction
xdx = dx/2;


% Sheet 1 EF on P
for n = 1:((X1*(10^-2))/dz)
    for m = 1:((Z1*(10^-2))/dz)
        
        
        % Define vector from P to location on sheet 1
        R = [(x - xdx) (y) (z - zdz)];
        Rx = R(1);
        Ry = R(2);
        Rz = R(3);
        
        % Unit Vector
        ar = R / sqrt((Rx^2) + (Ry^2) + (Rz^2));
        arx = ar(1);
        ary = ar(2);
        arz = ar(3);
        
        Ex = Ex + ((rho1 * (dx*dz))/(4*pi*e0*(Rx^2)))*arx;
        Ey = Ey + ((rho1 * (dx*dz))/(4*pi*e0*(Ry^2)))*ary;
        Ez = Ez + ((rho1 * (dx*dz))/(4*pi*e0*(Rz^2)))*arz;
        
        
        % Increment zdz
        zdz = zdz + dz;
        
    end
    
    zdz = dz/2;
    xdx = xdx + dx;
    
end

E1 = [Ex Ey Ez]

% Sheet 2 EF on P

Ex2 = 0;
Ey2 = 0;
Ez2 = 0;

SA2 = X1*(Z1*sqrt(1+(tan(theta))^2));
Z2 = SA2/X1;
Y2 = Y1 + Z1*tan(theta);

zdz = dz/2;          % initialize to first point in z direction
xdx = dx/2;          % initialize to first point in x direction
ydy = y;             % initialize y height

n = 0;
m = 0;

for n = 1:((X1*(10^-2))/dz)
    for m = 1:((Z1*(10^-2))/dz)
        
        
        % Define vector from P to location on sheet 1
        
        ydy = (Y1 + zdz*tan(theta))*(10^-2);
        
        R = [(x - xdx) (y - ydy) (z - zdz)];
        Rx = R(1);
        Ry = R(2);
        Rz = R(3);
        
        % Unit Vector
        ar = R / sqrt((Rx^2) + (Ry^2) + (Rz^2));
        arx = ar(1);
        ary = ar(2);
        arz = ar(3);
        
        Ex2 = Ex2 + (rho2 * (dx*dz))*arx/(4*pi*e0*(Rx^2));
        Ey2 = Ey2 + (rho2 * (dx*dz))*ary/(4*pi*e0*(Ry^2));
        Ez2 = Ez2 + (rho2 * (dx*dz))*arz/(4*pi*e0*(Rz^2));
        
        % Increment zdz
        zdz = zdz + dz;
        
    end
    zdz = dz/2;
    xdx = xdx + dx;
end

E2 = [Ex2 Ey2 Ez2]

E = E1 + E2

%% Section 3: Find E as function of x

clear;
clc;

% Initialize Values

X1 = 1;                                     % [cm]
Z1 = 23;                                    % [cm]
Y1 = ((X1+Z1)/2);                           % [cm]
birthYear = 1997;
theta = degtorad((birthYear - 1985) * 2);   % angle theta [rads]
rho1 = 0.1e2;                              % [C/m^2]
rho2 = 0.4e2;                              % [C/m^2]
e0 = 8.854e-12;                             % epsilon

EFy = zeros(1,10);                          % initialize vector for E in x dir.
EFx = zeros(1,10);                          % initialize vector for E in y dir.
EFz = zeros(1,10);                          % initialize vector for E in z dir.
xv = zeros(1,10);                           % initialize vector for x
% Calculations

% Calculate E field by summation of areas for sheet 1
dx = 0.1e-2;         % [m]
dz = 0.1e-2;         % [m]

n = 0;               % initialize counter for row
m = 0;               % initialize counter for column

Ex = 0;              % initialize running E-field sum in x direction
Ey = 0;              % initialize running E-field sum in y direction
Ez = 0;              % initialize running E-field sum in z direction

zdz = dz/2;          % initialize to first point in z direction
xdx = dx/2;

y = (Y1*(10^-2))/2;  % y constant
z = (Z1*(10^-2))/2;  % z constant
x = 0.2*(10^-2);   % initialize to first x value

% Sheet 1 EF on P
for q = 1:17
    for n = 1:((X1*(10^-2))/dz)
        for m = 1:((Z1*(10^-2))/dz)
            
            % Define vector from P to location on sheet 1
            R = [(x - xdx) (y) (z - zdz)];
            Rx = R(1);
            Ry = R(2);
            Rz = R(3);
            
            % Unit Vector
            ar = R / sqrt((Rx^2) + (Ry^2) + (Rz^2));
            arx = ar(1);
            ary = ar(2);
            arz = ar(3);
            
            Ex = Ex + (rho1 * (dx*dz))*arx/(4*pi*e0*(Rx^2));
            Ey = Ey + (rho1 * (dx*dz))*ary/(4*pi*e0*(Ry^2));
            Ez = Ez + (rho1 * (dx*dz))*arz/(4*pi*e0*(Rz^2));
            
            % Increment zdz
            zdz = zdz + dz;
            
        end
        
        zdz = dz/2;
        xdx = xdx + dx;
        
    end
    
    E1 = [Ex Ey Ez];
    
    % Sheet 2 EF on P
    
    Ex2 = 0;
    Ey2 = 0;
    Ez2 = 0;
    
    SA2 = X1*(Z1*sqrt(1+(tan(theta))^2));
    Z2 = SA2/X1;
    Y2 = Y1 + Z1*tan(theta);
    
    zdz = dz/2;          % initialize to first point in z direction
    xdx = dx/2;          % initialize to first point in x direction
    ydy = y;             % initialize y height
    
    n = 0;
    m = 0;
    
    for n = 1:((X1*(10^-2))/dz)
        for m = 1:((Z1*(10^-2))/dz)
            
            
            % Define vector from P to location on sheet 1
            
            ydy = (Y1 + zdz*tan(theta))*(10^-2);
            
            R = [(x - xdx) (y - ydy) (z - zdz)];
            Rx = R(1);
            Ry = R(2);
            Rz = R(3);
            
            % Unit Vector
            ar = R / sqrt((Rx^2) + (Ry^2) + (Rz^2));
            arx = ar(1);
            ary = ar(2);
            arz = ar(3);
            
            Ex2 = Ex2 + (rho2 * (dx*dz))*arx/(4*pi*e0*(Rx^2));
            Ey2 = Ey2 + (rho2 * (dx*dz))*ary/(4*pi*e0*(Ry^2));
            Ez2 = Ez2 + (rho2 * (dx*dz))*arz/(4*pi*e0*(Rz^2));
            
            % Increment zdz
            zdz = zdz + dz;
            
        end
        zdz = dz/2;
        xdx = xdx + dx;
    end
    E2 = [Ex2 Ey2 Ez2];
    
    E = E1 + E2;
    
    EFx(q) = E(1);
    EFy(q) = E(2);
    EFz(q) = E(3);
    xv(q) = x;
    
    x = x + ((0.7-0.2)/10)*(10-2); 
end   

% Plot resultant E in each direction
figure(1);
clf;

subplot(3,1,1);
plot(xv, EFx); %x plot
title('Resultant Field in X Direction');
xlabel('x [cm]');
ylabel('Electric Field [N/C]');
grid on; 

subplot(3,1,2);
plot(xv, EFy); %y plot
title('Resultant Field in Y Direction');
xlabel('y [cm]');
ylabel('Electric Field [N/C]');
grid on; 

subplot(3,1,3);
plot(xv, EFz); %z plot
title('Resultant Field in Z Direction');
xlabel('z [cm]');
ylabel('Electric Field [N/C]');
grid on; 


%% Section 4: Find E as function of 

clear;
clc;

% Initialize Values

X1 = 1;                                     % [cm]
Z1 = 23;                                    % [cm]
Y1 = ((X1+Z1)/2);                           % [cm]
birthYear = 1997;
theta = degtorad((birthYear - 1985) * 2);   % angle theta [rads]
rho1 = 0.1e-4;                              % [C/m^2]
rho2 = 0.4e-4;                              % [C/m^2]
e0 = 8.854e-12;                             % epsilon

EFy = zeros(1,10);                          % initialize vector for E in x dir.
EFx = zeros(1,10);                          % initialize vector for E in y dir.
EFz = zeros(1,10);                          % initialize vector for E in z dir.
xv = zeros(1,10);                           % initialize vector for x
% Calculations

% Calculate E field by summation of areas for sheet 1
dx = 0.1e-2;         % [m]
dz = 0.1e-2;         % [m]

n = 0;               % initialize counter for row
m = 0;               % initialize counter for column

Ex = 0;              % initialize running E-field sum in x direction
Ey = 0;              % initialize running E-field sum in y direction
Ez = 0;              % initialize running E-field sum in z direction

zdz = dz/2;          % initialize to first point in z direction
xdx = dx/2;

y = (Y1*(10^-2))/2;  % y constant
x = (X1*(10^-2))/2;   % initialize to first x value
z = (Z1*(10^-2)*0.2);  % z constant

% Sheet 1 EF on P
for q = 1:20
    for n = 1:((X1*(10^-2))/dz)
        for m = 1:((Z1*(10^-2))/dz)
            
            % Define vector from P to location on sheet 1
            R = [(x - xdx) (y) (z - zdz)];
            Rx = R(1);
            Ry = R(2);
            Rz = R(3);
            
            % Unit Vector
            ar = R / sqrt((Rx^2) + (Ry^2) + (Rz^2));
            arx = ar(1);
            ary = ar(2);
            arz = ar(3);
            
            Ex = Ex + (rho1 * (dx*dz))*arx/(4*pi*e0*(Rx^2));
            Ey = Ey + (rho1 * (dx*dz))*ary/(4*pi*e0*(Ry^2));
            Ez = Ez + (rho1 * (dx*dz))*arz/(4*pi*e0*(Rz^2));
            
            % Increment zdz
            zdz = zdz + dz;
            
        end
        
        zdz = dz/2;
        xdx = xdx + dx;
        
    end
    
    E1 = [Ex Ey Ez];
    
    % Sheet 2 EF on P
    
    Ex2 = 0;
    Ey2 = 0;
    Ez2 = 0;
    
    SA2 = X1*(Z1*sqrt(1+(tan(theta))^2));
    Z2 = SA2/X1;
    Y2 = Y1 + Z1*tan(theta);
    
    zdz = dz/2;          % initialize to first point in z direction
    xdx = dx/2;          % initialize to first point in x direction
    ydy = y;             % initialize y height
    
    n = 0;
    m = 0;
    
    for n = 1:((X1*(10^-2))/dz)
        for m = 1:((Z1*(10^-2))/dz)
            
            
            % Define vector from P to location on sheet 1
            
            ydy = (Y1 + zdz*tan(theta))*(10^-2);
            
            R = [(x - xdx) (y - ydy) (z - zdz)];
            Rx = R(1);
            Ry = R(2);
            Rz = R(3);
            
            % Unit Vector
            ar = R / sqrt((Rx^2) + (Ry^2) + (Rz^2));
            arx = ar(1);
            ary = ar(2);
            arz = ar(3);
            
            Ex2 = Ex2 + (rho2 * (dx*dz))*arx/(4*pi*e0*(Rx^2));
            Ey2 = Ey2 + (rho2 * (dx*dz))*ary/(4*pi*e0*(Ry^2));
            Ez2 = Ez2 + (rho2 * (dx*dz))*arz/(4*pi*e0*(Rz^2));
            
            % Increment zdz
            zdz = zdz + dz;
            
        end
        zdz = dz/2;
        xdx = xdx + dx;
    end
    E2 = [Ex2 Ey2 Ez2];
    
    E = E1 + E2;
    
    EFx(q) = E(1);
    EFy(q) = E(2);
    EFz(q) = E(3);
    xv(q) = x;
    
    x = x + ((0.7-0.2)/10)*(10-2); 
end   

% Plot resultant E in each direction
figure(2);
clf;

subplot(3,1,1);
plot(xv, EFx); %x plot
title('Resultant Field in X Direction');
xlabel('x [mm]');
ylabel('Electric Field [N/C]');
xlim([2,7])
grid on; 

subplot(3,1,2);
plot(xv, EFy); %y plot
title('Resultant Field in Y Direction');
xlabel('y [mm]');
ylabel('Electric Field [N/C]');
xlim([2,7])
grid on; 

subplot(3,1,3);
plot(xv, EFz); %z plot
title('Resultant Field in Z Direction');
xlabel('z [mm]');
ylabel('Electric Field [N/C]');
xlim([2,7])
grid on; 