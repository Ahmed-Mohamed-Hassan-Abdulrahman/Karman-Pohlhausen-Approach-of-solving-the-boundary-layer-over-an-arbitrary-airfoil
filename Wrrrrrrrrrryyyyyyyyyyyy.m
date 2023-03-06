function[vt,xc]=Wrrrrrrrrrryyyyyyyyyyyy(Airfoil)
%% Givens

% Flag to specify creating or loading airfoil
flagAirfoil.XFoilCreate = 1;                                                % Create specified NACA airfoil in XFOIL
flagAirfoil.XFoilLoad   = 0;                                                % Load Selig-format airfoil from directory

% User-defined knowns
Vinf = 50;                                                                   % Freestream velocity [] (just leave this at 1)
AoA  = 0;                                                                   % Angle of attack [deg]
NACA = Airfoil;                                                              % NACA airfoil to load [####(#)]

% Plotting flags
flagPlot = [1;          % Airfoil with panel normal vectors
            1;          % Geometry boundary pts, control pts, first panel, second panel
            1;          % Cp vectors at airfoil surface panels
            1;          % Pressure coefficient comparison (XFOIL vs. SPVP)
            1;          % Airfoil streamlines
            1];         % Pressure coefficient contour

%% XFOIL - CREATE/LOAD AIRFOIL

% PPAR menu options
PPAR.N  = '402';                                                            % "Number of panel nodes"
PPAR.P  = '4';                                                              % "Panel bunching parameter"
PPAR.T  = '1';                                                              % "TE/LE panel density ratios"
PPAR.R  = '1';                                                              % "Refined area/LE panel density ratio"
PPAR.XT = '1 1';                                                            % "Top side refined area x/c limits"
PPAR.XB = '1 1';                                                            % "Bottom side refined area x/c limits"

% Call XFOIL function to obtain the following:
% - Airfoil coordinates
% - Pressure coefficient along airfoil surface
% - Lift, drag, and moment coefficients
[xFoilResults,success] = XFOIL(NACA,PPAR,AoA,flagAirfoil);                  % Get XFOIL results for prescribed airfoil
if (success == 0)                                                           % If user canceled airfoil dialog box
    return;                                                                 % Exit the program
end

% Separate out results from XFOIL function results
afName  = xFoilResults.afName;                                              % Airfoil name
xFoilX  = xFoilResults.X;                                                   % X-coordinate for Cp result
xFoilY  = xFoilResults.Y;                                                   % Y-coordinate for Cp result
xFoilCP = xFoilResults.CP;                                                  % Pressure coefficient
XB      = xFoilResults.XB;                                                  % Boundary point X-coordinate
YB      = xFoilResults.YB;                                                  % Boundary point Y-coordinate
xFoilCL = xFoilResults.CL;                                                  % Lift coefficient
xFoilCD = xFoilResults.CD;                                                  % Drag coefficient
xFoilCM = xFoilResults.CM;                                                  % Moment coefficient

% Number of boundary points and panels
numPts = length(XB);                                                        % Number of boundary points
numPan = numPts - 1;                                                        % Number of panels (control points)

%% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

% Check for direction of points
edge = zeros(numPan,1);                                                     % Initialize edge value array
for i = 1:1:numPan                                                          % Loop over all panels
    edge(i) = (XB(i+1)-XB(i))*(YB(i+1)+YB(i));                              % Compute edge values
end
sumEdge = sum(edge);                                                        % Sum all edge values

% If panels are CCW, flip them (don't if CW)
if (sumEdge < 0)                                                            % If panels are CCW
    XB = flipud(XB);                                                        % Flip the X-data array
    YB = flipud(YB);                                                        % Flip the Y-data array
end

%% PANEL METHOD GEOMETRY - REF [1]

% Initialize variables
XC   = zeros(numPan,1);                                                     % Initialize control point X-coordinate array
YC   = zeros(numPan,1);                                                     % Initialize control point Y-coordinate array
S    = zeros(numPan,1);                                                     % Initialize panel length array
phiD = zeros(numPan,1);                                                     % Initialize panel orientation angle array [deg]

% Find geometric quantities of the airfoil
for i = 1:1:numPan                                                          % Loop over all panels
    XC(i)   = 0.5*(XB(i)+XB(i+1));                                          % X-value of control point
    YC(i)   = 0.5*(YB(i)+YB(i+1));                                          % Y-value of control point
    dx      = XB(i+1)-XB(i);                                                % Change in X between boundary points
    dy      = YB(i+1)-YB(i);                                                % Change in Y between boundary points
    S(i)    = (dx^2 + dy^2)^0.5;                                            % Length of the panel
    phiD(i) = atan2d(dy,dx);                                                % Angle of the panel (positive X-axis to inside face) [deg]
    if (phiD(i) < 0)                                                        % Make all panel angles positive [deg]
        phiD(i) = phiD(i) + 360;
    end
end

% Compute angle of panel normal w.r.t horizontal and include AoA
deltaD             = phiD + 90;                                             % Angle from positive X-axis to outward normal vector [deg]
betaD              = deltaD - AoA;                                          % Angle between freestream vector and outward normal vector [deg]
betaD(betaD > 360) = betaD(betaD > 360) - 360;                              % Make all panel angles between 0 and 360 [deg]

% Convert angles from [deg] to [rad]
phi  = phiD.*(pi/180);                                                      % Convert from [deg] to [rad]
beta = betaD.*(pi/180);                                                     % Convert from [deg] to [rad]

%% COMPUTE SOURCE AND VORTEX PANEL STRENGTHS - REF [10]

% Geometric integrals for SPM and VPM (normal [I,K] and tangential [J,L])
% - Refs [2], [3], [6], and [7]
[I,J] = COMPUTE_IJ_SPM(XC,YC,XB,YB,phi,S);                                  % Call COMPUTE_IJ_SPM function (Refs [2] and [3])
[K,L] = COMPUTE_KL_VPM(XC,YC,XB,YB,phi,S);                                  % Call COMPUTE_KL_VPM function (Refs [6] and [7])

% Populate A matrix

A = I + pi*eye(numPan,numPan);

% Right column of A matrix
for i = 1:1:numPan                                                          % Loop over all i panels (rows)
    A(i,numPan+1) = -sum(K(i,:));                                           % Add gamma term to right-most column of A matrix
end

% Bottom row of A matrix (Kutta condition)
for j = 1:1:numPan                                                          % Loop over all j panels (columns)
    A(numPan+1,j) = (J(1,j) + J(numPan,j));                                 % Source contribution of Kutta condition equation
end
A(numPan+1,numPan+1) = -sum(L(1,:) + L(numPan,:)) + 2*pi;                   % Vortex contribution of Kutta condition equation

% Populate b array

b = -Vinf*2*pi*cos(beta);

% Last element of b array (Kutta condition)
b(numPan+1) = -Vinf*2*pi*(sin(beta(1)) + sin(beta(numPan)));                % RHS of Kutta condition equation

% Compute result array
resArr = A\b;                                                               % Solve system of equations for all source strengths and single vortex strength

% Separate lambda and gamma values from result array
lambda = resArr(1:end-1);                                                   % All panel source strenths
gamma  = resArr(end);                                                       % Constant vortex strength

%% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS

% Compute velocities on each panel
Vt = zeros(numPan,1);                                                       % Initialize tangential velocity
Cp = zeros(numPan,1);                                                       % Initialize pressure coefficient
for i = 1:1:numPan
    term1 = Vinf*sin(beta(i));                                              % Uniform flow term
    term2 = (1/(2*pi))*sum(lambda.*J(i,:)');                                % Source panel terms when j is not equal to i
    term3 = gamma/2;                                                        % Vortex panel term when j is equal to i
    term4 = -(gamma/(2*pi))*sum(L(i,:));                                    % Vortex panel terms when j is not equal to i
    
    Vt(i) = term1 + term2 + term3 + term4;                                  % Compute tangential velocity on panel i
    Cp(i) = 1-(Vt(i)/Vinf)^2;                                               % Compute pressure coefficient on panel i
end
vt=Vt';
xc=XC';
end