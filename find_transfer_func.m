m = 17.0677e6;          % mass (kg)
Iz = 2.1732e10;         % yaw moment of inertia (kg m^2)
x_g = -3.7;              % CG x-ccordinate (m)
L = 161;                % length (m)
B = 21.8;               % beam (m)
T = 8.9;                % draft (m)
KT = 0.7;               % propeller coefficient (-)
Dia = 3.3;              % propeller diameter (m)
rho = 1025;             % density of water (m/s^3)
U_ref   = 7;            % desired surge speed (m/s)

% rudder coefficients
b = 2;
AR = 8;
CB = 0.8;

lambda = b^2 / AR;
tR = 0.45 - 0.28*CB;
CN = 6.13*lambda / (lambda + 2.25);
aH = 0.75;
xH = -0.4 * L;
xR = -0.5 * L;

% added mass matrix
Xudot = -8.9830e5;
Yvdot = -5.1996e6;
Yrdot =  9.3677e5;
Nvdot =  Yrdot;
Nrdot = -2.4283e10;
MA = -[Yvdot Yrdot; Nvdot Nrdot ];
T1 = 20; T2 = 20; T6 = 10;
Xu = -(m-Xudot)/T1;
Yv = -(m-Yvdot)/T2;
Nr = -(Iz-Nrdot)/T6;
D = -diag([Yv Nr]);

% rigid-body mass matrix
MRB = [ m    m*x_g
        m*x_g Iz ];
Minv = invQR(MRB + MA);

Y_delta = 0.25 * (1 + aH) * rho * AR * CN; 
N_delta = 0.25 * (xR + aH*xH) * rho * AR * CN;   

CRB = [0 m * U_ref
      0 m * x_g * U_ref];

CA = [0 (-Xudot * U_ref)
    (Xudot-Yvdot) * U_ref (- Yrdot * U_ref)];

N = CRB + CA + D;

b = 2 * U_ref * [- Y_delta, - N_delta]';
A = - Minv * N
B = Minv * b
C = [0 1]
Dd = 0;
[num, denom] = ss2tf(A,B,C,Dd);

s = tf(num, denom); 
