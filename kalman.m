% Nomoto parameters (use the values you identified)
K = -0.004883;      % example value
T = -99.47;       % example value
h = 0.1;      % sample time

% Continuous-time matrices
A = [0  1     0
     0 -1/T  -K/T
     0  0     0];

B = [0
     K/T
     0];
C = [1 0 0];
E = [0 0
     1 0
     0 1];

% Discrete-time system
[Ad, Bd] = c2d(A, B, h);
[~, Ed] = c2d(A, E, h);


function [x_hat, P_hat, x_prd, P_prd] = kalman_filter(x_prd, P_prd, y, u, Ad, Bd, Cd, Ed, Qd, Rd)
n = length(x_prd);

% Kalman gain
K = P_prd * Cd' / (Cd * P_prd * Cd' + Rd);

% Corrector (measurement update)
innov = y - Cd * x_prd;
% innov = mod(innov + pi, 2*pi) - pi;
x_hat = x_prd + K * innov;
IKC   = eye(n) - K * Cd;
P_hat = IKC * P_prd * IKC' + K * Rd * K';

% Predictor (time update)
x_prd = Ad * x_hat + Bd * u;
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';
end

