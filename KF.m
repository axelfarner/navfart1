function [x_hat, P_hat] = KF(x_prev, P_prev, y, u, Ad, Bd, Cd, Ed, Qd, Rd)
n = length(x_prev);

% Predictor (time update)
x_prd = Ad * x_prev + Bd * u;
P_prd = Ad * P_prev * Ad' + Ed * Qd * Ed';

% Kalman gain
K = P_prd * Cd' / (Cd * P_prd * Cd' + Rd);

% Corrector (measurement update)
innov = y - Cd * x_prd;
innov = ssa(y - Cd * x_prd);
x_hat = x_prd + K * innov;
IKC   = eye(n) - K * Cd;
P_hat = IKC * P_prd * IKC' + K * Rd * K';

end

