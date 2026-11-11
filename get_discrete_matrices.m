function [Ad, Bd, Ed, Cd] = get_discrete_matrices()
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
Cd = [1 0 0];
E = [0 0
     1 0
     0 1];

% Discrete-time system
[Ad, Bd] = c2d(A, B, h);
[~, Ed] = c2d(A, E, h);

end

