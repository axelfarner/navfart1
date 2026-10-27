function [delta_c] = PID_heading(e_psi,e_r,e_int)
   wb = 0.06;
   zeta = 1;
   wn = wb / sqrt(1 - 2*zeta + sqrt(4*zeta^4 - 4*zeta^2 + 2));
   T = 169.56;
   K = 0.0075;
   Kp = wn ^2 * T / K;
   Kd = (2 * zeta * wn * T - 1) / K;
   Ki = wn ^3 * T / (10 * K);
   delta_c = Kp * e_psi + Kd * e_r + Ki * e_int;
end