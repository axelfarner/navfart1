function [heading_desired,dy] = guidance_integrator(e_y, pi_h,h,y_int)
Delta = 1000; %look-ahead distance
Kp = 1 / Delta;
Ki = Kp;
k = 1;
dy = Delta * e_y / (Delta ^2 + (e_y + k * y_int)^2) * h;

heading_desired = pi_h - atan(e_y * Kp + Ki* (y_int + dy));

end