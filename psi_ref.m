function ref = psi_ref(i)
wref = -0.007;
h = 0.1;
t = h * i;
t0 = 400;

% if t > t0
%     ref = (t-t0) * wref;
% else
%     ref = 0;
% end
% 
% if t > t0
%     ref = sin(wref * i* h);
% else 
%     ref = 0;
% end
% for setpoint testing
if t > t0
    ref = 10 * pi/180;  % desired yaw angle (rad)
elseif t > 100
    ref = -20 * pi/180;
else 
    ref = 0;
end
end