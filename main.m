% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           Your name
% Study program:    Your study program

% Add folder for 3-D visualization files
addpath(genpath('flypath3d_v3'))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
T_final = 6800;	        % Final simulation time (s)
h = .1;                % Sampling time (s)

Rd = (.5 * pi / 180)^2; % measurement noise
qr = 1e-5;
qb = 1e-3;
Qd = diag([qr, qb]);
x_hat = [(-110 * pi / 180) 0 0]';
P_hat = diag([.1,.01,.1]);

U_ref   = 9;            % desired surge speed (m/s)
print = true;

% initial states
eta_0 = [0 0 (-110 * pi / 180)]';
nu_0  = [0 0 0]';
delta_0 = 0;
n_0 = 0;
Qm_0 = 0;
x = [nu_0' eta_0' delta_0 n_0, Qm_0]'; % The state vector can be extended with addional states here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ad, Bd, Ed, Cd] = get_discrete_matrices();

t = 0:h:T_final;                % Time vector
nTimeSteps = length(t);         % Number of time steps

simdata = zeros(nTimeSteps, 18); % Pre-allocate matrix for efficiency
rng(27); % seed the noise
noise = [normrnd(0,.5,nTimeSteps,1), normrnd(0,.1,nTimeSteps,1)];
noise = noise * (pi / 180); % convert noise to radians

e_int = 0; % integral term of error in yaw
y_int = 0; %integral term of error in course
xd = [x(6) x(3) 0]'; % initialize desired yaw and yawrate to be same as starting position

% compute Ki to use for anti integrator windup term
wb = 0.02;
zeta = 1;
wn = wb / sqrt(1 - 2*zeta + sqrt(4*zeta^4 - 4*zeta^2 + 2));
T = 169.56;
K = 0.0075;
Ki = wn ^3 * T / (10 * K);

% numerator and denomintar in transfer function for reference model
wn_refmodel = 0.03;
num = wn_refmodel^3;
den = [1 (wn_refmodel+2*zeta*wn_refmodel) (wn_refmodel^2 + 2*zeta*wn_refmodel^2) wn_refmodel^3];    
n = length(den) - 1; % system order

% Construct Companion Form as state-space representation of reference model:
A = [zeros(n-1, 1), eye(n-1);
     -fliplr(den(2:end))];
 
B = [zeros(n-1, 1); num(end)];

% ref model returns the first order approx of the increment xd, xd_dot
ref_model = @(xd, r) (A * xd + B * r) *h;
xk_prev = 0;
showPlot = false;
for i = 1:nTimeSteps
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 1a) Add current here 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vc = 1; % current speed (m/s)
    betaVc = 45 * pi/180; % current direction (rad)
    uc = Vc * cos(betaVc - x(6));
    vc = Vc * sin(betaVc - x(6));
    nu_c = [ uc vc 0 ]';
%     nu_c = [0 0 0]'; % remove current
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 1c) Add wind here 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta_Vw = 135 * pi/180; % wind direction (rad)
    Vw = 10; % wind speed (m/s)
    rho_air = 1.247; % air density (kg/m^3)
    c_y = 0.95;
    c_n = 0.15;
    L_oa = 161; % length overall (m)
    A_Lw = 10 * L_oa; % lateral projected area (m^2)

    if t(i) >= 200
        % wind velocity in body frame
        uw = Vw * cos(beta_Vw - x(6));
        vw = Vw * sin(beta_Vw - x(6));

        % ship velocity relative to the air
        u_rw = x(1) - uw;
        v_rw = x(2) - vw;
        V_rw = sqrt(u_rw^2 + v_rw^2); % relative wind speed (m/s)
        gamma_rw = -atan2(v_rw, u_rw); % wind angle of attack (rad)

        % wind coefficients for symmetrical ships, eqs. (10.22-10.23)
        Cy = c_y * sin(gamma_rw);
        Cn = c_n * sin(2*gamma_rw);

        Ywind = 0.5 * rho_air * V_rw^2 * Cy * A_Lw;
        Nwind = 0.5 * rho_air * V_rw^2 * Cn * A_Lw * L_oa;
    else
        Ywind = 0;
        Nwind = 0;
    end
    tau_wind = [0 Ywind Nwind]';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Guidance model and reference model 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ee
    % Guidance law
    [xk1,yk1,xk,yk,last] = WP_selector(x(4),x(5));

    % --- Reset ILOS integrator on new waypoint leg ---
    if xk ~= xk_prev
        y_int = 0;          % Reset integrator when entering new leg
        xk_prev = xk;
    end
    % -------------------------------------------------
    [e_y,pi_h] = crossTrackError(xk1,yk1,xk,yk,x(4),x(5),showPlot);
    

    %%%%%%%%%%%%%%%%%%%%%%
    % 2 c
    %%%%%%%%%%%%%%%%%%%%%%
    chi_d = guidance(e_y,pi_h); 
    crab_angle = atan2(x(2), x(1));      % crab angle
    if x(1) > 0.7 * U_ref % check that the boat is at about cruise speed to ensure crab angle well defined
        heading_d = chi_d - crab_angle;
    else
        heading_d = chi_d;
    end
    
    %%%%%%%%%%%%%%%%&&&&
    % 2 d
    %%%%%%%%%%%%%%%%%&&&
%     [chi_d, dy]= guidance_integrator(e_y,pi_h,h,y_int);
%     y_int = y_int + dy;
%     heading_d = chi_d;

    %%%%%%%%%%%%%%%%%%%%
    % Reference models
    %%%%%%%%%%%%%%%%%%%%

    xd_dot = ref_model(xd,heading_d);
    xd = xd + xd_dot;

    psi_d = xd(1);
    r_d = xd(2);
    u_d = U_ref;
    %%%%%%%%%%%%%%%%%%%%%%%
    % KALMAN FILTER
    %%%%%%%%%%%%%%%%%%%%%%%

    % predict the current state of x
    if i > 1
        [x_hat, P_hat] = KF(x_hat, P_hat, x(6) + noise(i,1),delta_c,Ad,Bd,Cd,Ed,Qd,Rd);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 2d) Add the heading controller here 
    % Define it as a function   
    
    % to prevent integrator windup, stop integrating when integral term alone puts you over rudder limitations
    integrator_increment = h * (psi_d- x(6)) * 0;
    if abs(1.1 * (e_int + integrator_increment) * Ki) < 40 * pi/180
        e_int = e_int + integrator_increment;
    end
    e_r = r_d - x_hat(2);
    e_psi = psi_d - x_hat(1);

    delta_c = PID_heading(e_psi,e_r,e_int);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part 2, 3e) Add open loop speed control here
    % Define it as a function
    %
    % The result should look like this:
    % n_c = open_loop_speed_control(U_ref);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_c = open_loop_speed_control(U_ref);                   % propeller speed [radians per second (rps)]
    % ship dynamics
    u = [delta_c n_c]';
    [xdot,tau_total] = ship(x,u,nu_c,tau_wind);
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d chi_d e_y x_hat'];     
 
    % Euler integration
    % x = euler2(xdot,x,h); 
    % Runge Kutta 4 integration
    x = rk4(@ship,h,x,u,nu_c,tau_wind);

    % --- Progress Update Logic ---
    if mod(i, floor(nTimeSteps / 10)) == 0 % Print an update every 10%
        progress = (i/nTimeSteps) * 100;
        fprintf('  %d%% complete\n', round(progress));
        
    elseif i == 1
        
        disp("  Simulation in progress:")
        fprintf('  %d%% complete\n', 0);
    end
    if last
        break
    end
end
simdata = simdata(1:i,:);
t = t(1:i);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u           = simdata(:,1);                 % m/s
v           = simdata(:,2);                 % m/s
r           = simdata(:,3);                 % rad/s
r_deg       = (180/pi) * r;                 % deg/s
x           = simdata(:,4);                 % m
y           = simdata(:,5);                 % m
psi         = simdata(:,6);                 % rad
psi_deg     = (180/pi) * psi;               % deg
delta_deg   = (180/pi) * simdata(:,7);      % deg
n           = (30/pi) * simdata(:,8);       % rpm
delta_c_deg = (180/pi) * simdata(:,9);      % deg
n_c         = (30/pi) * simdata(:,10);      % rpm
u_d         = simdata(:,11);                % m/s
psi_d       = simdata(:,12);                % rad
psi_d_deg   = (180/pi) * psi_d;             % deg
r_d         =  simdata(:,13);               % rad/s
r_d_deg     = (180/pi) * r_d;               % deg/s
pathplotter(x,y);
%%
figure(3)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions'); xlabel('(m)'); ylabel('(m)'); 
subplot(312)
plot(t,psi_deg,t,psi_d_deg,'linewidth',2);
title('Actual and desired yaw angle'); xlabel('Time (s)');  ylabel('Angle (deg)'); 
legend('actual yaw','desired yaw')
subplot(313)
plot(t,r_deg,t,r_d_deg,'linewidth',2);
title('Actual and desired yaw rates'); xlabel('Time (s)');  ylabel('Angle rate (deg/s)'); 
legend('actual yaw rate','desired yaw rate')

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocity'); xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('actual surge','desired surge')
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
title('Actual and commanded propeller speed'); xlabel('Time (s)'); ylabel('Motor speed (RPM)');
legend('actual RPM','commanded RPM')
subplot(313)
plot(t,delta_deg,t,delta_c_deg,'linewidth',2);
title('Actual and commanded rudder angle'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('actual rudder angle','commanded rudder angle')
%% Create objects for 3-D visualization 
% Since we only simulate 3-DOF we need to construct zero arrays for the 
% excluded dimensions, including height, roll and pitch
z = zeros(length(x),1);
phi = zeros(length(psi),1);
theta = zeros(length(psi),1);


% % create object 1: ship (ship1.mat)
new_object('flypath3d_v3/ship1.mat',[x,y,z,phi,theta,psi],...
'model','royalNavy2.mat','scale',(max(max(abs(x)),max(abs(y)))/1000),...
'edge',[0 0 0],'face',[0 0 0],'alpha',1,...
'path','on','pathcolor',[.89 .0 .27],'pathwidth',2);

% Plot trajectories 
flypath('flypath3d_v3/ship1.mat',...
'animate','on','step',500,...
'axis','on','axiscolor',[0 0 0],'color',[1 1 1],...
'font','Georgia','fontsize',12,...
'view',[-25 35],'window',[900 900],...
'xlim', [min(y)-0.1*max(abs(y)),max(y)+0.1*max(abs(y))],... 
'ylim', [min(x)-0.1*max(abs(x)),max(x)+0.1*max(abs(x))], ...
'zlim', [-max(max(abs(x)),max(abs(y)))/100,max(max(abs(x)),max(abs(y)))/20]); 
%% Crab and sideslip angles
% current in body frame (recomputed from logged heading)
uc_log = Vc * cos(betaVc - psi);
vc_log = Vc * sin(betaVc - psi);
u_r = u - uc_log;
v_r = v - vc_log;
U_r = sqrt(u_r.^2 + v_r.^2);

beta_c_deg = (180/pi) * atan2(v, u);      % crab angle
beta_deg   = (180/pi) * asin(v_r ./ U_r); % sideslip angle

figure(4)
figure(gcf)
plot(t(50:end),beta_c_deg(50:end),t(50:end),beta_deg(50:end),'linewidth',2);
title('Crab and sideslip angles'); xlabel('Time (s)'); ylabel('Angle (deg)');
legend('crab angle \beta_c','sideslip angle \beta')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2(a) / 2(b): Course, Desired Course, Heading, Crab, and Sideslip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chi_d_deg = (180/pi) * simdata(:,14);

% Course over ground (chi = psi + beta_c)
% Using velocity over ground in NED coordinates:
% x_dot_n = u*cos(psi) - v*sin(psi)
% y_dot_n = u*sin(psi) + v*cos(psi)
% chi = atan2(y_dot_n, x_dot_n)
x_dot_n = u .* cos(psi) - v .* sin(psi);
y_dot_n = u .* sin(psi) + v .* cos(psi);
chi_deg = (180/pi) * atan2(y_dot_n, x_dot_n);

% Wrap angles to [-180, 180] deg to avoid vertical wrap jumps across +/-180
wrap_deg = @(angle) mod(angle + 180, 360) - 180;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2: Course, Desired Course, Heading, Crab, Sideslip, and e_y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chi_d_deg = (180/pi) * simdata(:,14);
e_y_log   = simdata(:,15);                  % Cross-track error (m)

% Course over ground (chi = psi + beta_c)
x_dot_n = u .* cos(psi) - v .* sin(psi);
y_dot_n = u .* sin(psi) + v .* cos(psi);
chi_deg = (180/pi) * atan2(y_dot_n, x_dot_n);

% Wrap angles to [-180, 180] deg to avoid vertical wrap jumps across +/-180
wrap_deg = @(angle) mod(angle + 180, 360) - 180;

figure(5);
clf;

% Subplot 1: Course, Desired Course, and Heading/Yaw
subplot(3, 1, 1);
plot(t, wrap_deg(chi_deg), 'b', 'LineWidth', 1.5); hold on;
plot(t, wrap_deg(chi_d_deg), 'r--', 'LineWidth', 1.5);
plot(t, wrap_deg(psi_deg), 'k:', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Angle (deg)');
title('Course (\chi), Desired Course (\chi_d), and Heading/Yaw (\psi)');
legend('\chi (course)', '\chi_d (desired course)', '\psi (heading / yaw)', 'Location', 'best');

% Subplot 2: Crab and Sideslip Angles
subplot(3, 1, 2);
plot(t(100:end), beta_c_deg(100:end), 'LineWidth', 1.5); hold on;
plot(t(100:end), beta_deg(100:end), '--', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Angle (deg)');
title('Crab Angle (\beta_c) and Sideslip Angle (\beta)');
legend('\beta_c (crab)', '\beta (sideslip)', 'Location', 'best');

% Subplot 3: Cross-Track Error (e_y)
subplot(3, 1, 3);
plot(t, e_y_log, 'm', 'LineWidth', 1.5); hold on;
yline(0, 'k--', 'LineWidth', 1); % Zero-error reference line
grid on;
xlabel('Time (s)');
ylabel('Error (m)');
title('Cross-Track Error (e_y)');
legend('e_y', 'Zero Reference', 'Location', 'best');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 4(a): True vs. Noisy Yaw Angle and Yaw Rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncate pre-generated noise to match actual simulation length
noise_used = noise(1:length(t), :);

% Add zero-mean Gaussian noise (in deg and deg/s) to true states
psi_meas_deg = psi_deg + noise_used(:, 1) * 180 / pi;   % std = 0.5 deg
r_meas_deg   = r_deg   + noise_used(:, 2) * 180 / pi;   % std = 0.1 deg/s

figure(6);
clf;

% Subplot 1: Yaw Angle (True vs. Noisy Measurement)
subplot(2, 1, 1);
plot(t, psi_meas_deg, 'Color', [0.7 0.7 0.7], 'LineWidth', 1); hold on;
plot(t, psi_deg, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Angle (deg)');
title('Yaw Angle (\psi): True State vs. Noisy Measurement');
legend('\psi_{meas} (noisy)', '\psi (true)', 'Location', 'best');

% Subplot 2: Yaw Rate (True vs. Noisy Measurement)
subplot(2, 1, 2);
plot(t, r_meas_deg, 'Color', [0.7 0.7 0.7], 'LineWidth', 1); hold on;
plot(t, r_deg, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Yaw Rate (deg/s)');
title('Yaw Rate (r): True State vs. Noisy Measurement');
legend('r_{meas} (noisy)', 'r (true)', 'Location', 'best');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 4(b): Kalman Filter State Estimates vs. True States
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi_hat_deg  = (180/pi) * simdata(:, 16);   % Estimated yaw angle (deg)
r_hat_deg    = (180/pi) * simdata(:, 17);   % Estimated yaw rate (deg/s)
bias_hat_deg = (180/pi) * simdata(:, 18);   % Estimated rudder bias (deg)

figure(7);
clf;

% Subplot 1: True vs. Estimated Yaw Angle
subplot(3, 1, 1);
plot(t, psi_deg, 'b', 'LineWidth', 1.5); hold on;
plot(t, psi_hat_deg, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Angle (deg)');
title('Yaw Angle: True ($\psi$) vs. Kalman Filter Estimate ($\hat{\psi}$)', 'Interpreter', 'latex');
legend('$\psi$ (true)', '$\hat{\psi}$ (estimated)', 'Interpreter', 'latex', 'Location', 'best');

% Subplot 2: True vs. Estimated Yaw Rate
subplot(3, 1, 2);
plot(t, r_hat_deg, 'r', 'LineWidth', 1.5); hold on;
plot(t, r_deg, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Yaw Rate (deg/s)');
title('Yaw Rate: True ($r$) vs. Kalman Filter Estimate ($\hat{r}$)', 'Interpreter', 'latex');
legend('$\hat{r}$ (estimated)', '$r$ (true)', 'Interpreter', 'latex', 'Location', 'best');

% Subplot 3: Estimated Rudder Bias
subplot(3, 1, 3);
plot(t, bias_hat_deg, 'k', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Bias (deg)');
title('Estimated Rudder Bias ($\hat{b}$)', 'Interpreter', 'latex');
legend('$\hat{b}$ (estimated)', 'Interpreter', 'latex', 'Location', 'best');