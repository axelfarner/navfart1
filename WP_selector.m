function [xk1,yk1,xk,yk,last] = WP_selector(x,y)
waypoints = load("WP.mat").WP;
epsilon = 400;
last = false;
% find what waypoint we have reached
% find what is the next waypoint
persistent last_xk 
persistent last_yk 
persistent last_xk1
persistent last_yk1
xk1 = 0;
yk1 = 0;
xk = 0;
yk = 0;
for i = 1:size(waypoints,2)
    wp = waypoints(:,i);
    if abs(x - wp(1)) < epsilon && abs(y-wp(2)) < epsilon % check if waypoint reached
        if i == size(waypoints,2)
            last = true;
            break
        end
        last_xk = wp(1);
        last_yk = wp(2);
        last_xk1 = waypoints(1,i+1);
        last_yk1 = waypoints(2,i+1);  
    end
  xk1 = last_xk1;
  yk1 = last_yk1;
  xk = last_xk;
  yk = last_yk;
end