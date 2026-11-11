function [e_y,pi_h] = crossTrackError(xk1,yk1,xk,yk,x,y,showPlot)
pi_h = atan2(yk1 - yk, xk1-xk);
pi_h = pi_h;
if showPlot
    [~,~,e_y] = crosstrack(xk1, yk1, xk, yk, x,y,1);
else 
    [~,~,e_y] = crosstrack(xk1, yk1, xk, yk, x,y);
end
end