function [states,time]=car_simulation
clear all
close all
%2-D Environment
% Y
% ^
% |
% |_ _ _ > X
%% Dimensions ==> Only necessary for visualization
%vector of 41 evenly spaced values ranging from 0 to 5
x = linspace(0,5,41);
%vector of 40 evenly spaced values ranging from -0.75 to 0.75
y = linspace(-.75,.75,40);
%Matrices of the x,y values, sized to match the vector lengths
%(used for quiver plotting)
[x_mat y_mat]=meshgrid(x,y);
%% Time Values
t_step = .01; % [s] Commands changing at 100Hz
t_run = 2.7; % [s] Run for 2.7 seconds
time_in = 0; % [s] Start time at 0
time = 0; % [s] vector to collect all of the output times
%% Obstacle Cells (Matlab Tricks to keep things organized/ flexible in for-loops: Cells and Structures)
Ob{1}.x = 1.5; %x location
Ob{1}.y = 0; %y location
Ob{2}.x = 3.5;
Ob{2}.y = -0.3;
% Repulsive potential function parameters
tau = 20; %Repulsive maximum value sensitivity parameter
zeta = .08; %Repulsive region of influence parameter
% Repulsive potential function parameters
lam_rep.x = 0.5; %Relative weight of x component
lam_rep.y = 2.2;%Relative weight of y component
% Repulsive potential function parameters from the edge of the road
xi = 0.01;
kappa = 1;
%% Goal
G.x = 5;
G.y = 0;
% Attractive potential function parameters
lam_att.x = .5; %Relative weight of x component 
lam_att.y = 1.5;%Relative weight of y component

%% Initial Conditions
x_cg_0 = 0;
y_cg_0 = 0;
yaw_angle_0 = 0;
v_lat_0 = 0;
yaw_rate_0 = 0;
steer_angle_0 = 0;
%Initial Condition Vector used for ODE_45 input
x0 = [x_cg_0,y_cg_0,yaw_angle_0,v_lat_0,yaw_rate_0,steer_angle_0]';
% Matrix to collect all of the output states we care about
states = x0';
%% Sign determination of angle for rotational force
function [val_con]=signum(val)
if(val>=0)
val_con = 1;
else
val_con = -1;
end
end
%% Simulation Loops
X = x0;
for i = 1:t_step:t_run
x_cg = X(1);
y_cg = X(2);
yaw_angle = X(3);
v_lat = X(4);
yaw_rate = X(5);
steer_angle = X(6);
%Linear Attractive Potential Function Gradient structure
phi_att.x = -2*lam_att.x*(x_cg - G.x)/sqrt(lam_att.x*(x_cg-G.x)^2+lam_att.y*(y_cg-G.y)^2);
phi_att.y = -2*lam_att.y*(y_cg - G.y)/sqrt(lam_att.x*(x_cg- G.x)^2+lam_att.y*(y_cg-G.y)^2);
% Edge of road
map_pfg_edge.x = 0;
map_pfg_edge.y = 2*(y_cg-.75)*kappa/xi*exp(-((y_cg- .75)^2)/xi)+2*kappa/xi*(y_cg+.75)*exp(-((y_cg+.75)^2)/xi);
%Initialize the repulsive gradient components (to use in a for-loop sum)
phi_rep.x = 0;
phi_rep.y = 0;
phi_rot.x = 0;
phi_rot.y = 0;
%Sum up the repulsive gradient values from all of the obstacles
%(this is why I organized the obstacle values the way I did)
for k = 1:size(Ob,2)
phi_rep.x = phi_rep.x+2*lam_rep.x*tau/zeta*(x_cg-Ob{k}.x)*exp(-((lam_rep.x*(x_cg-Ob{k}.x)^2)+(lam_rep.y*(y_cg-Ob{k}.y)^2))/zeta);

phi_rep.y = phi_rep.y+2*lam_rep.y*tau/zeta*(y_cg-Ob{k}.y)*exp(-((lam_rep.x*(x_cg-Ob{k}.x)^2)+(lam_rep.y*(y_cg-Ob{k}.y)^2))/zeta);
alpha = atan2((Ob{k}.y-y_cg),(Ob{k}.x-x_cg));
%Rotational forces%
phi_rot.x = phi_rot.x-1*signum(alpha)*2*lam_rep.y*tau/zeta*(y_cg-Ob{k}.y)*exp(-((lam_rep.x*(x_cg-Ob{k}.x)^2)+(lam_rep.y*(y_cg-Ob{k}.y)^2))/zeta);
phi_rot.y = phi_rot.y+signum(alpha)*2*lam_rep.x*tau/zeta*(x_cg-Ob{k}.x)*exp(-((lam_rep.x*(x_cg-Ob{k}.x)^2)+(lam_rep.y*(y_cg-Ob{k}.y)^2))/zeta);
end
%Find the total gradient value for each component
phi_total.x = phi_rot.x+phi_att.x+map_pfg_edge.x;
phi_total.y = phi_rot.y+phi_att.y+map_pfg_edge.y;
%Find the magnitude of the vector made up of both components
phi_norm.magnitude = sqrt(phi_total.x^2+phi_total.y^2);
%Normalize the components with the total magnitude
phi_norm.x = phi_total.x/phi_norm.magnitude;
phi_norm.y = phi_total.y/phi_norm.magnitude;
%Calculate new commands based on the gradient values
%(We calculate the steering angle)
steer_angle = atan2(phi_norm.y,phi_norm.x)-yaw_angle;
%Prepare the input vector for ODE_45 based on the states from the
%previous run-loop, but with the new commands
x_in=[x_cg,y_cg,yaw_angle,v_lat,yaw_rate,steer_angle]';
%Run ODE_45 for the next time step
%[t,y] = ode45(yp,[t0,tf],y0)
% yp = the function
% t0,tf = initial and terminal values of t
% y0 = initial value of y at t0
%Basically a numerical integrator
[time_out,x_out]=ode45(@Dynamics,[time_in time_in+t_step],x_in);
%Append new output from ode_45 function onto previous
%(Just keeps adding the new output to the end)
time=[time;time_out];
states = [states;x_out];
%Collect the final values of the output to use in the next iteration.
X=x_out(end,:);
%This sets the next begin time
time_in = time_out(end,1);
end


%% PFG Map (for visualization only)
% No part of the trajectory calculations are based off of this
for i = 1:length(x)
for j = 1:length(y)
map_pfg_att.x(i,j) = 2*lam_att.x*(x(i)-G.x)/sqrt(lam_att.x*(x(i)-G.x)^2+lam_att.y*(y(j)-G.y)^2);
map_pfg_att.y(i,j) = 2*lam_att.y*(y(j)-G.y)/sqrt(lam_att.x*(x(i)- G.x)^2+lam_att.y*(y(j)-G.y)^2);
map_pfg_edge.x(i,j) = 0;
map_pfg_edge.y(i,j) = 2*(y(j)-.75)*kappa/xi*exp(-((y(j)- .75)^2)/xi)+2*(y(j)+.75)*kappa/xi*exp(-((y(j)+.75)^2)/xi);
map_pfg_rep.x(i,j) = 0;
map_pfg_rep.y(i,j) = 0;
map_pfg_rot.x(i,j) = 0;
map_pfg_rot.y(i,j) = 0;
for k = 1:size(Ob,2)
map_pfg_rep.x(i,j) = map_pfg_rep.x(i,j)+2*lam_rep.x*tau/zeta*(x(i)-Ob{k}.x)*exp(-(lam_rep.x*(x(i)-Ob{k}.x)^2+lam_rep.y*(y(j)-Ob{k}.y)^2)/zeta);
map_pfg_rep.y(i,j) = map_pfg_rep.y(i,j)+2*lam_rep.y*tau/zeta*(y(j)-Ob{k}.y)*exp(-(lam_rep.x*(x(i)-Ob{k}.x)^2+lam_rep.y*(y(j)-Ob{k}.y)^2)/zeta);
alpha = atan2((Ob{k}.y-y(j)),(Ob{k}.x-x(i)));
%Rotational Force%
map_pfg_rot.x(i,j) = map_pfg_rot.x(i,j)-1*signum(alpha)*2*lam_rep.y*tau/zeta*(y(j)-Ob{k}.y)*exp(-(lam_rep.x*(x(i)-Ob{k}.x)^2+lam_rep.y*(y(j)-Ob{k}.y)^2)/zeta);
map_pfg_rot.y(i,j) = map_pfg_rot.y(i,j)+signum(alpha)*2*lam_rep.x*tau/zeta*(x(i)-Ob{k}.x)*exp(-(lam_rep.x*(x(i)-Ob{k}.x)^2+lam_rep.y*(y(j)-Ob{k}.y)^2)/zeta);
end
map_pfg_total.x(i,j) = -map_pfg_att.x(i,j)+map_pfg_edge.x(i,j)+map_pfg_rot.x(i,j);
map_pfg_total.y(i,j) = -map_pfg_att.y(i,j)+map_pfg_edge.y(i,j)+map_pfg_rot.y(i,j);
map_pfg_norm.magnitude(i,j) = sqrt(map_pfg_total.x(i,j)^2+map_pfg_total.y(i,j)^2);
map_pfg_norm.x(i,j) = map_pfg_total.x(i,j)/map_pfg_norm.magnitude(i,j);
map_pfg_norm.y(i,j) = map_pfg_total.y(i,j)/map_pfg_norm.magnitude(i,j);
end
end

%Plot the PFG map
figure(1)
quiver(x_mat,y_mat,map_pfg_norm.x',map_pfg_norm.y')
hold on
xlabel('X Direction')
ylabel('Y Direction')
%Plot the obstacles
rectangle('Position',[1.3,-0.15,0.4,0.3],'LineWidth',3)
rectangle('Position',[3.3,-0.45,0.4,0.3],'LineWidth',3)
%% Trajectory output (overlayed on top of the pfg map)
figure(1)
plot(states(:,1),states(:,2),'r','LineWidth',2)
hold off
axis([0 5 -0.75 .75])
axis image
end