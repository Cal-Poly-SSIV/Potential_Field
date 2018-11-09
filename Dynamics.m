function [dX] = Dynamics(t,X)
%%Dynamics of bicycle model
%dv_x = 0 because we assumed that v_x is constant
%dv_y lateral acceleration
%C_front = sideslip coefficient of the front wheel [N/rad]
%C_rear = sideslip coefficient of the rear wheel [N/rad]
%a1 = distance between front wheel and CG [m]
%a2 = distance between rear wheel and CG [m]
%m = vehicle mass [kg]
%Iz = mass moment of inertia [kg*m²]
%v_long = longitudonal velocity [m/s](We assume that v_x is constant)
%v_lat = lateral velocity [m/s]
%x_cg = x-coordinate of center of gravity [m]
%y_cg = y-coordinate of center of gravity [m]
%Parameters example 402 from Vehicle Dynamics, Reza N. Jazar, page 637
C_front = 357.15;
C_rear = 357.15;
m = 3.117;
Iz = 0.0452;
a1 = 0.212;
a2 = 0.126;
v_long = 3;
x_cg = X(1);
y_cg = X(2);
yaw_angle = X(3);
v_lat = X(4);
yaw_rate = X(5);
steer_angle = X(6);
dv_lat = 1/(m*v_long)*(-a1*C_front*cos(steer_angle)+a2*C_rear)*yaw_rate - 1/(m*v_long)*(C_front*cos(steer_angle)+C_rear)*v_lat + 1/m*C_front*cos(steer_angle)*steer_angle-yaw_rate*v_long;
dyaw_rate = 1/(Iz*v_long)*(-a1^2*C_front*cos(steer_angle)-a2^2*C_rear)*yaw_rate - 1/(Iz*v_long)*(a1*C_front*cos(steer_angle)-a2*C_rear)*v_lat + 1/Iz*a1*C_front*cos(steer_angle)*steer_angle;
dsteer_angle = 0;
%%Calculate new coordinates of vehicle(bicycle model)
%v_x = longitudonal velocity [m/s](We assume that v_x is constant)
%v_y = lateral velocity [m/s]
%x_cg = x-coordinate of center of gravity [m]
%y_cg = y-coordinate of center of gravity [m]
dx_cg = v_long*cos(yaw_angle) - v_lat*sin(yaw_angle);
dy_cg = v_long*sin(yaw_angle) + v_lat*cos(yaw_angle);
dyaw_angle = yaw_rate;
dX = [dx_cg,dy_cg,dyaw_angle,dv_lat,dyaw_rate,dsteer_angle]';
end
