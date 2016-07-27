function [ int_time, ground_v,data_rate ] = whiskbroom( h,vel, del_y,Nx,Ny,locx,bits,Nband)
% Computes integration time, ground velocity for pushbroom sensors

% INPUTS:
% h = height or altitude of spacecraft. [m]
% v = velocity at one point in time [m/s]
% del_x = pixel size [m]
% FOV is field of view in [m]

% OUTPUTS:
% int_time = integration time [s]
% ground_v = ground velocity at one time [m/s]


% radius of the earth in [m]
Re = 6.371e6;

% compute ground velocty
ground_v = Re * vel / (Re+h);

% compute integration time
int_time = ( Ny*del_y ) / (2*ground_v*Nx);

data_rate = Nband*Nx*Ny*bits / (del_y/ground_v);

end

