function [ int_time, data_rate ] = whiskbroom( ground_v,del_y,Nx,Ny,bits,Nband)
% Computes integration time, ground velocity for pushbroom sensors

% INPUTS:
% h = height or altitude of spacecraft. [m]
% v = velocity at one point in time [m/s]
% del_x = pixel size [m]
% FOV is field of view in [m]

% OUTPUTS:
% int_time = integration time [s]
% ground_v = ground velocity at one time [m/s]




% compute integration time
int_time = ( Ny*del_y ) / (2*ground_v*Nx);

data_rate = Nband*Nx*Ny*bits / (del_y/ground_v);

end

