function [ int_time, data_rate ] = pushbroom( ground_v,del_x,Nx,Nband,bits)
% Computes integration time, ground velocity for pushbroom sensors

% INPUTS:
% ground_v = velocity at one point in time [m/s]
% del_x = pixel size [m]
% Nx is # crosstrack pixels
% Nband is number of bands
% 

% OUTPUTS:
% int_time = integration time [s]
% ground_v = ground velocity at one time [m/s]


% compute integration time
int_time = del_x / ground_v;

data_rate = Nband*Nx*bits / (del_x/ground_v);

end

