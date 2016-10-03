function [ int_time, ground_v,data_rate ] = conical_calc( h,vel, del_y,Nx,Ny,locx,bits,Nband)

Re = 6.371e6;

% compute ground velocty
ground_v = Re * vel / (Re+h);

% compute integration time
int_time = del_x / ground_v;

data_rate = Nband*Nx*bits / (del_x/ground_v);