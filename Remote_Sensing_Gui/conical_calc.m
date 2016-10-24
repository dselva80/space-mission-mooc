function [ int_time, data_rate ] = conical_calc( ground_v, del_y,Nx,Ny,locx,bits,Nband)

Re = 6.371e6;



% compute integration time
int_time = del_x / ground_v;

data_rate = Nband*Nx*bits / (del_x/ground_v);

