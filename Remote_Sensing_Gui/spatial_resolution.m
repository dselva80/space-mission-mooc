function [ del_x ] = spatial_resolution(h,lambda,D )
%Computes spatial resolution for circule diffraaction limited aperture

del_x = 1.22*h*lambda/D;

end

