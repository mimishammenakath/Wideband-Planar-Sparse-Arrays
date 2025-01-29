function [A ] = Steering_Matrix_CS_array_syn_Planar( r_x,r_y,u_vec,v_vec,U_num,V_num,k,MN)
%Creating the steering matrix of planarr array for multiple singal
%directions.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Software Used                   : MATLAB VERSION 2023(a)
% Platform Used                   : Intel(R) Xeon(R) Gold 6248 CPU @ 2.50GHz, 2494 Mhz, 20 Core(s), 40 Logical Processor(s)
%                                        : 128 GB RAM
% Operating System Used           : Windows 11 (64 bit)
% For Further Information Contact : Mimisha M Menakath
%                                   Indian Institute of Technology
%                                   Palakkad
%                        email id : 122014006@smail.iitpkd.ac.in
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%   Input arguments
%   r_x              : x position or locations of the sensors
%   r_x              : y position or locations of the sensors
%   u_vec            : signal azimuth  direction vector
%   v_vec            : signal elevation  direction vector
%   MN                : total number of sensors
%   k                : wavenumber
%   U_num                : length of u_vec
%   V_num              : length of v_vec
%   A                : Steering Matrix for multiple signal directions
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

k_x = k*u_vec;% x component of wavenumber 
k_y = k*v_vec;% y component of wavenumber 


K_X = repmat(k_x',MN,1); % repeating k_x for all the sensors
K_Y = repmat(k_y',MN,1); % repeating k_y for all the sensors

R_x = repmat(r_x,1,U_num);
R_y= repmat(r_y,1,V_num);
% creating a sensing matrix A
clear i;
 A = exp(i*(K_X.*R_x+K_Y.*R_y));
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
