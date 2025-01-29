function [b_u ] = BP_Planar_array(  r_x,r_y,u_vec,v_vec,U_num,V_num,k,MN,w,Mb,Nb)
%   This function will calculate the beam pattern of the planar array.
%   Input arguments
%     r_x     : x index of the sensor elements in the array
%     r_y     : y index of the sensor elements in the array 
%     u_vec : signal direction vector in u domain
%     v_vec : signal direction vector in v domain
%     k        : wavenumber
%     MN     : Number of sensor elements in the array
%     w       : weight vector of the sensor elements
%     U_num      : Number of signal direction in azimuth direction
%     V_num      : Number of signal direction in elevation direction
%  Output arguments
%    b_u      : Beampattern response of the planar array 
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
[A ] =  Steering_Matrix_CS_array_syn_Planar( r_x,r_y,u_vec,v_vec,U_num,V_num,k,MN);
% w = ones(M,1); %Creating a unit weight vector to obtained beam pattern

b_u = A'*w; %This gives BP of URA at different values of u_vec and v_vec.
b_u=reshape(b_u,Mb,Nb); % 3D beam pattern
% b_u = b_u/max(abs(b_u(:))); %Normalizing beampattern
% b_u = (M*sinc((M/2)*u_vec))./(sinc(u_vec/2));
end

