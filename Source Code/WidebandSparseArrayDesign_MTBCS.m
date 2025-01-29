%-- Script for for wideband sparse planar array design using MTBCS
%-- Authors: Mimisha M Menakath and Hareesh G
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified on 29 - January - 2025
%-------------------------------------------------------------------------%
clear all
close all;clearvars;
%% Input parameters

c = 1500;                                                              % speed of sound in water
fH=450e3;                                                              % Highest frequency in the bandwidth
fL=150e3;                                                              % Lowest frequency in the bandwidth
f0=300e3;                                                              % Center frequency in the bandwidth
lambda_ref=c/fH;
theta_str=0;                                                           % Azimuth steering angle
phi_str= 0;                                                            % Elevation steering angle
look_dir=(-90:5:90);                                                   % steeing direction in degree
indx=find(look_dir==theta_str);
look_dir_u =  sind(look_dir)-sind(theta_str);                          % steering direction in u domain
look_dir_v =  sind(look_dir)-sind(phi_str);                            % steering direction in v domain
[U,V]=meshgrid(look_dir_u,look_dir_v);
U=U(:);                                                                % steering vector in u domain
V=V(:);                                                                % steering vector in v domain
U_num=length(U);
V_num=length(V);
Mb= length(look_dir_u);                                                % total number of steering angles in azimuth direction
Nb= length(look_dir_v);                                                % total number of steering angles in elevation direction
M_ref=100;                                                             % number of elements in x direction of the reference array for maximum frequency
N_ref=100;                                                             % number of elements in y direction of the reference array for maximum frequency
MN_ref=M_ref*N_ref;                                                    % total number of elements in the reference array
L=(M_ref-1)*lambda_ref/2;                                              % Aperture size
delta_f=c/(sqrt(2)*L);                                                 % frequency resolution
f=fL:delta_f:fH;
lambda = c./f;                                                         % wavelength of sound in water
k = 2*pi./lambda;                                                      % wave number

%% Design of desired Beam Pattern

for i=1:length(f)
    d_ref = lambda_ref/2;                                              % element spacing in the reference array
    X = d_ref*[-(M_ref-1)/2:1:(M_ref-1)/2]';                           %  x component of position vector of sensors
    Y = d_ref*[-(N_ref-1)/2:1:(N_ref-1)/2]';                           % y component of position vector of sensors

    [r_x_ref,r_y_ref]=meshgrid(X,Y);
    r_x_ref=r_x_ref(:);                                                % x_index of sensors
    r_y_ref=r_y_ref(:);                                                % y_index of sensors
    w1=chebwin(M_ref,22);
    w2=chebwin(N_ref,22);
    w=w1*w2';
    weights=reshape(w,M_ref*N_ref,1);
    % weights = ones(M_ref*N_ref,1);                             % rectangular window
    [b_u ] = BP_Planar_array(r_x_ref,r_y_ref,U,V,U_num,V_num,k(i),MN_ref,weights,Mb,Nb);
    b_ref{i}=b_u(:);
    y_ref=cat(1,real(b_u(:)),imag(b_u(:)));
    B{i}=y_ref;
    
end


%% Initial Array definition (Candidate Locations)
M =200;                                                               % Number of elements in x direction of the array
N =200;                                                               % Number of elements in y direction of the array
MN=M*N;                                                               % total number of elements in the array
d_x = ((M_ref-1)*lambda_ref)/(2*(M-1));                               % interelement spacing in x direction
d_y =((N_ref-1)*lambda_ref)/(2*(N-1));                                % interelement spacing in y direction

% Creating the position vector of sensors along x-axis

r_xx = d_x*[-(M-1)/2:1:(M-1)/2]';

% Creating the position vector of sensors along y-axis

r_yy = d_y*[-(N-1)/2:1:(N-1)/2]';

[r_x,r_y]=meshgrid(r_xx,r_yy);
r_x=r_x(:);                                                          % x_index of sensors
r_y=r_y(:);                                                          % y_index of sensors


% Steering Matrix for compressive sensing model

for j=1:length(f)
    [A] =Steering_Matrix_CS_array_syn_Planar( r_x,r_y,U,V,U_num,V_num,k(j),MN);
    S1{j}=A';
    PHI_1=cat(1,real(S1{j}),imag(S1{j}));
    PHI_2=cat(1,-1*imag(S1{j}),real(S1{j}));
    PHI=cat(2,PHI_1,PHI_2);
    S{j}=PHI;
end

a = 1e3/0.1; b = 1;

tic
weights1 = mt_CS(S,B,a,b,1e-7);                                     % MTBCS Algorithm
toc



