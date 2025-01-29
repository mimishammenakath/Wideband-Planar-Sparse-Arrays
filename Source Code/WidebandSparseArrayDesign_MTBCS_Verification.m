%-- Script for finding MSE of the wideband sparse planar array design using MTBCS
%-- Authors: Mimisha M Menakath and Hareesh G
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified on 29 - January - 2025
%-------------------------------------------------------------------------%

c = 1500;                                                                % speed of sound in water
fH=450e3;
fL=150e3;
f0=300e3 ;
lambda_ref=c/fH;
theta_str=0;                                                             % Azimuth steering angle
phi_str= 0;                                                              % Elevation steering angle
look_dir=(-90:5:90);                                                     % steeing direction in degreee
indx=find(look_dir==theta_str);
look_dir_u =  sind(look_dir)-sind(theta_str);                            % steering direction in u domain
look_dir_v =  sind(look_dir)-sind(phi_str);                              % steering direction in v domain[U,V]=meshgrid(look_dir_u,look_dir_v);
[U,V]=meshgrid(look_dir_u,look_dir_v);
U=U(:);                                                                  % steering vector in u domain
V=V(:);                                                                  % steering vector in v domain
U_num=length(U);
V_num=length(V);
Mb= length(look_dir_u);                                                  % total number of steering angles in azimuth direction
Nb= length(look_dir_v);                                                  % total number of steering angles in elevation direction
M_ref=100;                                                               % number of elements in x direction of the reference array for maximum frequency
N_ref=100;                                                             % number of elements in y direction of the reference array for maximum frequency
MN_ref=M_ref*N_ref;                                                    % total number of elements in the reference array
L=(M_ref-1)*lambda_ref/2;                                              % Aperture size
delta_f=c/(sqrt(2)*L);                                                 % frequency resolution
f=fL:delta_f:fH;
lambda = c./f;                                                         % wavelength of sound in water
k = 2*pi./lambda;                                                      % wave number

%% Design of desired Beam Pattern

for i=1:length(f)
    d_ref =lambda_ref/2;                                               % element spacing in the reference array
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

%% Initial Array definition
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

R1=weights1(1:end/2,:);
I1=weights1(end/2+1:end,:);
Weight_new=R1+1i*I1;


for k=1:length(f)                                                          % frequency index
    W_new=Weight_new(:,k);
    Num=length(find(W_new~=0));                                            % Number of elements in the sparse array

    % creating Beam pattern of the sparse array

    B_u_CS = S1{k}*W_new;                                                  % sparse array pattern for kth frequency

    y=b_ref{k};                                                            % desired pattern for kth frequency,

    B_u_CS=reshape(B_u_CS,Mb,Nb);

    %Normalizing beampattern

    B_u_CS_norm =abs( B_u_CS)/max(abs(B_u_CS(:) ));

    % displying  3d beam pattern

    figure(1),
    mesh(look_dir_u,look_dir_v,20*log10(B_u_CS_norm));
    xlabel('Azimuth Angle in degrees');
    ylabel('Elevation Angle in degrees');
    zlabel(' Magnitude in dB')
    title('BeamPattern of Sparse Array')
    zlim([-60 0]), caxis([-60,0]);


    y=reshape(y,Mb,Nb);
    figure(2);
    mesh(look_dir_u,look_dir_v,20*log10(abs(y)./max(abs(y(:)))));
    xlabel('Azimuth Angle in degrees');
    ylabel('Elevation Angle in degrees');
    zlabel(' Magnitude in dB')
    title('BeamPattern of URA')
    zlim([-60 0]), caxis([-60,0]);

    % displying  2d beam pattern

    y1(:,k)=y(:,indx);
    y2(:,k)= B_u_CS(:,indx);
    % y22=y2(:,k);
    % plot(look_dir,20*log10(abs(y22)./max(y22(:))));ylim([-30,0]);
    % xlabel('Azimuth Angle in degrees');
    % ylabel('Magnitude in dB')
    % title('BeamPattern Comparison')
    % ylim([-60 0])

    W=reshape(W_new,M,N);
    figure(3),mesh(abs(W));
    figure(4),imagesc (abs(W));

    mse(k)=(sum(abs(B_u_CS(:)-y(:))).^2)./(sum(abs(y(:)).^2));
end

figure, plot(look_dir,20*log10(abs(y2)./max(y2(:))));ylim([-60,0]);
xlabel('Azimuth Angle (Degree)', 'FontSize',14);
ylabel('Magnitude (dB)', 'FontSize',14);
title('Frequency Invariant BeamPattern','FontSize',14);

figure, plot(look_dir,20*log10(abs(y1)./max(y1(:))));ylim([-60,0]);
xlabel('Azimuth Angle (Degree)', 'FontSize',14);
ylabel('Magnitude (dB)', 'FontSize',14);
title('Frequency Invariant BeamPattern','FontSize',14);

figure, surf(f,look_dir,20*log10(abs(y1)./max(abs(y1(:)))));zlim([-60,0]);
xlabel('Frequency (Hz)', 'FontSize',14);
ylabel('Azimuth Angle (Degree)', 'FontSize',14);
zlabel('Magnitude (dB)', 'FontSize',14);
title('Reference Wideband BeamPattern','FontSize',14); colorbar;caxis([-90,0])

figure,surf(f,look_dir,20*log10(abs(y2)./max(abs(y2(:)))));zlim([-60,0]);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Azimuth Angle (Degree)','FontSize',14);
zlabel('Magnitude (dB)','FontSize',14)
title('Synthesized Wideband BeamPattern','FontSize',14);colorbar;caxis([-90,0])

[locs] = find(W_new);
r_x_new=r_x(locs);
r_y_new=r_y(locs);
%Creating a new Weight vector
W_new2= zeros(size(W_new));
W_new2(locs)=W(locs);
figure,imagesc (r_x_new./lambda_ref,r_y_new./lambda_ref,abs(W));title('Weight Magnitude','FontSize',14);xlabel('Sensor Position in X direction (Lambda)','FontSize',14);
ylabel('Sensor Position in Y direction (Lambda)','FontSize',14);colorbar;
figure,scatter(r_x_new./lambda_ref,r_y_new./lambda_ref);title('Weight Locations','FontSize',14);xlabel('Sensor Position in X direction (Lambda)','FontSize',14);
ylabel('Sensor Position in Y direction (Lambda)','FontSize',14);

%finding minimum element spacing
D(:,1)=r_x_new./lambda_ref;
D(:,2)=r_y_new./lambda_ref;
d=pdist2(D,D);
mean_space=mean(d(:))
d(d==0)=nan;
[val,idx]=min(d(:));
val                                                                         % minimum spacing in lambda

mse_db=10*log10(mse);
figure, plot(f,mse_db);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('MSE (dB)','FontSize',14)
title('Error Plot','FontSize',14)

k=10
y11=y1(:,k);
figure, plot(look_dir,20*log10(abs(y11)./max(abs(y11(:)))));ylim([-30,0]);
hold on
y22=y2(:,k);
plot(look_dir,20*log10(abs(y22)./max(abs(y22(:)))));ylim([-30,0]);
hold on
y33=-3*ones(length(y11),1);
plot(look_dir,y33,'g');
xlabel('Azimuth Angle in degrees');
ylabel('Magnitude in dB')
title('BeamPattern Comparison')
ylim([-60 0])

x=1:Num;
Weight_new=Weight_new./max(Weight_new);
Weights_New=Weight_new(locs,:);
figure, mesh(f,x,abs(Weights_New));
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Sensor number','FontSize',14);
zlabel('Weight Magnitude','FontSize',14)
title('Weight Magnitude','FontSize',14);colorbar;
