% D2Q9 Lattice Boltzmann Bhatnagar–Gross–Krook (LBGK)
% Code developed by Tridip Das No slip flow in pipe
clear all;
close all;
clc;
%
disp ('D2Q9 LBGK model for pipe flow simulation');
% Define density & viscosity
rho=1; % gm/cc
viscosity = 1.0; % dynamic viscosity (cP)
viscosity = viscosity * 0.01; % g/(cm.s)
nu = viscosity/rho; % kinematic viscosity (cm^2/s)
del_x = 1; % grid spacing (cm)
del_t = 1; % time step (s)
b = 6; % vel of sound
D = 2;
d0 = 1/2;
tau = (1 + b*nu*(del_t/del_x)^2)/2;
w1=4/9; % weight factor for central component
w2=1/9; % weight factor for edge component
w3=1/36; % weight factor for diagonal component
c_squ=1/3;
%
nx=31; % grid points in x-direction
ny=31; % grid points in y-direction
F=repmat(rho/9,[nx ny 9]); % Initialize lattice
FEQ=F; % Initialize equilibriation point
meshsize=nx*ny; 
CI=[0:meshsize:meshsize*7]; % Repeated for each 8 components
BOUND(1:nx,[1 ny])=1; % Block the wall boundary
ON=find(BOUND); %matrix offset of each Occupied Node
TO_REFLECT=[ON+CI(1) ON+CI(2) ON+CI(3) ON+CI(4) ...
            ON+CI(5) ON+CI(6) ON+CI(7) ON+CI(8)];
REFLECTED= [ON+CI(5) ON+CI(6) ON+CI(7) ON+CI(8) ...
            ON+CI(1) ON+CI(2) ON+CI(3) ON+CI(4)];
%        
avg_u=1; 
prev_avg_u=1; 
time_steps=0; 
deltaU=1e-4; 
num_active_nodes=sum(sum(1-BOUND));
%% Time steps iteration loop
while (time_steps<50000 & 1e-5<abs((prev_avg_u-avg_u)/avg_u)) | time_steps<100
    %% Propagation of particles
    F(:,:,4)=F([2:nx 1],[ny 1:ny-1],4);
    F(:,:,3)=F(:,[ny 1:ny-1],3);
    F(:,:,2)=F([nx 1:nx-1],[ny 1:ny-1],2);
    F(:,:,5)=F([2:nx 1],:,5);
    F(:,:,1)=F([nx 1:nx-1],:,1);
    F(:,:,6)=F([2:nx 1],[2:ny 1],6);
    F(:,:,7)=F(:,[2:ny 1],7); 
    F(:,:,8)=F([nx 1:nx-1],[2:ny 1],8);
    BOUNCEDBACK=F(TO_REFLECT); %Densities bouncing back at next timestep
    density=sum(F,3);
    %% velocity calculations
    UX=(sum(F(:,:,[1 2 8]),3)-sum(F(:,:,[4 5 6]),3))./density;
    UY=(sum(F(:,:,[2 3 4]),3)-sum(F(:,:,[6 7 8]),3))./density;
    UX(1,1:ny)=UX(1,1:ny)+deltaU; %Increase inlet pressure
    UX(ON)=0; UY(ON)=0; density(ON)=0;
    U_SQU=UX.^2+UY.^2; 
    U_C2=UX+UY; 
    U_C4=-UX+UY; 
    U_C6=-U_C2; 
    U_C8=-U_C4;
    %% Calculate equilibrium distribution
    % stationary point
    FEQ(:,:,9)=w1*density.*(1-U_SQU/(2*c_squ));
    % nearest neighbours at the edge
    FEQ(:,:,1)=w2*density.*(1+UX/c_squ+0.5*(UX/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,3)=w2*density.*(1+UY/c_squ+0.5*(UY/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,5)=w2*density.*(1-UX/c_squ+0.5*(UX/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,7)=w2*density.*(1-UY/c_squ+0.5*(UY/c_squ).^2-U_SQU/(2*c_squ));
    % next nearest neighbours at the diagonal
    FEQ(:,:,2)=w3*density.*(1+U_C2/c_squ+0.5*(U_C2/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,4)=w3*density.*(1+U_C4/c_squ+0.5*(U_C4/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,6)=w3*density.*(1+U_C6/c_squ+0.5*(U_C6/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,8)=w3*density.*(1+U_C8/c_squ+0.5*(U_C8/c_squ).^2-U_SQU/(2*c_squ));
    F=F - (F - FEQ)/tau;
    F(REFLECTED)=BOUNCEDBACK;
    prev_avg_u=avg_u;
    avg_u=sum(sum(UX))/num_active_nodes; 
    time_steps=time_steps+1;
end
%% Flow field with arrow head plot
figure;
colormap(gray(2));image(2-BOUND');hold on;
quiver(1:nx,1:ny,UX(1:nx,:)',UY(1:nx,:)');
title(['Flow field after ',num2str(time_steps),' \Deltat']);xlabel('x cm');ylabel('y cm');
set(gcf,'Color',[1,1,1]);
%% surface plot
figure;
surfc(1:nx,1:ny,U_SQU(1:nx,:)','EdgeColor','none','LineWidth',2);
%Flow through random grid of cells
u_max = max(U_SQU(ny,:));
u_avg = u_max/2;
NRe = ny*u_avg/nu
uy(30) = zeros;
uy(2:31) = u_max*(1-((1:ny-1)-(ny-1)/2).^2/((ny-1)/2)^2);
%% Velocity profile
figure;
plot(1:ny,U_SQU(ny,:),1:ny,uy,'LineWidth',2)
xlabel('Pipe diameter (cm)')
ylabel('velocity (cm/s)')
legend('LBGK','Hagen–Poiseuille')
set(gcf,'Color',[1,1,1]);
%% End
