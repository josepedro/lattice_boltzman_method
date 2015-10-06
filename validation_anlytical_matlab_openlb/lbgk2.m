%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a simple implementation of the LBGK model.
%% By Andrey R. da Silva, September 2014 WorkS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, clc          % Clears everything from the memory buffer in Matlab
close all               % Close all previously open figures

%%%%%%%%
%%%%%%%%
%%%%%%%%        1) SETTING LATTICE AND PHYSICAL PARAMETERS
%%%%%%%%
%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nl = 300;                                         % Number of lines   (cells in the y direction)
Nc = 300;                                         % Number of columns (cells in the x direction)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Physical parameters (macro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 340;                                         % Sound velocity on the fluid [m/s]
rho_p = 1.2;                                      % physical density [kg/m^3]
rho_p = 1;                                        % Fluid density  [kg/m^3]
Lx = 1.5;                                           % Maximun dimenssion in the x direction [m]
Ly = 1.5;                                           % Maximun dimenssion on th y direction  [m]
Dx = Lx/Nc     ;                                   % Lattice space (pitch)
Dt = (1/sqrt(3))*Dx/c_p               ;             % lattice time step
freq = (1/sqrt(3))/20;
lamb = c_p/freq;                                 % Physical wavelength [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice parameters (micro - lattice unities)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega = 1.95;                                     % Relaxation frequency
tau = 1/omega;                                    % Relaxation time
rho_l = 1;                                        % avereged fluid density (latice density
cs = 1/sqrt(3);                                   % lattice speed of sound
cs2 = cs^2;                                      % Squared speed of sound cl^2
visc = cs2*(1/omega-0.5);                         % lattice viscosity
visc_phy = visc*(Dx^2)/Dt;                        % physical kinematic viscosity
ppw = lamb/Dx    ;                                % number of cells per wavelength

% defining a 2-cell-thick buffer frame around the lattice grid
% Nr and Mc are the effective size of the grid including the buffer

Nr = Nl+2;                                 
Mc = Nc+2;

steps = 250;                                       % number of time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice properties for the D2Q9 model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_c=9 ;                                             % number of directions
C_x=[1 0 -1  0 1 -1 -1  1 0];                       % velocity vectors in x
C_y=[0 1  0 -1 1  1 -1 -1 0];                       % velocity vectors in y
w0=16/36. ; w1=4/36. ; w2=1/36.;                    % lattice weights
f1=3.;
f2=4.5;
f3=1.5;                                             % coef. of the f equil.
f=zeros(Nr,Mc,N_c);                                 % allocating space for f
feq=zeros(Nr,Mc,N_c);                               % allocating space for feq

%%%%%%%%
%%%%%%%%
%%%%%%%%        2)  DEFINING INITIAL CONDITIONS
%%%%%%%%
%%%%%%%%

f(:,:,:)=rho_l/9;                                  % Filling the initial distribution function (at t=0) with initial values 
%f(Nl/2,Nc/2,9)=(1.001); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=0.001;

for ta = 1 : steps    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3) PROPAGATION (streaming)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    f(:,:,1) = [f(:,1:2,1) f(:,2:Mc-1,1)];
    f(:,:,2) = [f(1:2,:,2);f(2:Nr-1,:,2)];
    f(:,:,3) = [f(:,2:Mc-1,3) f(:,Mc-1:Mc,3)];    
    f(:,:,4) = [f(2:Nr-1,:,4);f(Nr-1:Nr,:,4)];
    f(:,:,5) = [f(:,1:2,5) f(:,2:Mc-1,5)];
    f(:,:,5) = [f(1:2,:,5);f(2:Nr-1,:,5)];
    f(:,:,6) = [f(:,2:Mc-1,6) f(:,Mc-1:Mc,6)];
    f(:,:,6) = [f(1:2,:,6);f(2:Nr-1,:,6)];
    f(:,:,7) = [f(:,2:Mc-1,7) f(:,Mc-1:Mc,7)];
    f(:,:,7) = [f(2:Nr-1,:,7);f(Nr-1:Nr,:,7)];
    f(:,:,8) = [f(:,1:2,8) f(:,2:Mc-1,8)];
    f(:,:,8) = [f(2:Nr-1,:,8);f(Nr-1:Nr,:,8)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4) CALCULATING NEW VALUES FOR rho and u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculating rho along the entire lattice
    rho=sum(f,3);   

    rt0= w0*rho;
    rt1= w1*rho;
    rt2= w2*rho;
     
    % Determining the velocities according to Eq.() (see slides)    
    ux = (C_x(1).*f(:,:,1)+C_x(2).*f(:,:,2)+C_x(3).*f(:,:,3)+C_x(4).*f(:,:,4)+C_x(5).*f(:,:,5)+C_x(6).*f(:,:,6)+C_x(7).*f(:,:,7)+C_x(8).*f(:,:,8))./rho ;
    uy = (C_y(1).*f(:,:,1)+C_y(2).*f(:,:,2)+C_y(3).*f(:,:,3)+C_y(4).*f(:,:,4)+C_y(5).*f(:,:,5)+C_y(6).*f(:,:,6)+C_y(7).*f(:,:,7)+C_y(8).*f(:,:,8))./rho ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 5) DETERMINING FEQ IN EACH DIRECTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    uxsq=ux.^2; 
    uysq=uy.^2; 
    usq=uxsq+uysq; 
       
    feq(:,:,1)= rt1 .*(1 +f1*ux +f2.*uxsq -f3*usq);
    feq(:,:,2)= rt1 .*(1 +f1*uy +f2*uysq -f3*usq);
    feq(:,:,3)= rt1 .*(1 -f1*ux +f2*uxsq -f3*usq);
    feq(:,:,4)= rt1 .*(1 -f1*uy +f2*uysq -f3*usq);
    feq(:,:,5)= rt2 .*(1 +f1*(+ux+uy) +f2*(+ux+uy).^2 -f3.*usq);
    feq(:,:,6)= rt2 .*(1 +f1*(-ux+uy) +f2*(-ux+uy).^2 -f3.*usq);
    feq(:,:,7)= rt2 .*(1 +f1*(-ux-uy) +f2*(-ux-uy).^2 -f3.*usq);
    feq(:,:,8)= rt2 .*(1 +f1*(+ux-uy) +f2*(+ux-uy).^2 -f3.*usq);
    feq(:,:,9)= rt0 .*(1 - f3*usq);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 6) PERFORMING COLLISION (RELAXATION)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    f= (1-omega)*f + omega*feq;
    %f(Nl/2,Nc/2,9)= f((Nl/2),Nc/2,9) + A*sin(2*pi*(cs/20)*ta);
    f(Nl/2,Nc/2,9)= f((Nl/2),Nc/2,9) + A*sin(2*pi*freq*ta);
    %F(Nli/2,Ncol/2,9)= F((Nli/2),Ncol/2,9) +Amp*sin(2*pi*(cl/20)*passo);

% % % % pltting the density disturbance rho_a over the lattice for each time step
 %surf(rho-1), view(2), shading flat, axis equal, caxis([-A,+A]), colorbar, colormap(hsv)
 %grid off
 %axis off
 %pause(.00001)
end %  End main time Evolution Loop

densities_vector(1:Nl/2) = 0;
for element_grid = 1:Nl/2
    densities_vector(element_grid) = rho(Nc/2,element_grid+149);
end
pressure_vector = (densities_vector - 1) * cs2;

%figure;
%plot(densities_vector);  


file = load('openlb_rho.mat');
openlb_rho = file.openlb_rho;             
openlb_pressure = (openlb_rho(1:length(openlb_rho)-1) - 1) * cs2;

[p pos]=cylin_wave((1/sqrt(3))/20,visc_phy,1/sqrt(3),A/20,1:Nc/2,pi/2);
figure;
%plot(pos,p,'b', [1:Nl/2], pressure_vector, 'r');
plot(pos,p,'b', [1:Nl/2], pressure_vector, 'r',...
[1:Nl/2], openlb_pressure, 'g');
legend('Analytical Solution','Matlab Simulation','OpenLB Simulation');

%load densi.mat densi
% figure;
% 
% plot(rho(Nl/2,1:Nc/2),'b');
% hold on
% 
% plot(densi(Nl/2,Nc/2:Nc),'r');
%save rho.mat rho
