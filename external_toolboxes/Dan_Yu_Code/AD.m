clear all;
close all;
%Source and receptors are not on the boundaries (except z = 0)
%10^5 system, Full State Measurements
%%
%SET PARAMETERS : Parameters can be changed
%Computational domain
xlim                            =   [ 0, 2000];  %[600m - 3000m]
ylim                            =   [-100, 400];
zlim                            =   [0, 50];
%Wind Profile parameters
Uwind                           =   4;    % wind speed (m/s) between 1-5 (m/s)
alp                             =   0/180 * pi;
beta                            =   0/180 * pi;
%Diffusion parameters
ay                              =   0.08;
by                              =   0.0001;
az                              =   0.06;
bz                              =   0.0015;
%discretization dimension
nx                              =   10;
ny                              =   10;
nz                              =   10;
dt                              =   1;

%%
%U = (U cos(alp)cos(beta), U cos(alpha) sin(beta), U sin(alp));
Ux = Uwind * cos(alp)*cos(beta);
Uy = Uwind * cos(alp)*sin(beta);
Uz = Uwind * sin(alp);
%Finite difference Scheme parameters
dx =( xlim(2) - xlim(1))/nx;
dy = (ylim(2) - ylim(1))/ny;
dz = (zlim(2) - zlim(1))/nz;
x0 = xlim(1) + [0 : nx]*dx;           % distance along wind direction (m)
y0 = ylim(1) + [0 : ny]*dy;           % cross-wind distance (m)
z0 = zlim(1) + [0 : nz]*dz;
[xmesh, ymesh] = meshgrid(x0,y0);
Ky = zeros(nx+1,1);
Kz = zeros(nx+1, 1);
for i = 1 : nx+1
    Ky(i,1) = Uwind * ay^2/2 * (2 * x0(i) + by *(x0(i))^2)/(1+by*x0(i))^2;
    Kz(i,1) = Uwind * az^2/2 * (2 * x0(i) + bz *(x0(i))^2)/(1+bz*x0(i))^2;
end

sy = Ky * dt/(dy^2);
sz = Kz * dt/(dz^2);
cx = Ux * dt/dx;
cy = Uy * dt/(dy);
cz = Uz * dt/(dz);
%Stable
if 1 - 2 * max(sy+sz)-cx^2 -cx -cy-cz < 0
    disp('Error: unstable system');
    
else
    NUM_SYS = (ny-1)*(nx-1)*(nz);
    %Primal system
    [SR, S, SC] = atmo3dlax2 (Uwind, alp, beta, nx, ny, nz, xlim, ylim, zlim, Ky, Kz, dt);
    %Adjoint system : Sad = S;
    SRad = SC;
    SCad = SR;
    %%
    %Source and receptor locations and discretization
    % Stack emission source data:
    source.n = 10;                         % # of sources
    source.x = [280, 300, 900, 1100, 500,1300,500,1700,1400,700];     % x-location (m)
    source.y = [ 75, 205, 25,  185,  250,230, 330,170, 240, 200];     % y-location (m)
    
    source.z = [ 15,  35,  15,   15, 10, 10, 10, 10, 15, 20];     % height (m)
    %source.label=[' S1'; ' S2'; ' S3'; ' S4'; 'S5';'S6';'S7';'S8';'S9';'S10'];
    %source.n = 1;
    %source.x = 1.1;
    %source.y = 0;
    %source.z = 0;
    
    tpy2kgps = 1.0 / 31536;               % conversion factor (tonne/yr to kg/s)
    source.Q = [35, 80, 5, 5, 10, 5, 5, 10, 5, 10] * tpy2kgps ; % emission rate (kg/s)
    %source.Q = 1;
    % Stack discrete source data: Transfer from continuous to discrete
    dissource.n = source.n;
    dissource.x = floor((source.x - xlim(1))/dx)+1;
    dissource.y = floor((source.y - ylim(1))/dy)+1;
    dissource.z = floor((source.z - zlim(1))/dz)+1;
    %%
    %USE IF NUM_SYS IS SMALL.
    A = zeros(NUM_SYS, NUM_SYS);
    for i = 1 : length(S)
        A(SR(i), SC(i)) = S(i);
    end
    %Receptors
    recept.n = 9;                                                 % # of receptors
    recept.x = [  60,  76, 267, 331, 514, 904, 1288, 1254, 972 ]; % x location (m)
    recept.y = [ 130,  70, -21, 308, 182,  75,  116,  383, 507 ]; % y location (m)
    recept.z = [   0,  10,  10,   1,  15,   2,    3,   12,  12 ]; % height (m)
    recept.label=[ ' R1 '; ' R2 '; ' R3 '; ' R4 '; ' R5 '; ' R6 '; ...
        ' R7 '; ' R8 '; ' R9 ' ];
    disrecept.n = recept.n;
    disrecept.x = floor((recept.x - xlim(1))/dx)+1;
    disrecept.y = floor((recept.y - ylim(1))/dy)+1;
    disrecept.z = floor((recept.z - zlim(1))/dz)+1;
    
    
    NUM_OUT                 = recept.n;
    NUM_IN                  = source.n;
    
    B = zeros(NUM_SYS, source.n);
    for ns = 1 : source.n
        B((nx-1)*(ny-1)*(dissource.z(1,ns)-1)+(nx-1)*(dissource.y(1,ns)-2)+dissource.x(1,ns)-1, ns) = 1;
    end
    Ft = (1 - (source.x - x0(dissource.x))/dx).*(1 -(source.y - y0(dissource.y))/dy) .* (1- (source.z - z0(dissource.z))/dz) .* source.Q *dt/(dx*dy*dz);
    C = zeros(recept.n, NUM_SYS);
    for ns = 1 : recept.n
        C(ns,(nx-1)*(ny-1) * (disrecept.z(1,ns)-1) + (nx-1) * (disrecept.y(1,ns) -2) + disrecept.x(1,ns)-1) = 1;
    end
end
