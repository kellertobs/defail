clear;                           %#ok<*NASGU> 

% set run parameters
runID    = 'Puc_4_Ty_4';         % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  20;                  % output frame plotted/saved every 'nop' time steps 
plot_op  =  0;                   % switch on (1) to display results
save_op  =  1;                   % switch on (1) to save output to file
plot_cv  =  0;                   % switch on (1) to live plot iterative convergence
bnchmrk  =  0;                   % switch on (1) to run solitary wave benchmark

% set model domain parameters
L        =  100;                 % domain dimension
N        =  200 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  L/(N-2);             % grid spacing
 
% set model timing parameters
M        =  1e4;                 % number of time steps to take
tend     =  10;                  % end time for simulation [s]
dt       =  1e-5;                % initial time step [s]

% set model liquid fraction parameters
f0       =  0.01;                % background liquid fraction
f1       =  0.10;                % amplitude of random noise
f2       =  0.00;                % amplitude of gaussian
smx      =  10/h^2;              % smoothness of initial random noise in x-direction
smz      =  10/h^2;              % smoothness of initial random noise in z-direction
wx       =  L/5;                 % horizontal half-width of gaussian
wz       =  L/5;                 % vertical half-width of initial gaussian
xpos     =  0;                   % x-position of initial gaussian (0 = middle of domain)
zpos     =  0;                   % z-position of initial gaussian (0 = middle of domain)
 
% set model rheology parameters
n        =  3;                   % permeability powerlaw
m        =  1;                   % compaction viscosity powerlaw
lmd      =  30;                  % liquid-weakening parameter
De       =  1e3;                 % visco-elastic Deborah number
Ty       =  1/4;                 % Griffith criterion tensile strength

% stress control parameters
Pu       =  1/4;                % ratio of pure-shear stress to buoyancy pressure
Si       =  0;                   % ratio of simple-shear stress to buoyancy pressure

% set numerical model parameters
nup      =  100;                 % nonlinear coefficients, residual norms updated every 'nup' iterations
CFL      =  1.0;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'FRM';               % advection scheme ('UPW2', 'UPW3', or 'FRM')
theta    =  0.50;                % time-stepping scheme selector (1=BE, 1/2=CN, 0=FE)
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-7;                % outer its absolute tolerance
maxit    =  1e4;                 % maximum outer its
alpha    =  0.95;                % inner its step size (multiple of stable step) [0,1]
beta     =  0.90;                % iterative damping parameter [0,1]
gamma    =  1/2;                 % rheology lag parameter [0,1]
delta    =  2.00;                % regularisation of eIIvp for failure [0,1]
etamin   =  1e-2;                % minimum viscosity for regularisation

% create output directory
if ~isfolder(['../out/',runID])
    mkdir(['../out/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 
    parfile = ['../out/',runID,'/',runID,'_par'];
    save(parfile);
end

% run code
addpath ../src
main
