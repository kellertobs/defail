clear;                           %#ok<*NASGU> 

% set run parameters
runID    = 'demo';               % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  10;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on (1) to display results
save_op  =  1;                   % switch on (1) to save output to file
plot_cv  =  1;                   % switch on (1) to live plot iterative convergence
bnchmrk  =  0;                   % switch on (1) to run solitary wave benchmark
demean   =  1;                   % remove mean from solution fields

% set model domain parameters
L        =  80;                  % domain dimension
N        =  200 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  L/(N-2);             % grid spacing

% set model timing parameters
M        =  1e3;                 % number of time steps to take
tend     =  1e3;                 % end time for simulation [s]
dt       =  1e-6;                % initial time step [s]

% set model liquid fraction parameters
f0       =  0.01;                % background liquid fraction
f1       =  0.50;                % amplitude of random noise
f2       =  0.00;                % amplitude of gaussian
smx      =  3/h^2;               % smoothness of initial random noise in x-direction
smz      =  3/h^2;               % smoothness of initial random noise in z-direction
wx       =  L/5;                 % horizontal half-width of gaussian
wz       =  L/5;                 % vertical half-width of initial gaussian
xpos     =  0;                   % x-position of initial gaussian (0 = middle of domain)
zpos     =  0;                   % z-position of initial gaussian (0 = middle of domain)

% set model rheology parameters
n        =  3;                   % permeability powerlaw
m        =  1;                   % compaction viscosity powerlaw
nn       =  3;                   % non-Newtonian shear viscosity powerlaw
lmd      =  30;                  % shear viscosity liquid-weakening parameter
Ty       =  1;                   % Griffith criterion tensile strength

% stress control parameters
Pu       =  0;                   % ratio of pure-shear stress to buoyancy pressure
Si       =  1;                   % ratio of simple-shear stress to buoyancy pressure

% set numerical model parameters
nup      =  50;                  % nonlinear coefficients, residual norms updated every 'nup' iterations
CFL      =  0.5;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'FRM';               % advection scheme ('UPW2', 'UPW3', or 'FRM')
theta    =  0.50;                % time-stepping scheme selector (1=BE, 1/2=CN, 0=FE)
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-6;                % outer its absolute tolerance
maxit    =  5e3;                 % maximum outer its
alpha    =  0.95;                % inner its step size (multiple of stable step) [0,1]
beta     =  alpha-15/N;          % iterative damping parameter [0,1]
gamma    =  0.50;                % rheology lag parameter [0,1]
delta    =  0.10;                % plastic strain-rate localisaton sharpening powerlaw (> 0)
eps      =  0.50;                % softening plastic stress cut-off (1 = harmonic sum, 0 = min() cut-off)
kappa    =  0.25;                % regularisation of eIIvp for failure [0,1]
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
