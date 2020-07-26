clear;                           %#ok<*NASGU>

% set parameter space to be tested
NN  = 200;
PU  = [1/8 ,1/8,1/4,1/4,1/2,1/2,1  ,1,2,2,4,4,8,8];
TY  = [1/16,1/8,1/8,1/4,1/4,1/2,1/2,1,1,2,2,4,4,8];

for i=1:length(TY)
    
    % set run parameters
    runID    = ['Pux',num2str(PU(i)),'_Ty',num2str(TY(i)),'_N',num2str(NN)]; runID(runID=='.') = '_'; % run identifier
    restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
    nop      =  10;                  % output frame plotted/saved every 'nop' time steps
    plot_op  =  1;                   % switch on (1) to display results
    save_op  =  1;                   % switch on (1) to save output to file
    plot_cv  =  0;                   % switch on (1) to live plot iterative convergence
    bnchmrk  =  0;                   % switch on (1) to run solitary wave benchmark
    
    % set model domain parameters
    L        =  100;                 % domain dimension
    N        =  NN + 2;              % number of grid points in z-direction (incl. 2 ghosts)
    h        =  L/(N-2);             % grid spacing
    
    % set model timing parameters
    M        =  1e3;                 % number of time steps to take
    tend     =  1e3;                 % end time for simulation [s]
    dt       =  4e-6;                % initial time step [s]
    
    % set model liquid fraction parameters
    f0       =  0.02;                % background liquid fraction
    f1       =  0.10;                % amplitude of random noise
    f2       =  0.00;                % amplitude of gaussian
    smx      =  5/h^2;               % smoothness of initial random noise in x-direction
    smz      =  5/h^2;               % smoothness of initial random noise in z-direction
    wx       =  L/5;                 % horizontal half-width of gaussian
    wz       =  L/5;                 % vertical half-width of initial gaussian
    xpos     =  0;                   % x-position of initial gaussian (0 = middle of domain)
    zpos     =  0;                   % z-position of initial gaussian (0 = middle of domain)
    
    % set model rheology parameters
    n        =  3;                   % permeability powerlaw
    m        =  1;                   % compaction viscosity powerlaw
    lmd      =  30;                  % liquid-weakening parameter
    De       =  1e4;                 % visco-elastic Deborah number
    Ty       =  TY(i);               % Griffith criterion tensile strength
    
    % stress control parameters
    Pu       = -PU(i);               % ratio of pure-shear stress to buoyancy pressure
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
    gamma    =  1/3;                 % rheology lag parameter [0,1]
    delta    =  1.00;                % regularisation of eIIvp for failure [0,1]
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
    
end
