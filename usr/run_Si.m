clear; close all;                %#ok<*NASGU>

% set parameter space to be tested
NN  =  500;
PU  =  [0,0,0,0,0];
SI  =  [1/100, 1/10, 1, 10, 100];
TY  =  [1/10 , 1/10, 1, 10, 10 ];

for i=1:length(TY)
    
    % set run parameters
    runID    = ['Si',num2str(abs(SI(i))),'_Ty',num2str(TY(i)),'_N',num2str(NN)]; runID(runID=='.') = '_'; % run identifier
    restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
    nop      =  10;                  % output frame plotted/saved every 'nop' time steps
    plot_op  =  1;                   % switch on (1) to display results
    save_op  =  1;                   % switch on (1) to save output to file
    plot_cv  =  1;                   % switch on (1) to live plot iterative convergence
    bnchmrk  =  0;                   % switch on (1) to run solitary wave benchmark
    demean   =  1;                   % remove mean from solution fields
    
    % set model domain parameters
    L        =  100;                 % domain dimension
    N        =  NN + 2;              % number of grid points in z-direction (incl. 2 ghosts)
    h        =  L/(N-2);             % grid spacing
    
    % set model timing parameters
    M        =  1e3;                 % number of time steps to take
    tend     =  100;                 % end time for simulation [s]
    
    % set model liquid fraction parameters
    f0       =  0.01;                % background liquid fraction
    f1       =  0.25;                % amplitude of random noise
    f2       =  0.00;                % amplitude of gaussian
    smx      =  L/20/h^2;            % smoothness of initial random noise in x-direction
    smz      =  L/20/h^2;            % smoothness of initial random noise in z-direction
    wx       =  L/5;                 % horizontal half-width of gaussian
    wz       =  L/5;                 % vertical half-width of initial gaussian
    xpos     =  0;                   % x-position of initial gaussian (0 = middle of domain)
    zpos     =  0;                   % z-position of initial gaussian (0 = middle of domain)
    
    % set model rheology parameters
    n        =  3;                   % permeability powerlaw
    m        =  1;                   % compaction viscosity powerlaw
    nn       =  3;                   % non-Newtonian shear viscosity powerlaw
    lmd      =  30;                  % shear viscosity liquid-weakening parameter
    
    % stress control parameters
    Ty       =  TY(i);               % Griffith criterion tensile strength
    Pu       =  PU(i);               % ratio of pure-shear stress to buoyancy pressure
    Si       =  SI(i);               % ratio of simple-shear stress to buoyancy pressure
    B        =  1;                   % relative magnitude of buoyancy force

    % set numerical model parameters
    nup      =  20;                  % nonlinear coefficients, residual norms updated every 'nup' iterations
    CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
    ADVN     =  'FRM';               % advection scheme ('UPW2', 'UPW3', or 'FRM')
    theta    =  0.50;                % time-stepping scheme selector (1=BE, 1/2=CN, 0=FE)
    rtol     =  1e-6;                % outer its relative tolerance
    atol     =  1e-6;                % outer its absolute tolerance
    minit    =  500;                 % maximum outer its
    maxit    =  3000;                % maximum outer its
    alpha    =  0.90;                % inner its step size (fraction of stable step) [0,1]
    beta     =  0.80;                % iterative damping parameter (fraction of previous step) [0,1]
    gamma    =  0.25;                % iterative relaxation for rheology updates [0,1]
    kappa    =  0.25;                % regularisation of eIIvp for failure [0,1]
    etamin   =  0.01;                % minimum viscosity for regularisation
    
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