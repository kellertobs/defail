clear; close all;                % #ok<*NASGU>

% set parameter space to be tested
NN  =  [50,100,200];
EE  =  0.*NN;

for ii = 1:length(NN)
    
    % set run parameters
    runID    = 'test';               % run identifier
    restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
    nop      =  1;                   % output frame plotted/saved every 'nop' time steps
    plot_op  =  1;                   % switch on (1) to display results
    save_op  =  0;                   % switch on (1) to save output to file
    plot_cv  =  1;                   % switch on (1) to live plot iterative convergence
    bnchmrk  =  0;                   % switch on (1) to run solitary wave benchmark
    demean   =  1;                   % remove mean from solution fields
    
    % set model domain parameters
    L        =  100;                 % domain dimension
    N        =  NN(ii) + 2;          % number of grid points in z-direction (incl. 2 ghosts)
    h        =  L/(N-2);             % grid spacing
    
    % set model timing parameters
    M        =  1;                   % number of time steps to take
    tend     =  1;                   % end time for simulation [s]
    
    % set model liquid fraction parameters
    f0       =  0.05;                % background liquid fraction
    f1       =  0.20;                % amplitude of random noise
    f2       =  0.00;                % amplitude of gaussian
    smx      =  (N/50)^2;            % smoothness of initial random noise in x-direction
    smz      =  (N/50)^2;            % smoothness of initial random noise in z-direction
    wx       =  L/5;                 % horizontal half-width of gaussian
    wz       =  L/5;                 % vertical half-width of initial gaussian
    xpos     =  0;                   % x-position of initial gaussian (0 = middle of domain)
    zpos     =  0;                   % z-position of initial gaussian (0 = middle of domain)
    
    % set model rheology parameters
    n        =  3;                   % permeability powerlaw
    m        =  1;                   % compaction viscosity powerlaw
    nn       =  1;                   % non-Newtonian shear viscosity powerlaw
    lmd      =  30;                  % shear viscosity liquid-weakening parameter
    
    % stress control parameters
    Pu       =  0;                   % ratio of pure-shear stress to buoyancy pressure
    Si       =  1;                   % ratio of simple-shear stress to buoyancy pressure
    B        =  1;                   % ratio of buoyancy force to
    
    % set numerical model parameters
    nup      =  20;                  % nonlinear coefficients, residual norms updated every 'nup' iterations
    CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
    ADVN     =  'FRM';               % advection scheme ('UPW2', 'UPW3', or 'FRM')
    theta    =  0.50;                % time-stepping scheme selector (1=BE, 1/2=CN, 0=FE)
    rtol     =  1e-12;               % outer its relative tolerance
    atol     =  1e-12;               % outer its absolute tolerance
    minit    =  500;                 % minimum solver iterations
    maxit    =  1e5;                 % maximum solver iterations
    alpha    =  0.98;                % inner its step size (fraction of stable step) [0,1]
    beta     =  0.9-4./(N-2);        % iterative damping parameter (fraction of previous step) [0,1]
    gamma    =  0.00;                % iterative relaxation for rheology updates [0,1]
    kappa    =  0.00;                % regularisation of eIIvp for failure [0,1]
    etamin   =  1e-16;               % minimum viscosity for regularisation
    
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
    
    % record and plot numerical error
    EE(ii) = err;
    
    figure(18);
    EL = EE(1).*(NN(1)./NN).^1;
    EQ = EE(1).*(NN(1)./NN).^2;
    EC = EE(1).*(NN(1)./NN).^3;
    HH = L./NN;
    loglog(HH,EE,'r*',HH,EL,'k--',HH,EQ,'k-',HH,EC,'k-.','LineWidth',1.5,'MarkerSize',10); hold on;
    legend('model results','first order','second order','third order','Interpreter','latex','FontSize',16,'Location','SouthEast','box','off');
    set(gca,'TickLabelInterpreter','latex','FontSize',13,'LineWidth',2);
    xlabel('grid step size','Interpreter','latex','FontSize',16);
    ylabel('numerical error','Interpreter','latex','FontSize',16);
    title('Benchmark to Manufactured Solution','Interpreter','latex','FontSize',20);
    drawnow;
    
end
