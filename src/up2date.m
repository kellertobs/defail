% update tensor magnitudes
eps(2:end-1,2:end-1) = (  (exx(2:end-1,2:end-1).^2 + ezz(2:end-1,2:end-1).^2 ...
                     + 2.*(exz(1:end-1,1:end-1).^2 + exz(2:end  ,1:end-1).^2 ...
                         + exz(1:end-1,2:end  ).^2 + exz(2:end  ,2:end  ).^2).*0.25)./2).^0.5 + 1e-16;
eps([1 end],:) = eps([end-1 2],:);                                         % periodic boundaries
eps(:,[1 end]) = eps(:,[end-1 2]);

tau(2:end-1,2:end-1) = (  (txx(2:end-1,2:end-1).^2 + tzz(2:end-1,2:end-1).^2 ...
                     + 2.*(txz(1:end-1,1:end-1).^2 + txz(2:end  ,1:end-1).^2 ...
                         + txz(1:end-1,2:end  ).^2 + txz(2:end  ,2:end  ).^2).*0.25)./2).^0.5 + 1e-16;
tau([1 end],:) = tau([end-1 2],:);                                         % periodic boundaries
tau(:,[1 end]) = tau(:,[end-1 2]);

% update yield criterion
tt   = 1-(f-1)*f0;                                                         % tensile strength reduced by porosity

fyt  = reshape(([p(:)  tau(:)] + [1 -1  ] - [-  tt(:) 0*tt(:)])/[1 1],N,N);% tensile yield factor
tyt = max( 1e-6,   0 + fyt );                                              % tensile yield stress
pyt = max( -tt , -tt + fyt );                                              % tensile yield pressure

tys = max( 1e-6,2*tt + p/2 );                                              % shear yield stress

ty  = min(tyt,tys) + bnchmrk*5;                                            % combined yield stress
py  = pyt + bnchmrk*5;                                                     % combined yield pressure

% update viscosities
etav  = exp(-lmd.*f0.*(f-1)) .* (1/2 + (eps./eps0/2).^-((1-n)./n)).^-1;    % shear viscosity
zetav = etav./max(-6,f0.*f);                                               % compaction viscosity

% update visco-plasticities
etay  = ty./max(1e-9,eps) + etamin;                                        % shear visco-plasticity
etay  = min(etav,etay);

zetay = -min(-1e-6,py)./max(1e-9,ups) + etamin./f/f0;                      % compaction visco-plasticity
zetay = min(zetav,zetay);

% apply regularisation
for k  = 1:ceil(kappa)                                                     
    kk = kappa/ceil(kappa);
    etay(2:end-1,2:end-1) = etay(2:end-1,2:end-1) + kk.*(diff(etay(2:end-1,:),2,2)+diff(etay(:,2:end-1),2,1))./8;
    etay([1 end],:) = etay([end-1 2],:);
    etay(:,[1 end]) = etay(:,[end-1 2]);

    zetay(2:end-1,2:end-1) = zetay(2:end-1,2:end-1) + kk.*(diff(zetay(2:end-1,:),2,2)+diff(zetay(:,2:end-1),2,1))./8;
    zetay([1 end],:) = zetay([end-1 2],:);
    zetay(:,[1 end]) = zetay(:,[end-1 2]);
end

% relax visco-plasticity update
 eta  =   etay.*(1-gamma) +  eta.*gamma;                                   % effective shear viscosity
zeta  =  zetay.*(1-gamma) + zeta.*gamma;                                   % effective compaction viscosity

etac = (eta(1:end-1,1:end-1)+eta(2:end,1:end-1) ...                        % evaluate in cell corners
     +  eta(1:end-1,2:end  )+eta(2:end,2:end  )).*0.25;

% get plastic strain rates
epsp = eps - tau./ etav;
upsp = ups +   p./zetav;

% update segregation coefficients
K    = f.^m + (epsp + upsp)/10;                                            % segregation coefficient

Kx   = (K(:,1:end-1)+K(:,2:end)).*0.5;                                     % evaluate on cell faces
Kz   = (K(1:end-1,:)+K(2:end,:)).*0.5;

% update iterative pseudo-time step sizes
dtW = (max( etav(1:end-1,:), etav(2:end,:))./(h/2)^2 ...
    +  max(zetav(1:end-1,:),zetav(2:end,:))./(h/2)^2).^-1;                 % W iterative step size
dtU = (max( etav(:,1:end-1), etav(:,2:end))./(h/2)^2 ...
    +  max(zetav(:,1:end-1),zetav(:,2:end))./(h/2)^2).^-1;                 % U iterative step size
dtP = (1./zetav + K./(h/2)^2).^-1;                                         % P iterative step size
