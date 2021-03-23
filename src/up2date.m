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

yieldt = max(1e-16,1 + p) + etamin.*eps + bnchmrk*5;

% update rheological parameters
etav  =  log10( exp(-lmd.*f0.*(f-1)) .* (eps./eps0).^((1-n)./n) );         % shear viscosity

etay  =  log10(yieldt)-log10(eps);                                         % shear visco-plasticity
etay  =  min(etav,etay);

zetay =  etay - max(-6,log10(f0.*f));                                      % cmpct viscosity

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    etay(2:end-1,2:end-1) = etay(2:end-1,2:end-1) + kk.*(diff(etay(2:end-1,:),2,2)+diff(etay(:,2:end-1),2,1))./8;
    etay([1 end],:) = etay([end-1 2],:);
    etay(:,[1 end]) = etay(:,[end-1 2]);
end

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    zetay(2:end-1,2:end-1) = zetay(2:end-1,2:end-1) + kk.*(diff(zetay(2:end-1,:),2,2)+diff(zetay(:,2:end-1),2,1))./8;
    zetay([1 end],:) = zetay([end-1 2],:);
    zetay(:,[1 end]) = zetay(:,[end-1 2]);
end

 eta  =   etay.*(1-gamma) +  eta.*gamma;                                   % effective shear viscosity
zeta  =  zetay.*(1-gamma) + zeta.*gamma;                                   % effective shear viscosity

etac = (eta(1:end-1,1:end-1)+eta(2:end,1:end-1) ...                        % evaluate in cell corners
     +  eta(1:end-1,2:end  )+eta(2:end,2:end  )).*0.25;

K    = f.^m;                                                               % segregation coefficient

% update iterative and physical time step sizes
dtW = (10.^(( eta(1:end-1,:)+ eta(2:end,:)).*0.5)./(h/2)^2 ...
    +  10.^((zeta(1:end-1,:)+zeta(2:end,:)).*0.5)./(h/2)^2).^-1;           % W iterative step size
dtU = (10.^(( eta(:,1:end-1)+ eta(:,2:end)).*0.5)./(h/2)^2 ...
    +  10.^((zeta(:,1:end-1)+zeta(:,2:end)).*0.5)./(h/2)^2).^-1;           % U iterative step size
dtP = (1./(10.^eta) + K./(h/2)^2).^-1;                                     % P iterative step size
