% update tensor magnitudes
eps(2:end-1,2:end-1) = (  (exx(2:end-1,2:end-1).^2 + ezz(2:end-1,2:end-1).^2 ...
                     + 2.*(exz(1:end-1,1:end-1).^2 + exz(2:end  ,1:end-1).^2 ...
                         + exz(1:end-1,2:end  ).^2 + exz(2:end  ,2:end  ).^2).*0.25)./2).^0.5;
eps([1 end],:) = eps([end-1 2],:);                                         % periodic boundaries
eps(:,[1 end]) = eps(:,[end-1 2]);

tau(2:end-1,2:end-1) = (  (txx(2:end-1,2:end-1).^2 + tzz(2:end-1,2:end-1).^2 ...
                     + 2.*(txz(1:end-1,1:end-1).^2 + txz(2:end  ,1:end-1).^2 ...
                         + txz(1:end-1,2:end  ).^2 + txz(2:end  ,2:end  ).^2).*0.25)./2).^0.5;
tau([1 end],:) = tau([end-1 2],:);                                         % periodic boundaries
tau(:,[1 end]) = tau(:,[end-1 2]);

yieldt = max(    1e-3,(p + tau + Ty)/2);                                   % yield pressure
yieldp = max(-Ty+1e-3,(p + tau - Ty)/2);                                   % yield pressure

% update rheological parameters
etav   =  exp(-lmd.*f0.*(f-1)) ...
      .* (eps./eps0./rmp).^((1-nn)./nn);                                   % shear viscosity

etay  =  yieldt./max(1e-16,eps);                                           % shear visco-plasticity

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    etay(2:end-1,2:end-1) = etay(2:end-1,2:end-1) + kk.*(diff(etay(2:end-1,:),2,2)+diff(etay(:,2:end-1),2,1))./8;
    etay([1 end],:) = etay([end-1 2],:);
    etay(:,[1 end]) = etay(:,[end-1 2]);
end

eta  =  (1-gamma).*(min(etav,etay) + etamin) + gamma.*eta;                 % effective shear viscosity

etac = (eta(1:end-1,1:end-1)+eta(2:end,1:end-1) ...                        % evaluate in cell corners
     +  eta(1:end-1,2:end  )+eta(2:end,2:end  )).*0.25;

zetav   = etav  .*max(1e-6,f0.*f).^-m;                                     % cmpct viscosity
zetamin = etamin.*max(1e-6,f0.*1).^-m;                                     % min cmpct viscosity

zetay   = (yieldp+Ty)./max(1e-16,ups+Ty./zetav);                           % cmpct visco-plasticity

for k  = 1:ceil(kappa)                                                     % regularisation
    kk = kappa/ceil(kappa);
    zetay(2:end-1,2:end-1) = zetay(2:end-1,2:end-1) + kk.*(diff(zetay(2:end-1,:),2,2)+diff(zetay(:,2:end-1),2,1))./8;
    zetay([1 end],:) = zetay([end-1 2],:);
    zetay(:,[1 end]) = zetay(:,[end-1 2]);
end

zeta = (1-gamma).*(min(zetav,zetay) + zetamin); + gamma.*zeta;             % effective cmpct viscosity

K    = f.^n;                                                               % segregation coefficient
   
% update iterative and physical time step sizes
dtW = ((eta(1:end-1,:)+ eta(2:end,:)).*0.5./(h/2)^2 ...
    + (zeta(1:end-1,:)+zeta(2:end,:)).*0.5./(h/2)^2).^-1;                  % W iterative step size
dtU = ((eta(:,1:end-1)+ eta(:,2:end)).*0.5./(h/2)^2 ...
    + (zeta(:,1:end-1)+zeta(:,2:end)).*0.5./(h/2)^2).^-1;                  % U iterative step size
dtP = (1./eta + K./(h/2)^2).^-1;                                           % P iterative step size

Vel = [U(:)+UBG(:);W(:)+WBG(:);u(:);w(:)];                                 % combine all velocity components
dt  = rmp.*CFL*min([h/2/max(abs(Vel)), 0.05./max(abs(Div_fV(:)))]);        % physical time step

% clean workspace
clear Vel k kk etav etay zetav zetamin zetay
