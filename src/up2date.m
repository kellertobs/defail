% update tensor magnitudes
eII(2:end-1,2:end-1) = (  (exx(2:end-1,2:end-1).^2 + ezz(2:end-1,2:end-1).^2 ...
                     + 2.*(exz(1:end-1,1:end-1).^2.* exz(2:end  ,1:end-1).^2 ...
                        .* exz(1:end-1,2:end  ).^2.* exz(2:end  ,2:end  ).^2).^0.25)./2).^0.5 + 1e-16;
eII([1 end],:) = eII([end-1 2],:);                                       % periodic boundaries
eII(:,[1 end]) = eII(:,[end-1 2]);

eIIvp = max(1e-6,eII  -(tII-tIIo)./(De*dt + etamin));                      % get visco-plastic strain rate
DvVvp = max(1e-6,Div_V-(p  -po  )./(De*dt + etamin));                      % get visco-plastic strain rate

tII(2:end-1,2:end-1) = (  (txx(2:end-1,2:end-1).^2 + tzz(2:end-1,2:end-1).^2 ...
                     + 2.*(txz(1:end-1,1:end-1).^2.* txz(2:end  ,1:end-1).^2 ...
                        .* txz(1:end-1,2:end  ).^2.* txz(2:end  ,2:end  ).^2).^0.25)./2).^0.5 + 1e-16;
tII([1 end],:) = tII([end-1 2],:);                                         % periodic boundaries
tII(:,[1 end]) = tII(:,[end-1 2]);

% update rheological parameters
etav  = exp(lmd.*f0.*(f-1));                                               % get shear viscosity

yield = max(5e-3,Ty + p);                                                  % get Griffith yield stress

etayi = etay;                                                             
etay  = min(etav,yield./eIIvp);                                            % get shear visco-plasticity

etay  = log(etay);
for d = 1:ceil(delta)                                                      % regularise visco-plasticity
    dd   = delta/ceil(delta);
    etay(2:end-1,2:end-1) = etay(2:end-1,2:end-1) + dd.*(diff(etay(2:end-1,:),2,2)+diff(etay(:,2:end-1),2,1))./8;
    etay([1 end],:) = etay([end-1 2],:);
    etay(:,[1 end]) = etay(:,[end-1 2]);
end
etay  =  exp(etay);

etay  = etay.*gamma + etayi.*(1-gamma);                                    % iterative relaxation
etayc = (etay(1:end-1,1:end-1)+etay(2:end,1:end-1) ...                     % evaluate in cell corners
       + etay(1:end-1,2:end  )+etay(2:end,2:end  ))./4;

eta   = (1./(etay  + etamin) + 1./(De*dt + etamin)).^-1;                   % get shear visco-elasto-plasticity
etac  = (1./(etayc + etamin) + 1./(De*dt + etamin)).^-1;
chi   = (1 + (De*dt + etamin)./(etay  + etamin)).^-1;                      % get shear visco-elastic evolution parameter
chic  = (1 + (De*dt + etamin)./(etayc + etamin)).^-1;

zeta  = eta./max(1e-6,(f0.*f).^m);                                         % get compaction visco-elasto-plasticity
xi    = chi;                                                               % get compaction visco-elastic evolution parameter

K     = f.^n;                                                              % get segregation coefficient

% update iterative and physical time step sizes
dtW = ((eta(1:end-1,:) +  eta(2:end,:))./2./(h/2)^2 ...
    + (zeta(1:end-1,:) + zeta(2:end,:))./2./(h/2)^2).^-1;                  % iterative step size
dtU = ((eta(:,1:end-1) +  eta(:,2:end))./2./(h/2)^2 ...
    + (zeta(:,1:end-1) + zeta(:,2:end))./2./(h/2)^2).^-1;                  % iterative step size
dtP = (1./eta + K./(h/2)^2).^-1;                                           % iterative step size
Vel = [U(:)+UBG(:);W(:)+WBG(:);u(:);w(:)];                                 % combine all velocity components
dt  = min(2*dto,CFL*min([h/max(abs(Vel))/2,0.01./max(abs(Div_fV(:)))]));   % physical time step

clear etayi eIIvp etav Vel d dd