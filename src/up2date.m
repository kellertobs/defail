% update tensor magnitudes
eII(2:end-1,2:end-1) = (  (exx(2:end-1,2:end-1).^2 + ezz(2:end-1,2:end-1).^2 ...
                     + 2.*(exz(1:end-1,1:end-1).^2 + exz(2:end  ,1:end-1).^2 ...
                         + exz(1:end-1,2:end  ).^2 + exz(2:end  ,2:end  ).^2).*0.25)./2).^0.5 + 1e-16;
eII([1 end],:) = eII([end-1 2],:);                                         % periodic boundaries
eII(:,[1 end]) = eII(:,[end-1 2]);

tII(2:end-1,2:end-1) = (  (txx(2:end-1,2:end-1).^2 + tzz(2:end-1,2:end-1).^2 ...
                     + 2.*(txz(1:end-1,1:end-1).^2 + txz(2:end  ,1:end-1).^2 ...
                         + txz(1:end-1,2:end  ).^2 + txz(2:end  ,2:end  ).^2).*0.25)./2).^0.5 + 1e-16;
tII([1 end],:) = tII([end-1 2],:);                                         % periodic boundaries
tII(:,[1 end]) = tII(:,[end-1 2]);

% yieldt = max(    1e-6,  p + Ty );                                          % get Griffith yield stress
% yieldp = max(-Ty+1e-6,-Ty + tII);                                          % get Griffith yield pressure
% ypcut  = max(-Ty+1e-6,-Ty + tII.*(1-delta));

yieldt = max(    1e-6,(p + tII + Ty)/2);                                     % get Griffith yield pressure
yieldp = max(-Ty+1e-6,(p + tII - Ty)/2);                                     % get Griffith yield pressure
ypcut  = max(-Ty+1e-6,(p + tII.*(1-delta) - Ty)/2);

% update rheological parameters
etav    = exp(-lmd.*f0.*(f-1)) ...
        .* (eII./eIIref).^((1-nn)./nn);                                      % get shear viscosity

t0 = etav.*eII;

etayi  =  etav.*(yieldt./max(1e-6,t0));
for k  = 1:ceil(kappa)                                                     % regularise visco-plasticity
    kk = delta/ceil(kappa);
    etayi(2:end-1,2:end-1) = etayi(2:end-1,2:end-1) + kk.*(diff(etayi(2:end-1,:),2,2)+diff(etayi(:,2:end-1),2,1))./8;
    etayi([1 end],:) = etayi([end-1 2],:);
    etayi(:,[1 end]) = etayi(:,[end-1 2]);
end
etay  =  etayi.*(1-gamma) +  etay.*gamma;                                  % iterative relaxation
eta   = (etav.^(-1/eps) + etay.^(-1/eps)).^(-eps) + etamin;
% eta   = min(etav,etay) + etamin;

zetav   = eta.*max(1e-6,(f0.*f)).^-m;                                      % get cmpct viscosity
zetamin = etamin.*max(1e-6,(f0.*f)).^-m;                                      % get cmpct viscosity

p0 = zetav.*Div_V;

zetayi  = zetav.*((yieldp + Ty)./max(1e-6,p0 + Ty));
for k  = 1:ceil(kappa)                                                     % regularise visco-plasticity
    kk = delta/ceil(kappa);
    zetayi(2:end-1,2:end-1) = zetayi(2:end-1,2:end-1) + kk.*(diff(zetayi(2:end-1,:),2,2)+diff(zetayi(:,2:end-1),2,1))./8;
    zetayi([1 end],:) = zetayi([end-1 2],:);
    zetayi(:,[1 end]) = zetayi(:,[end-1 2]);
end
zetay =  zetayi.*(1-gamma) +  zetay.*gamma;                                  % iterative relaxation

zeta  = (zetav.^(-1/eps) + zetay.^(-1/eps)).^(-eps) + zetamin;
% zeta = min(zetav,zetay) + zetamin;
% zeta = eta.*max(1e-6,(f0.*f)).^-m;

K    = f.^n;                                                               % get segregation coefficient

etac = (eta(1:end-1,1:end-1)+eta(2:end,1:end-1) ...                        % evaluate in cell corners
     +  eta(1:end-1,2:end  )+eta(2:end,2:end  )).*0.25;
   
% update iterative and physical time step sizes
dtW = ((eta(1:end-1,:)+ eta(2:end,:)).*0.5./(h/2)^2 ...
    + (zeta(1:end-1,:)+zeta(2:end,:)).*0.5./(h/2)^2).^-1;                  % iterative step size
dtU = ((eta(:,1:end-1)+ eta(:,2:end)).*0.5./(h/2)^2 ...
    + (zeta(:,1:end-1)+zeta(:,2:end)).*0.5./(h/2)^2).^-1;                  % iterative step size
dtP = (1./eta + K./(h/2)^2).^-1;                                           % iterative step size
Vel = [U(:)+UBG(:);W(:)+WBG(:);u(:);w(:)];                                 % combine all velocity components
if step > 1
    dt  = min(2*dto,CFL*min([h/2/max(abs(Vel)), ...
              0.05./max(abs(Div_fV(:)))]));                                % physical time step
end

clear etayi eIIvp Vel d dd