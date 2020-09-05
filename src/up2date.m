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

yieldt = max(    1e-6,  p + Ty );                                          % get Griffith yield stress
yieldp = max(-Ty+1e-6,-Ty + tII);                                          % get Griffith yield pressure

% update rheological parameters
etav   = exp(-lmd.*f0.*(f-1)) ...
       .* ((eII+abs(Div_V)./3)./eIIref).^((1-nn)./nn);                     % get shear viscosity

etayi  =  etav.*(yieldt./(etav.*eII)).^(1+delta);
for k  = 1:ceil(kappa)                                                     % regularise visco-plasticity
    kk = delta/ceil(kappa);
    etayi(2:end-1,2:end-1) = etayi(2:end-1,2:end-1) + kk.*(diff(etayi(2:end-1,:),2,2)+diff(etayi(:,2:end-1),2,1))./8;
    etayi([1 end],:) = etayi([end-1 2],:);
    etayi(:,[1 end]) = etayi(:,[end-1 2]);
end
etay  =  etayi.*(1-gamma) +  etay.*gamma;                                  % iterative relaxation

 eta = (etav.^(-1/eps) + etay.^(-1/eps)).^(-eps) + etamin;
zeta = eta.*max(1e-6,(f0.*f)).^-m;

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