% print start of run
fprintf('\n\n*****  defail  |  %s  |  %s  *****\n\n',runID,datetime);

% get smoothed random noise
rng(15);
a = randn(N,N);
for i = 1:max(smx,smz)
    a(2:end-1,2:end-1) = a(2:end-1,2:end-1) ...
                       + smz./max(smx,smz).*diff(a(:,2:end-1),2,1)./8 ...
                       + smx./max(smx,smz).*diff(a(2:end-1,:),2,2)./8;
    a([1 end],:) = (a(round(N/8)+[0 1],:)+a(round(N*7/8)+[0,1],:))/2;
    a(:,[1 end]) = (a(:,round(N/8)+[0 1])+a(:,round(N*7/8)+[0,1]))/2;
end
a([1 end],:) = a([end-1 2],:);
a(:,[1 end]) = a(:,[end-1 2]);
a = a./max(abs(a(:)));

% get coordinate arrays
x     = -L/2-h/2:h:L/2+h/2;
z     = -L/2-h/2:h:L/2+h/2;
[X,Z] = meshgrid(x,z);
xc    = (x(1:end-1)+x(2:end))./2;
zc    = (z(1:end-1)+z(2:end))./2;

% initialise solution and residual fields
csw  =  0;
if bnchmrk
    solwave;
else
    f  =  1 + f2.*exp(-(X+xpos).^2./wx^2).*exp(-(Z+zpos).^2./wz^2);
end
f      =  f + a.*f1;  fo = f;  fi = f;  res_f = 0.*f;  fmass0 = f0.*sum(f(:));
UP     = -(X(:,1:end-1)+X(:,2:end))/2/L*Pu*L;
US     = -(Z(:,1:end-1)+Z(:,2:end))/2/L*Si*L;  UBG = UP+US;
WP     =  (Z(1:end-1,:)+Z(2:end,:))/2/L*Pu*L;
WS     = -(X(1:end-1,:)+X(2:end,:))/2/L*Si*L;  WBG = WP+WS+csw;
clear WP WS UP US
U      =  0.*UBG;  Ui = U;  res_U = 0.*U;
W      =  0.*WBG;  Wi = W;  res_W = 0.*W;
P      =  0.*f;  Pi = P;  res_P = 0.*P;  
u      =  0.*U;
w      =  0.*W;
p      =  0.*P;  po = p;

% initialise parameter fields
Div_V  =  0.*P;  Div_v = 0.*P;  Div_fV = 0.*P;  Div_fVBG = 0.*P;
eIIref =  abs(Pu) + abs(Si) + 1e-6;  
exx    =  0.*P - Pu;  ezz = 0.*P + Pu;  exz = zeros(N-1,N-1) - Si;  eII = 0.*P + abs(Pu) + abs(Si);  
txx    =  0.*exx;  tzz = 0.*ezz;  txz = 0.*exz;  tII = 0.*eII;  tIIo = tII;
etay   =  ones(size(P));
zetay  =  ones(size(P))./(f0.*f).^m;
yieldp =  zeros(size(P));
yieldt =  Ty.*ones(size(P));
step   =  0;
time   =  0;
dto    =  dt;
it     =  0;  

% overwrite fields from file if restarting run
if     restart < 0  % restart from last continuation frame
    name = ['../out/',runID,'/',runID,'_par'];
    load(name);
    name = ['../out/',runID,'/',runID,'_cont'];
    load([name,'.mat']);
elseif restart > 0  % restart from specified continuation frame
    name = ['../out/',runID,'/',runID,'_',num2str(restart)];
    load([name,'.mat']);
    name = ['../out/',runID,'/',runID,'_par'];
    load(name);
end

load ocean;  % load custom colormap

while time < tend && step < M
    
    % plot results
    if ~mod(step,nop); up2date; output; end
    
    % increment time/step
    time = time+dt;
    step = step+1;
    
    fprintf(1,'\n\n*****  step %d;  dt = %4.4e;  time = %4.4e;\n\n',step,dt,time);
    tic;
    
    % store previous solution
    fo      = f;
    Div_fVo = Div_fV;  Div_fVBGo = Div_fVBG;
    txxo    = txx; tzzo = tzz; txzo = txz;  
    tIIo    = tII;
    po      = p;
    dto     = dt;
    
    % reset residuals and iteration count
    resnorm  = 1e3;
    resnorm0 = resnorm;
    it       = 0;
    
    Wi = W;  Ui = U;  Pi = P;  fi = f;
    
    % non-linear iteration loop
    while resnorm/resnorm0 >= rtol && resnorm >= atol && it <= maxit || it <= nup
                
        % store previous outer solution and residual fields   
        Wii = Wi;  Wi = W;
        Uii = Ui;  Ui = U;
        Pii = Pi;  Pi = P;
        fii = fi;  fi = f;
        
        % update fields
        if ~mod(it,nup); up2date; end

        % update liquid fraction
        if ~mod(it,nup)
            flxdiv;  flxdivBG;                                             % flux divergence for advection/compaction
            
            res_f = (f-fo)./dt - (theta.*Div_fV   + (1-theta).*Div_fVo  ) ...
                               - (theta.*Div_fVBG + (1-theta).*Div_fVBGo); % residual liquid evolution equation
            
            res_f([1 end],:) = res_f([end-1 2],:);                         % periodic boundaries
            res_f(:,[1 end]) = res_f(:,[end-1 2]);
            
            if demean; res_f = res_f - mean(res_f(:)); end

            f = fi - alpha.*res_f.*dt/10 + beta.*(fi-fii);                 % update solution
        end
        
        % update segregation velocities and compaction pressure
        w   = - (K(1:end-1,:)+K(2:end,:)).*0.5 .* (diff(P,1,1)./h + B);    % z-segregation velocity
        w([1 end],:) = [sum(w([1 end],:),1)./2; ...                        % periodic boundaries
                        sum(w([1 end],:),1)./2];
        w(:,[1 end]) = w(:,[end-1 2]);
        
        u   = - (K(:,1:end-1)+K(:,2:end)).*0.5 .* (diff(P,1,2)./h);        % x-segregation velocity
        u([1 end],:) = u([end-1 2],:);                                     % periodic boundaries
        u(:,[1 end]) = [sum(u(:,[1 end]),2)./2, ...
                        sum(u(:,[1 end]),2)./2];
                        
        p   = max(ypcut,-zeta .* Div_V);                                  % compaction pressure
%         p   = -zeta .* Div_V;                                  % compaction pressure

        p([1 end],:) = p([end-1 2],:);                                     % periodic boundaries
        p(:,[1 end]) = p(:,[end-1 2]);
        
        % update z-reference velocity
        Div_tz = diff(tzz(:,2:end-1),1,1)./h + diff(txz,1,2)./h;           % get z-stress divergence
        
        res_W(:,2:end-1) = - Div_tz + diff(P(:,2:end-1),1,1)./h ...        % residual z-momentum equation
                                    + diff(p(:,2:end-1),1,1)./h;  

        res_W([1 end],:) = [sum(res_W([1 end],:),1)./2; ...                % periodic boundaries
                            sum(res_W([1 end],:),1)./2];
        res_W(:,[1 end]) = res_W(:,[end-1 2]);
        
        if demean; res_W = res_W - mean(res_W(:)); end
        
        W = Wi - alpha.*res_W.*dtW + beta.*(Wi-Wii);                       % update solution

        % update x-reference velocity        
        Div_tx  = diff(txx(2:end-1,:),1,2)./h + diff(txz,1,1)./h;          % get x-stress divergence
        
        res_U(2:end-1,:) = - Div_tx + diff(P(2:end-1,:),1,2)./h ...
                                    + diff(p(2:end-1,:),1,2)./h;           % residual x-momentum equation
        
        res_U([1 end],:) = res_U([end-1 2],:);                             % periodic boundaries
        res_U(:,[1 end]) = [sum(res_U(:,[1 end]),2)./2, ...
                            sum(res_U(:,[1 end]),2)./2];

        if demean; res_U = res_U - mean(res_U(:)); end

        U = Ui - alpha.*res_U.*dtU + beta.*(Ui-Uii);                       % update solution
        
        % update velocity divergences
        Div_V(2:end-1,2:end-1) = diff(U(2:end-1,:),1,2)./h ...             % get velocity divergence
                               + diff(W(:,2:end-1),1,1)./h;
        Div_V([1 end],:) = Div_V([end-1 2],:);                             % periodic boundaries
        Div_V(:,[1 end]) = Div_V(:,[end-1 2]);
        
        Div_v(2:end-1,2:end-1) = diff(u(2:end-1,:),1,2)./h ...
                               + diff(w(:,2:end-1),1,1)./h;                % segregation velocity divergence
        Div_v([1 end],:) = Div_v([end-1 2],:);                             % periodic boundaries         
        Div_v(:,[1 end]) = Div_v(:,[end-1 2]);          
        
        % update strain rates
        exx(:,2:end-1)   = diff(U,1,2)./h - Div_V(:,2:end-1)./3 - Pu;      % get x-normal strain rate
        exx([1 end],:)   = exx([end-1 2],:);                               % periodic boundaries
        exx(:,[1 end])   = exx(:,[end-1 2]);               
        ezz(2:end-1,:)   = diff(W,1,1)./h - Div_V(2:end-1,:)./3 + Pu;      % get z-normal strain rate
        ezz([1 end],:)   = ezz([end-1 2],:);                               % periodic boundaries
        ezz(:,[1 end])   = ezz(:,[end-1 2]);          
        exz              = 1/2.*(diff(U,1,1)./h + diff(W,1,2)./h) - Si;    % get shear strain rate
        
        % update stresses
        txx = eta .* exx;                                                  % x-normal stress
        tzz = eta .* ezz;                                                  % z-normal stress
        txz = etac.* exz;                                                  % xz-shear stress  
        
        % update reference pressure
        res_P = Div_V + Div_v;                                             % residual mass equation

        res_P([1 end],:) = res_P([end-1 2],:);                             % periodic boundaries
        res_P(:,[1 end]) = res_P(:,[end-1 2]);
        
        if demean; res_P = res_P - mean(res_P(:)); end

        P = Pi - alpha.*res_P.*dtP + beta.*(Pi-Pii);                       % update solution
        
        % check and report convergence every nup iterations
        if ~mod(it,nup); report; end
        
        it = it+1;

    end
    
    clear Wi Wii Ui Uii Pi Pii wi wii ui uii pi pii fi fii fo Div_Vo txxo tzzo txzo Div_tz Div_tx 
    
    % print diagnostics
    fprintf(1,'\n         time to solution = %4.4f sec\n\n',toc);

    fprintf(1,'         min U   = %s%4.4f;   mean U   = %s%4.4f;   max U   = %s%4.4f;\n'  ,int8(min(  U(:))<0),min(  U(:)),int8(mean(  U(:))<0),mean(  U(:)),int8(max(  U(:))<0),max(  U(:)));
    fprintf(1,'         min W   = %s%4.4f;   mean W   = %s%4.4f;   max W   = %s%4.4f;\n'  ,int8(min( -W(:))<0),min( -W(:)),int8(mean( -W(:))<0),mean( -W(:)),int8(max( -W(:))<0),max( -W(:)));
    fprintf(1,'         min P   = %s%4.4f;   mean P   = %s%4.4f;   max P   = %s%4.4f;\n\n',int8(min(  P(:))<0),min(  P(:)),int8(mean(  P(:))<0),mean(  P(:)),int8(max(  P(:))<0),max(  P(:)));

    fprintf(1,'         min u   = %s%4.4f;   mean u   = %s%4.4f;   max u   = %s%4.4f;\n'  ,int8(min(  u(:))<0),min(  u(:)),int8(mean(  u(:))<0),mean(  u(:)),int8(max(  f(:))<0),max(  u(:)));
    fprintf(1,'         min w   = %s%4.4f;   mean w   = %s%4.4f;   max w   = %s%4.4f;\n'  ,int8(min( -w(:))<0),min( -w(:)),int8(mean( -w(:))<0),mean( -w(:)),int8(max(  f(:))<0),max( -w(:)));
    fprintf(1,'         min p   = %s%4.4f;   mean p   = %s%4.4f;   max p   = %s%4.4f;\n\n',int8(min(  p(:))<0),min(  p(:)),int8(mean(  p(:))<0),mean(  p(:)),int8(max(  f(:))<0),max(  p(:)));
    
    fprintf(1,'         min f   = %s%4.4f;   mean f   = %s%4.4f;   max f   = %s%4.4f;\n'  ,int8(min(  f(:))<0),min(  f(:)),int8(mean(  f(:))<0),mean(  f(:)),int8(max(  f(:))<0),max(  f(:)));
    fprintf(1,'         min K   = %s%4.4f;   mean K   = %s%4.4f;   max K   = %s%4.4f;\n\n',int8(min(  K(:))<0),min(  K(:)),int8(mean(  K(:))<0),mean(  K(:)),int8(max(  K(:))<0),max(  K(:)));
    
    fprintf(1,'         min eta = %s%4.4f;   mean eta = %s%4.4f;   max eta = %s%4.4f;\n'  ,int8(min(eta(:))<0),min(eta(:)),int8(mean(eta(:))<0),mean(eta(:)),int8(max(eta(:))<0),max(eta(:)));
    fprintf(1,'         min eII = %s%4.4f;   mean eII = %s%4.4f;   max eII = %s%4.4f;\n'  ,int8(min(eII(:))<0),min(eII(:)),int8(mean(eII(:))<0),mean(eII(:)),int8(max(eII(:))<0),max(eII(:)));
    fprintf(1,'         min tII = %s%4.4f;   mean tII = %s%4.4f;   max tII = %s%4.4f;\n\n',int8(min(tII(:))<0),min(tII(:)),int8(mean(tII(:))<0),mean(tII(:)),int8(max(tII(:))<0),max(tII(:)));

end

% plot results
output; up2date; 

diary off
