% print run header
fprintf('\n\n*****  defail  |  %s  |  %s  *****\n\n',runID,datetime);

% load custom colormap
load ocean;  

% produce smooth random perturbations
rng(5);
dr = randn(N,N);
for i = 1:max(smx,smz)
    dr(2:end-1,2:end-1) = dr(2:end-1,2:end-1) ...
                       + smz./max(smx,smz).*diff(dr(:,2:end-1),2,1)./8 ...
                       + smx./max(smx,smz).*diff(dr(2:end-1,:),2,2)./8;
    dr = dr - mean(dr(:));
    dr = dr./max(abs(dr(:)));
    dr([1 2 end-1 end],:) = 0;%dr([end-1 2],:);%dr([round(N/2)-1 round(N/2)],:);
    dr(:,[1 2 end-1 end]) = 0;%dr(:,[end-1 2]);%dr(:,[round(N/2)-1 round(N/2)]);
end

% get coordinate arrays
x     = -L/2-h/2:h:L/2+h/2;
z     = -L/2-h/2:h:L/2+h/2;
xc    = (x(1:end-1)+x(2:end))./2;
zc    = (z(1:end-1)+z(2:end))./2;
[X,Z] = meshgrid(x,z);

% initialise solution and residual fields
WP    =  (Z(1:end-1,:)+Z(2:end,:))/2/L*Pu*L;  % pure   shear WBG
WS    = -(X(1:end-1,:)+X(2:end,:))/2/L*Si*L;  % simple shear WBG
UP    = -(X(:,1:end-1)+X(:,2:end))/2/L*Pu*L;  % pure   shear UBG
US    = -(Z(:,1:end-1)+Z(:,2:end))/2/L*Si*L;  % simple shear UBG
WS    = 0.*WS; US = 2.*US; 
WBG   = WP+WS;
UBG   = UP+US;
clear  WP WS UP US

if bnchmrk
    mms;
    f = f_mms;
else
    f  =  1 + f1.*dr + f2.*exp(-(X+xpos).^2./wx^2).*exp(-(Z+zpos).^2./wz^2);
end
res_f = 0.*f;  df = 0.*f;

W      =  0.*WBG;  res_W = 0.*W;  dW = 0.*W;
U      =  0.*UBG;  res_U = 0.*U;  dU = 0.*U;
P      =  0.*f;    res_P = 0.*P;  dP = 0.*P;
u      =  0.*U;
w      =  0.*W;
p      =  0.*P;

% initialise parameter fields
ups    =  0.*P;  upss = 0.*P;  Div_fV = 0.*P;  %Div_fVBG = 0.*P;
eps0   =  abs(Pu) + abs(Si) + 1e-6;  
exx    =  0.*P - Pu;  ezz = 0.*P + Pu;  exz = zeros(N-1,N-1) - Si;  eps = 0.*P + (abs(Pu) + abs(Si));  
txx    =  0.*exx;  tzz = 0.*ezz;  txz = 0.*exz;  tau = 0.*eps;
eta    =  log10( exp(-lmd.*f0.*(f-1)) );
zeta   =  eta - max(-6,log10(f0.*f));
yieldt =  ones(size(P));

% initialise timing parameters
step   =  0;
time   =  0;
dt     =  0;

% print initial condition
if ~restart; up2date; output; end

% overwrite fields from file if restarting run
if     restart < 0  % restart from last continuation frame
    if isfile(['../out/',runID,'/',runID,'_cont.mat'])
        name = ['../out/',runID,'/',runID,'_par'];
        load(name);
        name = ['../out/',runID,'/',runID,'_cont.mat'];
        load(name);
    else
        restart = 0;
    end
elseif restart > 0  % restart from specified continuation frame
    name = ['../out/',runID,'/',runID,'_',num2str(restart),'.mat'];
    load(name);
    name = ['../out/',runID,'/',runID,'_par'];
    load(name);
end

% time stepping loop
while time < tend && step <= M
    
    fprintf(1,'\n\n*****  step %d;  dt = %4.4e;  time = %4.4e;\n\n',step,dt,time);
    tic;
    
    % store previous solution
    fo = f; Div_fVo = Div_fV;

    % reset residual norms and iteration count
    resnorm  = 1e3;
    resnorm0 = resnorm;
    it       = 0;
    
    % initialise iterative solution guess
    dWi = 0.*W;
    dUi = 0.*U;
    dPi = 0.*P;
    
    
    % non-linear iteration loop
    startup = 2*double(step<=0) + double(step>0);
    while resnorm/resnorm0 >= rtol/startup && resnorm >= atol && it <= maxit*startup || it <= minit
                
        % store previous iterative solution guess  
        Wi = W;
        Ui = U;
        Pi = P;
                
        
        % update fields
        if ~mod(it,nup); up2date; end

        
        % update liquid fraction
        if ~mod(it,nup) && step > 0
            flxdiv;                                                        % flux divergence for advection/compaction
            
            Vel = [U(:)+UBG(:);W(:)+WBG(:);u(:);w(:)];                     % combine all velocity components
            dt  = CFL*min([h/2/max(abs(Vel)), 0.05./max(abs(Div_fV(:)))]); % physical time step

            res_f = (f-fo)./dt - (theta.*Div_fV   + (1-theta).*Div_fVo  ); % residual liquid evolution equation
            
            df = - res_f.*dt/10;
            
            df([1 end],:) = df([end-1 2],:);                               % periodic boundaries
            df(:,[1 end]) = df(:,[end-1 2]);
            
            if demean; df = df - mean(df(:)); end                          % remove mean from update
            
            f = f + alpha*df;                                              % update solution
        end
        
        
        % update segregation velocities and compaction pressure
        w   = -(K(1:end-1,:)+K(2:end,:)).*0.5 .* (diff(P,1,1)./h + B);     % z-segregation velocity
        w([1 end],:) = [sum(w([1 end],:),1)./2; ...                        % periodic boundaries
                        sum(w([1 end],:),1)./2];
        w(:,[1 end]) = w(:,[end-1 2]);
        
        u   = - (K(:,1:end-1)+K(:,2:end)).*0.5 .* (diff(P,1,2)./h);        % x-segregation velocity
        u([1 end],:) = u([end-1 2],:);                                     % periodic boundaries
        u(:,[1 end]) = [sum(u(:,[1 end]),2)./2, ...
                        sum(u(:,[1 end]),2)./2];
                        
        p   = -10.^zeta .* ups;                                            % compaction pressure

        p([1 end],:) = p([end-1 2],:);                                     % periodic boundaries
        p(:,[1 end]) = p(:,[end-1 2]);
        
        
        % update strain rates
        exx(:,2:end-1)   = diff(U,1,2)./h - ups(:,2:end-1)./3 - Pu;        % x-normal strain rate
        exx([1 end],:)   = exx([end-1 2],:);                               % periodic boundaries
        exx(:,[1 end])   = exx(:,[end-1 2]);               
        ezz(2:end-1,:)   = diff(W,1,1)./h - ups(2:end-1,:)./3 + Pu;        % z-normal strain rate
        ezz([1 end],:)   = ezz([end-1 2],:);                               % periodic boundaries
        ezz(:,[1 end])   = ezz(:,[end-1 2]);          
        exz              = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h) - Si;      % shear strain rate
        
        % update stresses
        txx = 10.^eta .* exx;                                              % x-normal stress
        tzz = 10.^eta .* ezz;                                              % z-normal stress
        txz = 10.^etac.* exz;                                              % xz-shear stress  
        
        
        % update z-reference velocity
        Div_tz = diff(tzz(:,2:end-1),1,1)./h + diff(txz,1,2)./h;           % z-stress divergence
        
        res_W(:,2:end-1) = - Div_tz + diff(P(:,2:end-1),1,1)./h ...        % residual z-momentum equation
                                    + diff(p(:,2:end-1),1,1)./h;
                                
        if bnchmrk; res_W = res_W - src_W_mms; end
        
        dW = -res_W.*dtW;
        
        dW([1 end],:) = [sum(dW([1 end],:),1)./2; ...                      % periodic boundaries
                         sum(dW([1 end],:),1)./2];
        dW(:,[1 end]) = dW(:,[end-1 2]);
                
        if demean; dW = dW-mean(dW(:)); end                                % remove mean from update

        W = Wi + alpha*dW + beta*dWi;                                      % update z-velocity solution

        dWi = W - Wi;                                                      % store update step
        

        % update x-reference velocity        
        Div_tx  = diff(txx(2:end-1,:),1,2)./h + diff(txz,1,1)./h;          % x-stress divergence
        
        res_U(2:end-1,:) = - Div_tx + diff(P(2:end-1,:),1,2)./h ...        % residual x-momentum equation
                                    + diff(p(2:end-1,:),1,2)./h;
                                
        if bnchmrk; res_U = res_U - src_U_mms; end
        
        dU = -res_U.*dtU;
        
        dU([1 end],:) = dU([end-1 2],:);                                   % periodic boundaries
        dU(:,[1 end]) = [sum(dU(:,[1 end]),2)./2, ...
                         sum(dU(:,[1 end]),2)./2];

        if demean; dU = dU-mean(dU(:)); end % remove mean from update
        
        dU = dU-mean(dU(:));                                               % remove mean from update
        
        U = Ui + alpha*dU + beta*dUi;                                      % update x-velocity solution
      
        dUi = U - Ui;                                                      % store update step
        
                
        % update velocity divergences
        ups(2:end-1,2:end-1) = diff(U(2:end-1,:),1,2)./h ...               % velocity divergence
                             + diff(W(:,2:end-1),1,1)./h;
        ups([1 end],:) = ups([end-1 2],:);                                 % periodic boundaries
        ups(:,[1 end]) = ups(:,[end-1 2]);
        
        upss(2:end-1,2:end-1) = diff(u(2:end-1,:),1,2)./h ...
                              + diff(w(:,2:end-1),1,1)./h;                 % segregation velocity divergence
        upss([1 end],:) = upss([end-1 2],:);                               % periodic boundaries         
        upss(:,[1 end]) = upss(:,[end-1 2]);          

        
        % update reference pressure
        res_P = ups + upss;                                                % residual mass equation
        
        if bnchmrk; res_P = res_P - src_P_mms; end
        
        dP = -res_P.*dtP;
        
        dP([1 end],:) = dP([end-1 2],:);                                   % periodic boundaries
        dP(:,[1 end]) = dP(:,[end-1 2]);
        
        if demean; dP = dP-mean(dP(:)); end                                % remove mean from update
            
        P = Pi + alpha*dP + beta*dPi;                                      % update pressure solution

        dPi = P - Pi;                                                      % store update step

        
        % check and report convergence every max(100,nup) iterations
        if ~mod(it,max(100,nup)); report; end
        
        it = it+1;

    end
    
    % print diagnostics
    fprintf(1,'\n         time to solution = %4.4f sec\n\n',toc);

    fprintf(1,'         min U    = %s%4.4f;   mean U    = %s%4.4f;   max U    = %s%4.4f;\n'  ,int8(min(  U(:))<0),min(  U(:)),int8(mean(  U(:))<0),mean(  U(:)),int8(max(  U(:))<0),max(  U(:)));
    fprintf(1,'         min W    = %s%4.4f;   mean W    = %s%4.4f;   max W    = %s%4.4f;\n'  ,int8(min( -W(:))<0),min( -W(:)),int8(mean( -W(:))<0),mean( -W(:)),int8(max( -W(:))<0),max( -W(:)));
    fprintf(1,'         min P    = %s%4.4f;   mean P    = %s%4.4f;   max P    = %s%4.4f;\n'  ,int8(min(  P(:))<0),min(  P(:)),int8(mean(  P(:))<0),mean(  P(:)),int8(max(  P(:))<0),max(  P(:)));
    fprintf(1,'         min f    = %s%4.4f;   mean f    = %s%4.4f;   max f    = %s%4.4f;\n\n',int8(min(  f(:))<0),min(  f(:)),int8(mean(  f(:))<0),mean(  f(:)),int8(max(  f(:))<0),max(  f(:)));

    fprintf(1,'         min u    = %s%4.4f;   mean u    = %s%4.4f;   max u    = %s%4.4f;\n'  ,int8(min(  u(:))<0),min(  u(:)),int8(mean(  u(:))<0),mean(  u(:)),int8(max(  f(:))<0),max(  u(:)));
    fprintf(1,'         min w    = %s%4.4f;   mean w    = %s%4.4f;   max w    = %s%4.4f;\n'  ,int8(min( -w(:))<0),min( -w(:)),int8(mean( -w(:))<0),mean( -w(:)),int8(max(  f(:))<0),max( -w(:)));
    fprintf(1,'         min p    = %s%4.4f;   mean p    = %s%4.4f;   max p    = %s%4.4f;\n\n',int8(min(  p(:))<0),min(  p(:)),int8(mean(  p(:))<0),mean(  p(:)),int8(max(  f(:))<0),max(  p(:)));
    
    
    fprintf(1,'         min  eta = %s%4.4f;   mean  eta = %s%4.4f;   max  eta = %s%4.4f;\n'  ,int8(min( eta(:))<0),min( eta(:)),int8(mean( eta(:))<0),mean( eta(:)),int8(max( eta(:))<0),max( eta(:)));
    fprintf(1,'         min zeta = %s%4.4f;   mean zeta = %s%4.4f;   max zeta = %s%4.4f;\n'  ,int8(min(zeta(:))<0),min(zeta(:)),int8(mean(zeta(:))<0),mean(zeta(:)),int8(max(zeta(:))<0),max(zeta(:)));
    fprintf(1,'         min K    = %s%4.4f;   mean K    = %s%4.4f;   max K    = %s%4.4f;\n\n',int8(min(  K(:))<0),min(  K(:)),int8(mean(  K(:))<0),mean(  K(:)),int8(max(  K(:))<0),max(  K(:)));

    fprintf(1,'         min ups  = %s%4.4f;   mean ups  = %s%4.4f;   max ups  = %s%4.4f;\n'  ,int8(min(ups(:))<0),min(ups(:)),int8(mean(ups(:))<0),mean(ups(:)),int8(max(ups(:))<0),max(ups(:)));
    fprintf(1,'         min eps  = %s%4.4f;   mean eps  = %s%4.4f;   max eps  = %s%4.4f;\n'  ,int8(min(eps(:))<0),min(eps(:)),int8(mean(eps(:))<0),mean(eps(:)),int8(max(eps(:))<0),max(eps(:)));
    fprintf(1,'         min tau  = %s%4.4f;   mean tau  = %s%4.4f;   max tau  = %s%4.4f;\n\n',int8(min(tau(:))<0),min(tau(:)),int8(mean(tau(:))<0),mean(tau(:)),int8(max(tau(:))<0),max(tau(:)));

    
    % update parameters and plot results
    if ~mod(step,nop); up2date; output; end
    
    
    % increment time step
    time = time+dt;
    step = step+1;
    
end

% plot results
up2date; output; 

diary off
