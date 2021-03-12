clear;
close all;

runID = 'Si2_B0_1_N250';
start = 5;
step  = 5;
stop  = 20;
sat   = 5;

rng(5);
Fc = randn(4,3).*2;

for frame = start:step:stop
    
    load(['../out/',runID,'/',runID,'_',int2str(frame)]);
    
    % process solution variables:
    %   remove boundaries, interpolate to centre nodes, get magnitude
    f  =  f(2:end-1,2:end-1);
    W  = -(W(1:end-1,2:end-1)+W(2:end,2:end-1))./2;
    U  =  (U(2:end-1,1:end-1)+U(2:end-1,2:end))./2;
    V  =  sqrt(W.^2 + U.^2);
    P  =  P(2:end-1,2:end-1);
    
    w  = -(w(1:end-1,2:end-1)+w(2:end,2:end-1))./2;
    u  =  (u(2:end-1,1:end-1)+u(2:end-1,2:end))./2;
    v  =  sqrt(w.^2 + u.^2);
    p  =  p(2:end-1,2:end-1);
    
    eps  =  eps(2:end-1,2:end-1);
    ups  =  ups(2:end-1,2:end-1);
    tau  =  tau(2:end-1,2:end-1);
    eta  =  eta(2:end-1,2:end-1);
   zeta  = zeta(2:end-1,2:end-1);

    vn   =  v  ./std(v  (:));  vn = vn-mean(vn(:));
    en   =  eps./std(eps(:));  en = en-mean(en(:));
    un   =  ups./std(ups(:));  un = un-mean(un(:));
    pn   =  p  ./std(p  (:));  pn = pn-mean(pn(:));
    tn   =  tau./std(tau(:));  tn = tn-mean(tn(:));
    
    X    =  [f(:),v(:),p(:),eps(:),ups(:),tau(:)];
    [PC_C, PC_A, PC_V] = pca(X,'Algorithm','svd','Centered','on','VariableWeights','variance');
    PC = PC_C.';
    PA = PC_A;
    PV = PC_V;
    CV = cumsum(PC_V(1:end-1)./sum(PC_V(1:end-1)));
    
    figure();
    subplot(2,3,1); plot(X(:,1),X(:,1),'k.'); hold on;
    subplot(2,3,2); plot(X(:,1),X(:,2),'k.'); hold on;
    subplot(2,3,3); plot(X(:,1),X(:,3),'k.'); hold on;
    subplot(2,3,4); plot(X(:,1),X(:,4),'k.'); hold on;
    subplot(2,3,5); plot(X(:,1),X(:,5),'k.'); hold on;
    subplot(2,3,6); plot(X(:,1),X(:,6),'k.'); hold on;

    Fp = PC(1:3,:);
    Ap = PA(:,1:3);
    Xp = Ap*Fp + mean(X);
    
    subplot(2,3,1); plot(Xp(:,1),Xp(:,1),'b.'); hold on;
    subplot(2,3,2); plot(Xp(:,1),Xp(:,2),'b.'); hold on;
    subplot(2,3,3); plot(Xp(:,1),Xp(:,3),'b.'); hold on;
    subplot(2,3,4); plot(Xp(:,1),Xp(:,4),'b.'); hold on;
    subplot(2,3,5); plot(Xp(:,1),Xp(:,5),'b.'); hold on;
    subplot(2,3,6); plot(Xp(:,1),Xp(:,6),'b.'); hold on;
    
    RGB = [Ap(:,1),Ap(:,2),Ap(:,3)];
    RGB = (RGB - (mean(RGB)-10/sat*std(RGB))) ./ ((mean(RGB)+10/sat*std(RGB)) - (mean(RGB)-10/sat*std(RGB)));
    
    figure();
    image(reshape(RGB,N-2,N-2,3)); axis ij equal tight; box on; drawnow;
    
    [Ic,Fc] = kmeans(Ap,4,'Start',Fc,'MaxIter',200);
    
    figure();
    scatter3(Ap(:,1),Ap(:,2),Ap(:,3),10,Ic,'filled'); axis equal tight; drawnow;
    
    figure();
    imagesc(reshape(Ic,N-2,N-2));  axis ij equal tight; box on; colormap(ocean); drawnow;
        
end