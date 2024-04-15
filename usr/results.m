clear;

runID = 'demo';
start = 0;
step  = 1;
stop  = 0;
sat   = 5;
smth  = 0.1;
c     = 3;
cc    = 3;

rng(24);
Fc = randn(cc,c).*2;

% load refspct;

for frame = start:step:stop
    
    load(['../out/',runID,'/',runID,'_',int2str(frame)]);
    Xax = x(2:end-1).*sqrt(f0);
    
    % process solution variables:
    % remove boundaries, interpolate to centre nodes, get magnitude
    f  =  f(2:end-1,2:end-1);
    W  = -(W(1:end-1,2:end-1)+W(2:end,2:end-1))./2;
    U  =  (U(2:end-1,1:end-1)+U(2:end-1,2:end))./2;
    V  =  sqrt(W.^2 + U.^2);
    P  =  P(2:end-1,2:end-1);
    
    w  = -(w(1:end-1,2:end-1)+w(2:end,2:end-1))./2;
    u  =  (u(2:end-1,1:end-1)+u(2:end-1,2:end))./2;
    v  =  sqrt(w.^2 + u.^2);
    p  =  p(2:end-1,2:end-1);
    
    eps  = imnlmfilt( eps(2:end-1,2:end-1),'DegreeOfSmoothing',smth);
    ups  = imnlmfilt( ups(2:end-1,2:end-1),'DegreeOfSmoothing',smth);
    tau  = imnlmfilt( tau(2:end-1,2:end-1),'DegreeOfSmoothing',smth);
    eta  = imnlmfilt( eta(2:end-1,2:end-1),'DegreeOfSmoothing',smth);
   zeta  = imnlmfilt(zeta(2:end-1,2:end-1),'DegreeOfSmoothing',smth);
    K    = imnlmfilt(   K(2:end-1,2:end-1),'DegreeOfSmoothing',smth);
    
    X    =  [f(:),v(:),p(:),tau(:),eps(:),ups(:)];
    k    =  size(X,2);
    
    [PC_C, PC_A, PC_V] = pca(X,'Algorithm','svd','Centered','on','VariableWeights','variance');
    PC = PC_C.';
    PA = PC_A;
    PV = PC_V;
    CV = cumsum(PC_V(1:end-1)./sum(PC_V(1:end-1)));
    
    Fp = PC(1:c,:);
    Ap = PA(:,1:c);
    Xp = Ap*Fp + mean(X);
    
    figure(1); clf;
    ik = 1;
    for kk=1:k
        for kkk = 1:k
            subplot(k,k,ik);
            plot(X(1:10:end,kkk)./1e3,X(1:10:end,kk)./1e3,'k.','MarkerSize',10); axis ij tight; hold on;
            ik = ik+1;
        end
    end
    drawnow;
    
    ik = 1;
    for kk=1:k
        for kkk = 1:k
            subplot(k,k,ik);
            plot(Xp(1:10:end,kkk)./1e3,Xp(1:10:end,kk)./1e3,'b.','MarkerSize',10); axis ij tight; hold off;
            ik = ik+1;
        end
    end
    drawnow;
    
    RGB = [Ap(:,1),Ap(:,2),Ap(:,3)];
    RGB = (RGB - (mean(RGB)-10/sat*std(RGB))) ./ ((mean(RGB)+10/sat*std(RGB)) - (mean(RGB)-10/sat*std(RGB)));
    
    figure(2); clf;
    image(Xax,Xax,reshape(RGB,N-2,N-2,3)); axis ij equal tight; box on; drawnow;
    
    [Ic,Fc] = kmeans(Ap,cc,'Start',Fc,'MaxIter',200);
    
    figure(3); clf;
    scatter3(Ap(:,1),Ap(:,2),Ap(:,3),10,Ic,'filled'); axis equal tight; colormap(copper); colorbar; drawnow;
    
    figure(4); clf;
    imagesc(Xax,Xax,reshape(Ic,N-2,N-2));  axis ij equal tight; box on; colormap(copper); colorbar; drawnow;
    
    maxsum = zeros(180,cc);
    for th = 1:1:181
        for ic = 1:cc
            maxsum(th,ic) = max(sum(imrotate(reshape(Ic==ic,N-2,N-2),th-1),1));        
        end

    end
    figure(5); clf;
    for ic = 1:cc
        subplot(cc,1,ic);
        plot(    [0,0],[-4,4],'k-'); hold on;
        plot(   -[45,45],[-4,4],'k--')
        plot(    [45,45],[-4,4],'k--')
        plot(   -[20,20],[-4,4],'k:')
        plot(    [20,20],[-4,4],'k:')
        plot(-90+[20,20],[-4,4],'k:')
        plot(+90-[20,20],[-4,4],'k:')
        plot(-90:90,(maxsum(:,ic)-mean(maxsum(:,ic)))./std(maxsum(:,ic)),'k'); axis tight; box on;
    end
    drawnow;

    Fs  = 1/(h*sqrt(f0));             % grid spacing in compaction lengths

    F  = reshape(sum(Ap.^2./c,2).^0.5,N-2,N-2);  % Melt fraction reshaped to 2D array
    
    figure(6); clf; 
    imagesc(Xax,Xax,F); axis equal tight; colorbar; colormap(copper); drawnow;

    n = 2^nextpow2(N-2);             % next higher power of two
    S = fft(F,n,2);                  % raw FFT spectrum
    
    S2 = abs(S/(N-2));               % double-sided spectrum
    S1 = 2*S2(:,1:n/2+1);            % single-sided spectrum
    
    figure(7); clf; 
    subplot(2,1,1);
    plot(1./((2*Fs/n):(Fs/n):(Fs/2-Fs/n)),sum(S1(:,3:n/2))); axis tight;
    title('horizontal FFT-spectrum of |PC|','Interpreter','latex','FontSize',18);
    xlabel('Period [$\delta_0$]','Interpreter','latex','FontSize',16);
    ylabel('Spectral Energy','Interpreter','latex','FontSize',16);
    
    S = fft(rot90(F),n,2);           % raw FFT spectrum
    
    S2 = abs(S/(N-2));               % double-sided spectrum
    S1 = 2*S2(:,1:n/2+1);            % single-sided spectrum
    
    figure(7);
    subplot(2,1,2)
    plot(1./((2*Fs/n):(Fs/n):(Fs/2-Fs/n)),sum(S1(:,3:n/2))); axis tight;
    title('vertical FFT-spectrum of |PC|','Interpreter','latex','FontSize',18);
    xlabel('Period [$\delta_0$]','Interpreter','latex','FontSize',16);
    ylabel('Spectral Energy','Interpreter','latex','FontSize',16);
    drawnow;
end