clear;

runID = 'demo';        % run identifier
start = 10;             % start frame
step  = 1;             % step frame
stop  = 10;             % stop frame
smth  = 0.2;           % smoothing level for filtering
c     = 5;             % number of principal components retained
cc    = 4;             % number of clusters
sat   = 5;             % colour saturation for false colour plots

figno = 1;

for frame = start:step:stop
    
    % load output frame
    load(['../out/',runID,'/',runID,'_',int2str(frame)]);
    Xax = x(2:end-1).*sqrt(f0);
    
    % process output
    % remove boundaries, interpolate to centre nodes, get magnitude
    f  =  f(2:end-1,2:end-1);
    W  =  abs(W(1:end-1,2:end-1)+W(2:end,2:end-1))./2;
    U  =  abs(U(2:end-1,1:end-1)+U(2:end-1,2:end))./2;
    V  =  sqrt(W.^2 + U.^2);
    P  =  P(2:end-1,2:end-1);
    
    w  =  abs(w(1:end-1,2:end-1)+w(2:end,2:end-1))./2;
    u  =  abs(u(2:end-1,1:end-1)+u(2:end-1,2:end))./2;
    v  =  sqrt(w.^2 + u.^2);
    p  =  p(2:end-1,2:end-1);
    
    % filter output by non-local means filter
    eps  = imnlmfilt( eps(2:end-1,2:end-1),'DegreeOfSmoothing',smth);
    ups  = imnlmfilt( ups(2:end-1,2:end-1),'DegreeOfSmoothing',smth);
    tau  = imnlmfilt( tau(2:end-1,2:end-1),'DegreeOfSmoothing',smth);
    eta  = imnlmfilt( eta(2:end-1,2:end-1),'DegreeOfSmoothing',smth);
   zeta  = imnlmfilt(zeta(2:end-1,2:end-1),'DegreeOfSmoothing',smth);
    K    = imnlmfilt(   K(2:end-1,2:end-1),'DegreeOfSmoothing',smth);
    
    % assemble output vector
    X    =  [f(:),V(:),P(:),v(:),p(:),tau(:),eps(:),ups(:),eta(:),zeta(:)];
    k    =  size(X,2);
    
    % principal component analysis
    [PC_C, PC_A, PC_V] = pca(X,'Algorithm','svd','Centered','on','VariableWeights','variance');
    PC = PC_C.';
    PA = PC_A;
    PV = PC_V;
    CV = cumsum(PC_V./sum(PC_V));
    
    figure(1); clf;
    plot(1:k,CV,'k-o','MarkerSize',8,'LineWidth',2); axis tight; box on; hold on;
    plot(c,CV(c),'ro','MarkerSize',8,'LineWidth',2);

    % reduce data dimensionality by discarding higher principal components
    Fp = PC(1:c,:);
    Ap = PA(:,1:c);
    Xp = Ap*Fp + mean(X);
    
    % plot correlation matrix of unreduced vs reduced data
    figure(figno); clf; figno=figno+1;
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
    
    % plot false colour images of PCs
    figure(figno); clf; figno=figno+1;
    subplot(1,2,1)
    RGB = [Ap(:,1),Ap(:,2),Ap(:,3)];
    RGB = (RGB - (mean(RGB)-10/sat*std(RGB))) ./ ((mean(RGB)+10/sat*std(RGB)) - (mean(RGB)-10/sat*std(RGB)));
    image(Xax,Xax,reshape(RGB,N-2,N-2,3)); axis ij equal tight; box on; drawnow;
    subplot(1,2,2)
    RGB = [Ap(:,end-2),Ap(:,end-1),Ap(:,end)];
    RGB = (RGB - (mean(RGB)-10/sat*std(RGB))) ./ ((mean(RGB)+10/sat*std(RGB)) - (mean(RGB)-10/sat*std(RGB)));
    image(Xax,Xax,reshape(RGB,N-2,N-2,3)); axis ij equal tight; box on; drawnow;    

    % clustering analysis
    rng(24);                  % random number seed
    Fc = randn(cc,c,10).*2;   % random starting guess
    [Ic,Fc] = kmeans(Ap,cc,'Start',Fc,'MaxIter',200,'Replicates',10);
    
    % extract mean properties of clusters
    Xc = Fc*Fp + mean(X);
    [~,isort] = sort(Xc(:,1),'ascend');
    Xc = Xc(isort,:);
    Ics = 0*Ic;
    for ic=1:cc
        Ics(Ic==ic) = isort(ic);
    end    

    % plot clusters in data space of first three PCs
    figure(figno); clf; figno=figno+1;
    scatter3(Ap(:,1),Ap(:,2),Ap(:,3),10,Ic,'filled'); axis equal tight; colormap(copper); colorbar; drawnow;
    
    % plot clusters as colour map
    figure(figno); clf; figno=figno+1;
    imagesc(Xax,Xax,reshape(Ic,N-2,N-2));  axis ij equal tight; box on; colormap(copper); colorbar; drawnow;

    % process correlation angles of cluster patterns
    maxsum = zeros(180,cc);
    for th = 1:1:181
        for ic = 1:cc
            maxsum(th,ic) = max(sum(imrotate(reshape(Ic==ic,N-2,N-2),th-1),1));        
        end

    end

    % process correlation angles of cluster patterns
    figure(figno); clf; figno=figno+1;
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

    % prepare FFT to process wavelengths of localisation patterns
    Fs  = 1/h;              % frequency of grid spacing
    n   = 2^nextpow2(N-2);  % next higher power of two

    F  = reshape(sum(Ap.^2./c,2).^1,N-2,N-2);  % get intensity of PC weights
%     F  = reshape(Ap(:,1:3),N-2,N-2);  % get intensity of PC1 weights

    figure(figno); clf; figno=figno+1;
    imagesc(Xax,Xax,log10(F)); axis equal tight; colorbar; colormap(copper); drawnow;

    % get horizontal FFT-spectrum
    S = fft(F,n,2);                  % raw FFT spectrum
    
    S2 = abs(S/(N-2));               % double-sided spectrum
    S1 = 2*S2(:,1:n/2+1);            % single-sided spectrum
    
    % plot horizontal FFT-spectrum
    figure(figno); clf;
    subplot(2,1,1);
    plot(1./((2*Fs/n):(Fs/n):(Fs/2-Fs/n)),sum(S1(:,3:n/2))); axis tight;
    title('horizontal FFT-spectrum of |PC|','Interpreter','latex','FontSize',18);
    xlabel('Period [$\delta_0$]','Interpreter','latex','FontSize',16);
    ylabel('Spectral Energy','Interpreter','latex','FontSize',16);
    
    % get vertical FFT-spectrum
    S = fft(rot90(F),n,2);           % raw FFT spectrum
    
    S2 = abs(S/(N-2));               % double-sided spectrum
    S1 = 2*S2(:,1:n/2+1);            % single-sided spectrum
    
    % plot vertical FFT-spectrum
    figure(figno);
    subplot(2,1,2)
    plot(1./((2*Fs/n):(Fs/n):(Fs/2-Fs/n)),sum(S1(:,3:n/2))); axis tight;
    title('vertical FFT-spectrum of |PC|','Interpreter','latex','FontSize',18);
    xlabel('Period [$\delta_0$]','Interpreter','latex','FontSize',16);
    ylabel('Spectral Energy','Interpreter','latex','FontSize',16);
    drawnow;
end