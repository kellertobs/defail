% get frame number and start preparing output
frame = floor(step/nop);
fprintf('\n*****  preparing output frame %d for %s \n',frame,runID);

% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};
UN = {'Units','Centimeters'};

axh = 6.00; axw = 7.50;
ahs = 1.00; avs = 1.00;
axb = 0.75; axt = 0.90;
axl = 1.75; axr = 0.90;
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 3*axw + 2*ahs + axr;

set(0,'DefaultFigureVisible',plot_op)

% prepare and plot figure for solution fields
fh1 = figure(1); clf; colormap(ocean);
set(fh1,UN{:},'Position',[5 5 fw fh]);
set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh1,'Color','w','InvertHardcopy','off');
set(fh1,'Resize','off');
ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(3) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
ax(4) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(5) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(6) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

set(fh1, 'CurrentAxes', ax(1))
imagesc(xc,zc,-W(:      ,2:end-1)-0.*WBG(:      ,2:end-1)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['matrix z-velocity $W$'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
text(0,0.4*L,['time ',num2str(time,4)],TX{:},FS{:},'HorizontalAlignment','center','VerticalAlignment','middle')
set(fh1, 'CurrentAxes', ax(2))
imagesc(xc,zc, U(2:end-1,:      )+0.*UBG(2:end-1,:      )); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['matrix x-velocity $U$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(fh1, 'CurrentAxes', ax(3))
imagesc(xc,zc, P(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['pore pressure $P$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(fh1, 'CurrentAxes', ax(4))
imagesc(xc,zc,-w(:      ,2:end-1)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['segr. z-velocity $w$'],TX{:},FS{:});
set(fh1, 'CurrentAxes', ax(5))
imagesc(xc,zc, u(2:end-1,:      )); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['segr. x-velocity $u$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
set(fh1, 'CurrentAxes', ax(6))
imagesc(xc,zc, p(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['comp. pressure $p$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
drawnow;

% prepare and plot figure for material properties
fh2 = figure(2); clf; colormap(ocean);
set(fh2,UN{:},'Position',[10 10 fw fh]);
set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh2,'Color','w','InvertHardcopy','off');
set(fh2,'Resize','off');
ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(3) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
ax(4) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(5) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(6) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

set(fh2, 'CurrentAxes', ax(1))
imagesc(xc,zc,         f(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['liquid fraction $\phi$'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
text(0,0.4*L,['time ',num2str(time,4)],TX{:},FS{:},'HorizontalAlignment','center','VerticalAlignment','middle')
set(fh2, 'CurrentAxes', ax(2))
imagesc(xc,zc,     ups(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['decompaction rate $\dot{\upsilon}$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(fh2, 'CurrentAxes', ax(3))
imagesc(xc,zc,log10(eps(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ strain rate $\dot{\varepsilon}$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(fh2, 'CurrentAxes', ax(4))
imagesc(xc,zc,(eta(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ shear viscosity $\eta$'],TX{:},FS{:});
set(fh2, 'CurrentAxes', ax(5))
imagesc(xc,zc,(zeta(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ comp. viscosity $\zeta$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
set(fh2, 'CurrentAxes', ax(6))
imagesc(xc,zc,log10(tau(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ shear stress $\tau$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
drawnow;

% prepare and plot figure for solution residuals
if plot_cv
    fh3 = figure(3); clf; colormap(ocean);
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh3,UN{:},'Position',[15 15 fw fh]);
    set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh3,'Color','w','InvertHardcopy','off');
    set(fh3,'Resize','off');
    ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(3) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(4) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

    if bnchmrk
        set(fh3, 'CurrentAxes', ax(1))
        imagesc(xc,zc,-(W-W_mms)./(1e-16+norm(W_mms(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['err. matrix z-velocity $W$'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
        set(fh3, 'CurrentAxes', ax(2))
        imagesc(xc,zc, (U-U_mms)./(1e-16+norm(U_mms(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['err. matrix x-velocity $U$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        set(fh3, 'CurrentAxes', ax(3))
        imagesc(xc,zc, (P-P_mms)./(1e-16+norm(P_mms(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['err. pore pressure $P$'],TX{:},FS{:});
        set(fh3, 'CurrentAxes', ax(4))
        imagesc(xc,zc, (f-f_mms)./(1e-16+norm(f_mms(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['err. liquid fraction $\phi$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
        drawnow;
    else
        set(fh3, 'CurrentAxes', ax(1))
        imagesc(xc,zc,(dW./(1e-16+norm(W(:)+WBG(:),2)./N))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. matrix z-velocity $W$'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
        set(fh3, 'CurrentAxes', ax(2))
        imagesc(xc,zc,(dU./(1e-16+norm(U(:)+UBG(:),2)./N))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. matrix x-velocity $U$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        set(fh3, 'CurrentAxes', ax(3))
        imagesc(xc,zc,(dP./(1e-16+norm(P(:),2)       ./N))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. pore pressure $P$'],TX{:},FS{:});
        set(fh3, 'CurrentAxes', ax(4))
        imagesc(xc,zc,(df./(1e-16+norm(f(:),2)       ./N))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. liquid fraction $\phi$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
        drawnow;
    end
    
    figure(4); clf;
    p0 = linspace(-1,max(p(:)),1e3);
    plot(p0,eps0.*ones(size(p0)),'k',p0,1+p0+bnchmrk*5,'r',p(:),yieldt(:),'r.','LineWidth',2); axis equal tight; box on; hold on;
    scatter(p(:),tau(:),20,(eta(:)),'filled'); colorbar; colormap(ocean);
    title('Failure Criterion');
    drawnow;
end

% save output frame
if save_op
    name = ['../out/',runID,'/',runID,'_sol_',num2str(frame)];
    print(fh1,name,'-dpng','-r300');
    name = ['../out/',runID,'/',runID,'_mat_',num2str(frame)];
    print(fh2,name,'-dpng','-r300');
    
    clear fh1 fh2 fh3 ax cb fw axw axh avs ahs axl axr axt axb fh fw TX TL FS TS UN p0 dWi dWii dUi dUii dPi dPii
    clear Wi Ui Pi dtWi dtUi dtPi dtWii dtUii dtPii fo Div_fVo Div_fVBGo Div_tz Div_tx dtW dtU dtP etac k kk

    name = ['../out/',runID,'/',runID,'_cont'];
    save([name,'.mat']);
    name = ['../out/',runID,'/',runID,'_',num2str(frame)];
    save([name,'.mat']);
    
    if step == 0
        logfile = ['../out/',runID,'/',runID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end