close all;

frame = floor(step/nop);

fprintf('\n*****  saving output frame %d to out/%s \n',frame,runID);

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
fh1 = figure; colormap(ocean);
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
imagesc(xc,zc,-W(:      ,2:end-1)-WBG(:      ,2:end-1)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['matrix z-velocity $W$'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
set(fh1, 'CurrentAxes', ax(2))
imagesc(xc,zc, U(2:end-1,:      )+UBG(2:end-1,:      )); axis ij equal tight; box on; cb = colorbar;
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
fh2 = figure; colormap(ocean);
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
set(fh2, 'CurrentAxes', ax(2))
imagesc(xc,zc,log10(   K(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ segr. coefficient $K$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(fh2, 'CurrentAxes', ax(3))
imagesc(xc,zc,log10(eII(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ strain rate $\varepsilon_{II}$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(fh2, 'CurrentAxes', ax(4))
imagesc(xc,zc,log10(eta(2:end-1,2:end-1)+etamin)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ eff. viscosity $\eta$'],TX{:},FS{:});
set(fh2, 'CurrentAxes', ax(5))
imagesc(xc,zc,log10(chi(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ elastic par. $\chi$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
set(fh2, 'CurrentAxes', ax(6))
imagesc(xc,zc,log10( tII(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ shear stress $\tau_{II}$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
drawnow;

% prepare and plot figure for solution residuals
if plot_cv
    fh3 = figure; colormap(ocean);
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh3,UN{:},'Position',[15 15 fw fh]);
    set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh3,'Color','w','InvertHardcopy','off');
    set(fh3,'Resize','off');
    ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(3) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(4) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

    set(fh3, 'CurrentAxes', ax(1))
    imagesc(xc,zc,log10(abs(res_W.*dtW  ./(1e-16+norm(W(:)+WBG(:),2)./N)))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. matrix z-velocity $W$'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
    set(fh3, 'CurrentAxes', ax(2))
    imagesc(xc,zc,log10(abs(res_U.*dtU  ./(1e-16+norm(U(:)+UBG(:),2)./N)))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. matrix x-velocity $U$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    set(fh3, 'CurrentAxes', ax(3))
    imagesc(xc,zc,log10(abs(res_P.*dtP  ./(1e-16+norm(P(:),2)       ./N)))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. pore pressure $P$'],TX{:},FS{:});
    set(fh3, 'CurrentAxes', ax(4))
    imagesc(xc,zc,log10(abs(res_f.*dt/10./(1e-16+norm(f(:),2)       ./N)))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. liquid fraction $\phi$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    drawnow;
end

if save_op
    name = ['../out/',runID,'/',runID,'_sol_',num2str(frame)];
    print(fh1,name,'-dpng','-r300');
    name = ['../out/',runID,'/',runID,'_mat_',num2str(frame)];
    print(fh2,name,'-dpng','-r300');
    
    name = ['../out/',runID,'/',runID,'_cont'];
    save([name,'.mat'],'U','W','P','u','w','p','f','etay','eta','zeta','chi','K','exx','ezz','exz','txx','tzz','txz','eII','tII','Div_V','Div_v','time','step','fmass0');
    name = ['../out/',runID,'/',runID,'_',num2str(frame)];
    save([name,'.mat'],'U','W','P','u','w','p','f','etay','eta','zeta','chi','K','exx','ezz','exz','txx','tzz','txz','eII','tII','Div_V','Div_v','time','step','fmass0');
    
    if step == 1
        logfile = ['../out/',runID,'/',runID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end
    