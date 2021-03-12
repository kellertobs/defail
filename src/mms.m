% create manufactured solution
clear x z;
syms U_mms(x,z) W_mms(x,z) P_mms(x,z) u_mms(x,z) w_mms(x,z) p_mms(x,z) f_mms(x,z)

fprintf(1,'  ***  compose manufactured solution\n\n');

% compose manufactured solution variables
W_mms(x,z) =   0.50.*(cos(4*(x-z)*pi/L).*sin(2*(x+z)*pi/L));
U_mms(x,z) =   0.25.*(sin(4*(x-z)*pi/L).*cos(2*(x+z)*pi/L));
P_mms(x,z) =   1.00.*(cos(4*(x-z)*pi/L).*cos(2*(x+z)*pi/L));
f_mms(x,z) = 1+0.25.*(sin(4*(x-z)*pi/L).*sin(2*(x+z)*pi/L));

fprintf(1,'       W = %s \n',char(W_mms));
fprintf(1,'       U = %s \n',char(U_mms));
fprintf(1,'       P = %s \n',char(P_mms));
fprintf(1,'       f = %s \n',char(f_mms));
fprintf(1,'       . ');

% update strain rates
exx_mms(x,z) = diff(U_mms,x) - (diff(W_mms,z) + diff(U_mms,x))./3 - Pu;    % x-normal strain rate
ezz_mms(x,z) = diff(W_mms,z) - (diff(W_mms,z) + diff(U_mms,x))./3 + Pu;    % z-normal strain rate
exz_mms(x,z) = 1/2.*(diff(U_mms,z)+diff(W_mms,x))                 - Si;    % shear strain rate
fprintf(1,' . ');

% update coefficients
 eps_mms(x,z) = sqrt(1/2*(exx_mms.^2 + ezz_mms.^2 + 2*exz_mms.^2));
 eta_mms(x,z) = exp(-lmd*f0.*(f_mms-1)) .* (eps_mms./(abs(Pu)+abs(Si))).^((1-n)./n);
zeta_mms(x,z) = eta_mms ./ (f0.*f_mms);
K_mms(x,z)    = f_mms.^m;
fprintf(1,' . ');

% update stresses
txx_mms(x,z) = eta_mms .* exx_mms;                                         % x-normal stress
tzz_mms(x,z) = eta_mms .* ezz_mms;                                         % z-normal stress
txz_mms(x,z) = eta_mms .* exz_mms;                                         % xz-shear stress
fprintf(1,' . ');

% update segregation/compaction variables
w_mms(x,z) = -  K_mms   .* (diff(P_mms,z) + B);
u_mms(x,z) = -  K_mms   .* (diff(P_mms,x));
p_mms(x,z) = - zeta_mms .* (diff(W_mms,z) + diff(U_mms,x));
fprintf(1,' . ');

% manufactured solution residuals
res_W_mms = - (diff(tzz_mms,z) + diff(txz_mms,x)) + diff(P_mms,z) + diff(p_mms,z);
res_U_mms = - (diff(txx_mms,x) + diff(txz_mms,z)) + diff(P_mms,x) + diff(p_mms,x);
res_P_mms =   (diff(W_mms,z) + diff(U_mms,x)) + (diff(w_mms,z) + diff(u_mms,x));
fprintf(1,' . ');

% plot manufactured solution
figure(15);
colormap(ocean);
subplot(2,4,1); fcontour( -W_mms,[-L/2,L/2]); axis ij equal tight; colorbar; box on; title('manufactured $W$','Interpreter','latex');
subplot(2,4,2); fcontour(  U_mms,[-L/2,L/2]); axis ij equal tight; colorbar; box on; title('manufactured $U$','Interpreter','latex');
subplot(2,4,3); fcontour(  P_mms,[-L/2,L/2]); axis ij equal tight; colorbar; box on; title('manufactured $P$','Interpreter','latex');
subplot(2,4,4); fcontour(  f_mms,[-L/2,L/2]); axis ij equal tight; colorbar; box on; title('manufactured $f$','Interpreter','latex');
fprintf(1,' . ');
subplot(2,4,5); fcontour( -w_mms,[-L/2,L/2]); axis ij equal tight; colorbar; box on; title('manufactured $w$','Interpreter','latex');
subplot(2,4,6); fcontour(  u_mms,[-L/2,L/2]); axis ij equal tight; colorbar; box on; title('manufactured $u$','Interpreter','latex');
subplot(2,4,7); fcontour(  p_mms,[-L/2,L/2]); axis ij equal tight; colorbar; box on; title('manufactured $p$','Interpreter','latex');
subplot(2,4,8); fcontour(eta_mms,[-L/2,L/2]); axis ij equal tight; colorbar; box on; title('manufactured $\eta$','Interpreter','latex');
drawnow;
fprintf(1,' . \n');

% evaluate mms source terms on appropriate coordinate grids
fprintf(1,'\n  ***  evaluate manufactured solution\n\n');
x_mms  = -L/2-h/2:h:L/2+h/2;
z_mms  = -L/2-h/2:h:L/2+h/2;
xu_mms = (x_mms(1:end-1)+x_mms(2:end))./2;
zw_mms = (z_mms(1:end-1)+z_mms(2:end))./2;

fprintf(1,'       Patience, my young Padawan!\n');
fprintf(1,'       . ');

[x,z] = meshgrid(x_mms,zw_mms);
src_W_mms = double(subs(res_W_mms)); fprintf(1,' . ');
[x,z] = meshgrid(xu_mms,z_mms);
src_U_mms = double(subs(res_U_mms)); fprintf(1,' . ');
[x,z] = meshgrid(x_mms,z_mms);
src_P_mms = double(subs(res_P_mms)); fprintf(1,' . ');

% plot manufactured residuals and evaluated source terms
figure(16);
colormap(ocean);
subplot(2,3,1); fcontour(-res_W_mms,[-L/2,L/2]); axis ij equal tight; colorbar; box on; title('manufactured $W$-res','Interpreter','latex');
subplot(2,3,2); fcontour(-res_U_mms,[-L/2,L/2]); axis ij equal tight; colorbar; box on; title('manufactured $U$-res','Interpreter','latex');
subplot(2,3,3); fcontour(-res_P_mms,[-L/2,L/2]); axis ij equal tight; colorbar; box on; title('manufactured $P$-res','Interpreter','latex');
subplot(2,3,4); imagesc(x_mms,zw_mms,-src_W_mms); axis ij equal tight; colorbar; box on; title('evaluated $W$-res','Interpreter','latex');
subplot(2,3,5); imagesc(xu_mms,z_mms,-src_U_mms); axis ij equal tight; colorbar; box on; title('evaluated $U$-res','Interpreter','latex');
subplot(2,3,6); imagesc(x_mms ,z_mms,-src_P_mms); axis ij equal tight; colorbar; box on; title('evaluated $P$-res','Interpreter','latex');
drawnow;

% evaluate analytical solution on appropriate coordinate grids
[x,z] = meshgrid(x_mms,zw_mms);
W_mms = double(subs(W_mms)); fprintf(1,' . ');
[x,z] = meshgrid(xu_mms,z_mms);
U_mms = double(subs(U_mms)); fprintf(1,' . ');
[x,z] = meshgrid(x_mms,z_mms);
P_mms = double(subs(P_mms)); fprintf(1,' . ');
f_mms = double(subs(f_mms)); fprintf(1,' . \n');

M = 1;  % only initial step for benchmarking
