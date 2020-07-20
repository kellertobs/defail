% generate solitary wave for benchmarking

A        =  f2+1;
fsw      =  linspace(1,A,1e6);

zsw      =  sqrt(A+0.5) .* (-2.*sqrt(A-fsw) + 1/sqrt(A-1) .* log( (sqrt(A-1)-sqrt(A-fsw))./(sqrt(A-1)+sqrt(A-fsw))));
zsw(1)   =  -20;
zsw      =  [zsw,-zsw(1:end-1)];
[zsw,i]  =  sort(zsw);
fsw      =  [fsw,fsw(1:end-1)];
fsw      =  fsw(i);

zsw      =  zsw.*sqrt(1+2/3);

f1D      =  interp1(zsw,fsw,z,'linear',1);

f        =  repmat(f1D.',1,length(z));
fsw      =  f;
csw      =  (2*A + 1)./f0;