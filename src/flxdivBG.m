vz = WBG;
vx = UBG;
wp = vz(2:end  ,2:end-1);
wm = vz(1:end-1,2:end-1);
up = vx(2:end-1,2:end  );
um = vx(2:end-1,1:end-1);

a  = 1/f0 - f;

agh                    = zeros(N+2,N+2);
agh(2:end-1,2:end-1)   = a;

Xbnd = round(X(:,end-1)./h);  Xbnd(Xbnd<2) = Xbnd(Xbnd<2)+(N-3);  Xbnd(Xbnd>N-2) = Xbnd(Xbnd>N-2)-(N-3);
Zbnd = round(Z(end-1,:)./h);  Zbnd(Zbnd<2) = Zbnd(Zbnd<2)+(N-3);  Zbnd(Zbnd>N-2) = Zbnd(Zbnd>N-2)-(N-3);

idx = sub2ind(size(dr), repmat((1:N).',1,4), Xbnd+[-1 0 1 2]);
agh(2:end-1,[1 2 end-1 end]) = 1/f0 - 1 - f1.*dr(idx);

idx = sub2ind(size(dr), Zbnd+[-1 0 1 2].',repmat((1:N),4,1));
agh([1 2 end-1 end],2:end-1) = 1/f0 - 1 - f1.*dr(idx);

acc = agh(3:end-2,3:end-2);
ajp = agh(4:end-1,3:end-2);  ajpp = agh(5:end-0,3:end-2);
ajm = agh(2:end-3,3:end-2);  ajmm = agh(1:end-4,3:end-2);
aip = agh(3:end-2,4:end-1);  aipp = agh(3:end-2,5:end-0);
aim = agh(3:end-2,2:end-3);  aimm = agh(3:end-2,1:end-4);

Div_fVBG(2:end-1,2:end-1)   =     up .*(-(aipp-aip)./h./8 + (aip + acc)./h./2 + (acc-aim )./h./8) ...
                            - abs(up).*(-(aipp-aip)./h./8 + (aip - acc)./h./4 - (acc-aim )./h./8) ...
                            -     um .*(-(aip -acc)./h./8 + (acc + aim)./h./2 + (aim-aimm)./h./8) ...
                            + abs(um).*(-(aip -acc)./h./8 + (acc - aim)./h./4 - (aim-aimm)./h./8) ...
                            +     wp .*(-(ajpp-ajp)./h./8 + (ajp + acc)./h./2 + (acc-ajm )./h./8) ...
                            - abs(wp).*(-(ajpp-ajp)./h./8 + (ajp - acc)./h./4 - (acc-ajm )./h./8) ...
                            -     wm .*(-(ajp -acc)./h./8 + (acc + ajm)./h./2 + (ajm-ajmm)./h./8) ...
                            + abs(wm).*(-(ajp -acc)./h./8 + (acc - ajm)./h./4 - (ajm-ajmm)./h./8);

clear vz vx wp wm up um a agh acc ajp ajpp ajm ajmm aip aipp aim aimm