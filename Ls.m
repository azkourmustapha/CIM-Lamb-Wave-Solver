function Lst = Ls(k,F)
%Case cs & cl are known
%cs = 2100;
%cl = 2148;
%Case cs & cl are unknown
E = 70*1.0e+9;
v = 0.32;
D = 2700;
cl = sqrt((E*(1-v))/(D*(1+v)*(1-2*v)));  
cs = sqrt(E/(2*D*(1+v)));
ep = 2;
d = ep*1.0e-3; %thickness
coeff = 1.0e+6/ep;
p = abs(sqrt(1-((2*pi*F*coeff/k)/cs).^2));
q = abs(sqrt(1-((2*pi*F*coeff/k)/cl).^2));
if F <= cs*k/(2*pi*coeff)
    Lst = (-4*p*q*sinh(k*d*q/2)*cosh(k*d*p/2)+(1+p^2)^2*sinh(k*d*p/2)*cosh(k*d*q/2));
elseif F > cs*k/(2*pi*coeff) && F <= k*cl/(2*pi*coeff)
    Lst = (-4*p*q*sinh(k*d*q/2)*cos(k*d*p/2)+(1-p^2)^2*sin(k*d*p/2)*cosh(k*d*q/2));
else
    Lst = (4*p*q*sin(k*d*q/2)*cos(k*d*p/2)+(1-p^2)^2*sin(k*d*p/2)*cos(k*d*q/2));
end
end

