tic 
fti               = 0;
ftf               = 10; 
dft               = 1e-2;                 %___initial step frequency*thickness in (MHz.mm).
RunFT             = ftf;                  %___limit of frequency*thickness in (MHz.mm).
n                 = round((ftf-fti)/dft);
epsilon           = 1e-5;
nom               = 11;
ki                = 0;
kf                = 10995; 
nsub              = 600;
dk                = kf/nsub;
k                 = ki:dk:kf;
%______________________________Symmetric modes_____________________________
fts               = zeros(length(k),nom);
Cps               = zeros(length(k),nom);
Ws                = zeros(length(k),nom);
Cgs               = zeros(length(k),nom);
ERS               = zeros(length(k),nom);
VERS              = zeros(length(k),1);
IMS               = zeros(length(k),nom);
RATIOS            = zeros(length(k),nom); 
thick             = 2; 
kilo              = 1e3;
mega              = 1e6;
nbs               = 0;
j                 = 0;
for m             = 1:length(k)
    for i         = 1:n+1
        ai        = fti+(i-1)*dft;
        while ai < ftf
            aj    = ai+dft;
            di    = aj-ai;
            m0    = ai+di/2;
            m1    = ai+di/3;
            m2    = ai+2*di/3;
            Lsai  = Ls(k(m),ai); %__ai = freq*thick, freq: frequency, thick: plate thickness, k(m): wavenumber k, Ls: symmetric Lamb modes function.
            Lsaj  = Ls(k(m),aj);
            Lsm0  = Ls(k(m),m0);
            Lsm1  = Ls(k(m),m1);
            Lsm2  = Ls(k(m),m2);
            if Lsai == 0
                nbs        = nbs+1;
                fts(m,nbs) = ai;
                Cps(m,nbs) = (2*pi*ai*mega)/(k(m)*thick);
                Ws(m,nbs)  = Cps(m,nbs)/(ai*mega/thick);
                if m >= 2
                    if fts(m,nbs)*fts(m-1,nbs) ~= 0
                        Cgs(m,nbs) = 2*pi*(fts(m,nbs)-fts(m-1,nbs))*mega/(thick*(k(m)-k(m-1)));
                    else
                        break
                    end
                end
            end
            if Lsaj == 0
                nbs        = nbs+1;
                fts(m,nbs) = aj;
                Cps(m,nbs) = (2*pi*aj*mega)/(k(m)*thick);
                Ws(m,nbs)  = Cps(m,nbs)/(aj*mega/thick);
                if m >= 2
                    if fts(m,nbs)*fts(m-1,nbs) ~= 0
                        Cgs(m,nbs) = 2*pi*(fts(m,nbs)-fts(m-1,nbs))*mega/(thick*(k(m)-k(m-1)));
                    else
                        break
                    end
                end
            end
            if Lsm0 == 0
                nbs        = nbs+1;
                fts(m,nbs) = m0;
                Cps(m,nbs) = (2*pi*m0*mega)/(k(m)*thick);
                Ws(m,nbs)  = Cps(m,nbs)/(m0*mega/thick);
                if m >= 2
                    if fts(m,nbs)*fts(m-1,nbs) ~= 0
                        Cgs(m,nbs) = 2*pi*(fts(m,nbs)-fts(m-1,nbs))*mega/(thick*(k(m)-k(m-1)));
                    else
                        break
                    end
                end
            end
            if Lsm1 == 0
                nbs        = nbs+1;
                fts(m,nbs) = m1;
                Cps(m,nbs) = (2*pi*m1*mega)/(k(m)*thick);
                Ws(m,nbs)  = Cps(m,nbs)/(m1*mega/thick);
                if m >= 2
                    if fts(m,nbs)*fts(m-1,nbs) ~= 0
                        Cgs(m,nbs) = 2*pi*(fts(m,nbs)-fts(m-1,nbs))*mega/(thick*(k(m)-k(m-1)));
                    else
                        break
                    end
                end
            end
            if Lsm2 == 0
                nbs        = nbs+1;
                fts(m,nbs) = m2;
                Cps(m,nbs) = (2*pi*m2*mega)/(k(m)*thick);
                Ws(m,nbs)  = Cps(m,nbs)/(m2*mega/thick);
                if m >= 2
                    if fts(m,nbs)*fts(m-1,nbs) ~= 0
                        Cgs(m,nbs) = 2*pi*(fts(m,nbs)-fts(m-1,nbs))*mega/(thick*(k(m)-k(m-1)));
                    else
                        break
                    end
                end
            end
            if sign(Lsai/Lsaj) == 1
                break
            end
            if sign(Lsai/Lsaj) == -1
                err = max(max(max(abs(Lsai),abs(Lsaj)),max(abs(Lsm1),abs(Lsm2))),abs(Lsm0));
                while err/epsilon > 1
                    di   = aj-ai;
                    m0   = ai+di/2;
                    m1   = ai+di/3;
                    m2   = ai+2*di/3;
                    Lsai = Ls(k(m),ai);
                    Lsaj = Ls(k(m),aj);
                    Lsm0 = Ls(k(m),m0);
                    Lsm1 = Ls(k(m),m1);
                    Lsm2 = Ls(k(m),m2);
                    %______________________________________________________
                    lambda1 = (Lsm1*ai)/(m1^2*(ai-m1));
                    lambda2 = (Lsai*m1)/(ai^2*(ai-m1));
                    lambda3 = 1/ai+1/m1;
                    lambda4 = 1/ai^2+1/m1^2+1/(ai*m1);
                    T1      = m2*(m2^2/(ai*m1)-m2/m1+1);
                    T2      = (m2^3/ai)*(lambda4-1/(ai^2))-lambda4*m2^2+1;
                    T3      = (m2^3/ai)*(Lsai/(ai^2)-(lambda1-lambda2))+(lambda1-lambda2)*m2^2;
                    Z1      = (aj^3)*((1/(ai*m1))*(lambda3-T2/T1)+1/ai^2)+(aj^2)*(lambda3*(T2/T1)-lambda4)-(T2/T1)*aj+1;
                    Z2      = (aj^3)*(-Lsai/(ai^3)+(lambda1-lambda2)/ai+(1/(ai^2)-lambda3/ai)*((Lsm2/T1)-(T3/T1)))+(aj^2)*(lambda3*((Lsm2/T1)-(T3/T1))-(lambda1-lambda2))-((Lsm2/T1)-(T3/T1))*aj;
                    cof4    = (Lsaj+Z2)/Z1;
                    cof3    = Lsm2/T1-T2*cof4/T1-T3/T1;
                    cof2    = lambda1-lambda2-lambda3*cof3-lambda4*cof4;
                    cof1    = Lsai/(ai^3)-(cof2/ai)-cof3/(ai^2)-cof4/(ai^3);
                    %______________________________________________________
                    if ai == 0 || ai-m1 == 0 || T1 == 0 || Z1 == 0
                        ai = aj;
                        aj = ai+dft;
                    end
                    %_Calculation of the solutions of the cubic polynomial:
                    if cof1 ~= 0
                        A     = cof2/cof1;
                        B     = cof3/cof1;
                        C     = cof4/cof1;
                        p     = (B-A^2/3);
                        q     = (2*A^3)/27-A*B/3+C;
                        delta = sqrt(4*p^3+27*q^2);
                        dis   = 18*cof1*cof2*cof3*cof4-4*cof2^3*cof4+cof2^2*cof3^2-4*cof1*cof3^3-27*cof1^2*cof4^2;
                        if p == 0 && q == 0
                            s = -A/3;
                            if s > ai && s < aj
                                j      = j+1;
                                cis(j) = s; % cis: cubic interpolation solution.
                            end
                            if s < ai || s > aj
                                j      = j+1;
                                cis(j) = m0;
                            end
                        end
                        if p == 0 && q ~= 0
                            s = -A/3+(-q)^(1/3);
                            if s > ai && s < aj
                                j      = j+1;
                                cis(j) = s;
                            end
                            if s < ai || s > aj
                                j      = j+1;
                                cis(j) = m0;
                            end
                        end
                        if p > 0 && q == 0
                            s = -A/3;
                            if s > ai && s < aj
                                j      = j+1;
                                cis(j) = s;
                            end
                            if s < ai || s > aj
                                j      = j+1;
                                cis(j) = m0;
                            end
                        end
                        if p < 0 && q == 0
                            s1    = -A/3;
                            s2    = -A/3+sqrt(-p);
                            s3    = -A/3-sqrt(-p);
                            Lss1  = Ls(k(m),s1);
                            Lss2  = Ls(k(m),s2);
                            Lss3  = Ls(k(m),s3);
                            mn    = min(s1,min(s2,s3));
                            mx    = max(s1,max(s2,s3));
                            Lsmin = min(abs(Lss1),min(abs(Lss2),abs(Lss3)));
                            if s1-s2 == 0 && s2-s3 == 0
                                s = s1;
                                if s > ai && s < aj
                                    j      = j+1;
                                    cis(j) = s;
                                end
                                if s < ai || s > aj
                                    j      = j+1;
                                    cis(j) = m0;
                                end
                            end
                            if (s1-s2)*(s1-s3)*(s2-s3) ~= 0
                                if (s1-mn)*(s1-mx) ~= 0
                                    mid = s1;
                                end
                                if (s2-mn)*(s2-mx) ~= 0
                                    mid = s2;
                                end
                                if (s3-mn)*(s3-mx) ~= 0
                                    mid = s3;
                                end
                                Lsmid = Ls(k(m),mid);
                                if mx < ai || mn > aj
                                    j      = j+1;
                                    cis(j) = m0;
                                end
                                if mn > ai && mx < aj
                                    if Lsmin == abs(Lsmn)
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                    if Lsmin == abs(Lsmid)
                                        j      = j+1;
                                        cis(j) = mid;
                                    end
                                    if Lsmin == abs(Lsmx)
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                                if mn < ai && mx > aj
                                    if mid > ai && mid < aj
                                        j      = j+1;
                                        cis(j) = mid;
                                    end
                                    if mid < ai || mid > aj
                                        j      = j+1;
                                        cis(j) = m0;
                                    end
                                end
                                if mn > ai && mx > aj
                                    if mid < aj
                                        mn    = min(mn,mid);
                                        mx    = max(mn,mid);
                                        Lsmn  = Ls(k(m),mn);
                                        Lsmx  = Ls(k(m),mx);
                                        Lsmin = min(abs(Lsmn),abs(Lsmx));
                                        if Lsmin == abs(Lsmn)
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                        if Lsmin == abs(Lsmid)
                                            j      = j+1;
                                            cis(j) = mid;
                                        end
                                    end
                                    if mid > aj
                                        if mn < aj
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                        if mn > aj
                                            j      = j+1;
                                            cis(j) = m0;
                                        end
                                    end
                                end
                                if mn < ai && mx < aj
                                    if mx < ai
                                        j      = j+1;
                                        cis(j) = m0;
                                    end
                                    if mx > ai
                                        if mid > ai
                                            mn    = min(mid,mx);
                                            mx    = max(mid,mx);
                                            Lsmn  = Ls(k(m),mn);
                                            Lsmx  = Ls(k(m),mx);
                                            Lsmin = min(abs(Lsmn),abs(Lsmx));
                                            if Lsmin == abs(Lsmid)
                                                j      = j+1;
                                                cis(j) = mid;
                                            end
                                            if Lsmin == abs(Lsmx)
                                                j      = j+1;
                                                cis(j) = mx;
                                            end
                                        end
                                        if mid < ai
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                end
                            end
                            if s1-s2 == 0 && s1-s3 ~= 0
                                mn    = min(s1,s3);
                                mx    = max(s1,s3);
                                Lsmn  = Ls(k(m),mn);
                                Lsmx  = Ls(k(m),mx);
                                Lsmin = min(abs(Lsmn),abs(Lsmx));
                                if mn > aj || mx < ai
                                    j      = j+1;
                                    cis(j) = m0;
                                end
                                if mn > ai && mx < aj
                                    if Lsmin == abs(Lsmn)
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                    if Lsmin == abs(Lsmx)
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                                if mn > ai && mn < aj
                                    if mx > aj
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                end
                                if mx > ai && mx < aj
                                    if mn < ai
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                            end
                            if s1-s3 == 0 && s1-s2 ~= 0
                                mn    = min(s1,s2);
                                mx    = max(s1,s2);
                                Lsmn  = Ls(k(m),mn);
                                Lsmx  = Ls(k(m),mx);
                                Lsmin = min(abs(Lsmn),abs(Lsmx));
                                if mn > aj || mx < ai
                                    j      = j+1;
                                    cis(j) = m0;
                                end
                                if mn > ai && mx < aj
                                    if Lsmin == abs(Lsmn)
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                    if Lsmin == abs(Lsmx)
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                                if mn > ai && mn < aj
                                    if mx > aj
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                end
                                if mx > ai && mx < aj
                                    if mn < ai
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                            end
                            if s2-s3 == 0 && s2-s1 ~= 0
                                mn    = min(s1,s2);
                                mx    = max(s1,s2);
                                Lsmn  = Ls(k(m),mn);
                                Lsmx  = Ls(k(m),mx);
                                Lsmin = min(abs(Lsmn),abs(Lsmx));
                                if mn > aj || mx < ai
                                    j      = j+1;
                                    cis(j) = m0;
                                end
                                if mn > ai && mx < aj
                                    if Lsmin == abs(Lsmn)
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                    if Lsmin == abs(Lsmx)
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                                if mn > ai && mn < aj
                                    if mx > aj
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                end
                                if mx > ai && mx < aj
                                    if mn < ai
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                            end
                        end
                        if p ~= 0 && q ~= 0
                            z1 = -A/3+((sqrt(3)*delta-9*q)^(1/3))/(2^(1/3)*3^(2/3))-(2/3)^(1/3)*p/((sqrt(3)*delta-9*q)^(1/3));
                            z2 = -A/3+(1+1i*sqrt(3))*p/(2^(2/3)*3^(1/3)*((sqrt(3)*delta-9*q)^(1/3)))-(1-1i*sqrt(3))*((sqrt(3)*delta-9*q)^(1/3))/(2*2^(1/3)*3^(2/3));
                            z3 = -A/3+(1-1i*sqrt(3))*p/(2^(2/3)*3^(1/3)*((sqrt(3)*delta-9*q)^(1/3)))-(1+1i*sqrt(3))*((sqrt(3)*delta-9*q)^(1/3))/(2*2^(1/3)*3^(2/3));
                            if dis >= 0
                                s1   = real(z1);
                                s2   = real(z2);
                                s3   = real(z3);
                                Lss1 = Ls(k(m),s1);
                                Lss2 = Ls(k(m),s2);
                                Lss3 = Ls(k(m),s3);
                                mn = min(s1,min(s2,s3));
                                mx = max(s1,max(s2,s3));
                                Lsmin = min(abs(Lss1),min(abs(Lss2),abs(Lss3)));
                                if s1-s2 == 0 && s2-s3 == 0
                                    s = s1;
                                    if s > ai && s < aj
                                        j      = j+1;
                                        cis(j) = s;
                                    end
                                    if s < ai || s > aj
                                        j      = j+1;
                                        cis(j) = m0;
                                    end
                                end
                                if (s1-s2)*(s1-s3)*(s2-s3) ~= 0
                                    if (s1-mn)*(s1-mx) ~= 0
                                        mid = s1;
                                    end
                                    if (s2-mn)*(s2-mx) ~= 0
                                        mid = s2;
                                    end
                                    if (s3-mn)*(s3-mx) ~= 0
                                        mid = s3;
                                    end
                                    Lsmid = Ls(k(m),mid);
                                    if mx < ai || mn > aj
                                        j      = j+1;
                                        cis(j) = m0;
                                    end
                                    if mn > ai && mx < aj
                                        if Lsmin == abs(Lsmn)
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                        if Lsmin == abs(Lsmid)
                                            j      = j+1;
                                            cis(j) = mid;
                                        end
                                        if Lsmin == abs(Lsmx)
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                    if mn < ai && mx > aj
                                        if mid > ai && mid < aj
                                            j      = j+1;
                                            cis(j) = mid;
                                        end
                                        if mid < ai || mid > aj
                                            j      = j+1;
                                            cis(j) = m0;
                                        end
                                    end
                                    if mn > ai && mx > aj
                                        if mid < aj
                                            mn    = min(mn,mid);
                                            mx    = max(mn,mid);
                                            Lsmn  = Ls(k(m),mn);
                                            Lsmx  = Ls(k(m),mx);
                                            Lsmin = min(abs(Lsmn),abs(Lsmx));
                                            if Lsmin == abs(Lsmn)
                                                j      = j+1;
                                                cis(j) = mn;
                                            end
                                            if Lsmin == abs(Lsmid)
                                                j      = j+1;
                                                cis(j) = mid;
                                            end
                                        end
                                        if mid > aj
                                            if mn < aj
                                                j      = j+1;
                                                cis(j) = mn;
                                            end
                                            if mn > aj
                                                j      = j+1;
                                                cis(j) = m0;
                                            end
                                        end
                                    end
                                    if mn < ai && mx < aj
                                        if mx < ai
                                            j      = j+1;
                                            cis(j) = m0;
                                        end
                                        if mx > ai
                                            if mid > ai
                                                mn    = min(mid,mx);
                                                mx    = max(mid,mx);
                                                Lsmn  = Ls(k(m),mn);
                                                Lsmx  = Ls(k(m),mx);
                                                Lsmin = min(abs(Lsmn),abs(Lsmx));
                                                if Lsmin == abs(Lsmid)
                                                    j      = j+1;
                                                    cis(j) = mid;
                                                end
                                                if Lsmin == abs(Lsmx)
                                                    j      = j+1;
                                                    cis(j) = mx;
                                                end
                                            end
                                            if mid < ai
                                                j      = j+1;
                                                cis(j) = mx;
                                            end
                                        end
                                    end
                                end
                                if s1-s2 == 0 && s1-s3 ~= 0
                                    mn    = min(s1,s3);
                                    mx    = max(s1,s3);
                                    Lsmn  = Ls(k(m),mn);
                                    Lsmx  = Ls(k(m),mx);
                                    Lsmin = min(abs(Lsmn),abs(Lsmx));
                                    if mn > aj || mx < ai
                                        j      = j+1;
                                        cis(j) = m0;
                                    end
                                    if mn > ai && mx < aj
                                        if Lsmin == abs(Lsmn)
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                        if Lsmin == abs(Lsmx)
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                    if mn > ai && mn < aj
                                        if mx > aj
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                    end
                                end
                                if s1-s3 == 0 && s1-s2 ~= 0
                                    mn    = min(s1,s2);
                                    mx    = max(s1,s2);
                                    Lsmn  = Ls(k(m),mn);
                                    Lsmx  = Ls(k(m),mx);
                                    Lsmin = min(abs(Lsmn),abs(Lsmx));
                                    if mn > aj || mx < ai
                                        j = j+1;
                                        cis(j) = m0;
                                    end
                                    if mn > ai && mx < aj
                                        if Lsmin == abs(Lsmn)
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                        if Lsmin == abs(Lsmx)
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                    if mn > ai && mn < aj
                                        if mx > aj
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                    end
                                    if mx > ai && mx < aj
                                        if mn < ai
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                end
                                if s2-s3 == 0 && s2-s1 ~= 0
                                    mn    = min(s1,s2);
                                    mx    = max(s1,s2);
                                    Lsmn  = Ls(k(m),mn);
                                    Lsmx  = Ls(k(m),mx);
                                    Lsmin = min(abs(Lsmn),abs(Lsmx));
                                    if mn > aj || mx < ai
                                        j      = j+1;
                                        cis(j) = m0;
                                    end
                                    if mn > ai && mx < aj
                                        if Lsmin == abs(Lsmn)
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                        if Lsmin == abs(Lsmx)
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                    if mn > ai && mn < aj
                                        if mx > aj
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                    end
                                    if mx > ai && mx < aj
                                        if mn < ai
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                end
                            end
                            if dis < 0
                                if z1 > ai && z1 < aj
                                    j      = j+1;
                                    cis(j) = z1;
                                end
                                if z1 < ai || z1 > aj
                                    j      = j+1;
                                    cis(j) = m0;
                                end
                            end
                        end
                    end
                    if cof1 == 0
                        if cof2 ~= 0
                            discr = cof3^2-4*cof2*cof4;
                            s1    = (-cof3-sqrt(discr))/(2*cof2);
                            s2    = (-cof3+sqrt(discr))/(2*cof2);
                            mn    = min(s1,s2);
                            mx    = max(s1,s2);
                            Lsmn  = Ls(k(m),mn);
                            Lsmx  = Ls(k(m),mx);
                            Lsmin = min(abs(Lsmn),abs(Lsmx));
                            if mn > aj || mx < ai
                                j      = j+1;
                                cis(j) = m0;
                            end
                            if mn > ai && mx < aj
                                if Lsmin == abs(Lsmn)
                                    j      = j+1;
                                    cis(j) = mn;
                                end
                                if Lsmin == abs(Lsmx)
                                    j      = j+1;
                                    cis(j) = mx;
                                end
                            end
                            if mn > ai && mn < aj
                                if mx > aj
                                    j      = j+1;
                                    cis(j) = mn;
                                end
                            end
                        end
                        if cof2 == 0
                            s = -cof4/cof3;
                            if s > ai && s < aj
                                j      = j+1;
                                cis(j) = s;
                            end
                            if s < ai || s > aj
                                j      = j+1;
                                cis(j) = m0;
                            end
                        end
                    end
                    Lscis   = Ls(k(m),cis(j));
                    xi      = 1e-10;
                    LsCISxi = Ls(k(m),cis(j)+xi);
                    dLsCIS  = (LsCISxi-Lscis)/xi;
                    if abs(Lscis)/epsilon < 1
                        nbs = nbs+1;
                        fts(m,nbs)    = cis(j);
                        ERS(m,nbs)    = abs(Lscis/dLsCIS);
                        Cps(m,nbs)    = (2*pi*cis(j)*mega)/(k(m)*thick);
                        Ws(m,nbs)     = Cps(m,nbs)/(cis(j)*mega/thick);
                        IMS(m,nbs)    = abs(Lscis);
                        RATIOS(m,nbs) = abs(Lscis/Lsm0);
                        if m >= 2
                            if fts(m,nbs)*fts(m-1,nbs) ~= 0
                                Cgs(m,nbs) = 2*pi*(fts(m,nbs)-fts(m-1,nbs))*mega/(thick*(k(m)-k(m-1)));
                            else
                                break
                            end
                        end
                        ai = aj;
                        err = abs(Lscis);
                        break
                    end
                    if cis(j) <= m1
                        P1 = Lsai*Lscis;
                        P2 = Lscis*Lsm1;
                        P3 = Lsm1*Lsm2;
                        P4 = Lsm2*Lsaj;
                        if P1 <= 0
                            ai   = ai;
                            aj   = cis(j);
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Lsai = Ls(k(m),ai);
                            Lsaj = Ls(k(m),aj);
                            Lsm1 = Ls(k(m),m1);
                            Lsm2 = Ls(k(m),m2);
                            Lsm0 = Ls(k(m),m0);
                        end
                        if P2 <= 0
                            ai   = cis(j);
                            aj   = m1;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Lsai = Ls(k(m),ai);
                            Lsaj = Ls(k(m),aj);
                            Lsm1 = Ls(k(m),m1);
                            Lsm2 = Ls(k(m),m2);
                            Lsm0 = Ls(k(m),m0);
                        end
                        if P3 <= 0
                            ai   = m1;
                            aj   = m2;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Lsai = Ls(k(m),ai);
                            Lsaj = Ls(k(m),aj);
                            Lsm1 = Ls(k(m),m1);
                            Lsm2 = Ls(k(m),m2);
                            Lsm0 = Ls(k(m),m0);
                        end
                        if P4 <= 0
                            ai   = m2;
                            aj   = aj;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Lsai = Ls(k(m),ai);
                            Lsaj = Ls(k(m),aj);
                            Lsm1 = Ls(k(m),m1);
                            Lsm2 = Ls(k(m),m2);
                            Lsm0 = Ls(k(m),m0);
                        end
                    end
                    if cis(j) > m1 && cis(j) <= m2
                        P1 = Lsai*Lsm1;
                        P2 = Lsm1*Lscis;
                        P3 = Lscis*Lsm2;
                        P4 = Lsm2*Lsaj;
                        if P1 <= 0
                            ai   = ai;
                            aj   = m1;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = ai+di/2;
                            Lsai = Ls(k(m),ai);
                            Lsaj = Ls(k(m),aj);
                            Lsm1 = Ls(k(m),m1);
                            Lsm2 = Ls(k(m),m2);
                            Lsm0 = Ls(k(m),m0);
                        end
                        if P2 <= 0
                            ai   = m1;
                            aj   = cis(j);
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = ai+di/2;
                            Lsai = Ls(k(m),ai);
                            Lsaj = Ls(k(m),aj);
                            Lsm1 = Ls(k(m),m1);
                            Lsm2 = Ls(k(m),m2);
                        end
                        if P3 <= 0
                            ai   = cis(j);
                            aj   = m2;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = ai+di/2;
                            Lsai = Ls(k(m),ai);
                            Lsaj = Ls(k(m),aj);
                            Lsm1 = Ls(k(m),m1);
                            Lsm2 = Ls(k(m),m2);
                        end
                        if P4 <= 0
                            ai   = m2;
                            aj   = aj;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = ai+di/2;
                            Lsai = Ls(k(m),ai);
                            Lsaj = Ls(k(m),aj);
                            Lsm1 = Ls(k(m),m1);
                            Lsm2 = Ls(k(m),m2);
                        end
                    end
                    if cis(j) > m2
                        P1 = Lsai*Lsm1;
                        P2 = Lsm1*Lsm2;
                        P3 = Lsm2*Lscis;
                        P4 = Lscis*Lsaj;
                        if P1 <= 0
                            ai   = ai;
                            aj   = m1;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = ai+di/2;
                            Lsai = Ls(k(m),ai);
                            Lsaj = Ls(k(m),aj);
                            Lsm1 = Ls(k(m),m1);
                            Lsm2 = Ls(k(m),m2);
                        end
                        if P2 <= 0
                            ai   = m1;
                            aj   = m2;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = ai+di/2;
                            Lsai = Ls(k(m),ai);
                            Lsaj = Ls(k(m),aj);
                            Lsm1 = Ls(k(m),m1);
                            Lsm2 = Ls(k(m),m2);
                        end
                        if P3 <= 0
                            ai   = m2;
                            aj   = cis(j);
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = ai+di/2;
                            Lsai = Ls(k(m),ai);
                            Lsaj = Ls(k(m),aj);
                            Lsm1 = Ls(k(m),m1);
                            Lsm2 = Ls(k(m),m2);
                        end
                        if P4 <= 0
                            ai   = cis(j);
                            aj   = aj;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = ai+di/2;
                            Lsai = Ls(k(m),ai);
                            Lsaj = Ls(k(m),aj);
                            Lsm1 = Ls(k(m),m1);
                            Lsm2 = Ls(k(m),m2);
                        end
                    end
                    err = abs(Lscis);
                end
            end
            ai = aj;
        end
    end
    VERS(m,1) = max(ERS(m,:));
    nbs       = 0;
end
%___________________________ Anti-symmetric modes _________________________
fta              = zeros(length(k),nom);
Cpa              = zeros(length(k),nom);
Wa               = zeros(length(k),nom);
Cga              = zeros(length(k),nom);
ERA              = zeros(length(k),nom);
VERA             = zeros(length(k),1);
IMA              = zeros(length(k),nom);
RATIOA           = zeros(length(k),nom); 
nbs              = 0;
j                = 0;
for m            = 1:length(k)
    for i        = 1:n+1
        ai       = fti+(i-1)*dft;
        while ai < ftf
            aj   = ai+dft;
            di   = aj-ai;
            m0   = ai+di/2;
            m1   = ai+di/3;
            m2   = ai+2*di/3;
            Laai = La(k(m),ai); %__ai = freq*thick, freq: frequency, thick: plate thickness, k(m): wavenumber k, La: anti-symmetric Lamb modes function.
            Laaj = La(k(m),aj);
            Lam0 = La(k(m),m0);
            Lam1 = La(k(m),m1);
            Lam2 = La(k(m),m2);
            if Laai == 0 && ai ~= fti
                nbs = nbs;
                break
            end
            if Laai == 0 && ai == fti
                nbs        = nbs+1;
                fta(m,nbs) = fti;
                Cpa(m,nbs) = (2*pi*fti*mega)/(k(m)*thick);
                Wa(m,nbs)  = Cpa(m,nbs)/(fti*mega/thick);
                if m >= 2
                    if fta(m,nbs)*fta(m-1,nbs) ~= 0
                        Cga(m,nbs) = 2*pi*(fta(m,nbs)-fta(m-1,nbs))*mega/(thick*(k(m)-k(m-1)));
                    else
                        break
                    end
                end
            end
            if Laaj == 0
                nbs        = nbs+1;
                fta(m,nbs) = aj;
                Cpa(m,nbs) = (2*pi*aj*mega)/(k(m)*thick);
                Wa(m,nbs)  = Cpa(m,nbs)/(aj*mega/thick);
                if m >= 2
                    if fta(m,nbs)*fta(m-1,nbs) ~= 0
                        Cga(m,nbs) = 2*pi*(fta(m,nbs)-fta(m-1,nbs))*mega/(thick*(k(m)-k(m-1)));
                    else
                        break
                    end
                end
            end
            if Lam0 == 0
                nbs        = nbs+1;
                fta(m,nbs) = m0;
                Cpa(m,nbs) = (2*pi*m0*mega)/(k(m)*thick);
                Wa(m,nbs)  = Cpa(m,nbs)/(m0*mega/thick);
                if m >= 2
                    if fta(m,nbs)*fta(m-1,nbs) ~= 0
                        Cga(m,nbs) = 2*pi*(fta(m,nbs)-fta(m-1,nbs))*mega/(thick*(k(m)-k(m-1)));
                    else
                        break
                    end
                end
            end
            if Lam1 == 0
                nbs        = nbs+1;
                fta(m,nbs) = m1;
                Cpa(m,nbs) = (2*pi*m1*mega)/(k(m)*thick);
                Wa(m,nbs)  = Cpa(m,nbs)/(m1*mega/thick);
                if m >= 2
                    if fta(m,nbs)*fta(m-1,nbs) ~= 0
                        Cga(m,nbs) = 2*pi*(fta(m,nbs)-fta(m-1,nbs))*mega/(thick*(k(m)-k(m-1)));
                    else
                        break
                    end
                end
            end
            if Lam2 == 0
                nbs        = nbs+1;
                fta(m,nbs) = m2;
                Cpa(m,nbs) = (2*pi*m2*mega)/(k(m)*thick);
                Wa(m,nbs)  = Cpa(m,nbs)/(m2*mega/thick);
                if m >= 2
                    if fta(m,nbs)*fta(m-1,nbs) ~= 0
                        Cga(m,nbs) = 2*pi*(fta(m,nbs)-fta(m-1,nbs))*mega/(thick*(k(m)-k(m-1)));
                    else
                        break
                    end
                end
            end
            if sign(Laai/Laaj) == 1
                break
            end
            if sign(Laai/Laaj) == -1
                err = max(max(max(abs(Laai),abs(Laaj)),max(abs(Lam1),abs(Lam2))),abs(Lam0));
                while err/epsilon > 1
                    di   = aj-ai;
                    m0   = ai+di/2;
                    m1   = ai+di/3;
                    m2   = ai+2*di/3;
                    Laai = La(k(m),ai);
                    Laaj = La(k(m),aj);
                    Lam0 = La(k(m),m0);
                    Lam1 = La(k(m),m1);
                    Lam2 = La(k(m),m2);
                    %______________________________________________________
                    lambda1 = (Lam1*ai)/(m1^2*(ai-m1));
                    lambda2 = (Laai*m1)/(ai^2*(ai-m1));
                    lambda3 = 1/ai+1/m1;
                    lambda4 = 1/ai^2+1/m1^2+1/(ai*m1);
                    T1      = m2*(m2^2/(ai*m1)-m2/m1+1);
                    T2      = (m2^3/ai)*(lambda4-1/(ai^2))-lambda4*m2^2+1;
                    T3      = (m2^3/ai)*(Laai/(ai^2)-(lambda1-lambda2))+(lambda1-lambda2)*m2^2;
                    Z1      = (aj^3)*((1/(ai*m1))*(lambda3-T2/T1)+1/ai^2)+(aj^2)*(lambda3*(T2/T1)-lambda4)-(T2/T1)*aj+1;
                    Z2      = (aj^3)*(-Laai/(ai^3)+(lambda1-lambda2)/ai+(1/(ai^2)-lambda3/ai)*((Lam2/T1)-(T3/T1)))+(aj^2)*(lambda3*((Lam2/T1)-(T3/T1))-(lambda1-lambda2))-((Lam2/T1)-(T3/T1))*aj;
                    cof4    = (Laaj+Z2)/Z1;
                    cof3    = Lam2/T1-T2*cof4/T1-T3/T1;
                    cof2    = lambda1-lambda2-lambda3*cof3-lambda4*cof4;
                    cof1    = Laai/(ai^3)-(cof2/ai)-cof3/(ai^2)-cof4/(ai^3);
                    %______________________________________________________
                    if ai == 0 || ai-m1 == 0 || T1 == 0 || Z1 == 0
                        ai = aj;
                        aj = ai+dft;
                        break
                    end
                    %_Calculation of the solutions of the cubic polynomial:
                    if cof1 ~= 0
                        A     = cof2/cof1;
                        B     = cof3/cof1;
                        C     = cof4/cof1;
                        p     = (B-A^2/3);
                        q     = (2*A^3)/27-A*B/3+C;
                        delta = sqrt(4*p^3+27*q^2);
                        dis   = 18*cof1*cof2*cof3*cof4-4*cof2^3*cof4+cof2^2*cof3^2-4*cof1*cof3^3-27*cof1^2*cof4^2;
                        if p == 0 && q == 0
                            s = -A/3;
                            if s > ai && s < aj
                                j      = j+1;
                                cis(j) = s;
                            end
                            if s < ai || s > aj
                                j      = j+1;
                                cis(j) = m0;
                            end
                        end
                        if p == 0 && q ~= 0
                            s = -A/3+(-q)^(1/3);
                            if s > ai && s < aj
                                j      = j+1;
                                cis(j) = s;
                            end
                            if s < ai || s > aj
                                j      = j+1;
                                cis(j) = m0;
                            end
                        end
                        if p > 0 && q == 0
                            s = -A/3;
                            if s > ai && s < aj
                                j      = j+1;
                                cis(j) = s;
                            end
                            if s < ai || s > aj
                                j      = j+1;
                                cis(j) = m0;
                            end
                        end
                        if p < 0 && q == 0
                            s1    = -A/3;
                            s2    = -A/3+sqrt(-p);
                            s3    = -A/3-sqrt(-p);
                            Las1  = La(k(m),s1);
                            Las2  = La(k(m),s2);
                            Las3  = La(k(m),s3);
                            mn    = min(s1,min(s2,s3));
                            mx    = max(s1,max(s2,s3));
                            Lamin = min(abs(Las1),min(abs(Las2),abs(Las3)));
                            if s1-s2 == 0 && s2-s3 == 0
                                s = s1;
                                if s > ai && s < aj
                                    j      = j+1;
                                    cis(j) = s;
                                end
                                if s < ai || s > aj
                                    j      = j+1;
                                    cis(j) = m0;
                                end
                            end
                            if (s1-s2)*(s1-s3)*(s2-s3) ~= 0
                                if (s1-mn)*(s1-mx) ~= 0
                                    mid = s1;
                                end
                                if (s2-mn)*(s2-mx) ~= 0
                                    mid = s2;
                                end
                                if (s3-mn)*(s3-mx) ~= 0
                                    mid = s3;
                                end
                                Lamid = La(k(m),mid);
                                if mx < ai || mn > aj
                                    j      = j+1;
                                    cis(j) = m0;
                                end
                                if mn > ai && mx < aj
                                    if Lamin == abs(Lamn)
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                    if Lamin == abs(Lamid)
                                        j      = j+1;
                                        cis(j) = mid;
                                    end
                                    if Lamin == abs(Lamx)
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                                if mn < ai && mx > aj
                                    if mid > ai && mid < aj
                                        j      = j+1;
                                        cis(j) = mid;
                                    end
                                    if mid < ai || mid > aj
                                        j      = j+1;
                                        cis(j) = m0;
                                    end
                                end
                                if mn > ai && mx > aj
                                    if mid < aj
                                        mn    = min(mn,mid);
                                        mx    = max(mn,mid);
                                        Lamn  = La(k(m),mn);
                                        Lamx  = La(k(m),mx);
                                        Lamin = min(abs(Lamn),abs(Lamx));
                                        if Lamin == abs(Lamn)
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                        if Lamin == abs(Lamid)
                                            j      = j+1;
                                            cis(j) = mid;
                                        end
                                    end
                                    if mid > aj
                                        if mn < aj
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                        if mn > aj
                                            j      = j+1;
                                            cis(j) = m0;
                                        end
                                    end
                                end
                                if mn < ai && mx < aj
                                    if mx < ai
                                        j      = j+1;
                                        cis(j) = m0;
                                    end
                                    if mx > ai
                                        if mid > ai
                                            mn    = min(mid,mx);
                                            mx    = max(mid,mx);
                                            Lamn  = La(k(m),mn);
                                            Lamx  = La(k(m),mx);
                                            Lamin = min(abs(Lamn),abs(Lamx));
                                            if Lamin == abs(Lamid)
                                                j      = j+1;
                                                cis(j) = mid;
                                            end
                                            if Lamin == abs(Lamx)
                                                j      = j+1;
                                                cis(j) = mx;
                                            end
                                        end
                                        if mid < ai
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                end
                            end
                            if s1-s2 == 0 && s1-s3 ~= 0
                                mn    = min(s1,s3);
                                mx    = max(s1,s3);
                                Lamn  = La(k(m),mn);
                                Lamx  = La(k(m),mx);
                                Lamin = min(abs(Lamn),abs(Lamx));
                                if mn > aj || mx < ai
                                    j      = j+1;
                                    cis(j) = m0;
                                end
                                if mn > ai && mx < aj
                                    if Lamin == abs(Lamn)
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                    if Lamin == abs(Lamx)
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                                if mn > ai && mn < aj
                                    if mx > aj
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                end
                                if mx > ai && mx < aj
                                    if mn < ai
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                            end
                            if s1-s3 == 0 && s1-s2 ~= 0
                                mn    = min(s1,s2);
                                mx    = max(s1,s2);
                                Lamn  = La(k(m),mn);
                                Lamx  = La(k(m),mx);
                                Lamin = min(abs(Lamn),abs(Lamx));
                                if mn > aj || mx < ai
                                    j      = j+1;
                                    cis(j) = m0;
                                end
                                if mn > ai && mx < aj
                                    if Lamin == abs(Lamn)
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                    if Lamin == abs(Lamx)
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                                if mn > ai && mn < aj
                                    if mx > aj
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                end
                                if mx > ai && mx < aj
                                    if mn < ai
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                            end
                            if s2-s3 == 0 && s2-s1 ~= 0
                                mn    = min(s1,s2);
                                mx    = max(s1,s2);
                                Lamn  = La(k(m),mn);
                                Lamx  = La(k(m),mx);
                                Lamin = min(abs(Lamn),abs(Lamx));
                                if mn > aj || mx < ai
                                    j      = j+1;
                                    cis(j) = m0;
                                end
                                if mn > ai && mx < aj
                                    if Lamin == abs(Lamn)
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                    if Lamin == abs(Lamx)
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                                if mn > ai && mn < aj
                                    if mx > aj
                                        j      = j+1;
                                        cis(j) = mn;
                                    end
                                end
                                if mx > ai && mx < aj
                                    if mn < ai
                                        j      = j+1;
                                        cis(j) = mx;
                                    end
                                end
                            end
                        end
                        if p ~= 0 && q ~= 0
                            z1 = -A/3+((sqrt(3)*delta-9*q)^(1/3))/(2^(1/3)*3^(2/3))-(2/3)^(1/3)*p/((sqrt(3)*delta-9*q)^(1/3));
                            z2 = -A/3+(1+1i*sqrt(3))*p/(2^(2/3)*3^(1/3)*((sqrt(3)*delta-9*q)^(1/3)))-(1-1i*sqrt(3))*((sqrt(3)*delta-9*q)^(1/3))/(2*2^(1/3)*3^(2/3));
                            z3 = -A/3+(1-1i*sqrt(3))*p/(2^(2/3)*3^(1/3)*((sqrt(3)*delta-9*q)^(1/3)))-(1+1i*sqrt(3))*((sqrt(3)*delta-9*q)^(1/3))/(2*2^(1/3)*3^(2/3));
                            if dis >= 0
                                s1    = real(z1);
                                s2    = real(z2);
                                s3    = real(z3);
                                Las1  = La(k(m),s1);
                                Las2  = La(k(m),s2);
                                Las3  = La(k(m),s3);
                                mn    = min(s1,min(s2,s3));
                                mx    = max(s1,max(s2,s3));
                                Lamin = min(abs(Las1),min(abs(Las2),abs(Las3)));
                                if s1-s2 == 0 && s2-s3 == 0
                                    s = s1;
                                    if s > ai && s < aj
                                        j      = j+1;
                                        cis(j) = s;
                                    end
                                    if s < ai || s > aj
                                        j      = j+1;
                                        cis(j) = m0;
                                    end
                                end
                                if (s1-s2)*(s1-s3)*(s2-s3) ~= 0
                                    if (s1-mn)*(s1-mx) ~= 0
                                        mid = s1;
                                    end
                                    if (s2-mn)*(s2-mx) ~= 0
                                        mid = s2;
                                    end
                                    if (s3-mn)*(s3-mx) ~= 0
                                        mid = s3;
                                    end
                                    Lamid = La(k(m),mid);
                                    if mx < ai || mn > aj
                                        j      = j+1;
                                        cis(j) = m0;
                                    end
                                    if mn > ai && mx < aj
                                        if Lamin == abs(Lamn)
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                        if Lamin == abs(Lamid)
                                            j = j+1;
                                            cis(j) = mid;
                                        end
                                        if Lamin == abs(Lamx)
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                    if mn < ai && mx > aj
                                        if mid > ai && mid < aj
                                            j      = j+1;
                                            cis(j) = mid;
                                        end
                                        if mid < ai || mid > aj
                                            j      = j+1;
                                            cis(j) = m0;
                                        end
                                    end
                                    if mn > ai && mx > aj
                                        if mid < aj
                                            mn    = min(mn,mid);
                                            mx    = max(mn,mid);
                                            Lamn  = La(k(m),mn);
                                            Lamx  = La(k(m),mx);
                                            Lamin = min(abs(Lamn),abs(Lamx));
                                            if Lamin == abs(Lamn)
                                                j      = j+1;
                                                cis(j) = mn;
                                            end
                                            if Lamin == abs(Lamid)
                                                j      = j+1;
                                                cis(j) = mid;
                                            end
                                        end
                                        if mid > aj
                                            if mn < aj
                                                j      = j+1;
                                                cis(j) = mn;
                                            end
                                            if mn > aj
                                                j      = j+1;
                                                cis(j) = m0;
                                            end
                                        end
                                    end
                                    if mn < ai && mx < aj
                                        if mx < ai
                                            j      = j+1;
                                            cis(j) = m0;
                                        end
                                        if mx > ai
                                            if mid > ai
                                                mn    = min(mid,mx);
                                                mx    = max(mid,mx);
                                                Lamn  = La(k(m),mn);
                                                Lamx  = La(k(m),mx);
                                                Lamin = min(abs(Lamn),abs(Lamx));
                                                if Lamin == abs(Lamid)
                                                    j      = j+1;
                                                    cis(j) = mid;
                                                end
                                                if Lamin == abs(Lamx)
                                                    j      = j+1;
                                                    cis(j) = mx;
                                                end
                                            end
                                            if mid < ai
                                                j      = j+1;
                                                cis(j) = mx;
                                            end
                                        end
                                    end
                                end
                                if s1-s2 == 0 && s1-s3 ~= 0
                                    mn    = min(s1,s3);
                                    mx    = max(s1,s3);
                                    Lamn  = La(k(m),mn);
                                    Lamx  = La(k(m),mx);
                                    Lamin = min(abs(Lamn),abs(Lamx));
                                    if mn > aj || mx < ai
                                        j = j+1;
                                        cis(j) = m0;
                                    end
                                    if mn > ai && mx < aj
                                        if Lamin == abs(Lamn)
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                        if Lamin == abs(Lamx)
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                    if mn > ai && mn < aj
                                        if mx > aj
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                    end
                                end
                                if s1-s3 == 0 && s1-s2 ~= 0
                                    mn    = min(s1,s2);
                                    mx    = max(s1,s2);
                                    Lamn  = La(k(m),mn);
                                    Lamx  = La(k(m),mx);
                                    Lamin = min(abs(Lamn),abs(Lamx));
                                    if mn > aj || mx < ai
                                        j      = j+1;
                                        cis(j) = m0;
                                    end
                                    if mn > ai && mx < aj
                                        if Lamin == abs(Lamn)
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                        if Lamin == abs(Lamx)
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                    if mn > ai && mn < aj
                                        if mx > aj
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                    end
                                    if mx > ai && mx < aj
                                        if mn < ai
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                end
                                if s2-s3 == 0 && s2-s1 ~= 0
                                    mn    = min(s1,s2);
                                    mx    = max(s1,s2);
                                    Lamn  = La(k(m),mn);
                                    Lamx  = La(k(m),mx);
                                    Lamin = min(abs(Lamn),abs(Lamx));
                                    if mn > aj || mx < ai
                                        j      = j+1;
                                        cis(j) = m0;
                                    end
                                    if mn > ai && mx < aj
                                        if Lamin == abs(Lamn)
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                        if Lamin == abs(Lamx)
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                    if mn > ai && mn < aj
                                        if mx > aj
                                            j      = j+1;
                                            cis(j) = mn;
                                        end
                                    end
                                    if mx > ai && mx < aj
                                        if mn < ai
                                            j      = j+1;
                                            cis(j) = mx;
                                        end
                                    end
                                end
                            end
                            if dis < 0
                                if z1 > ai && z1 < aj
                                    j      = j+1;
                                    cis(j) = z1;
                                end
                                if z1 < ai || z1 > aj
                                    j      = j+1;
                                    cis(j) = m0;
                                end
                            end
                        end
                    end
                    if cof1 == 0
                        if cof2 ~= 0
                            discr = cof3^2-4*cof2*cof4;
                            s1    = (-cof3-sqrt(discr))/(2*cof2);
                            s2    = (-cof3+sqrt(discr))/(2*cof2);
                            mn    = min(s1,s2);
                            mx    = max(s1,s2);
                            Lamn  = La(k(m),mn);
                            Lamx  = La(k(m),mx);
                            Lamin = min(abs(Lamn),abs(Lamx));
                            if mn > aj || mx < ai
                                j      = j+1;
                                cis(j) = m0;
                            end
                            if mn > ai && mx < aj
                                if Lamin == abs(Lamn)
                                    j      = j+1;
                                    cis(j) = mn;
                                end
                                if Lamin == abs(Lamx)
                                    j      = j+1;
                                    cis(j) = mx;
                                end
                            end
                            if mn > ai && mn < aj
                                if mx > aj
                                    j      = j+1;
                                    cis(j) = mn;
                                end
                            end
                            if mx > ai && mx < aj
                                if mn < ai
                                    j      = j+1;
                                    cis(j) = mx;
                                end
                            end
                        end
                        if cof2 == 0
                            s = -cof4/cof3;
                            if s > ai && s < aj
                                j      = j+1;
                                cis(j) = s;
                            end
                            if s < ai || s > aj
                                j      = j+1;
                                cis(j) = m0;
                            end
                        end
                    end
                    Lacis   = La(k(m),cis(j));
                    xi      = 1e-10;
                    LaCISxi = La(k(m),cis(j)+xi);
                    dLaCIS  = (LaCISxi-Lacis)/xi;
                    if abs(Lacis)/epsilon < 1
                        nbs           = nbs+1;
                        fta(m,nbs)    = cis(j);
                        ERA(m,nbs)    = abs(Lacis/dLaCIS);
                        Cpa(m,nbs)    = (2*pi*cis(j)*mega)/(k(m)*thick);
                        Wa(m,nbs)     = Cpa(m,nbs)/(cis(j)*mega/thick);
                        IMA(m,nbs)    = abs(Lacis);
                        RATIOA(m,nbs) = abs(Lacis/Lam0);
                        if m >= 2
                            if fta(m,nbs)*fta(m-1,nbs) ~= 0
                                Cga(m,nbs) = 2*pi*(fta(m,nbs)-fta(m-1,nbs))*mega/(thick*(k(m)-k(m-1)));
                            else
                                break
                            end
                        end
                        ai  = aj;
                        err = abs(Lacis);
                        break
                    end
                    if cis(j) <= m1
                        P1 = Laai*Lacis;
                        P2 = Lacis*Lam1;
                        P3 = Lam1*Lam2;
                        P4 = Lam2*Laaj;
                        if P1 <= 0
                            ai   = ai;
                            aj   = cis(j);
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = ai+di/2;
                            Laai = La(k(m),ai);
                            Laaj = La(k(m),aj);
                            Lam1 = La(k(m),m1);
                            Lam2 = La(k(m),m2);
                            Lam0 = La(k(m),m0);
                        end
                        if P2 <= 0
                            ai   = cis(j);
                            aj   = m1;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = ai+di/2;
                            Laai = La(k(m),ai);
                            Laaj = La(k(m),aj);
                            Lam1 = La(k(m),m1);
                            Lam2 = La(k(m),m2);
                            Lam0 = La(k(m),m0);
                        end
                        if P3 <= 0
                            ai   = m1;
                            aj   = m2;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Laai = La(k(m),ai);
                            Laaj = La(k(m),aj);
                            Lam1 = La(k(m),m1);
                            Lam2 = La(k(m),m2);
                            Lam0 = La(k(m),m0);
                        end
                        if P4 <= 0
                            ai   = m2;
                            aj   = aj;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Laai = La(k(m),ai);
                            Laaj = La(k(m),aj);
                            Lam1 = La(k(m),m1);
                            Lam2 = La(k(m),m2);
                            Lam0 = La(k(m),m0);
                        end
                    end
                    if cis(j) > m1 && cis(j) <= m2
                        P1 = Laai*Lam1;
                        P2 = Lam1*Lacis;
                        P3 = Lacis*Lam2;
                        P4 = Lam2*Laaj;
                        if P1 <= 0
                            ai   = ai;
                            aj   = m1;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Laai = La(k(m),ai);
                            Laaj = La(k(m),aj);
                            Lam1 = La(k(m),m1);
                            Lam2 = La(k(m),m2);
                            Lam0 = La(k(m),m0);
                        end
                        if P2 <= 0
                            ai   = m1;
                            aj   = cis(j);
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Laai = La(k(m),ai);
                            Laaj = La(k(m),aj);
                            Lam1 = La(k(m),m1);
                            Lam2 = La(k(m),m2);
                        end
                        if P3 <= 0
                            ai   = cis(j);
                            aj   = m2;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Laai = La(k(m),ai);
                            Laaj = La(k(m),aj);
                            Lam1 = La(k(m),m1);
                            Lam2 = La(k(m),m2);
                        end
                        if P4 <= 0
                            ai   = m2;
                            aj   = aj;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Laai = La(k(m),ai);
                            Laaj = La(k(m),aj);
                            Lam1 = La(k(m),m1);
                            Lam2 = La(k(m),m2);
                        end
                    end
                    if cis(j) > m2
                        P1 = Laai*Lam1;
                        P2 = Lam1*Lam2;
                        P3 = Lam2*Lacis;
                        P4 = Lacis*Laaj;
                        if P1 <= 0
                            ai   = ai;
                            aj   = m1;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Laai = La(k(m),ai);
                            Laaj = La(k(m),aj);
                            Lam1 = La(k(m),m1);
                            Lam2 = La(k(m),m2);
                        end
                        if P2 <= 0
                            ai   = m1;
                            aj   = m2;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Laai = La(k(m),ai);
                            Laaj = La(k(m),aj);
                            Lam1 = La(k(m),m1);
                            Lam2 = La(k(m),m2);
                        end
                        if P3 <= 0
                            ai   = m2;
                            aj   = cis(j);
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Laai = La(k(m),ai);
                            Laaj = La(k(m),aj);
                            Lam1 = La(k(m),m1);
                            Lam2 = La(k(m),m2);
                        end
                        if P4 <= 0
                            ai   = cis(j);
                            aj   = aj;
                            di   = aj-ai;
                            m1   = ai+di/3;
                            m2   = ai+2*di/3;
                            m0   = (ai+aj)/2;
                            Laai = La(k(m),ai);
                            Laaj = La(k(m),aj);
                            Lam1 = La(k(m),m1);
                            Lam2 = La(k(m),m2);
                        end
                    end
                    err = abs(Lacis);
                end
            end
            ai = aj;
        end
    end
    VERA(m,1) = max(ERA(m,:));
    nbs       = 0;
end
m               = 1:length(k);
r               = 1:nom;
xrs             = fts(m,r);
xra             = fta(m,r);
yrs             = Cgs(m,r);
yra             = Cga(m,r);
yprs            = Cps(m,r);
ypra            = Cpa(m,r);
ws              = Ws(m,r);
wa              = Wa(m,r);
vers            = VERS(m,1);
vera            = VERA(m,1);
mxxs            = max(VERS);
mnns            = min(VERS);
mxxa            = max(VERA);
mnna            = min(VERA);
mxx             = max(mxxs,mxxa);
xrs(xrs == 0)   = NaN;
xra(xra == 0)   = NaN;
yprs(yprs == 0) = NaN;
ypra(ypra == 0) = NaN;
ws(ws == 0)     = NaN;
wa(wa == 0)     = NaN;
yrs(yrs == 0)   = NaN;
yra(yra == 0)   = NaN;
ers             = ERS(m,r);
era             = ERA(m,r);
mxs             = max(ERS);
mns             = min(ERS);
mxa             = max(ERA);
mna             = min(ERA);
mx              = max(mxs,mxa);
ims             = IMS(m,r);
ima             = IMA(m,r);
ratios          = RATIOS(m,r);
ratioa          = RATIOA(m,r);
d               = thick/(2*pi*kilo);
mil             = 1e-3;
Cpi             = 0;
Cpf             = 12;
Cgi             = -2;
Cgf             = 6;
Wi              = 0;
Wf              = 10;
dF              = RunFT/50;
F               = 0:dF:RunFT;
j               = 1:length(F);
s               = 2;
%_________________________________Figures__________________________________
%**************************************************************************
figure(22)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(k(m)*mil,vera,'-b','linewidth',s);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
[~, objh] = legend({'$A$'}, 'Interpreter','latex', 'location', 'Northeast', 'Fontsize', 25);
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 26); 
objhl = findobj(objh, 'type', 'patch'); 
set(objhl, 'Markersize', 26); 
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Wavenumber [rad/mm]');
ylabel('Absolute maximal error','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman');
set(gca,'FontSize',20);
xlim([0 kf*mil]);
ylim([0 mxxa]);
%__________________________________________________________________________
figure(21)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(k(m)*mil,vers,'-r','linewidth',s);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
[~, objh] = legend({'$S$'}, 'Interpreter','latex', 'location', 'Northeast', 'Fontsize', 25);
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 26); 
objhl = findobj(objh, 'type', 'patch'); 
set(objhl, 'Markersize', 26); 
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Wavenumber [rad/mm]');
ylabel('Absolute maximal error','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman');
set(gca,'FontSize',20);
xlim([0 kf*mil]);
ylim([0 mxxs]);
%__________________________________________________________________________
figure(20)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot3(xra,k(m)*mil,ima,'^b','linewidth',s);
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
[~, objh] = legend({'$A$'}, 'Interpreter','latex', 'location', 'Northeast', 'Fontsize', 25);
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 26);
objhl = findobj(objh, 'type', 'patch'); 
set(objhl, 'Markersize', 26); 
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
g = xlabel('Frequency.thickness [MHz.mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [3.5, -2.5, 0]);
h = ylabel('Wavenumber [rad/mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [-2.0, 0.35, 0]);
set(g,'Rotation',+9.25);
set(h,'Rotation',-16.25);
zlabel('Absolute value of $La(\mathrm{cis})$','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([0 kf*mil]);
%__________________________________________________________________________
figure(19)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot3(xrs,k(m)*mil,ims,'^r','linewidth',s);
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
[~, objh] = legend({'$S$'}, 'Interpreter','latex', 'location', 'Northeast', 'Fontsize', 25);
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 26); 
objhl = findobj(objh, 'type', 'patch');
set(objhl, 'Markersize', 26); 
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
g = xlabel('Frequency.thickness [MHz.mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [3.5, -2.5, 0]);
h = ylabel('Wavenumber [rad/mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [-2.0, 0.35, 0]);
set(g,'Rotation',+9.25);
set(h,'Rotation',-16.25);
zlabel('Absolute value of $Ls(\mathrm{cis})$','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([0 kf*mil]);
%__________________________________________________________________________
figure(18)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot3(xra,k(m)*mil,ratioa,'^b','linewidth',s);
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
[~, objh] = legend({'$A$'}, 'Interpreter','latex', 'location', 'Northeast', 'Fontsize', 25);
objhl = findobj(objh, 'type', 'line');
set(objhl, 'Markersize', 26); 
objhl = findobj(objh, 'type', 'patch'); 
set(objhl, 'Markersize', 26); 
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
g = xlabel('Frequency.thickness [MHz.mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [3.5, -2.5, 0]);
h = ylabel('Wavenumber [rad/mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [-2.0, 0.35, 0]);
set(g,'Rotation',+9.25);
set(h,'Rotation',-16.25);
zlabel('$\left\vert La(\mathrm{cis})/La(\mathrm{m0})\right\vert$','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([0 kf*mil]);
%__________________________________________________________________________
figure(17)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot3(xrs,k(m)*mil,ratios,'^r','linewidth',s);
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
[~, objh] = legend({'$S$', '$A$'}, 'Interpreter','latex', 'location', 'Northeast', 'Fontsize', 25);
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 26); 
objhl = findobj(objh, 'type', 'patch'); 
set(objhl, 'Markersize', 26); 
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
g = xlabel('Frequency.thickness [MHz.mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [3.5, -2.5, 0]);
h = ylabel('Wavenumber [rad/mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [-2.0, 0.35, 0]);
set(g,'Rotation',+9.25);
set(h,'Rotation',-16.25);
zlabel('\textbf{$\left\vert Ls(\mathrm{cis})/Ls(\mathrm{m0})\right\vert$}','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([0 kf*mil]);
%__________________________________________________________________________
figure(16)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot3(xrs,k(m)*mil,ims,'^r','linewidth',s);
grid on
hold on
plot3(xra,k(m)*mil,ima,'^b','linewidth',s);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
[~, objh] = legend({'$S$', '$A$'}, 'Interpreter','latex', 'location', 'Northeast', 'Fontsize', 25);
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 26); 
objhl = findobj(objh, 'type', 'patch'); 
set(objhl, 'Markersize', 26); 
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
g = xlabel('Frequency.thickness [MHz.mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [3.5, -2.5, 0]);
h = ylabel('Wavenumber [rad/mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [-2.0, 0.35, 0]);
set(g,'Rotation',+9.25);
set(h,'Rotation',-16.25);
zlabel('Absolute values of $Ls(\mathrm{cis})$ , $La(\mathrm{cis})$','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([0 kf*mil]);
%__________________________________________________________________________
figure(15)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot3(xrs,k(m)*mil,ratios,'^r','linewidth',s);
grid on
hold on
plot3(xra,k(m)*mil,ratioa,'^b','linewidth',s);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
[~, objh] = legend({'$S$', '$A$'}, 'Interpreter','latex', 'location', 'Northeast', 'Fontsize', 25);
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 26); 
objhl = findobj(objh, 'type', 'patch'); 
set(objhl, 'Markersize', 26); 
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
g = xlabel('Frequency.thickness [MHz.mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [3.5, -2.5, 0]);
h = ylabel('Wavenumber [rad/mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [-2.0, 0.35, 0]);
set(g,'Rotation',+9.25);
set(h,'Rotation',-16.25);
zlabel('\textbf{$\left\vert Ls(\mathrm{cis})/Ls(\mathrm{m0})\right\vert$ , $\left\vert La(\mathrm{cis})/La(\mathrm{m0})\right\vert$}','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([0 kf*mil]);
%__________________________Absolute maximal error__________________________
figure(14)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(k(m)*mil,vers,'-r','linewidth',s);
hold on
plot(k(m)*mil,vera,'-b','linewidth',s);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
[~, objh] = legend({'$S$', '$A$'}, 'Interpreter','latex', 'location', 'Northeast', 'Fontsize', 25);
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 26); 
objhl = findobj(objh, 'type', 'patch'); 
set(objhl, 'Markersize', 26); 
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Wavenumber [rad/mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman');
ylabel('Absolute maximal error','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman');
set(gca,'FontSize',20);
xlim([0 kf*mil]);
ylim([0 mxx]);
%______________________________Absolute error______________________________
figure(13)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot3(xrs,k(m)*mil,ers,'^r','linewidth',s);
grid on
hold on
plot3(xra,k(m)*mil,era,'^b','linewidth',s);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
[~, objh] = legend({'$S$', '$A$'}, 'Interpreter','latex', 'location', 'Northeast', 'Fontsize', 25);
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 26); 
objhl = findobj(objh, 'type', 'patch'); 
set(objhl, 'Markersize', 26); 
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
g = xlabel('Frequency.thickness [MHz.mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [3.5, -2.5, 0]);
h = ylabel('Wavenumber [rad/mm]','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman','Position', [-2.0, 0.35, 0]);
set(g,'Rotation',+9.25);
set(h,'Rotation',-16.25);
zlabel('Absolute error','interpreter','latex','FontWeight','bold','FontSize',23,'FontName','Times New Roman');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([0 kf*mil]);
%________________________Anti-symmetric wavelength_________________________
figure(12)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(xra,kilo*wa,'.b','linewidth',s);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
set(0,'DefaultTextFontName', 'Times New Roman');
s2 = 18;
text(0.47,9,'A0','Color','blue','FontSize',s2);
text(2.11,9,'A1','Color','blue','FontSize',s2);
text(4.79,9,'A2','Color','blue','FontSize',s2);
text(5.89,9,'A3','Color','blue','FontSize',s2);
text(7.5,9,'A4','Color','blue','FontSize',s2);
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Frequency.thickness [MHz.mm]');
ylabel('Wavelength [mm]');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([Wi Wf]);
%___________________________Symmetric wavelength___________________________
figure(11)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(xrs,kilo*ws,'.r','linewidth',s);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
set(0,'DefaultTextFontName', 'Times New Roman');
s2 = 18;
text(1.27,9,'S0','Color','red','FontSize',s2);
text(2.88,9,'S1','Color','red','FontSize',s2);
text(3.72,9,'S2','Color','red','FontSize',s2);
text(6.4,9,'S3','Color','red','FontSize',s2);
text(8.83,9,'S4','Color','red','FontSize',s2);
text(9.6,9,'S5','Color','red','FontSize',s2);
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Frequency.thickness [MHz.mm]');
ylabel('Wavelength [mm]');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([Wi Wf]);
%_____________________Anti-symmetric group velocity________________________
figure(10)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(xra,mil*yra,'.b','linewidth',s);
hold on
plot(F(j),0*F(j),'--black','linewidth',0.25);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
set(0,'DefaultTextFontName', 'Times New Roman');
s2 = 18;
text(0.9,3.45,'A0','Color','blue','FontSize',s2);
text(2.54,4,'A1','Color','blue','FontSize',s2);
text(6.14,3.8,'A2','Color','blue','FontSize',s2);
text(5.75,0.5,'A3','Color','blue','FontSize',s2);
text(7.5,0.5,'A4','Color','blue','FontSize',s2);
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Frequency.thickness [MHz.mm]');
ylabel('Group velocity [m/ms]');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([Cgi Cgf]);
%________________________Symmetric group velocity__________________________
figure(9)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(xrs,mil*yrs,'.r','linewidth',s);
hold on
plot(F(j),0*F(j),'--black','linewidth',0.25);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
set(0,'DefaultTextFontName', 'Times New Roman');
s2 = 18;
text(0.9,5.45,'S0','Color','red','FontSize',s2);
text(3.93,4.95,'S1','Color','red','FontSize',s2);
text(6.93,4.85,'S2','Color','red','FontSize',s2);
text(8.96,4.95,'S3','Color','red','FontSize',s2);
text(8.88,0.5,'S4','Color','red','FontSize',s2);
text(9.53,0.5,'S5','Color','red','FontSize',s2);
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Frequency.thickness [MHz.mm]');
ylabel('Group velocity [m/ms]');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([Cgi Cgf]);
%_____________________Anti-symmetric phase velocity________________________
figure(8)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(xra,mil*ypra,'.b','linewidth',s);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
set(0,'DefaultTextFontName', 'Times New Roman');
s2 = 18;
text(0.6,2.6,'A0','Color','blue','FontSize',s2);
text(1.5,10.5,'A1','Color','blue','FontSize',s2);
text(4.45,10.5,'A2','Color','blue','FontSize',s2);
text(6.8,10.5,'A3','Color','blue','FontSize',s2);
text(8.2,10.5,'A4','Color','blue','FontSize',s2);
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Frequency.thickness [MHz.mm]');
ylabel('Phase velocity [m/ms]');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([Cpi Cpf]);
%________________________Symmetric phase velocity__________________________
figure(7)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(xrs,mil*yprs,'.r','linewidth',s);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
set(0,'DefaultTextFontName', 'Times New Roman');
s2 = 18;
text(0.6,5.78,'S0','Color','red','FontSize',s2);
text(2.5,10.5,'S1','Color','red','FontSize',s2);
text(3.82,10.5,'S2','Color','red','FontSize',s2);
text(6.34,10.5,'S3','Color','red','FontSize',s2);
text(9.25,10.5,'S4','Color','red','FontSize',s2);
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Frequency.thickness [MHz.mm]');
ylabel('Phase velocity [rad/mm]');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([Cpi Cpf]);
%________________________Anti-symmetric wavenumber_________________________
figure(6)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(xra,k(m)*mil,'.b','linewidth',s);
hold on
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
set(0,'DefaultTextFontName', 'Times New Roman');
s2 = 18;
text(0.38,1.25,'A0','Color','blue','FontSize',s2);
text(2,0.5,'A1','Color','blue','FontSize',s2);
text(4.38,0.5,'A2','Color','blue','FontSize',s2);
text(5.8,0.5,'A3','Color','blue','FontSize',s2);
text(7.5,0.5,'A4','Color','blue','FontSize',s2);
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Frequency.thickness [MHz.mm]');
ylabel('Wavenumber [rad/mm]');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([0 kf*mil]);
%___________________________Symmetric wavenumber___________________________
figure(5)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(xrs,k(m)*mil,'.r','linewidth',s);
hold on
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
set(0,'DefaultTextFontName', 'Times New Roman');
s2 = 18;
text(1.38,1.25,'S0','Color','red','FontSize',s2);
text(2.55,0.5,'S1','Color','red','FontSize',s2);
text(3.67,0.5,'S2','Color','red','FontSize',s2);
text(6.35,0.5,'S3','Color','red','FontSize',s2);
text(8.85,0.5,'S4','Color','red','FontSize',s2);
text(9.55,0.5,'S5','Color','red','FontSize',s2);
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Frequency.thickness [MHz.mm]');
ylabel('Wavenumber [rad/mm]');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([0 kf*mil]);
%________________________________Wavelength________________________________
figure(4)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(xrs,kilo*ws,'.r','linewidth',s);
hold on
plot(xra,kilo*wa,'.b','linewidth',s);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
set(0,'DefaultTextFontName', 'Times New Roman');
s2 = 18;
text(0.47,9,'A0','Color','blue','FontSize',s2);
text(2.11,9,'A1','Color','blue','FontSize',s2);
text(4.79,9,'A2','Color','blue','FontSize',s2);
text(5.89,9,'A3','Color','blue','FontSize',s2);
text(7.5,9,'A4','Color','blue','FontSize',s2);
text(1.27,9,'S0','Color','red','FontSize',s2);
text(2.88,9,'S1','Color','red','FontSize',s2);
text(3.72,9,'S2','Color','red','FontSize',s2);
text(6.4,9,'S3','Color','red','FontSize',s2);
text(8.83,9,'S4','Color','red','FontSize',s2);
text(9.6,9,'S5','Color','red','FontSize',s2);
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Frequency.thickness [MHz.mm]');
ylabel('Wavelength [mm]');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([Wi Wf]);
%______________________________Group velocity______________________________
figure(3)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(xrs,mil*yrs,'.r','linewidth',2);
hold on
plot(xra,mil*yra,'.b','linewidth',2);
hold on
plot(F(j),0*F(j),'--black','linewidth',0.25);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
set(0,'DefaultTextFontName', 'Times New Roman');
s2 = 18;
text(0.9,3.45,'A0','Color','blue','FontSize',s2);
text(2.54,4,'A1','Color','blue','FontSize',s2);
text(6.14,3.8,'A2','Color','blue','FontSize',s2);
text(5.75,0.5,'A3','Color','blue','FontSize',s2);
text(7.5,0.5,'A4','Color','blue','FontSize',s2);
text(0.9,5.45,'S0','Color','red','FontSize',s2);
text(3.93,4.95,'S1','Color','red','FontSize',s2);
text(6.93,4.85,'S2','Color','red','FontSize',s2);
text(8.96,4.95,'S3','Color','red','FontSize',s2);
text(8.88,0.5,'S4','Color','red','FontSize',s2);
text(9.53,0.5,'S5','Color','red','FontSize',s2);
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Frequency.thickness [MHz.mm]');
ylabel('Group velocity [m/ms]');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([Cgi Cgf]);
%______________________________Phase velocity______________________________
figure(2)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(xrs,mil*yprs,'.r','linewidth',2);
hold on
plot(xra,mil*ypra,'.b','linewidth',2);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
set(0,'DefaultTextFontName', 'Times New Roman');
s2 = 18;
text(0.6,2.6,'A0','Color','blue','FontSize',s2);
text(1.5,10.5,'A1','Color','blue','FontSize',s2);
text(4.45,10.5,'A2','Color','blue','FontSize',s2);
text(6.8,10.5,'A3','Color','blue','FontSize',s2);
text(8.2,10.5,'A4','Color','blue','FontSize',s2);
text(0.6,5.78,'S0','Color','red','FontSize',s2);
text(2.5,10.5,'S1','Color','red','FontSize',s2);
text(3.82,10.5,'S2','Color','red','FontSize',s2);
text(6.34,10.5,'S3','Color','red','FontSize',s2);
text(9.25,10.5,'S4','Color','red','FontSize',s2);
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Frequency.thickness [MHz.mm]');
ylabel('Phase velocity [m/ms]');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([Cpi Cpf]);
%_________________________________Wavenumber_______________________________
figure(1)
set(gcf, 'WindowStyle', 'normal', 'WindowState', 'maximized');
plot(xrs,k(m)*mil,'.r','linewidth',2);
hold on
plot(xra,k(m)*mil,'.b','linewidth',2);
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
s1 = 25;
set(findall(gcf,'type','axes'),'fontsize',s1);
set(0,'DefaultTextFontName', 'Times New Roman');
s2 = 18;
text(0.38,1.25,'A0','Color','blue','FontSize',s2);
text(2,0.5,'A1','Color','blue','FontSize',s2);
text(4.38,0.5,'A2','Color','blue','FontSize',s2);
text(5.8,0.5,'A3','Color','blue','FontSize',s2);
text(7.5,0.5,'A4','Color','blue','FontSize',s2);
text(1.38,1.25,'S0','Color','red','FontSize',s2);
text(2.55,0.5,'S1','Color','red','FontSize',s2);
text(3.67,0.5,'S2','Color','red','FontSize',s2);
text(6.35,0.5,'S3','Color','red','FontSize',s2);
text(8.85,0.5,'S4','Color','red','FontSize',s2);
text(9.55,0.5,'S5','Color','red','FontSize',s2);
set(gca, 'color', [1 1 1])
set(gcf, 'color', [1 1 1])
xlabel('Frequency.thickness [MHz.mm]');
ylabel('Wavenumber [rad/mm]');
set(gca,'FontSize',20);
xlim([0 RunFT]);
ylim([0 kf*mil]);
toc

