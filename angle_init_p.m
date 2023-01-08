function [phi0] = angle_init_p(phii,rci,s1,s2,pi1)
    ai = 2*phii;
%     pi1 = ones(size(pi1));
%     pi1(pi1<0.99)=0;
    M  = pi1.*(s2^2-s1^2);
    Ni = pi1.*(s1^2+s2^2 - 2*rci.^2);
    %%
    
    As2 = sum(M.*sin(ai).*cos(ai));
%     Ac2 = -As2;
    Asc = sum(M.*cos(ai).^2) - sum(M.*sin(ai).^2);
    As  = sum(Ni.*cos(ai));
    Ac  = -sum(Ni.*sin(ai));

    %%
    
    plnm = [(4*As2^2 + Asc^2) (-4*As2*Ac + 2*Asc*As)  (Ac^2 + As^2 - Asc^2 - 4*As2^2) (2*As2*Ac - 2*Asc*As)  (As2^2 - As^2)];
    
    xr  = roots(plnm);
%     xr  = double(roots_4(vpa(plnm)));
    yr1 = sqrt(1-xr.^2);
    yr2 = -sqrt(1-xr.^2);
    
    ind = abs(imag(xr))<1e-6;
    xr = real(xr(ind));
    yr1 = real(yr1(ind));
    yr2 = real(yr2(ind));
    
    ind1 = abs(As2*yr1.^2 - As2*xr.^2 + Asc*xr.*yr1 + As*yr1 + Ac*xr)<1e-6;
    ind2 = abs(As2*yr2.^2 - As2*xr.^2 + Asc*xr.*yr2 + As*yr2 + Ac*xr)<1e-6 ;
    
    xr1 = xr(ind1);
    yr1 = yr1(ind1);
    
    xr2 = xr(ind2);
    yr2 = yr2(ind2);
    
    %%
    
    a1 = atan2(yr1,xr1);
    a2 = atan2(yr2,xr2);
    
    phi1 = a1/2;
    phi2 = a2/2;
    
    phi00s =[phi1;phi2];
    
    s0 = zeros(size(phi00s));
    for i = 1:length(phi00s)
        phi00 = phi00s(i);
        s0(i) = sum( pi1.*(s1.^2*sin(phi00-phii).^2 + s2.^2.*cos(phi00-phii).^2 - rci.^2).^2);
    end
    [~,ind] = min(s0);
    phi0 = phi00s(ind);

end