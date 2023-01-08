 sigma_profile = zeros(N,N_comp);

for i = 1:N_comp
    sigma_profile(:,i) = sqrt(s1(i)^2.*sin(phi-phi0(i)).^2+s2(i).^2.*cos(phi-phi0(i)).^2);

    pd(:,i) = sum(pd(:,i))/N*1./sqrt(2.*pi.*sigma_profile(:,i).^2).*exp(- ( r- (-mue(1,i)*sin(phi) + mue(2,i)*cos(phi)) ).^2./2./sigma_profile(:,i).^2 );
%     pd(:,i) = 1./sqrt(2.*pi.*sigma_profile(:,i).^2).*exp(- ( r- (-mue(1,i)*sin(phi) + mue(2,i)*cos(phi)) ).^2./2./sigma_profile(:,i).^2 );
end
p0 = sum(pd,2);
pd = pd./p0;

for i = 1:N_comp
    A = [-sum(pd(:,i).*sin(phi).^2)             sum(pd(:,i).*sin(phi).*cos(phi));
         -sum(pd(:,i).*sin(phi).*cos(phi))      sum(pd(:,i).*cos(phi).^2)];
    b = [sum(pd(:,i).*r.*sin(phi)); 
         sum(pd(:,i).*r.*cos(phi))];
     
    sol = A\b;
    mue(1,i) = sol(1);
    mue(2,i) = sol(2);
end
Rc = sqrt(mue(1,:).^2 + mue(2,:).^2);
phic = atan2(mue(2,:),mue(1,:));    