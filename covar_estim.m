d = zeros(N,N_comp);
rc = zeros(N,N_comp);
for i = 1:N_comp
d(:,i)  = -mue(1,i)*sin(phi) + mue(2,i)*cos(phi);
rc(:,i) = r-d(:,i); 
end

%%


for i = 1:N_comp
    % Calculate the moments
    M2 = sum(pd(:,i).*(rc(:,i)).^2)/(sum(pd(:,i)));
    M4 = sum(pd(:,i).*(rc(:,i)).^4)/(sum(pd(:,i)));
    
    % solving system
    
    if M4/3-M2^2 <0
        disp('error type1')
    end
    if M2 + sqrt(2)*sqrt(M4/3-M2^2) < 0
        disp('error type2')
    end

    y(1)=sqrt(abs(M2 + sqrt(2)*sqrt(abs(M4/3-M2^2))));
    y(2)=sqrt(abs(M2 - sqrt(2)*sqrt(abs(M4/3-M2^2))));

    s1(i) = y(1);
    s2(i) = y(2);
    
    
    % Profile of variance changes with rotation with f(phi)=  sqrt(a^2.*sin(phi-phishi).^2+b.^2*cos(phi-phishi).^2)))
    % Find the best phishi for our data.

     
    phi0(i)=angle_init_p(phi,rc(:,i),s1(i),s2(i),pd(:,i));

    % 
    if kor==1
        for iter_cor = 1:N_iter_cor
            A = [sum( pd(:,i).*sin(phi-phi0(i)).^4 )                                sum(pd(:,i).*sin(phi-phi0(i)).^2.*cos(phi-phi0(i)).^2);
                sum(pd(:,i).*sin(phi-phi0(i)).^2.*cos(phi-phi0(i)).^2)             sum(pd(:,i).*cos(phi-phi0(i)).^4)];
            b = [sum( pd(:,i).*rc(:,i).^2.*sin(phi-phi0(i)).^2) ; 
                 sum( pd(:,i).*rc(:,i).^2.*cos(phi-phi0(i)).^2) ];
            
            y = sqrt(abs(A\b));
        
            if y(1)<y(2)
                s1(i) = y(2);
                s2(i) = y(1);
            else    
                s1(i) = y(1);
                s2(i) = y(2);
            end
                
            phi0(i)=angle_init_p(phi,rc(:,i),s1(i),s2(i),pd(:,i));
        end
    end
    
    
    
    % eigendecomposition and covar matrix
    V = [cos(phi0(i)) cos(phi0(i)+pi/2) ; 
         sin(phi0(i)) sin(phi0(i)+pi/2)  ];
    D = diag([s1(i) s2(i)]);
    
    Sigmae(:,:,i)=V*D.^2*V'; % This is covar
end