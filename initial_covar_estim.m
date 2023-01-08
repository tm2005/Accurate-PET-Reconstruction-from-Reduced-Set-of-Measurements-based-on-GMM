Sigmae = zeros(2,2,N_comp);
% centering
d = zeros(N,N_comp);
rc = zeros(N,N_comp);
for i = 1:N_comp
d(:,i)  = -mue(1,i)*sin(phi) + mue(2,i)*cos(phi);
rc(:,i) = r-d(:,i); 
end

%%

for i = 1:N_comp
    % Calculate the moments
    M2 = sum(rc(init_ind(:,i),i).^2)/(sum(init_ind(:,i)));
    M4 = sum(rc(init_ind(:,i),i).^4)/(sum(init_ind(:,i)));
    
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

     
    phi0(i)=angle_init(phi(init_ind(:,i)),rc(init_ind(:,i),i),s1(i),s2(i));

    % 
    if kor==1
        for iter_cor = 1:N_iter_cor
            A = [sum( sin(phi(init_ind(:,i))-phi0(i)).^4)                                sum(sin(phi(init_ind(:,i))-phi0(i)).^2.*cos(phi(init_ind(:,i))-phi0(i)).^2);
                sum(sin(phi(init_ind(:,i))-phi0(i)).^2.*cos(phi(init_ind(:,i))-phi0(i)).^2)     sum(cos(phi(init_ind(:,i))-phi0(i)).^4)];
            b = [sum( rc(init_ind(:,i),i).^2.*sin(phi(init_ind(:,i))-phi0(i)).^2) ; sum(rc(init_ind(:,i),i).^2.*cos(phi(init_ind(:,i))-phi0(i)).^2) ];
            
            y = sqrt(abs(A\b));
        
            if y(1)<y(2) 
                s1(i) = y(2);
                s2(i) = y(1);
            else    
                s1(i) = y(1);
                s2(i) = y(2);
            end
        

%             phi0(i)=angle_init(phi(init_ind(:,i)),rc(init_ind(:,i),i),s1(i),s2(i));
        end
    end
    
    
    
    % eigendecomposition and covar matrix
    V = [cos(phi0(i)) cos(phi0(i)+pi/2) ; 
         sin(phi0(i)) sin(phi0(i)+pi/2)  ];
    D = diag([s1(i) s2(i)]);
    
    Sigmae(:,:,i)=V*D.^2*V'; % This is covar
end

Sigmaeo = Sigmae;