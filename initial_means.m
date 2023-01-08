%% initial step
init_ind = zeros(N,N_comp,'logical');
mue = zeros(2,N_comp);
for i = 1:N_comp
    init_ind(i:N_comp:N,i) = 1; % random 

    A = [-sum(sin(phi(init_ind(:,i))).^2)                        sum(sin(phi(init_ind(:,i))).*cos(phi(init_ind(:,i))));
        -sum(sin(phi(init_ind(:,i))).*cos(phi(init_ind(:,i))))   sum(cos(phi(init_ind(:,i))).^2)];
    b = [sum(r(init_ind(:,i)).*sin(phi(init_ind(:,i)))); 
         sum(r(init_ind(:,i)).*cos(phi(init_ind(:,i))))];
    
    sol = A\b;
    mue(1,i) = sol(1);
    mue(2,i) = sol(2);
end
Rc = sqrt(mue(1,:).^2 + mue(2,:).^2);
phic = atan2(mue(2,:),mue(1,:));

d = zeros(N,N_comp);
dc = zeros(N,N_comp);


for iter = 1:N_iter_init_means
    for i = 1:N_comp
    d(:,i)  = -mue(1,i)*sin(phi) + mue(2,i)*cos(phi);
    dc(:,i) = r-d(:,i); 
    end
    [~,ind_min]=min(abs(dc),[],2);
    for i = 1:N_comp
        init_ind(:,i) = (ind_min==i);

        A = [-sum(sin(phi(init_ind(:,i))).^2)                        sum(sin(phi(init_ind(:,i))).*cos(phi(init_ind(:,i))));
            -sum(sin(phi(init_ind(:,i))).*cos(phi(init_ind(:,i))))   sum(cos(phi(init_ind(:,i))).^2)];
        b = [sum(r(init_ind(:,i)).*sin(phi(init_ind(:,i)))); 
             sum(r(init_ind(:,i)).*cos(phi(init_ind(:,i))))];
        
        sol = A\b;
        mue(1,i) = sol(1);
        mue(2,i) = sol(2);
    end
    Rc = sqrt(mue(1,:).^2 + mue(2,:).^2);
    phic = atan2(mue(2,:),mue(1,:));    
end
mueo = mue;


