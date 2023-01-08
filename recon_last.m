close all
clear all
clc

for nsim=1:1
    nsim
simulation_many_Gaussians
[Xg,Yg] = meshgrid(-Rs:Rs/360:Rs,-Rs:Rs/360:Rs);
dx = Rs/360;

Z = zeros(size(Xg(:)));
for i = 1:N_comp
    Z = Z + Ns(i)/N*mvnpdf([Xg(:) Yg(:)],mus(:,i)',Sigmas(:,:,i));
end
Z = reshape(Z, size(Xg));
figure(10), surfl(Xg, Yg, Z), shading interp, colormap copper, title('Ground truth');
Pe = Z;
% N_comp = 2;
%%
N_iter_init_means = 500;

kor = 1;
N_iter_cor = 1 ;

initial_means
initial_covar_estim


Z = zeros(size(Xg(:)));
for i = 1:N_comp
    Z = Z + sum(init_ind(:,i))/N*mvnpdf([Xg(:) Yg(:)],mue(:,i)',Sigmae(:,:,i));
end
Z = reshape(Z, size(Xg));
% figure, contour(Xg, Yg, Z,10), shading interp, colormap hot, title('Apr'), ax = axis; hold on,...
%      scatter(X(1:1:end,1), X(1:1:end,2), 14, [0.5 0.5 0.5], '+','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1);
%%
pd = double(init_ind);
for i=1:N_comp
    Ne(i) = sum(pd(:,i));
end
Neo = Ne;
%%
N_iter_init_means = 200;
N_erroro = inf;
for count = 1:100
% count
    assign_by_prob_and_update_means
    covar_estim
    % 
    % Z = zeros(size(Xg(:)));
    % for i = 1:N_comp
    %     Z = Z + sum(init_ind(:,i))/N*mvnpdf([Xg(:) Yg(:)],mue(:,i)',Sigmae(:,:,i));
    % end
    % Z = reshape(Z, size(Xg));
    % figure(333), contour(Xg, Yg, Z,10), shading interp, colormap hot, title('Apr'), ax = axis; hold on,...
    %      scatter(X(1:1:end,1), X(1:1:end,2), 14, [0.5 0.5 0.5], '+','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1),hold off;
    % count
    
    for i=1:N_comp
        Ne(i) = sum(pd(:,i));
    end
    
    N_error = norm(Neo-Ne,1);
    
    abs(Neo-Ne);
    if (N_error > N_erroro ) || abs(N_error) <25
%         mueo=mue;
%         Sigmaeo = Sigmae;
%         Neo = Ne;
%         N_erroro = N_error;
        break;
    else
        mueo=mue;
        Sigmaeo = Sigmae;
        Neo = Ne;
        N_erroro = N_error;
    end

end

% mue = mueo;
% Sigmae = Sigmaeo;
% Neo = Ne;

% errori
% errori_idejni
% %%
% Z = zeros(size(Xg(:)));
% for i = 1:N_comp
%     Z = Z + Ne(i)/N*mvnpdf([Xg(:) Yg(:)],mue(:,i)',Sigmae(:,:,i));
% end
% Z = reshape(Z, size(Xg));
% P2 = Z;

% figure(109), contour(Xg, Yg, Z,10), shading interp, colormap hot, title('Apr'), ax = axis; hold on,...
%      scatter(X(1:1:end,1), X(1:1:end,2), 14, [0.5 0.5 0.5], '+','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1),hold off;


% %%
Z = zeros(size(Xg(:)));
for i = 1:N_comp
    Z = Z + Ne(i)/N*mvnpdf([Xg(:) Yg(:)],mue(:,i)',Sigmae(:,:,i));
end
Z = reshape(Z, size(Xg));
figure(11), surfl(Xg, Yg, Z), shading interp, colormap copper, title('Estim');
Z = reshape(Z, size(Xg));
% figure(18), imagesc(Z), shading interp, colormap copper, title('Estim');
end
%%
% Ne
% figure, stem(m1error)
% figure, stem(m2error)
% figure, stem(m3error)
% 
% figure, stem(s1error)
% figure, stem(s2error)
% figure, stem(s3error)
% 
% figure, stem(tau1error*N)
% figure, stem(tau2error*N)
% figure, stem(tau3error*N)