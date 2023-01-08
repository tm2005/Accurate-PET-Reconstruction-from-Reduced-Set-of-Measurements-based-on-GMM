% close all
% clear all
% clc
clearvars x y k A B C D r phi
%%

Rs = 3;

N_comp = 3;

sigmas1 = [2.5;   2;  2; 3; 2]/10;
sigmas2 = [2.5;   3;  1; 3; 2]/10;

p =[0; 0.5; 0.3 ; -0.15; 0.1];

mus = [ 0  0-0.4      1.25    -0.6    0.3;     %.........xs
        0  0-0.4       -1       1       0.8 ];   %.........ys 
 
Ns = [3500; 2500; 1000; 1500; 1800];

Ns = round(Ns);
%%

X=[];


Sigmas = zeros(2,2,N_comp);
for i = 1:N_comp
    Sigmas(:,:,i) = [sigmas1(i)^2               sigmas1(i)*sigmas2(i)*p(i);
                    sigmas1(i)*sigmas2(i)*p(i)  sigmas2(i)^2];
    Xt = mvnrnd(mus(:,i),Sigmas(:,:,i),Ns(i));
    X = [X;Xt];
end
N = length(X);
mus = mus(:,1:N_comp);

phi = pi*rand(N,1)-pi/2;
k = tan(phi);

D = X(:,2)-k.*X(:,1);

A = 1+k.^2;
B = 2.*k.*D;
C = D.^2-Rs^2;

x(1,:) = -B./2./A + sqrt(B.^2-4.*A.*C)./2./A;
x(2,:) = -B./2./A - sqrt(B.^2-4.*A.*C)./2./A;

y(1,:) = k'.*x(1,:)+D';
y(2,:) = k'.*x(2,:)+D';


r = D./sqrt(1.^2+k.^2);
phi1 = atan( (y(2,:)-y(1,:))./(x(2,:)-x(1,:)) )';

p = randperm(length(r)); % ovo dodaje još slučajnosti

r = r(p);
phi = phi1(p);

% 
figure, scatter(phi(1:1:end), r(1:1:end), 15, [0.5 0.5 0.5], '+'), title('Random samples');
% 
% 
% figure, scatter(X(1:1:end,1), X(1:1:end,2), 12, [0.6 0.5 0.2], '+'), title('Random samples'), axis equal, xlim([-Rs Rs]), ylim([-Rs Rs]);
% 
% figure, scatter(phi(1:1:end), r(1:1:end), 15, [0.5 0.5 0.5], '+'), title('Random samples'), axis equal, xlim([-Rs Rs]), ylim([-Rs Rs]);

% save('simulation_two_comp','Rs','r','phi','N','Ns','mus','Sigmas','X')
