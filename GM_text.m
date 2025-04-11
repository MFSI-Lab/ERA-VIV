%% LTR controller

clear all
clc
format long;

ngm   = 1;
fid   = fopen('../Common/param.dat','rt');
num   = fscanf(fid,'%d',1)  %% number of time steps in simulation
dt    = fscanf(fid,'%f',1)  %% time step
stp   = fscanf(fid,'%d',1)  %% number of time steps between two snapshots
p     = fscanf(fid,'%d',1)  %% total number of computed bpod modes
nstab = fscanf(fid,'%d',1)  %% number of bpod modes in rom
fclose(fid);
nt =  p;       %% total size of stored state space model
ns = 2*ngm + nstab;   %% actual size of state-space model

%% define matrix
Num = 101;
GM_positive = zeros(Num,Num);
GM_negative = zeros(Num,Num);
PM  = zeros(Num,Num);
rho = zeros(Num,Num);
rates = zeros(Num,Num);
%% Read reduced state-space model
file = fopen('../ROM/ArRe80m20U12.txt','r+');
dat  = fscanf(file,'%g',[nt,nt]);
fclose(file);
A    = dat(1:ns,1:ns);

file = fopen('../ROM/BrRe80m20U12.txt','r+');
dat  = fscanf(file,'%g',[nt,1]);
fclose(file);
B    = dat(1:ns,1);

file = fopen('../ROM/CrRe80m20U12.txt','r+');
dat  = fscanf(file,'%g',[1,nt]);
fclose(file);
C    = dat(1,1:ns);

log((eig(A)))/0.02
%% end of read ROM matrix
%% control parameter
R = 1e-11;
G = 1e-11;
for m = 1:Num
    R = 1e-10*10^(0.2*(m-1));
    for n = 1:Num
        G = 1e-10*10^(0.2*(n-1));
%% solve Riccati equation
        H       = [A + B*(R^-1)*B'*(A^-1)'*C'*C  -B*(R^-1)*B'*(A^-1)';
                        -(A^-1)'*C'*C                  (A^-1)'];
        [U,H]   = schur(H,'complex');
        [Us,Hs] = ordschur(U,H,'udi');
        X       = Us(ns+1:2*ns,1:ns)/Us(1:ns,1:ns);
        X       = real(X);
        K       = -(R + B'*X*B)^-1*B'*X*A;
%% solve the Riccati equation for noise
        H       = [A' + C'*(G^-1)*C*((A')^-1)'*B*B'  -C'*(G^-1)*C*((A')^-1)';
                        -((A')^-1)'*B*B'                 ((A')^-1)'];
        [U,H]   = schur(H,'complex');
        [Us,Hs] = ordschur(U,H,'udi');
        X       = Us(ns+1:2*ns,1:ns)/Us(1:ns,1:ns);
        X       = real(X);
        L       = -A*X*C'*(G + C*X*C')^-1;
        
        J  = A + L*C + B*K;
%% Performance and robusteness evaluations, how to do ???
        B1 = [zeros(ns,1) B]; 
        Bd = B; 
        B2 = B; 
        C1 = [C; zeros(1,ns)]; 
        Cd = [zeros(1,ns)]; 
        C2 = C; 
        D  = [0 0 0 0; 0 0 0 1; 0 1 0 1; 1 0 0 0];
        sys1 = ss(A,[B1 Bd B2],[C1; Cd; C2],D,dt,'InputName',{'g','w','wdu','u'},'OutputName',{'zy','zu','zd','y'});
        sys2 = ss(J,-L,K,0,dt,'InputName','y','OutputName','u');
        H = feedback(sys1,sys2,4,4,1);

% GM+
        sys = ss(A,B2,C2,0,dt);

        gmpmaxi = 10;
        f = @(x) max(abs(zero(1-x*sys*sys2)))-1;
        for i = 1 : 25
            if f(1)*f(gmpmaxi)<=0 
                gmp = fzero(f,[1 gmpmaxi],optimset('TolX',1.e-6));
                break
            else
                gmpmaxi=gmpmaxi*10;
            end
        end

% GM-
        if f(0)*f(1)>=0
            gmm = 0;
        else
            gmm = fzero(f,[0 1]);
        end

% PM    
        g = @(x) max(abs(zero(1-exp(-1i*x)*sys*sys2))) - 1;
        if g(0)*g(pi/2)>=0
            pm=pi/2;
        else
            pm=fzero(g,[0 pi/2]);
        end
% save in matrix
        GM_positive(n,m) = 20*log10(gmp);
        GM_negative(n,m) = 20*log10(gmm);
        PM(n,m)  = pm*180/pi;
        rho(n,m) = 1./norm(H(3,3),inf);
        rates(n,m) = max(real(log(eig(J))/dt));
        m*n
    end
end
sumGM_negative = sum(sum(GM_negative))/121
sumGM_positive = sum(sum(GM_positive))/121
sumGM_PM       = sum(sum(PM))/121
sumGM_rho      = sum(sum(rho))/121
contourf(GM_positive,20);
colorbar
% figure
det = GM_positive - GM_negative;
% 
% x = -10:0.5:10;
% y = -10:0.5:10;
% surf(x,y,det);
% colormap('jet');               %设置颜色
% colorbar
% shading interp;   
% material shiny;
% set(gcf,'Color',[1 1 1]);
% zlabel('GM','FontName','Times New Roman','FontSize', 20);
% xlabel('R','FontName','Times New Roman','FontSize', 20);
% ylabel('G','FontName','Times New Roman','FontSize', 20);
