%% this code using Re = 80 control strategy for higer Reynolds number (Re = 100);
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
Num = 11;
GM_positive = zeros(Num,Num);
GM_negative = zeros(Num,Num);
PM  = zeros(Num,Num);
rho = zeros(Num,Num);
rates = zeros(Num,Num);
%% Read ROM at Re = 100
file = fopen('../ROM/ArRe90m20U12.txt','r+');
dat  = fscanf(file,'%g',[nt,nt]);
fclose(file);
A100    = dat(1:ns,1:ns);

file = fopen('../ROM/BrRe90m20U12.txt','r+');
dat  = fscanf(file,'%g',[nt,1]);
fclose(file);
B100    = dat(1:ns,1);

file = fopen('../ROM/CrRe90m20U12.txt','r+');
dat  = fscanf(file,'%g',[1,nt]);
fclose(file);
C100    = dat(1,1:ns);
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
    R = 1e-10*10^(2*(m-1));
    for n = 1:Num
        G = 1e-10*10^(2*(n-1));
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
        J100  = A100 + L*C100 + B100*K;
        real(log(eig(J))/dt)
%% save in matrix
        rates(n,m) = max(real(log(eig(J))/dt));
        m*n
    end
end
contourf(rates,20);
colorbar
% figure
% det = GM_positive - GM_negative;
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
