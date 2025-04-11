%% LTR controller

clear all
clc
format long;
%% read the control parameters, where is it from??
R = 1e-3;
G = 1e2;

ngm   = 1;
fid   = fopen('../Common/param.dat','rt');
num   = fscanf(fid,'%d',1)  %% number of time steps in simulation
dt    = fscanf(fid,'%f',1)  %% time step
stp   = fscanf(fid,'%d',1)  %% number of time steps between two snapshots
p     = fscanf(fid,'%d',1)  %% total number of computed bpod modes
nstab = fscanf(fid,'%d',1)  %% number of bpod modes in rom
fclose(fid);

p = 20;
nt =  p;       %% total size of stored state space model
ns =  p;   %% actual size of state-space model

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
%% The optimization problen for control yields algebraic matrix Riccati equation for an auxiliary variable R;
%% Q = C'*C;
% solve Riccati equation
H       = [A + B*(R^-1)*B'*(A^-1)'*C'*C  -B*(R^-1)*B'*(A^-1)';
                 -(A^-1)'*C'*C                  (A^-1)'];
[U,H]   = schur(H,'complex');
[Us,Hs] = ordschur(U,H,'udi');
X       = Us(ns+1:2*ns,1:ns)/Us(1:ns,1:ns);
X       = real(X);
K       = -(R + B'*X*B)^-1*B'*X*A;

log((eig(A+B*K)))/0.02

%% solve the Riccati equation for noise
H       = [A' + C'*(G^-1)*C*((A')^-1)'*B*B'  -C'*(G^-1)*C*((A')^-1)';
                 -((A')^-1)'*B*B'                 ((A')^-1)'];
[U,H]   = schur(H,'complex');
[Us,Hs] = ordschur(U,H,'udi');
X       = Us(ns+1:2*ns,1:ns)/Us(1:ns,1:ns);
X       = real(X);
L       = -A*X*C'*(G + C*X*C')^-1;

log((eig(A+B*K+L*C)))/0.02
%% Write controller in state-space form to file
J  = A + L*C + B*K;
JJ = A + B*K;
%% Performance
dt  = 0.02;
sys0 = ss(A,B,C,0,dt);
sys1 = ss(J,B,C,0,dt);
sys2 = ss(JJ,B,C,0,dt);

X   = 1:1:3001;
t   = dt*X;
[y0,t,x0] = initial(sys0,B,t); 
[y1,t,x1] = initial(sys1,B,t); 
[y2,t,x2] = initial(sys2,B,t); 

% plot(t,y1,'r');
% hold on
% plot(t,y2,'b');

%% performance format 2
w(:,1) = B; u1(1) = 0;
for i = 1:3001
   w(:,i+1) = J*w(:,i);
   y(i) = C*w(:,i);  
   u1(i+1) = K*w(:,i+1);
end
hold on
plot(u1);

AAfile = fopen('J.txt','w');
for i = 1:size(J,1)
  for jj = 1:size(J,1)
    fprintf(AAfile,'%.15g\n',J(jj,i));
  end
end
fclose(AAfile);

BBfile = fopen('L.txt','w');
CCfile = fopen('K.txt','w');
for i = 1:size(J,1)
  fprintf(BBfile,'%.15g\n',L(i));
  fprintf(CCfile,'%.15g\n',K(i));
end
fclose(BBfile);
fclose(CCfile);
DDfile = fopen('M.txt','w');
fprintf(DDfile,'%.15g\n',0.);
fclose(DDfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% transfer function, method 1
% sysc = ss(J,-L,K,0,dt);
% [mag,phase,ww] = bode(sys0);
% figure
% plot(ww(20:90,1),(mag(1,20:90)),'r+-');
% grid
%% transfer function, method 2
s = 0:0.01:3;
z = exp(s*1i*dt);
I = eye(size(A));
for i = 1 : 301
    TF(i) = C*(z(i)*I - A)^(-1)*B;
end
figure
plot(s,abs(log(TF)/dt));
%% transfer function, method 3
[b,a] = ss2tf(A,B,C,0);
% figure
% plot(a,b);  %% a,b are coefficient of the transfer function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        gmp = fzero(f,[1 gmpmaxi],optimset('TolX',1.e-6))
        break
    else
        gmpmaxi=gmpmaxi*10;
    end
end

% GM-
if f(0)*f(1)>=0
   gmm = 0
else
   gmm = fzero(f,[0 1])
end

% PM    
g = @(x) max(abs(zero(1-exp(-1i*x)*sys*sys2))) - 1;
if g(0)*g(pi/2)>=0
   pm=pi/2
else
   pm=fzero(g,[0 pi/2])
end
%% Write to file
file=fopen(['perfo-lqg.txt'],'wt');
fprintf(file,'%s\n','VARIABLES= "N2zg",  "N2zw", "N2ug", "N2uw", "GM+", "GM-", "PM", "rho"');
fprintf(file,'%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n',...
               norm(H(1,1),2),norm(H(1,2),2),norm(H(2,1),2),norm(H(2,2),2),20*log10(gmp),20*log10(gmm),pm*180/pi,1./norm(H(3,3),inf));
fclose(file);

