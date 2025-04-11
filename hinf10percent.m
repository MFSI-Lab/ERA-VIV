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
p  = 18;
nt = 2*ngm + p;       %% total size of stored state space model
ns = 2*ngm + p;   %% actual size of state-space model

%% Read reduced state-space model
file = fopen('../ROM/Ar.txt','r+');
dat  = fscanf(file,'%g',[nt,nt]);
fclose(file);
A    = dat(1:ns,1:ns);

file = fopen('../ROM/Br.txt','r+');
dat  = fscanf(file,'%g',[nt,1]);
fclose(file);
B    = dat(1:ns,1);

file = fopen('../ROM/Cr.txt','r+');
dat  = fscanf(file,'%g',[1,nt]);
fclose(file);
C    = dat(1,1:ns);
%%
% Input-multiplicative perturbations (matrices Bd and Cd)
G    = 1.e-2;	% is G/W
%%
D    = [0 0 1; 0 sqrt(G) 0];
sys1 = ss(A,[B zeros(ns,1) B],[zeros(1,ns); C],D,dt,'InputName',{'wd','g','u'},'OutputName',{'zd','y'});

%% Initialize algorithm with H2 solution
% Solve Riccati equation
H       = [      A       -B*B'*(A^-1)';
             zeros(ns,ns)    (A^-1)'];
[U,H]   = schur(H,'complex');
[Us,Hs] = ordschur(U,H,'udi');
Xc      = Us(ns+1:2*ns,1:ns)/Us(1:ns,1:ns);
Xc      = real(Xc);

% Solve Riccati equation
H       = [A' + C'*(G^-1)*C*((A')^-1)'*B*B'  -C'*(G^-1)*C*((A')^-1)';
                 -((A')^-1)'*B*B'                 ((A')^-1)'];
[U,H]   = schur(H,'complex');
[Us,Hs] = ordschur(U,H,'udi');
Xe      = Us(ns+1:2*ns,1:ns)/Us(1:ns,1:ns);
Xe      = real(Xe);

% Controller in state-space form
Cc   = -B'*Xc;
Bc   = Xe*C'/G;
Ac   = A+B*Cc-Bc*C;
sys2 = ss(Ac,Bc,Cc,0,dt,'InputName','y','OutputName','u');
% Evaluation of closed-loop infinity norm
H    = feedback(sys1,sys2,3,2,1);
normold = norm(H(1,1:2),inf)
%% Loop on gamma to obtain Hinf solution
gammaold = 1.e30;

jjinit = log10(norm(H(1,1:2),inf)*2);
for jj = jjinit:-0.01:-100
    gamma = 10^jj;
% solve Riccati equation
    H       = [       A       -B*B'*(A^-1)'*(1.-1./gamma^2);
                 zeros(ns,ns)            (A^-1)'];
    [U,H]   = schur(H,'complex');
    [Us,Hs] = ordschur(U,H,'udi');
    Xc      = Us(ns+1:2*ns,1:ns)/Us(1:ns,1:ns);
    Xc      = real(Xc);
    
% solve Riccati equation
    H       = [A' + C'*(G^-1)*C*((A')^-1)'*B*B'  -C'*(G^-1)*C*((A')^-1)';
                    -((A')^-1)'*B*B'                 ((A')^-1)'];
    [U,H]   = schur(H,'complex');
    [Us,Hs] = ordschur(U,H,'udi');
    Xe      = Us(ns+1:2*ns,1:ns)/Us(1:ns,1:ns);
    Xe      = real(Xe);
    
    % exit criterion 
    if min(real(eig(eye(ns)-Xe*Xc/gamma^2))) <= 0
        break
    end
    % update controller 
    Cc   = -B'*Xc;  
	Bc   = (eye(ns)-Xe*Xc/gamma^2)\(Xe*C'/G);
	Ac   = A+(1.-1./gamma^2)*B*Cc-Bc*C;
    sys2 = ss(Ac,Bc,Cc,0,dt,'InputName','y','OutputName','u');
    % update controller and norm
    H    = feedback(sys1,sys2,3,2,1);
    normnew = norm(H(1,1:2),inf);
    if(normnew>gamma)
        break;
    end
    normold  = normnew
    gammaold = gamma;
    Acold = Ac;
    Bcold = Bc;
    Ccold = Cc;
end
%% Final controller 10 percent above minimum
gamma = gammaold*1.1;

H       = [       A       -B*B'*(A^-1)'*(1.-1./gamma^2);
             zeros(ns,ns)            (A^-1)'];
[U,H]   = schur(H,'complex');
[Us,Hs] = ordschur(U,H,'udi');
Xc      = Us(ns+1:2*ns,1:ns)/Us(1:ns,1:ns);
Xc      = real(Xc);

% solve Riccati equation

H       = [A' + C'*(G^-1)*C*((A')^-1)'*B*B'  -C'*(G^-1)*C*((A')^-1)';
                -((A')^-1)'*B*B'                 ((A')^-1)'];
[U,H]   = schur(H,'complex');
[Us,Hs] = ordschur(U,H,'udi');
Xe      = Us(ns+1:2*ns,1:ns)/Us(1:ns,1:ns);
Xe      = real(Xe);
% update controller 
Ccold = -B'*Xc;  
Bcold = (eye(ns)-Xe*Xc/gamma^2)\(Xe*C'/G);
Acold = A + (1.-1./gamma^2)*B*Ccold-Bcold*C;
%% Write controller in state-space form to file
AAfile=fopen('J.txt','w');
for i=1:size(Acold,1)
  for jj=1:size(Acold,1)
   fprintf(AAfile,'%.15g\n',Acold(jj,i));
  end
end
fclose(AAfile);
BBfile=fopen('L.txt','w');
CCfile=fopen('K.txt','w');
for i=1:size(Acold,1)
 fprintf(BBfile,'%.15g\n',Bcold(i));
 fprintf(CCfile,'%.15g\n',Ccold(i));
end
fclose(BBfile);
fclose(CCfile);
DDfile=fopen('M.txt','w');
fprintf(DDfile,'%.15g\n',0.);
fclose(DDfile);

%% performance format 2
w(:,1) = B; u1(1) = 0;
for i = 1:3001
   w(:,i+1) = Acold*w(:,i);
   y(i) = C*w(:,i);  
   u1(i+1) = Ccold*w(:,i+1);
end
hold on
plot(u1);
%% Performance and robusteness evaluations, how to do ???
B1 = [zeros(ns,1) B]; 
Bd = B; 
B2 = B; 
C1 = [C; zeros(1,ns)]; 
Cd = [zeros(1,ns)]; 
C2 = C; 
D  = [0 0 0 0; 0 0 0 1; 0 1 0 1; 1 0 0 0];
sys1 = ss(A,[B1 Bd B2],[C1; Cd; C2],D,dt,'InputName',{'g','w','wdu','u'},'OutputName',{'zy','zu','zd','y'});
sys2=ss(Acold,Bcold,Ccold,0,dt,'InputName','y','OutputName','u');
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
file=fopen(['perfo-im-10percent.txt'],'wt');
fprintf(file,'%s\n','VARIABLES= "N2zg",  "N2zw", "N2ug", "N2uw", "GM+", "GM-", "PM", "rho"');
fprintf(file,'%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n',...
               norm(H(1,1),2),norm(H(1,2),2),norm(H(2,1),2),norm(H(2,2),2),20*log10(gmp),20*log10(gmm),pm*180/pi,1./norm(H(3,3),inf));
fclose(file);
