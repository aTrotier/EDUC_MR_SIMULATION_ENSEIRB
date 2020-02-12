%% script spoiling

T1 = 1000;
T2 =  50;
flip = pi/3;
TE = 3;
TR = 1000;
dfreq=0;
%% Question 1 : spoiler naturel

TR = 5 % TR << T2
Msig_NoSpoil = abs(sssignal(flip,T1,T2,TE,TR,dfreq))
Msig_spoiled = abs(srsignal(flip,T1,T2,TE,TR,dfreq))

% Résultats différents

TR = 1000 % TR >> T2
Msig_NoSpoil = abs(sssignal(flip,T1,T2,TE,TR,dfreq))
Msig_PerfSpoiled = abs(srsignal(flip,T1,T2,TE,TR,dfreq))

% Résultats identique

%% Question 2 : gradient spoiling
flip = pi/3;
T1 = 600;
T2 = 100;
TE = 2;
TR = 10;
dfreq = 0;
phi= pi/2;

% add one spoiler
[Msig,Mss] = gssignal(flip,T1,T2,TE,TR,dfreq,phi)

% gresignal(flip,T1,T2,TE,TR,dfreq) permet de calculer l'effet avec
% plusieurs isochromats dans un voxel dephasé différement selon leur
% position
%%
% no spoil
Msig_NoSpoil = abs(sssignal(flip,T1,T2,TE,TR,dfreq))

% grad spoil
Msig_gradSpoil = abs(gresignal(flip,T1,T2,TE,TR,dfreq))

% comparison with perfect spoil :
Msig_PerfSpoiled = abs(srsignal(flip,T1,T2,TE,TR,dfreq))

% Mxy signal is higher with grad_spoil than perfect_spoil

%% Question 3
inc = 117;
phi(1)=0;
for n = 2:5
    phi(n)=(n-1)*inc+phi(n-1)
end

%%
clear Msig
df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 2;		% ms.
TR = 10;
flip = pi/6;	% radians.
inc = 117/180*pi;

Nex = 100;

Nf = 200;	% Simulate 200 different gradient-spoiled spins.
phi = [1:Nf]/Nf*2*pi;

M=zeros(3,Nf,Nex+1);
Msig = zeros(Nex,1);



[Ate,Bte] = freeprecess(TE,T1,T2,0);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,0);

M = [zeros(2,Nf);ones(1,Nf)];


Rfph = 0;
Rfinc = inc;

for n=1:Nex
    
    A = Ate * throt(flip,Rfph);
    B = Bte;
    M = A*M+B;
    
    Msig(n) = mean( squeeze(M(1,:)+1j*M(2,:)) ) * exp(-1j*Rfph);
    
    M=Atr*M+Btr;
    
    for k=1:Nf
        M(:,k) = zrot(phi(k))*M(:,k);
    end
    
    Rfph = Rfph+Rfinc;
    Rfinc = Rfinc+inc;
end


% ===== Plot the Results ======
figure;
time = [0:Nex-1]*TR+TE;
subplot(2,1,1);
plot(time,abs(Msig));
xlabel('Time (ms)');
ylabel('Magnitude');
grid on;


subplot(2,1,2);
plot(time,angle(Msig));
xlabel('Time (ms)');
ylabel('Phase');
grid on;

% compare perfect spoiled and this one
Msig_PerfSpoiled = abs(srsignal(flip,T1,T2,TE,TR,dfreq))
abs(Msig(end))

%% Autre possibilité d'écriture : 
 clear Msig
df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 2;		% ms.
TR = 10;
flip = pi/6;	% radians.
inc = 117/180*pi;

Nex = 100;

Nf = 200;	% Simulate 200 different gradient-spoiled spins.
phi = [1:Nf]/Nf*2*pi;

M=zeros(3,Nf,Nex+1);
Msig = zeros(Nex);



[Ate,Bte] = freeprecess(TE,T1,T2,0);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,0);

M = [zeros(2,Nf);ones(1,Nf)];

for k=1:Nf
    Rfph = 0;
    Rfinc = inc;
    M=[0;0;1];
    for n=1:Nex
        
        A = Ate * throt(flip,Rfph);
        B = Bte;
        M = A*M+B;
        
        Msig(n,k) = mean( squeeze(M(1)+1j*M(2)) ) * exp(-1j*Rfph);
        
        M=Atr*M+Btr;
        
        M = zrot(phi(k))*M;
        
        Rfph = Rfph+Rfinc;
        Rfinc = Rfinc+inc;
    end
end
Msig=mean(Msig,2);

% ===== Plot the Results ======
figure;
time = [0:Nex-1]*TR+TE;
subplot(2,1,1);
plot(time,abs(Msig));
xlabel('Time (ms)');
ylabel('Magnitude');
grid on;


subplot(2,1,2);
plot(time,angle(Msig));
xlabel('Time (ms)');
ylabel('Phase');
grid on;

% compare perfect spoiled and this one
Msig_PerfSpoiled = abs(srsignal(flip,T1,T2,TE,TR,dfreq))
abs(Msig(end))

%% compare across flip angle perfect spoiled and rf spoil

df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 2;		% ms.
TR = 10;
flip = [0:0.01:0.5 ]*pi;	% radians.
inc = 117/180*pi;
Nex = 100;


for k=1:length(flip)
	sig1(k)=spgrsignal(flip(k),T1,T2,TE,TR,df,Nex,inc);
    
	[Msig,M]=srsignal(flip(k),T1,T2,TE,TR,df);
	sig2(k)=M(1)+1i*M(2);
    
    [Msig,M]=gresignal(flip(k),T1,T2,TE,TR,dfreq);
    sig3(k)=M(1)+1i*M(2);
end


% ===== Plot the Results ======

plot(flip*180/pi,abs(sig1),'r-',flip*180/pi,abs(sig2),'b--',flip*180/pi,abs(sig3),'g--');
xlabel('Flip (rad)');
ylabel('Signal Magnitude');
grid on;
title('Signal vs Flip for RF-spoiled GRE and SR');
legend('RF-spoiled GRE','Saturation-Recovery Approximation','Gradient Spoiled');

