%% Free Precession

% b) plot the evolution of Mx,My and Mz vs time from M=[1,0,0]' with a time
% step of 1 ms and plot the response for 1000 ms.

dT = 1;		% 1ms delta-time.
T = 1000;	% total duration
N = ceil(T/dT)+1; % number of time steps.
df = 10;	% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.

% ===== Get the Propagation Matrix ======

[A,B] = freeprecess(dT,T1,T2,df);


% ===== Simulate the Decay ======

M = zeros(3,N);	% Keep track of magnetization at all time points.
M(:,1)=[1;0;0];	% Starting magnetization.

for k=2:N
	M(:,k) = A*M(:,k-1)+B;
end


% ===== Plot the Results ======

time = [0:N-1]*dT;
plot(time,M(1,:),'b-',time,M(2,:),'r--',time,M(3,:),'g-.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;


phi = 2*pi*df*dT/1000

zrot(phi)

%% 1) Gradient Echo
% The sequence simply consists of 60-degree excitation pulses about the y-axis, spaced TR apart.

%% a)Magnetization after first excitation at TE
% initial parameters :
M=[0,0,1]';  %Initial value along axis [x,y,z]

df = 0;	% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 1;     %ms
TR = 500;     %ms
alpha = 60; %deg

% What is the magnetization 1 ms after the first excitation????

%step 1 : Magnetization after flip;
Rflip=yrot(60,'deg');
M=Rflip*M;

%step 2 : Relaxation process
[Ate,Bte]=freeprecess(TE,T1,T2,df);

M=Ate*M+Bte

%% b) M after the second excitation ?
% first step : magnetization at TR
M0=[0,0,1]';
[Atr,Btr]=freeprecess(TR,T1,T2,df);
M0=Rflip*M0;
MTR=Atr*M0+Btr;

MTR=Rflip*MTR;
MTE=Ate*MTR+Bte

%% c) Magnetization over the 10 excitations (Mx/y/z vs time)


df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 1;		% ms.
dT = 1;
TR = 500;	% ms.
flip = pi/3;	% radians.
Ntr = round(TR/dT);
Nex = 10;	% 20 excitations.

M = [0;0;1];
Rflip = yrot(flip);
[A1,B1] = freeprecess(dT,T1,T2,df);


M(1,Nex*Ntr)=0;	%	Allocate to record all M's.
		% 	Not necessary, but makes program faster.

Mcount=1;
for n=1:Nex
	M(:,Mcount) = Rflip*M(:,Mcount);	

	for k=1:Ntr
		Mcount=Mcount+1;
		M(:,Mcount)=A1*M(:,Mcount-1)+B1;
	end;
end;

time = [0:Mcount-1]*dT;
plot(time,M(1,:),'b-',time,M(2,:),'r--',time,M(3,:),'g-.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;

%% d) After 4 excitations we reach steady states. We can calculate this Steady state directly :

df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 1;		% ms.
TR = 500;	% ms.
flip = pi/3;	% radians.


M = [0;0;1];
Rflip = yrot(flip);
[Atr,Btr] = freeprecess(TR,T1,T2,df);

% A*Rz=Rz*A so
%
% M1 = Atr * Rflip * M + Btr.
%
% But M1=M in steady state, so
%
% 	Mss = Atr*Rflip * Mss + Btr.
%	(I-Atr*Rflip)*Mss = Btr.
%   (I-Atr*Rflip)^(-1)*(I-Atr*Rflip)*Mss = (I-Atr*Rflip)^(-1) * Btr 

Mss = inv(eye(3)-Atr*Rflip)*Btr

%% e) Steady state at time TE ??

df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 1;		% ms.
TR = 500;	% ms.
flip = pi/3;	% radians.
M = [0;0;1];

Rflip = yrot(flip);
[Atr,Btr] = freeprecess(TR,T1,T2,df);
[Ate,Bte] = freeprecess(TE,T1,T2,df);
[Atetr,Btetr] = freeprecess(TR-TE,T1,T2,df);

%	Calculation using B-1d first:

Mss = inv(eye(3)-Atr*Rflip)*Btr;
Mte1 = Ate*Rflip*Mss+Bte


% 	Direct calculation at TE

% 	Starting at TE, M=M1
%	At TR, M=M2, and M2=Atetr*M1+Btetr.
%	At TE, M=M3, and M3=Ate*Rflip*M2+Bte.
%			M3=Ate*Rflip*(Atetr*M1+Btetr)+Bte.
%
%	But M3=M1=Mte2 in steady state:

Mte2 = inv(eye(3)-Ate*Rflip*Atetr)* (Ate*Rflip*Btetr+Bte)


% direct calculation :
[mSig,Mss]=sssignal(flip,T1,T2,TE,TR,df) % with second method (starting at TE)

%% B-1f 

% check the value if we put magnetization along x to 0

df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 1;		% ms.
TR = 500;	% ms.
flip = pi/3;	% radians.


M = [0;0;1];
Rflip = yrot(flip);
[Atr,Btr] = freeprecess(TR,T1,T2,df);
Mss = inv(eye(3)-Atr*Rflip)*Btr

spoiler = [0 0 0;0 0 0;0 0 1];
Mss = inv(eye(3)-spoiler*Atr*Rflip)*Btr