%% Script part 3

%% Simulation Echo de spin sur un TR
dT = 1;		% 1ms delta-time.
T = 1000;	% total duration
df = 10;	% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 50;	% ms.
TR = 500;	% ms.

N1 = round(TE/2/dT);
N2 = round((TR-TE/2)/dT);


% ===== Get the Propagation Matrix ======

[A,B] = freeprecess(dT,T1,T2,df);

% ===== Simulate the Decay ======

M = zeros(3,N1+N2);	% Keep track of magnetization at all time points.
M(:,1)=[0;0;1];	% Starting magnetization.

Rflip = yrot(pi/2);
Rrefoc = xrot(pi);

M(:,2)=A*Rflip*M(:,1)+B;
for k=3:(N1+1)
	M(:,k) = A*M(:,k-1)+B;
end;

M(:,N1+2)=A*Rrefoc*M(:,N1+1)+B;

for k=2:N2-1
	M(:,k+N1+1) = A*M(:,k+N1)+B;
end;

% ===== Plot the Results ======
figure;
time = [0:N1+N2-1]*dT;
plot(time,M(1,:),'b-',time,M(2,:),'r--',time,M(3,:),'g-.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;

%% Bloch Equation Simulation : Plusieurs spins
% -----------------------------------------
% 

dT = 1;		% 1ms delta-time.
T = 1000;	% total duration
Nf = 100;		% Number of frequencies.
df = 100*(rand(1,Nf)-.5);	% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 50;	% ms.
TR = 500;	% ms.

N1 = round(TE/2/dT);
N2 = round((TR-TE/2)/dT);

% ===== Get the Propagation Matrix ======

clear Msig M;

for f=1:length(df)
  [A,B] = freeprecess(dT,T1,T2,df(f));

  % -------- This section taken from b2a ----------------
  M = zeros(3,N1+N2);	% Keep track of magnetization at all time points.
  M(:,1)=[0;0;1];	% Starting magnetization.
  
  Rflip = yrot(pi/2);
  Rrefoc = xrot(pi);
  
  M(:,2)=A*Rflip*M(:,1)+B;
  for k=3:(N1+1)
	  M(:,k) = A*M(:,k-1)+B;
  end;
  
  M(:,N1+2)=A*Rrefoc*M(:,N1+1)+B;
  
  for k=2:N2-1
	  M(:,k+N1+1) = A*M(:,k+N1)+B;
  end;
  % -----------------------------------------------------
  Msig(f,:)=M(1,:)+1i*M(2,:);	% Just keep the signal component.

end;

% ===== Plot the Results ======

figure(1);
time = [0:N1+N2-1]*dT;
subplot(2,1,1);
plot(time,abs(Msig.'));
xlabel('Time (ms)');
ylabel('Magnitude');
axis([min(time) max(time) -1 1]);
grid on;

subplot(2,1,2);
plot(time,angle(Msig.'));
xlabel('Time (ms)');
ylabel('Phase (radians)');
axis([min(time) max(time) -pi pi]);
grid on;


figure(2);
time = [0:N1+N2-1]*dT;
plot(time,abs(mean(Msig)));
xlabel('Time (ms)');
ylabel('Net Magnitude');
axis([min(time) max(time) -1 1]);
grid on;

%% Steady-state

df = 0;	% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 5;	% ms.
TR = 50;	% ms.

[Msig,Mss] = sesignal(T1,T2,TE,TR,df)

%% Fast Spin Echo

df = 0;	% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 50;	% ms.
TR = 1000;	% ms.
ETL = 8;

[Msig,Mss] = fsesignal(T1,T2,TE,TR,df,ETL);
Mss

ETL = 1;
[Msig,Mss] = fsesignal(T1,T2,TE,TR,df,ETL);
Mss