%% Script part 4 Steady-state

%% BSSFP TE EFFECT

T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 0:2.5:10;	% ms.
TR = 10;	% ms.
flip = pi/3;

df = [-100:100]; 	%Hz

Sig = zeros(length(df),length(TE));

for n=1:length(TE)

	for k=1:length(df)
		[Msig,Mss] = sssignal(flip,T1,T2,TE(n),TR,df(k));
		Sig(k,n)=Mss(1)+i*Mss(2);
	end;

end;
	% ===== Plot the Results ======

subplot(2,1,1);
plot(df,abs(Sig));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
legend('TE=0', 'TE=2.5', 'TE=5.0', 'TE=7.5', 'TE=10');

subplot(2,1,2);
plot(df,angle(Sig));
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
axis([min(df) max(df) -pi pi]);
grid on;
legend('TE=0', 'TE=2.5', 'TE=5.0', 'TE=7.5', 'TE=10');

%% BSSFP TR effect

TE = 0;
T1 = 600;	% ms.
T2 = 100;	% ms.
TR = [2,6,10];	% ms.
flip = pi/3;

df = [-500:500]; 	%Hz

Sig = zeros(length(df),length(TE));

for n=1:length(TR)

	for k=1:length(df)
		[Msig,Mss] = sssignal(flip,T1,T2,TR(n)/2,TR(n),df(k));
		Sig(k,n)=Mss(1)+i*Mss(2);
	end;

end;
	% ===== Plot the Results ======
figure;
subplot(2,1,1);
plot(df,abs(Sig));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
legend( 'TR=2',  'TR=6', 'TR=10');

subplot(2,1,2);
plot(df,angle(Sig));
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
axis([min(df) max(df) -pi pi]);
grid on;
legend( 'TR=2',  'TR=6', 'TR=10');


%% Shifting of zero signal 
% Apply a RF phase increment use another function called ssfp

T1 = 600;	% ms.
T2 = 100;	% ms.
TE = 2.5;	% ms.
TR = 5;	% ms.
flip = pi/3;
phi=[0 pi/2 pi 1.5*pi];

df = [-500:10:500]; 	%Hz

Sig = zeros(length(df),length(TE));

for n=1:length(phi)

	for k=1:length(df)
		Mss = ssfp(flip,T1,T2,TE,TR,df(k),phi(n));
		Sig(k,n)=Mss;
	end;

end;
	% ===== Plot the Results ======
figure;
subplot(2,1,1);
plot(df,abs(Sig));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
legend( '\phi=0', '\phi=\pi/2', '\phi=\pi', '\phi=3\pi/2');

subplot(2,1,2);
plot(df,angle(Sig));
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
axis([min(df) max(df) -pi pi]);
grid on;
legend( '\phi=0', '\phi=\pi/2', '\phi=\pi', '\phi=3\pi/2');

