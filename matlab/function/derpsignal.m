% 
%	function [Msig,Mss] = derpsignal(T1,T2,TE,TR,df)
% 
%	Calculate the steady state signal at TE for a multi-echo spin-echo
%	sequence, given T1,T2,TR,TE in ms.  Force the
%	transverse magnetization to zero before each excitation.
%	df is the resonant frequency in Hz.  flip is in radians.
%

% Sequence :
% 90y--TE/2--180x--TE--180x---90y---TR

function [Msig,Mss] = derpsignal(T1,T2,TE,TR,df)

% ===== Get the Propagation Matrix ======

[Ate2,Bte2]=freeprecess(TE/2,T1,T2,df);
[Ate,Bte]=freeprecess(TE,T1,T2,df);
[Atr,Btr]=freeprecess(TR-2*TE,T1,T2,df);

% Propagate A and B
Aeq=(yrot(pi/2)*Ate2)*(xrot(pi)*Ate)*(xrot(pi)*Ate2)*(yrot(-pi/2)*Atr);

Beq=(xrot(pi)*Ate)*(xrot(pi)*Ate2)*(yrot(-pi/2)*Atr)*Bte2+...
    Ate2*xrot(pi)*yrot(-pi/2)*Atr*Bte+...
    yrot(-pi/2)*Atr*Bte2+Btr;


Mss = inv(eye(3)-Aeq)*Beq;

% Calculate signal at echo.
M = Ate2*xrot(pi)*Ate2*yrot(pi/2)*Mss + Bte2+Ate2*xrot(pi)*Bte2;
Msig=M(1)+1i*M(2);

end




