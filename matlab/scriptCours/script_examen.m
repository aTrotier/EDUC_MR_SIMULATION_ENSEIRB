% 4 a)-------------------------------------------------------------------------------
% Voir papier Van Uijen adn Den Boef ! Driven-equelibrium Radiofrequency
% pulses in NMR Imaging, 1984

%t1=te/2
%t2=te
%t3=te/2
%t4=tr-2*te
% Mz(90-)=M0*(1-exp(-(t1+t2+t3)/T1)  et Mx(90-)=M0*exp(-(t1+t2+t3)/T2)
% Mz(90+)=Mx(90-)
% Mz(TR)=M0-(M0-M(TR-1)*exp(-(t1+t2+t3)/T2)*exp(-t4/T1)

% => Mzz=M0*(1-exp(-t4/T1))/(1-exp(-(t1+t2+t3)/T2)*exp(-t4/T1))
%   => S(TE)=Mzz*exp(-TE/T2)

% 4) b) colormap of signal -------------------------------------------------------------------------------


matrixSig=zeros(length(0:1:3000),length(0:1:1000));
idxT1=0;
TE=15;
TR=400;
for T1=0:1:3000
    idxT1=idxT1+1;
    idxT2=0;
    for T2=0:1:1000
        idxT2=idxT2+1;
        
        Mzz=(1-exp(-(TR-2*TE)/T1))/(1-exp(-(2*TE)/T2)*exp(-(TR-2*TE)/T1));
        matrixSig(idxT1,idxT2)=Mzz*exp(-TE/T2);
    end
end

figure;imagesc(matrixSig);colorbar();


%4)c) it is a T2 contrast (most of the variation came from T2 variation)
figure;contour(matrixSig,'ShowText','on');

%% par simulation 

df = 0;	% Hz off-resonance.
T1 = 600;	% ms.
T2 = 300;	% ms.

TE = 15;	% ms.
TR = 400;	% ms.

% ===== Get the Propagation Matrix ======

[Ate2,Bte2]=freeprecess(TE/2,T1,T2,df);
[Ate,Bte]=freeprecess(TE,T1,T2,df);
[Atr,Btr]=freeprecess(TR-2*TE,T1,T2,df);



Aeq=(yrot(pi/2)*Ate2)*(xrot(pi)*Ate)*(xrot(pi)*Ate2)*(yrot(-pi/2)*Atr);

Beq=(xrot(pi)*Ate)*(xrot(pi)*Ate2)*(yrot(-pi/2)*Atr)*Bte2+...
    Ate2*xrot(pi)*yrot(-pi/2)*Atr*Bte+...
    yrot(-pi/2)*Atr*Bte2+Btr;

Mss = inv(eye(3)-Aeq)*Beq

% comparaison avec m?thode analytique
Mzz=(1-exp(-(TR-2*TE)/T1))/(1-exp(-(2*TE)/T2)*exp(-(TR-2*TE)/T1))
     
% Encapsulation :
[Msig,Mss2] = derpsignal(T1,T2,TE,TR,0)
%%
clear matrixSig Msig
idxT1=0;
TE=15;
TR=400;
for T1=0:100:3000
    idxT1=idxT1+1
    idxT2=0;
    for T2=0:10:1000
        idxT2=idxT2+1;
        
        Mzz=(1-exp(-(TR-2*TE)/T1))/(1-exp(-(2*TE)/T2)*exp(-(TR-2*TE)/T1));
        matrixSig(idxT1,idxT2)=Mzz*exp(-TE/T2);
        [Msig(idxT1,idxT2),Mss2] = derpsignal(T1,T2,TE,TR,0);
    end
end

% figure;imagesc(flipud(matrixSig));colorbar();
% ylabels = cellstr(num2str([3000:-300:0]'));  %time labels
% set(gca, 'YTick', 0:3:31,'YTickLabel', ylabels)
% xlabels = cellstr(num2str([0:100:1000]'));  %time labels
% set(gca, 'XTick', 0:10:101,'XTickLabel', xlabels)


figure;imagesc(abs(Msig));colorbar();
ylabels = cellstr(num2str([3000:-300:0]'));  %time labels
set(gca, 'YTick', 0:3:31,'YTickLabel', ylabels)
xlabels = cellstr(num2str([0:100:1000]'));  %time labels
set(gca, 'XTick', 0:10:101,'XTickLabel', xlabels)

%4)c) it is a T2 contrast (most of the variation came from T2 variation)
figure;contour((matrixSig),'ShowText','on')
ylabels = cellstr(num2str([0:300:3000]'));  %time labels
set(gca, 'YTick', 0:3:31,'YTickLabel', ylabels)
xlabels = cellstr(num2str([0:100:1000]'));  %time labels
set(gca, 'XTick', 0:10:101,'XTickLabel', xlabels)