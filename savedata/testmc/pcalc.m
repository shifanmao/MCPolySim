close all
clear

L=4;
LTOT=4;
LP=1;

NRUN=1;
NP=3068;
%home='../AJS-testmc-2-6-14C/savedata/';
home='../../data/';
if NRUN==1
    R=load([home,'r0']);    
else
    R=load([home,'data1/r0']);
end
N=length(R)/NP;

I0=70;
IF=100;
SKIP=1;

NH=50;
DX=1.5/NH;
P=zeros(NH,1);
X=transpose([(DX/2):DX:(1.5-DX/2)]);
COUNT=0;

IL=round(L/LTOT*N);
L=(IL-1)/(N-1)*LTOT;
I1=([1:1:NP]-1)*N+1;
I2=([1:1:NP]-1)*N+IL;

for IR=1:NRUN
    for IS=I0:SKIP:IF
        if NRUN == 1
            R=load([home,'r',int2str(IS)]);            
        else
            R=load([home,'data',int2str(IR),'/r',int2str(IS)]);
        end
        REND=sqrt(sum(power(R(I2,1:3)-R(I1,1:3),2),2))/L;
        P=P+transpose(hist(REND,X));
        COUNT=COUNT+NP;
    end
    IR
end

P=P/(COUNT*DX);

FNUM=round(100*L/(2*LP));
G=load(['data/out',int2str(FNUM),'.txt']);
XG=transpose(linspace(0,1,5001));

figure(1)
plot(XG,G.*power(XG,2)*(4*pi),'b-','LineWidth',2)
hold on
plot(X,P,'k.-','LineWidth',2)
axis([0 1 0 3])
set(gca,'FontSize',14)
xlabel('R/L','FontSize',18)
ylabel('P(R)','FontSize',18)
saveas(gcf,'figures/ptest.eps','epsc')

figure(2)
semilogy(XG,G.*power(XG,2)*(4*pi),'b-','LineWidth',2)
hold on
semilogy(X,P,'k.-','LineWidth',2)

figure(3)
plot(XG,G*(4*pi),'b-','LineWidth',2)
hold on
plot(X,P./power(X,2),'k.-','LineWidth',2)

figure(4)
semilogy(XG,G*(4*pi),'b-','LineWidth',2)
hold on
semilogy(X,P./power(X,2),'k.-','LineWidth',2)
