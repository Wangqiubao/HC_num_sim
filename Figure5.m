%%%Figure 5. Effect of different time delays on periodic oscillation when k = 0.1%%
clear 
 clc
global mu  gamma1 delta rho varphi  alpha omega beta0 betag N eta gamma2 tau epsilon  Stepsize K k1
mu=0.007;
gamma1=0.00015;
delta=0.000017;
rho=0.25;
varphi=0.02;
alpha=0.8;
omega=75000;
beta0=0.9;
betag=0.81;
N=10000000;
eta=0.9;
gamma2=1/14;
epsilon=0.01;
tau=20;
k1=0.1;
K=0;

randn('state',100);
Stepsize=0.01;
FinalTime=4500;
NumberOfSample=100;
InitialValue=1000;
DiffusionCoef=0.1;
DriftCoef=1;
dimension=2;
m=ceil(tau/Stepsize);
TimeGrid = 0 : Stepsize : FinalTime;
NumberOfSubinterval = length(TimeGrid)-1;
dW = sqrt(Stepsize) * randn(NumberOfSample,NumberOfSubinterval+1); 
Solution = ones(NumberOfSample,NumberOfSubinterval+1).*552000;
Esolution=ones(NumberOfSample,NumberOfSubinterval+1).*1405470;
Asolution = ones(NumberOfSample,NumberOfSubinterval+1).*87494;
Qsolution=ones(NumberOfSample,NumberOfSubinterval+1).*71666;
Vsolution = ones(NumberOfSample,NumberOfSubinterval+1).*4741153;
betat=ones(NumberOfSample,NumberOfSubinterval+1).*beta0;
Rsolution= ones(NumberOfSample,NumberOfSubinterval+1).*3142180;
for k=1:m
    
Solution(:,k) =463863;
Esolution(:,k) =201495;
Asolution(:,k) =52454;
Qsolution(:,k)=49096;
Vsolution(:,k) =3636852;
Rsolution(:,k) =3142000;
 end

for i = m: NumberOfSubinterval
   betat(:,i)=beta0./(1+(k1+K.*sqrt(epsilon).*dW(:,i)).*atan(Qsolution(:,i)-Qsolution(:,i-m+1)));
   for j=1:NumberOfSample
   if betat(j,i)>1
       betat(j,i)=1;
   elseif betat(j,i)<0
       betat(j,i)=0;
   end
   end
   Solution(:,i+1)=Solution(:,i)+Stepsize.*((1-eta).*omega+gamma1.*Vsolution(:,i)-betat(:,i).*Asolution(:,i).*Solution(:,i)./N-mu.*Solution(:,i));
   Qsolution(:,i+1)=Qsolution(:,i)+Stepsize.*((1-alpha).*varphi.*Esolution(:,i)-(gamma2+mu+delta).*Qsolution(:,i));
   Asolution(:,i+1)=Asolution(:,i)+Stepsize.*(alpha.*varphi.*Esolution(:,i)-(rho+mu+delta).*Asolution(:,i));
   Vsolution(:,i+1)=Vsolution(:,i)  +Stepsize.*(eta.*omega-betag.*Asolution(:,i).*Vsolution(:,i)./N-(gamma1+mu).*Vsolution(:,i));
   Esolution(:,i+1)=Esolution(:,i)+Stepsize.*(betat(:,i).*Asolution(:,i).*Solution(:,i)./N +betag.*Asolution(:,i).*Vsolution(:,i)./N-(varphi+mu).*Esolution(:,i));
 Rsolution(:,i+1)=Rsolution(:,i)+Stepsize.*(rho.*Asolution(:,i)+gamma2.*Qsolution(:,i)-mu.*Rsolution(:,i));
end
figure(3)
        plot(TimeGrid,mean(Qsolution)/N,'LineWidth',1);
              h1=xlabel('$t$','FontSize',15);
      h2=ylabel('$Q/N$','FontSize',15);
        %h3=legend('$k$=0.2','$k$=0.3','$k$=0.4','$k$=0.5');
        h3=legend('$\tau$=10','$\tau$=15','$\tau$=20');
      set(h3,'FontName','Times New Roman','FontSize',11,'box','off','Interpreter','latex');
      set(h1,'Interpreter','latex');
      set(h2,'Interpreter','latex');
      set(gca,'FontSize',17,'Fontname','Time New Roman');
      shading interp;
      hold on
     axis([4000,4500,0.00714,0.00717])
    