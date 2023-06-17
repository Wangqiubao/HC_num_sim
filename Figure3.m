%%%%Figure3. The time series diagram of S,E,A,V,Q,R with the parameter tau=10 <tau_0%%%
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
tau=10;
k1=0.003;
K=0;

%The explicit Euler Maruyama method for linear problem with a small%%%%%%%%%%%%%
%number  of samples %
randn('state',100);
Stepsize=0.01;
FinalTime=5000;
NumberOfSample=10;
InitialValue=1000;
DiffusionCoef=0.1;
DriftCoef=1;
dimension=2;
m=ceil(tau/Stepsize);
TimeGrid = 0 : Stepsize : FinalTime;
NumberOfSubinterval = length(TimeGrid)-1;
dW = sqrt(Stepsize) * randn(NumberOfSample,NumberOfSubinterval+1); 
Solution = ones(NumberOfSample,NumberOfSubinterval+1).*552030.54;
Esolution=ones(NumberOfSample,NumberOfSubinterval+1).*1405470.79;
Asolution = ones(NumberOfSample,NumberOfSubinterval+1).*87494.34;
Qsolution=ones(NumberOfSample,NumberOfSubinterval+1).*71666.04;
Vsolution = ones(NumberOfSample,NumberOfSubinterval+1).*4741153.57;
betat=ones(NumberOfSample,NumberOfSubinterval+1).*beta0;
Rsolution= ones(NumberOfSample,NumberOfSubinterval+1).*3142180;
for k=1:m

Solution(:,k) =552030.54*rand(1);
Esolution(:,k) =1405470.79*rand(1);
Asolution(:,k) =87494.34*rand(1);
Qsolution(:,k)=71666.046*rand(1);
Vsolution(:,k) =4741153.57*rand(1);
Rsolution(:,k) =3142180*rand(1);
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
figure(1)  
plot(TimeGrid,mean(Solution)/N,'b-','LineWidth',1) 
      h1=xlabel('$t$','FontSize',15);
      h2=ylabel('$S/N$','FontSize',15);
      set(h1,'Interpreter','latex');
      set(h2,'Interpreter','latex');
      set(gca,'FontSize',17,'Fontname','Time New Roman');
      shading interp;

figure(2)
        plot(TimeGrid,mean(Asolution)/N,'b-','LineWidth',1)
              h1=xlabel('$t$','FontSize',15);
      h2=ylabel('$A/N$','FontSize',15);
      set(h1,'Interpreter','latex');
      set(h2,'Interpreter','latex');
      set(gca,'FontSize',17,'Fontname','Time New Roman');
      shading interp;

figure(3)
        plot(TimeGrid,mean(Qsolution)/N,'b-','LineWidth',1);
              h1=xlabel('$t$','FontSize',15);
      h2=ylabel('$Q/N$','FontSize',15);
      set(h1,'Interpreter','latex');
      set(h2,'Interpreter','latex');
      set(gca,'FontSize',17,'Fontname','Time New Roman');
      shading interp;

figure(4)
        plot(TimeGrid,mean(Esolution)/N,'b-','LineWidth',1)
              h1=xlabel('$t$','FontSize',15);
      h2=ylabel('$E/N$','FontSize',15);
      set(h1,'Interpreter','latex');
      set(h2,'Interpreter','latex');
      set(gca,'FontSize',17,'Fontname','Time New Roman');
      shading interp;

figure(5)
        plot(TimeGrid,mean(Vsolution)/N,'b-','LineWidth',1)
              h1=xlabel('$t$','FontSize',15);
      h2=ylabel('$V/N$','FontSize',15);
      set(h1,'Interpreter','latex');
      set(h2,'Interpreter','latex');
      set(gca,'FontSize',17,'Fontname','Time New Roman');
      shading interp;

   figure(6) 
   plot(TimeGrid,mean(Rsolution)/N,'b-','LineWidth',1)
              h1=xlabel('$t$','FontSize',15);
      h2=ylabel('$R/N$','FontSize',15);
      set(h1,'Interpreter','latex');
      set(h2,'Interpreter','latex');
      set(gca,'FontSize',17,'Fontname','Time New Roman');
      shading interp;
