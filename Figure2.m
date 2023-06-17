%%%Figure2. P(D) bifurcation diagram of the system%%%
clear clc
global mu1 mu2 mu3 c
mu1=-0.2;mu2=-1;mu3=1;%%(a)
%mu1=0.999;mu2=-1;mu3=1;%%(b)
%mu1=2;mu2=-1;mu3=1;%%(c)
zz=10;
z=-zz:0.01:zz;
[y1,y2]=meshgrid(z);
A=sqrt(y1.^2+y2.^2);
c=1;
p=c./mu3.*A.^(2.*mu1./mu3-2).*exp(mu2.*A.^2./mu3);
pst=real(p);
 meshz(y1,y2,pst)
    subplot(3,1,[1,2])
    meshz(y1,y2,pst)
    view(-26,20)
    h1=xlabel('$y_1$','FontSize',20);
  h2=ylabel('$y_2$','FontSize',20);
 h3=zlabel('$p_{st}$','FontSize',20);
 set(gca,'FontSize',20,'Fontname', 'Times New Roman');
  set(h1,'Interpreter','latex');
  set(h2,'Interpreter','latex');
  set(h3,'Interpreter','latex');
    subplot(3,1,3)
  contour(y1,y2,pst,30)
    view(-26,90)
