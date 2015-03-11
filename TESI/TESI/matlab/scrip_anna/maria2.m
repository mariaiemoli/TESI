close all
clear all 
clc

%punto di incontro

b=[0 0];

%frattura 1
a1=[-1 0];                          % punto finale
N1=10;                              % numero suddivisioni
t1=[0:1/N1:1]';                     
x1=t1*(b(1)-a1(1)) + a1(1);         % mesh da -1 a 0, tipo mesh velocità
y1=t1*(b(2)-a1(2)) + a1(2);         
xx1=0.5*(x1(1:end-1)+x1(2:end));    % mesh da -1+h/2 a 0-h/2
yy1=0.5*(y1(1:end-1)+y1(2:end));
h1=1/N1*((a1(1)-b(1))^2 +(a1(2)-b(2))^2 )^0.5;  % lunghezza elemento mesh
tau1=[b(1)-a1(1), b(2)-a1(2)];      % tangente
n1=[-b(2)+a1(2), b(1)-a1(1)];       
n1=n1/(n1(1)^2+n1(2)^2)^0.5;        % normale
d1=0.1;                             % spessore
U1=1;                               % velocità

%frattura 2
a2=[1 1];
N2=15;
t2=[0:1/N2:1]';
x2=t2*(b(1)-a2(1)) + a2(1);
y2=t2*(b(2)-a2(2)) + a2(2);
xx2=0.5*(x2(1:end-1)+x2(2:end));
yy2=0.5*(y2(1:end-1)+y2(2:end));

h2=1/N2*((a2(1)-b(1))^2 +(a2(2)-b(2))^2 )^0.5;
tau2=[b(1)-a2(1), b(2)-a2(2)];
n2=[-b(2)+a2(2), b(1)-a2(1)];
n2=n2/(n2(1)^2+n2(2)^2)^0.5;
d2=0.2;
U2=-0.9;

%frattura 3
a3=[1 -1];
N3=12;
t3=[0:1/N3:1]';
x3=t3*(b(1)-a3(1)) + a3(1);
y3=t3*(b(2)-a3(2)) + a3(2);
xx3=0.5*(x3(1:end-1)+x3(2:end));
yy3=0.5*(y3(1:end-1)+y3(2:end));

h3=1/N3*((a3(1)-b(1))^2 +(a3(2)-b(2))^2 )^0.5;
tau3=[b(1)-a3(1), b(2)-a3(2)];
n3=[-b(2)+a3(2), b(1)-a3(1)];
n3=n3/(n3(1)^2+n3(2)^2)^0.5;
d3=0.3;
U3=-0.1;

figure(1)
plot(x1,y1,'bo-')
hold on
plot(x2,y2,'ro-')
plot(x3,y3,'go-')



figure(2)
%calcolo punti del triangolo
tt=[tau1(1) -tau3(1);tau1(2) -tau3(2)]\(a3' + n3'*d3/2-a1'+n1'*d1/2);
P1=tt(1)*tau1 + a1-n1*d1/2;

tt=[tau1(1) -tau2(1);tau1(2) -tau2(2)]\(a2' - n2'*d2/2-a1'-n1'*d1/2);
P2=tt(1)*tau1 + a1+n1*d1/2;
tt=[tau2(1) -tau3(1);tau2(2) -tau3(2)]\(a3' - n3'*d3/2-a2'-n2'*d2/2);
P3=tt(1)*tau2 + a2+n2*d2/2;
AA=norm(cross([P3-P2, 0], [P3-P1,0]))/2;        % area triangolo
AA1=norm(cross([b-P2, 0], [b-P1,0]))/2;
AA2=norm(cross([b-P2, 0], [b-P3,0]))/2;
AA3=norm(cross([b-P3, 0], [b-P1,0]))/2;


line([P1(1), P2(1)],[P1(2) P2(2)])
hold on
line([P1(1), P3(1)],[P1(2) P3(2)])
line([P3(1), P2(1)],[P3(2) P2(2)])

line([P1(1), b(1)],[P1(2) b(2)],'color','r')
line([P2(1), b(1)],[P2(2) b(2)],'color','r')
line([P3(1), b(1)],[P3(2) b(2)],'color','r')

PM1=(P1+b)/2;
PM2=(P2+b)/2;
PM3=(P3+b)/2;
Upm1=(-U1*(PM1-P3)-U2*(PM1-P1)-U3*(PM1-P2))/2/AA;
Nm1C=cross([b-P1,0], [0 0 1]);
Nm1C=Nm1C/norm(Nm1C);
Nm1B=cross(-[b-P2,0], [0 0 1]);
Nm1B=Nm1B/norm(Nm1B);

quiver(PM1(1), PM1(2),Upm1(1),Upm1(2),0.1)

Upm2=(-U1*(PM2-P3)-U2*(PM2-P1)-U3*(PM2-P2))/2/AA;
Nm2A=-Nm1B;
Nm2C=cross([P3-b,0], [0 0 1]);
Nm2C=Nm2C/norm(Nm2C);
quiver(PM2(1), PM2(2),Upm2(1),Upm2(2),0.1)
Upm3=(-U1*(PM3-P3)-U2*(PM3-P1)-U3*(PM3-P2))/2/AA;
Nm3A=-Nm1C;
Nm3B=-Nm2C;
quiver(PM3(1), PM3(2),Upm3(1),Upm3(2),0.1)



%------------------------fine della geometria-------------------------------

%------------------tempo-----------------------------------------

Nstep = 51; %50;
%dt=0.2*min([h1*d1,h2*d2,h3*d3]./abs([U1,U2,U3]))
dt=0.2*min([h1,h2,h3]./abs([U1,U2,U3]))

%--------------------cc--------------------

BC1=1;
BC2=1;
BC3=0; 


%---------condizioni iniziali---------------------------------

s1=zeros(length(x1)-1,Nstep);
s2=zeros(length(x2)-1,Nstep);
s3=zeros(length(x3)-1,Nstep);

sA=zeros(1,Nstep);  %confina con B attraverso PM2, con C attraverso PM1
sB=zeros(1,Nstep);  %confina con C attraverso PM3
sC=zeros(1,Nstep);


for tt=1:Nstep-1

   %la uno
   sL=BC1;%
   ff=[];
   for j=1:N1
       sR=s1(j,tt);
       [f]=U1*(sL*(U1>0)+sR*(U1<0));
       ff=[ff; f];  % aggiunge ogni volta un elemento
       sL=sR;
   end  
   sR=sA(tt);
   [f]=U1*(sL*(U1>0)+sR*(U1<0));
   ff=[ff; f];
   Fluxes=ff(2:end)-ff(1:end-1);
   s1(:,tt+1)=s1(:,tt)-Fluxes*dt/h1; %/d1;
     f1A=-f;

  
   %la due
   sL=BC2;%
   ff=[];
   for j=1:N2
       sR=s2(j,tt);
       [f]=U2*(sL*(U2>0)+sR*(U2<0));
       ff=[ff; f];
       sL=sR;
   end  
   sR=sA(tt);
   [f]=U2*(sL*(U2>0)+sR*(U2<0));
   ff=[ff; f];
   Fluxes=ff(2:end)-ff(1:end-1);
   s2(:,tt+1)=s2(:,tt)-Fluxes*dt/h2; %/d2;
     f2A=-f;
    %la tre
   sL=BC3;%
   ff=[];
   for j=1:N3
       sR=s3(j,tt);
       [f]=U3*(sL*(U3>0)+sR*(U3<0));
       ff=[ff; f];
       sL=sR;
   end  
   sR=sA(tt);
   [f]=U3*(sL*(U3>0)+sR*(U3<0));
   ff=[ff; f];
   Fluxes=ff(2:end)-ff(1:end-1);
   s3(:,tt+1)=s3(:,tt)-Fluxes*dt/h3; %/d3;
        f3A=-f;
     %bilancio su A
  
  fnet=f1A+f2A+f3A;
  sA(tt+1)=sA(tt)-dt*fnet/AA;

  

  
end

figure(3)
plot3(xx1,yy1,s1(:,end))
hold on
plot3(x1(end),y1(end),sA(end),'*')
%plot3(x2(end),y2(end),sB(end),'r*')
%plot3(x3(end),y3(end),sC(end),'g*')
plot3(xx2,yy2,s2(:,end))
plot3(xx3,yy3,s3(:,end))

sA(tt+1);

figure
subplot(3,1,1)
plot(0.5*(t1(2:end)+t1(1:end-1)), s1(:,end))
subplot(3,1,2)
plot(0.5*(t2(2:end)+t2(1:end-1)), s2(:,end))
subplot(3,1,3)
plot(0.5*(t3(2:end)+t3(1:end-1)), s3(:,end))
