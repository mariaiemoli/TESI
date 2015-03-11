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

figure(1)
plot(x1,y1,'bo-')
hold on
plot(x2,y2,'ro-')


%%

%------------------------fine della geometria-------------------------------

%------------------tempo-----------------------------------------

Nstep = 50;
%dt=0.2*min([h1*d1,h2*d2,h3*d3]./abs([U1,U2,U3]));
dt = 0.2*h1*d1./abs(U1);

%--------------------cc--------------------

BC1=1;
BC2=1;




%---------condizioni iniziali---------------------------------

s1=zeros(length(x1)-1,Nstep);
s2=zeros(length(x2)-1,Nstep);

sA=ones(1,Nstep);  %confina con B attraverso PM2, con C attraverso PM1


for tt=1:Nstep-1

   %la due
   sL=BC2;%
   ff=[];
   for j=1:N2
       sR=s2(j,tt);
       [f]=U2*(sL*(U2>0)+sR*(U2<0));
       ff=[ff; f];
       sL=sR;
   end  
   sR=sA(1);
   [f]=U2*(sL*(U2>0)+sR*(U2<0));
   ff=[ff; f];
   Fluxes=ff(2:end)-ff(1:end-1);
   s2(:,tt+1)=s2(:,tt)-Fluxes*dt/h2/d2;
     f2A=-f;


  
end
%%
figure(3)
plot3(xx1,yy1,s1(:,end))
hold on
plot3(x1(end),y1(end),sA(end),'*')
%plot3(x2(end),y2(end),sB(end),'r*')
%plot3(x3(end),y3(end),sC(end),'g*')

sA(tt+1)

figure

plot(0.5*(t1(2:end)+t1(1:end-1)), s1(:,end))

figure

plot(0.5*(t2(2:end)+t2(1:end-1)), s2(:,end))