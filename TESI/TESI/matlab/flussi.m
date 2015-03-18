clear all
close all
clc

x=linspace(0,1);
flux=inline('10.*(x.^2)./(10.*x.^2+20.*(1-x).^2 )');
flux2=inline('50.*(x.^2)./(50.*x.^2+5.*(1-x).^2 )');

plot(x,flux(x))
hold on
plot(x,flux2(x),'g')

%%

x=linspace(0,1); 
flux3=inline('(x)./( x+0.5*(1-x) )');
flux4=inline('(x.^2)./(x.^2+0.5*(1-x).^2 )');
flux5=inline('4.*(x).*(1-x)./(x.^2+2.*(1-x).^2)');

figure()
subplot(1,2,1)
plot(x,flux4(x))
subplot(1,2,2)
plot(x,flux5(x))

%%


x=linspace(0,1);
flux5=inline('10.*(x.^2).*20.*(1-x).^2./(10.*x.^2+20.*(1-x).^2 )');
flux6=inline('50.*(x.^2).*5.*(1-x).^2 ./(50.*x.^2+5.*(1-x).^2 )');

plot(x,flux5(x))
hold on
%plot(x,flux6(x),'g')

%%

x=linspace(0,1);
flux=inline('(x.^2).*(1-x).^2./(x.^2+(1-x).^2 )');
flux2=inline('-2.*(x).*(2.*x.^4-5.*x.^3+6.*x.^2-4.*x+1)./(2.*x.^2-2.*x+1).^2');

subplot(1,2,1)
plot(x,flux(x))
subplot(1,2,2)
plot(x,flux2(x))
