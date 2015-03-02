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

figure()
plot(x,flux3(x))
hold on
plot(x,flux4(x),'g')

%%


x=linspace(0,1);
flux5=inline('10.*(x.^2).*20.*(1-x).^2./(10.*x.^2+20.*(1-x).^2 )');
flux6=inline('50.*(x.^2).*5.*(1-x).^2 ./(50.*x.^2+5.*(1-x).^2 )');

plot(x,flux5(x))
hold on
%plot(x,flux6(x),'g')