%% RISULTATI KAASSCHIETER
clear all
close all
clc

nx= 80;

A=-0.4;
B=0.4;

h = (B-A)/nx;

x=linspace(A+h/2,B-h/2,nx-1);

yy = linspace (0,1,nx);

%% 1a T=1

figure()

% cfl = 1.4
subplot(2,2,1)

v1 = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.8603e-250, 3.9824e-96, 3.08268e-37, 1.0618e-14, 4.49497e-06, 0.00876535, 0.132462, 0.283396, 0.342814, 0.368552, 0.382969, 0.392187, 0.398306, 0.402246, 0.404581, 0.405797, 0.406332, 0.406523, 0.406576, 0.406587, 0.406588, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25 ];
plot (x, v1,'Linewidth',2)
hold on
plot (0,yy,'--', 'Linewidth',2)

axis([-0.4 0.4 0 1])

%% 1b

subplot(2,2,2)

v1 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.49966, 0.471618, 0.373809, 0.317627, 0.302541, 0.29906, 0.298304, 0.298153, 0.298128, 0.298125, 0.298125, 0.298125, 0.298125, 0.298125, 0.298125, 0.298125, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 ];

plot (x, v1,'Linewidth',2)
hold on
plot (0,yy,'--', 'Linewidth',2)
axis([-0.4 0.4 0 1])

%% 2a

subplot(2,2,3)

v1 = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.05027e-194, 1.5924e-74, 7.08797e-29, 1.82718e-11, 8.01607e-05, 0.0259373, 0.186441, 0.310245, 0.358137, 0.38249, 0.398654, 0.411177, 0.421776, 0.431234, 0.439959, 0.448195, 0.4561, 0.463799, 0.471409, 0.479083, 0.487158, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734111, 0.734193, 0.735652, 0.760315, 0.916598, 0.99923, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ];

plot (x, v1,'Linewidth',2)
hold on
plot (0,yy,'--', 'Linewidth',2)
axis([-0.4 0.4 0 1])

%% 2b

subplot(2,2,4)

v1 = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.66008e-200, 7.83563e-77, 9.35283e-30, 8.44326e-12, 5.97336e-05, 0.0232844, 0.180643, 0.307921, 0.35715, 0.381927, 0.398269, 0.410887, 0.421546, 0.431046, 0.439804, 0.448066, 0.455995, 0.463715, 0.471343, 0.479036, 0.487129, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734107, 0.734106, 0.734085, 0.729459, 0.507816, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4 ];

plot (x, v1,'Linewidth',2)
hold on
plot (0,yy,'--', 'Linewidth',2)
axis([-0.4 0.4 0 1])


%% 3a

subplot(2,2,1)

v1 = [ 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.861025, 0.861025, 0.861025, 0.861025, 0.860986, 0.842727, 0.298048, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25 ];

plot (x, v1,'Linewidth',2)
hold on
plot (0,yy,'--', 'Linewidth',2)
axis([-0.4 0.4 0 1])

%% 3b
subplot(2,2,2)

v1 = [ 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.580262, 0.296315, 0.224091, 0.217499, 0.216829, 0.216763, 0.216757, 0.216756, 0.216756, 0.216756, 0.216756, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15 ];

plot (x, v1,'Linewidth',2)
hold on
plot (0,yy,'--', 'Linewidth',2)
axis([-0.4 0.4 0 1])


%% 4a

subplot(2,2,3)

v1 = [ 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752228, 0.752337, 0.754395, 0.788326, 0.943756, 0.999729, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ];

plot (x, v1,'Linewidth',2)
hold on
plot (0,yy,'--', 'Linewidth',2)
axis([-0.4 0.4 0 1])


%% 4b

subplot(2,2,4)

v1 = [ 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.752222, 0.7522, 0.750137, 0.588007, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4 ];

plot (x, v1,'Linewidth',2)
hold on
plot (0,yy,'--', 'Linewidth',2)
axis([-0.4 0.4 0 1])

