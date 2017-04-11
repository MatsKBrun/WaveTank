close all; clearvars;
% Script to solve BBM with DirichletBC 
% using 6th order SBP-SAT
%
% u_t + a*u_x + b*u*u_x + c*u_txx = 0
% u(x=0) = g(t), u(x=L) = 0
%
% Data:
% h0 = 0.36m
% dt = 0.004s
% 
% Posisjon til målere:
% wg1 = 6.1m    
% wg2 = 10.6m   
% wg3 = 10.9m 
% wg4 = 11.2m 
%
% Data i tabell:
% time [s] | wg01 [m] | wg02 [m] | wg03 [m] | wg04 [m]
%
% Mats.Brun@uib.no

% import data
Data = readtable('f0p8-a1p1-int.dat');
Data = Data{:,{'x_Time','wg01', 'wg02','wg03','wg04'}};
Data(1:1000,:) = [];

% parameters
h0 = 0.36;                        % depth
c0 = sqrt(9.825*h0);              % char. wavespeed
L = 20;                           % length of domain
m = 321;                          % number of grid points
x = linspace(0,L,m)';             % spatial grid
h = x(2) - x(1);                  % grid size
N = 5000;                         % number of time steps
t = Data(:,1)-Data(1,1);          % time vector from Data
dt = 0.0001;                       % step size
tt = (0:dt/2:N*dt)';                % time vector for interpolation
alpha = 0.18; % or greater

% spline interpolation
wg1 = spline(t,Data(:,2),tt);
wg2 = spline(t,Data(:,3),tt);
wg3 = spline(t,Data(:,4),tt);
wg4 = spline(t,Data(:,5),tt);

% lin interpolation
%wg1 = interp1(t,Data(:,2),tt);
%wg2 = interp1(t,Data(:,3),tt);
%wg3 = interp1(t,Data(:,4),tt);
%wg4 = interp1(t,Data(:,5),tt);

% data (no interpolation)
%wg1 = Data(:,2); wg2 = Data(:,3); wg3 = Data(:,4); wg4 = Data(:,5);

% Import SBP operators
[ H,D1,D2,D3,S2_1,S2_m,S_1,S_m,e_1,e_m ] = SBP6(m,h);

% initial cond.
Hmax = abs(max(Data(:,2)));
u = zeros(m,1);

% initiate plot
fig1 = figure(1);
p = plot(x,u);
axis([0,L/2,-Hmax,Hmax]);
title('KdV vs data')
hold on
yi0 = wg1(1); xi0 = 0;
yi1 = wg2(1); xi1 = 4.5;
yi2 = wg3(1); xi2 = 4.8;
yi3 = wg4(1); xi3 = 5.1;
p1 = plot(xi0,yi0,'r*','YDataSource','yi0','XDataSource','xi0');
p2 = plot(xi1,yi1,'r*','YDataSource','yi1','XDataSource','xi1');
p3 = plot(xi2,yi2,'r*','YDataSource','yi2','XDataSource','xi2');
p4 = plot(xi3,yi3,'r*','YDataSource','yi3','XDataSource','xi3');
legend('KdV','Data')
drawnow;

% calculate derivative of boundary data
dg = @(t) (-1/60*wg1(t-3) + 3/20*wg1(t-2) - 3/4*wg1(t-1) + 3/4*wg1(t+1) - 3/20*wg1(t+2) + 1/60*wg1(t+3))/dt;

% coeffs
a = c0; b = 3/2*c0/h0; c = 1/6*h0^2;
ap = a + b*max(u(1),0);
an = a + b*min(u(end),0);

% operators
II = eye(m);

A = II - c*D2 + c*H\(S_1'*e_1') + c^2/(2*alpha*h)*H\(e_1*e_1') - c*H\(S_m'*e_m') + c^2/(2*alpha*h)*H\(e_m*e_m');

B = @(ap,an) -a*D1 - ap/2*H\(e_1*e_1') + an/2*H\(e_m*e_m');

C = -b/2*D1;

G = @(t,ap,an) ap/2*H\e_1*wg1(t) + c*H\S_1'*dg(t) + c^2/(2*alpha*h)*H\e_1*dg(t);

% first time step and iteration
for i = 6:2:N
    k1 = B(ap,an)*u + C*u.^2 + G(i,ap,an);
    utemp = u + k1*dt/2;
    
    k2 = B(ap,an)*utemp + C*utemp.^2 + G(i+1,ap,an);
    utemp = u + k2*dt/2;
    
    k3 = B(ap,an)*utemp + C*utemp.^2 + G(i+1,ap,an);
    utemp = u + k3*dt;
    
    k4 = B(ap,an)*utemp + C*utemp.^2 + G(i+2,ap,an);
    
    u = u + dt*A\(k1 + 2*k2 + 2*k3 + k4)/6;
    
    % update plot
    yi0 = wg1(i);
    yi1 = wg2(i);
    yi2 = wg3(i);
    yi3 = wg4(i);
    set(p,'YData',u);
    refreshdata(p1);
    refreshdata(p2);
    refreshdata(p3);
    refreshdata(p4);
    drawnow;
    
end

  





