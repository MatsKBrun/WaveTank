close all; clear; clc;
% Script to solve BBM with Dirichlet BC 
% using 6th order SBP-SAT
%
% u_t + a*u_x + b*u*u_x + c*u_xxx = 0
% u(x=0) = g1(t), u(x=L) = g2(t), u_x(x=L) = g3(t)
% Solitary wave as initial data
%
% Mats.Brun@uib.no

% parameters
h0 = 0.36;                         % depth
c0 = sqrt(9.825*h0);               % char. wavespeed
Hmax = 0.08;                       % wave height
L = 10;                            % length of domain
m = 101;                           % number of grid points
x = linspace(-L,L,m)';             % spatial grid
h = x(2) - x(1);                   % grid size
N = 4000;                          % number of time steps       
dt = 0.0005;                        % step size
t = 0;                             % initial time
alpha = 0.18;

% coeffs
a = c0; b = 3/2*c0/h0; d = 1/6*c0*h0^2;
c = b*Hmax/3 + a;                  % intrinsic wave speed

% Import SBP operators
[ H,D1,D2,D3,S2_1,S2_m,S_1,S_m,e_1,e_m ] = SBP6(m,h);

% Solitary wave
sol = @(x,t) Hmax*sech(1/2*sqrt((c-a)/d/c)*(x - c*t)).^2;
dtsol = @(x,t) Hmax*c*sqrt((c-a)/d/c)*sech(1/2*sqrt((c-a)/d/c)*(x - c*t)).^2*tanh(1/2*sqrt((c-a)/d/c)*(x - c*t));

% boundary and initial cond.
u = sol(x,t);
g1 = @(t) sol(x(1),t);
g2 = @(t) sol(x(end),t);
dg1 = @(t) dtsol(x(1),t);
dg2 = @(t) dtsol(x(end),t);

% initiate plot
fig1 = figure(1);
p = plot(x,u,'o');
axis([-L, L, 0, Hmax]);
title('Solitary wave')
hold on
p2 = plot(x,sol(x,t),'r');
legend('Numerical','Exact');
drawnow;

% variable coeffs
ap = a + b*max(u(1),0);
an = a + b*min(u(end),0);


% operators
II = eye(m);

A = II - d*D2 + b*H\(S_1'*e_1') + -b^2/(2*alpha*h)*H\(e_1*e_1') - b*H\(S_m'*e_m') + b^2/(2*alpha*h)*H\(e_m*e_m');

B = @(ap,an) -a*D1 - ap/2*H\(e_1*e_1') + an/2*H\(e_m*e_m');

C = -b/2*D1;

G = @(t,ap,an) ap/2*H\e_1*g1(t) - an/2*H\e_m*g2(t) + b*H\S_1'*dg1(t) + b^2/(2*alpha*h)*H\e_1*dg1(t) -...
    b*H\S_m'*dg2(t) + b^2/(2*alpha*h)*H\e_m*dg2(t);

for i = 1:N
   t1 = t;
   k1 = B(ap,an)*u + C*u.^2 + G(t1,ap,an);
   utemp = u + k1*dt/2;
   
   t2 = t + dt/2;
   k2 = B(ap,an)*utemp + C*utemp.^2 + G(t2,ap,an);
   utemp = u + k2*dt/2;
   
   t3 = t + dt/2;
   k3 = B(ap,an)*utemp + C*utemp.^2 + G(t3,ap,an);
   utemp = u + k3*dt;
   
   t4 = t + dt;
   k4 = B(ap,an)*utemp + C*utemp.^2 + G(t4,ap,an);
   
   u = u + dt*A\(k1 + 2*k2 + 2*k3 + k4)/6;
   

   
   % update plot
   set(p,'YData',u);
   set(p2,'YData',sol(x,t));
   drawnow;
    
   ap = a + b*max(u(1),0);
   an = a + b*min(u(end),0);
end


l2_error = sqrt(sum((sol(x,t-dt) - u).^2)/m);
l2_error
T_final = t-dt


