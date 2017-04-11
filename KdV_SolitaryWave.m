close all; clear; clc;
% Script to solve KdV with Dirichlet/Neumann BC 
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
c = (Hmax + 2*h0)*c0/2/h0;         % wave speed
L1 = 3;                            % left boundary 
L2 = 13;                           % right boundary
m = 201;                           % number of grid points
x = linspace(L1,L2,m)';            % spatial grid
h = x(2) - x(1);                   % grid size
N = 4000;                          % number of time steps       
dt = 0.001;                        % step size
t = 0;                             % initial time

% Import SBP operators
[ H,D1,D2,D3,S2_1,S2_m,S_1,S_m,e_1,e_m ] = SBP6(m,h);

% Solitary wave
sol = @(x,t) Hmax*sech(sqrt(3*Hmax/h0)/2/h0*(x - c*t)).^2;
dsol = @(x,t) -Hmax*sqrt(3*Hmax/h0)/h0*sech(sqrt(3*Hmax/h0)/2/h0*(x - c*t)).^2*...
    tanh(sqrt(3*Hmax/h0)/2/h0*(x - c*t));

% boundary and initial cond.
u = sol(x,t);
g1 = @(t) sol(x(1),t);
g2 = @(t) sol(x(end),t);
g3 = @(t) dsol(x(end),t);

% initiate plot
fig1 = figure(1);
p = plot(x,u,'o');
axis([L1, L2, 0, Hmax]);
title('Solitary wave')
hold on
p2 = plot(x,sol(x,t),'r');
legend('Numerical','Exact');
drawnow;

% coeffs
a = c0; b = 3/2*c0/h0; c = 1/6*c0*h0^2;
ap = a + b*max(u(1),0);
an = a + b*min(u(end),0);

% operators
A = - c*D3 - H\(2*c*(S2_1'*e_1') -...
              2*c*(S2_m'*e_m') + (S_m'*S_m)*c)/2;
          
G = @(t,ap,an) H\((ap*e_1 + 2*c*S2_1')*g1(t) -...
               (an*e_m + 2*c*S2_m')*g2(t) +...
               S_m'*g3(t)*c)/2;
B = -b/2*D1;

E = @(ap,an) -a*D1 - H\(ap*(e_1*e_1') -...
              an*(e_m*e_m'))/2;
          
II = eye(m);

%first time step
t = dt;
utemp = u;
u = (II - dt/2*A)\( (II + dt/2*A)*u + dt*B*u.^2 + dt*E(ap,an)*u + dt*G(t,ap,an) );
t = t+dt;
Etemp = E(ap,an);
ap = a + b*max(u(1),0);
an = a + b*min(u(end),0);
for i = 3:N+1
   utemp2 = u;
   En = E(ap,an);
   u = (II - dt/2*A)\( (II + dt/2*A)*u + dt*(3*(B*u.^2 + En*u) - (B*utemp.^2 + Etemp*utemp))/2 + dt*G(t,ap,an));
   utemp = utemp2;
   Etemp = En; 
   
   % update plot
   set(p,'YData',u);
   set(p2,'YData',sol(x,t));
   drawnow;
    
   t = t + dt;
   ap = a + b*max(u(1),0);
   an = a + b*min(u(end),0);
end


l2_error = sqrt(sum((sol(x,t-dt) - u).^2)/m);
l2_error
T_final = t-dt









