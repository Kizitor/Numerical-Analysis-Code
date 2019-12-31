
%NAme: Kizitor Chukwuma Date: 12/6/19
function z = Lab6()
f = @(t, y) (sin(y)/(1 + t));
div = @(t, y) (sin(y)/((t + 1)*cos(y)));
[a,a1] = forwEuler(f, 0, 20, 500, 1);
%disp(a1);
[b, b1] = RK4(f, 0, 20, 500, 1);
%disp(b1);
[c, c1] = backEuler(f, 0, 20, 500, 1);
%disp(c1);
sol = @(t) (2*atan((t + 1)/cot(1/2)));
diff = 0;
N = 500;
T= 20;
t0 = 0;
h = (T-t0)/N;
for i = 1:N
    t = t0 + i*h;
    if(abs(sol(t) - a1(i)) > diff)
        diff = abs(sol(t) - a1(i));
    end
end
disp(diff);

diff = 0;
N = 500;
T= 20;
t0 = 0;
h = (T-t0)/N;
for i = 1:N
    t = t0 + i*h;
    if(abs(sol(t) - b1( i)) > diff)
        diff = abs(sol(t) - b1(i));
    end
end
disp(diff);

diff = 0;
N = 500;
T= 20;
t0 = 0;
h = (T-t0)/N;
for i = 1:N
    t = t0 + i*h;
    if(abs(sol(t) - c1( i)) > diff)
        diff = abs(sol(t) - c1(i));
    end
end
disp(diff);


[a,a1] = forwEuler(f, 0, 20, 1000, 1);
%disp(a1);
[b, b1] = RK4(f, 0, 20, 1000, 1);
%disp(b1);
[c, c1] = backEuler(f, 0, 20, 1000, 1);
%disp(c1);

diff = 0;
N = 1000;
T= 20;
t0 = 0;
h = (T-t0)/N;
h1 = h;
for i = 1:N
    t = t0 + i*h;
    if(abs(sol(t) - a1( i)) > diff)
        diff = abs(sol(t) - a1(i));
    end
end
e1 = diff;

diff = 0;
N = 1000;
T= 20;
t0 = 0;
h = (T-t0)/N;
h11 = h;
for i = 1:N
    t = t0 + i*h;
    if(abs(sol(t) - b1( i)) > diff)
        diff = abs(sol(t) - b1(i));
    end
end
e11 = diff;

diff = 0;
N = 1000;
T= 20;
t0 = 0;
h = (T-t0)/N;
h111 = h;
for i = 1:N
    t = t0 + i*h;
    if(abs(sol(t) - c1( i)) > diff)
        diff = abs(sol(t) - c1(i));
    end
end
e111 = diff;

[a,a1] = forwEuler(f, 0, 20, 2000, 1);
%disp(a1);
[b, b1] = RK4(f, 0, 20, 2000, 1);
%disp(b1);
[c, c1] = backEuler(f, 0, 20, 2000, 1);
%disp(c1);

diff = 0;
N = 2000;
T= 20;
t0 = 0;
h = (T-t0)/N;
h2 = h;
for i = 1:N
    t = t0 + i*h;
    if(abs(sol(t) - a1( i)) > diff)
        diff = abs(sol(t) - a1(i));
    end
end
e2 = diff;

diff = 0;
N = 2000;
T= 20;
t0 = 0;
h = (T-t0)/N;
h22 = h;
for i = 1:N
    t = t0 + i*h;
    if(abs(sol(t) - b1( i)) > diff)
        diff = abs(sol(t) - b1(i));
    end
end
e22 = diff;

diff = 0;
N = 2000;
T= 20;
t0 = 0;
h = (T-t0)/N;
h222 = h;
for i = 1:N
    t = t0 + i*h;
    if(abs(sol(t) - c1( i)) > diff)
        diff = abs(sol(t) - c1(i));
    end
end
e222 = diff;


r1= log(e1/e2)/log(h1/h2);
disp(r1);
r2 = log(e11/e22)/log(h11/h22);
disp(r2);
r3 = log(e111/e222)/log(h111/h222);
disp(r3);
% All the rate match the theoretical onvergence rate except for backward
% euler due to my implementation.


end
function [time, approx] = forwEuler(ODE, t0, T, N, y0)
time = zeros(N+1, 1);
approx = zeros(N+1, 1);
h = (T-t0)/N;
t= t0;
w = y0;

time(1) = t;
approx(1) = w;

for i = 1:N
    w = w + h*ODE(t, w);
    t = t0 + i*h;
    
    time(i + 1) = t;
    approx(i + 1) = w;
end

end

function [time, approx] = RK4(ODE, t0, T, N, y0)
time = zeros(N+1, 1);
approx = zeros(N+1, 1);
h = (T-t0)/N;
t= t0;
w = y0;

time(1) = t;
approx(1) = w;

for i = 1:N
    k1 = h*ODE(t, w);
    k2 = h*ODE(t + 0.5*h, w+0.5*k1);
    k3 = h*ODE(t + 0.5*h, w+0.5*k2);
    k4 = h*ODE(t + h, w + k3);
    
    w = w + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
    t = t + h;
    time(i + 1) = t;
    approx(i + 1) = w;
end
end

function [time, approx] = backEuler(ODE, t0, T, N, y0)
time = zeros(N+1, 1);
approx = zeros(N+1, 1);
h = (T-t0)/N;
t= t0;
w = y0;
ODE = @(t, y) (sin(y)/(1 + t));
time(1) = t;
approx(1) = w;

for i = 1:N
    t = t0 + i*h;
    w = fixedpoint(@(Y) (w+ h*(sin(Y)/(1 + t)) - Y), [0 20], w, 10^(-10), 100);
    time(i + 1) = t;
    approx(i + 1) = w;
end
end


function [ x ] = fixedpoint(g,I,y,tol,m)
a=I(1);b=I(2);
if(y<a | y>b)
    error('The starting iteration does not lie in I.')
end
x=y;
gx=g(y);
while(abs(x-gx)>tol & m>0)
    if(gx<a | gx>b)
        error('The point g(x) does not lie in I.')
    end
    y=x;
    x=g(y);
    m=m-1;
end
end