function f = lab3()
m = dlmread('lab3mat.txt');
b = dlmread('lab3b.txt');
%m = [29 0 8 -2; 0 18 0 0; 8 0 23 10; -2 0 10 20];
%b = [-5; 36; 61; 32];


x= zeros(1,1500);
x = transpose(x);

%Second time is with modded mymult
% Total number of iterations: 15000 Time: 66.4253 seconds/4.7270 seconds
[s, n, i] = steepestDesc(x, m, b, .00000000000001, 15000);

% Total number of iterations: 1500 Time: 3.3732 seconds/ 0.0429 seconds
[s2, n2, i2]  = conjGrad(x, m, b, .00000000000001, 1500);
%disp(i);
%disp(i2);
%disp(s);
%tic
%[Q2,R2] = modifiedGS(m);

%y2 = Q2' * b;
%g2 = backward(R2, y2); %total is 45.7899 seconds
%t = toc;
%disp(t);
end

function f = mymult(A,x)
N = 1500;
f= 2.01*[x(1:N)] - [x(2:N);x(1)] - [x(N);x(1:N-1)];
end

function [sol, norm, iter]  = steepestDesc(x, A, b, tol, maxit)

tic
k = 0;
while( k < maxit )
    k = k + 1;
    r = b - mymult(A, x);
    d = r;
    alpha = dot(r,r)/dot(r, mymult(A, r));
    if (alpha <= tol)
        break;
    end
    x = x + alpha * r;
end
sol = x;
norm = vecnorm(r);
iter = k;
t = toc;
disp(t);
end

function [sol, norm, iter] = conjGrad(x, A, b, tol, maxit)

tic
r = b - mymult(A, x);
d = r;

k = 0;
while( k < maxit )
    k = k + 1;
    alpha = dot(r,r)/dot(d, mymult(A, d));
    if (alpha <= tol)
        break;
    end
    x = x + alpha*d;
    rNew = r - alpha*mymult(A, d);
    
    s = dot(rNew, rNew)/dot(r, r);
    d = r + s*d;
    
    r = rNew;
    k = k + 1;
end
sol = x;
norm = vecnorm(r);
iter = k;
t= toc;
disp(t)
end

function [y,p] = modifiedGS(n)
Name = 'Modified GS';
disp(Name);
r = zeros(size(n,2), size(n,2));
v = zeros(size(n,2), size(n,2));
q = zeros(size(n,2), size(n,2));
for j = 1:size(n,2)
    v(:,j) = n(:,j);
    for i = 1:j-1
        r(i,j)= dot(q(i,:),n(:,j));
        v(:,j) = v(:,j) - dot(v(:,j),v(:,j-1))*v(:,j-1);
        v(:,j) = v(:,j) - transpose((r(i,j)*q(i,:)));
    end
    r(j,j) = dot(v(:,j), v(:,j));
    q(:,j) = v(:,j)/r(j,j);
end

y = q;


p = r;

end

function[x]=backward(U,b)
disp('Backsub: ');
S=size(U);
m=S(1);
x=zeros(1,m);
x(1,m)=b(end)./U(m,m);
%bacward substitution
for k=m-1:-1:1
   
  
        x1=1/U(k,k).*(b(k)-sum(U(k,k+1:end).*x(k+1:end)));
        x(k)=x1;
end
x=x';
end
