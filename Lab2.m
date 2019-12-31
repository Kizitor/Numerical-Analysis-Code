function f = lab2()

z = [1 3 2 1; 4 2 1 2; 2 1 2 3; 1 2 4 1];
disp(z);

b = [-2;2;1;-1];

[Q,R] = classicalGS(z);

y = Q' *b;

g = backward(R, y);
disp('Ax=b, x:');
disp(g);

disp(y);
[Q2,R2] = modifiedGS(z);

y2 = Q2' * b;
g2 = backward(R2, y2);
disp('Ax=b, x:');
disp(g2);

m = dlmread('lab2mat.txt');
classicalGS(m);
modifiedGS(z);
end

function [y,p] = classicalGS(n)
Name = 'Classical GS';
disp(Name);
r = zeros(size(n,2), size(n,2));
v = zeros(size(n,2), size(n,2));
q = zeros(size(n,2), size(n,2));
for j = 1:size(n,2)
    v(:,j) = n(:,j);
    for i = 1:j-1
        r(i,j)= dot(q(i,:),n(:,j));
        
        v(:,j) = v(:,j) - (r(i,j)*q(i,:))';
    end
    r(j,j) = dot(v(:,j), v(:,j));
    q(:,j) = v(:,j)/r(j,j);
end

y = q;
S = 'Orthogonal Matrix';
disp(S);
disp(q);

p = r;
S = 'Upper Triangular Matrix';
disp(S);
disp(r);
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
        v(:,j) = v(:,j) - (r(i,j)*q(i,:))';
    end
    r(j,j) = dot(v(:,j), v(:,j));
    q(:,j) = v(:,j)/r(j,j);
end

y = q;
S = 'Orthogonal Matrix';
disp(S);
disp(q);

p = r;
S = 'Upper Triangular Matrix';
disp(S);
disp(r);
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