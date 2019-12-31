function p = Lab4()
    y= compTrap((@(x) x*sin(x)), 0, pi/2,18);
    
    vecT = zeros(10,1);
    vecS = zeros(10,1);
    for i = 1:10
        vecT(i) = compTrap((@(x) exp(1/2 * x)), 0, 1, (2^(i + 2)));
        vecS(i) = compSimp((@(x) exp(1/2 * x)), 0, 1, (2^(i + 2)));
    end
    act = 2 * (exp(1)^(1/2)) - 2;
    errT = zeros(10,1);
    errS = zeros(10,1);
    
    for i = 1:10
        errT(i) = act - vecT(i);
        errS(i) = act - vecS(i);
    end
    
    rT = zeros(10,1);
    rS = zeros(10,1);
    
    for i = 2:10
        rT(i) = log(errT(i)/errT(i - 1))/log(1/2);
        rS(i) = log(errS(i)/errS(i - 1))/log(1/2);
    end
    for i = 2:10
        fprintf(1,'i=%6f --Trap: Error=%6.15f Rate=%6.15f -- Simp: Error=%6.15f Rate=%6.15f\n', i, errT(i), rT(i), errS(i), rS(i));
        %The rate of convergence starts to significantly decrease with
        %small h.
    end
    
    
    
    for i = 1:10
        vecT(i) = compTrap((@(x) x^(9/4)), 0, 1, (2^(i + 2)));
        vecS(i) = compSimp((@(x) x^(9/4)), 0, 1, (2^(i + 2)));
    end
    act = 4/13;
    errT = zeros(10,1);
    errS = zeros(10,1);
    
    for i = 1:10
        errT(i) = act - vecT(i);
        errS(i) = act - vecS(i);
    end
    
    rT = zeros(10,1);
    rS = zeros(10,1);
    
    for i = 2:10
        rT(i) = log(errT(i)/errT(i - 1))/log(1/2);
        rS(i) = log(errS(i)/errS(i - 1))/log(1/2);
    end
    for i = 2:10
        fprintf(1,'i=%6f --Trap: Error=%6.15f Rate=%6.15f -- Simp: Error=%6.15f Rate=%6.15f\n', i, errT(i), rT(i), errS(i), rS(i));
        %The rate of convergence starts to significantly decrease with
        %small h.
    end
    
    x = compSimp((@(x) exp(6/7 * x) * sin(4*x^3 + 2*x^2)), -1, 0,2000);
    disp(x);
end






function f = compTrap(fun, a, b, N)
    h = (b-a)/N;
    
    f = h/2 * (fun(a) + fun(b));
    
    sum = 0;
    x0 = a;
    for i = 1:N
        x = x0 + i*h;
        sum = sum + fun(x);
    end
    
    sum = sum *h;
    
    f = f + sum;
end

function g = compSimp(fun, a, b, N)
    if (mod(N, 2) == 0)
       h = (b-a)/N; 
       g = fun(a) + fun(b);
       for i=1:2:N-1
           g = g + 4*fun(a+i*h);
       end
       for i = 2:2:N-2
           g = g +2*fun(a+i*h);
       end
       g = h/3 *g;
    end
end
