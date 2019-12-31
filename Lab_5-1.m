function f = Lab_5()
    f = myBisect((@(x) log(x) - 2), 6, 8, 10^-7, 200);
    disp(f); % Took 25 iterations
    f = myNewton((@(x) log(x) - 2), (@(x) 1/x), 6, 10^-7, 200);
    disp(f); % Took 5 iterations
    f = myNewton((@(x)tan(x) - x), (@(x) (sec(x))^2), 100, 10^-7, 20);
    disp(f);
    f = myNewton((@(x)tan(x) - x), (@(x) (sec(x))^2), 101, 10^-7, 20);
    disp(f);
    f = myNewton((@(x)tan(x) - x), (@(x) (sec(x))^2), 102, 10^-7, 20);
    disp(f);
    %Max iterations are being schiebved
    %k = -2
end

function f = myBisect(fun, a, b, tol, maxit)
    it = 0;
    yo = -1;
    f = -1;
    if (fun(a)*fun(b) >= 0)
        f = -1;
        fprintf("Bad end points");
    end
        while ((b-a) >= tol && it <= maxit || yo ~= 0)
            m = (a+b)/2;
            yo = fun(m);
            if (fun(m) == 0)
                f = m;
            elseif (fun(a)*fun(m) < 0)
                b = m;
                it = it + 1;
            else
                a = m;
                it = it + 1;
            end
        end
        
        if (f ~= m)
            fprintf(1, 'Expected end reached\n');
             err = b-a;
            if ((b-a) < tol)
                fprintf(1, 'Max tolerance reached: %6.15f Iterations: %6f\n', err, it);
            elseif (it > maxit)
                fprintf(1, 'Max iterations reached: %6.15f Tolerance: %6f\n', it, err);
            end
        end
end

function n = myNewton(fun, div, x0, tol, maxit)
    x1 = -1;
    it = 1;
    x1 = x0 - ((fun(x0)/div(x0)));
    yo = fun(x1);
    
    while (abs(x1 - x0) >= tol && it <= maxit)
        x0 = x1;
        x1 = x0 - ((fun(x0)/div(x0)));
        yo = fun(x1);
        it = it + 1;
    end
    n = x1;
end