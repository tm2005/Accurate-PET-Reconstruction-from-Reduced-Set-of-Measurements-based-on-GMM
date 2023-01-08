function [x] = roots_3(p)
    p =p/p(1);

    if abs(p(4))<1e-12
        x(1) = -p(2)/2 + sqrt(p(2)^2-4*p(3))/2;
        x(2) = -p(2)/2 - sqrt(p(2)^2-4*p(3))/2;
        x(3) = 0;
    else
        a1 = 2*p(2)^3 - 9*p(2)*p(3) + 27* p(4);
        b1 = (p(2)^2 -3*p(3))^3;
        
        if abs(imag(-a1/2 + sqrt(a1^2-4*b1)/2)) < 1e-15
            x0 = (-p(2)  + nthroot(-a1/2 + sqrt(a1^2-4*b1)/2,3) + nthroot(-a1/2 - sqrt(a1^2-4*b1)/2,3) )/3;
        else
            x0 = (-p(2)  + (-a1/2 + sqrt(a1^2-4*b1)/2)^(1/3) + (-a1/2 - sqrt(a1^2-4*b1)/2)^(1/3) )/3;
        end
    
        a2 = 1*x0+p(2);
        b2 = a2*x0 + p(3);
        x(1) = -a2/2 + sqrt(a2^2-4*b2)/2;
        x(2) = -a2/2 - sqrt(a2^2-4*b2)/2;
        x(3) = x0;
    end
    x=x';
end