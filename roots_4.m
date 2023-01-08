function [x] = roots_4(p0)
    p0 = p0/p0(1);

    % x^4 + ax^3 + bx^2 + cx +d = 0
    a = p0(2);
    b = p0(3);
    c = p0(4);
    d = p0(5);
    
    % x^4 + px^2 + qx + r =0
    p = (b - (3*a^2)/8);
    q = (a^3/8 - (b*a)/2 + c);
    r = - (3*a^4)/256 + (b*a^2)/16 - (c*a)/4 + d;
    
    rez = roots_3([1 p p^2/4-r -q^2/8]);
    rez(abs(rez)<1e-15)=[];
    
    if isempty(rez) == 1
        x = [0 0 0 0];
    else
        alpha0 = rez(end);
        x(1) =  sqrt(2*alpha0)/2 + sqrt(-2*alpha0-2*p-q*2/sqrt(2*alpha0))/2;
        x(2) =  sqrt(2*alpha0)/2 - sqrt(-2*alpha0-2*p-q*2/sqrt(2*alpha0))/2;
        x(3) = -sqrt(2*alpha0)/2 + sqrt(-2*alpha0-2*p+q*2/sqrt(2*alpha0))/2;
        x(4) = -sqrt(2*alpha0)/2 - sqrt(-2*alpha0-2*p+q*2/sqrt(2*alpha0))/2;
    end
    
    x=x'- a/4;
end