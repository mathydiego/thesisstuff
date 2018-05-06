function eS = skewsymexpm(S)
% SKEWSYMEXPM(S) computes the exponential of a real skew-symmetric matrix 
% up to size 9 using the formulas from Diego Avalos's thesis.
% See also expm.
[n, m] = size(S);
if n ~= m
    msg1 = 'Input must be a square real skew-symmetric matrix.';
    error(msg1);
    quit;
elseif S+S' ~= zeros(n)
    msg2 = 'Input must be a real skew-symmetric matrix.';
    error(msg2);
    quit;
elseif n > 9
    msg3 = 'Input must be a real skew-symmetric matrix of size 9 or less.';
    error(msg3);
    quit;
end

coeff = charpoly(S);

if S == zeros(n)
    eS = eye(n);
    
elseif n == 2
    a = coeff(3); th = sqrt(a);
    c0 = cos(th); c1 = sin(th)/th;
    eS = c0*eye(2) + c1*S;
    
elseif n == 3
    a = coeff(3); th = sqrt(a);
    c1 = sin(th)/th; c2 = (1-cos(th))/th^2;
    eS = eye(3) + c1*S + c2*S^2;
    
elseif n == 4
    a = coeff(3); b = coeff(5);
    k = unique(abs(roots([1, a, b])));
    if length(k) == 1
        th = sqrt(a/2); 
        c0 = cos(th); c1 = sin(th)/th;
        eS = c0*eye(4) + c1*S;
    elseif length(k) == 2 & b == 0
        th = sqrt(a);
        c1 = sin(th)/th; c2 = (1-cos(th))/th^2;
        eS = eye(4) + c1*S + c2*S^2;
    else
        th1 = sqrt((a-sqrt(a^2-4*b))/2); th2 = sqrt((a+sqrt(a^2-4*b))/2);       
        c0 = (th2^2*cos(th1)-th1^2*cos(th2))/(th2^2-th1^2);
        c1 = (th2^3*sin(th1)-th1^3*sin(th2))/(th1*th2*(th2^2-th1^2));
        c2 = (cos(th1)-cos(th2))/(th2^2-th1^2);
        c3 = (th2*sin(th1)-th1*sin(th2))/(th1*th2*(th2^2-th1^2));
        eS = c0*eye(4) + c1*S + c2*S^2 + c3*S^3;
    end
    
elseif n == 5
    a = coeff(3); b = coeff(5);
    k = unique(abs(roots([1, a, b]))); k = k(k~=0); e = imag(eig(S));
    if length(k) == 1
        m = length(e(e~=0))/2;
        th = sqrt(a/m);
        c1 = sin(th)/th; c2 = (1-cos(th))/th^2;
        eS = eye(5) + c1*S + c2*S^2;
    else 
        th1 = sqrt((a-sqrt(a^2-4*b))/2); th2 = sqrt((a+sqrt(a^2-4*b))/2);         
        c1 = (th2^3*sin(th1)-th1^3*sin(th2))/(th1*th2*(th2^2-th1^2));
        c2 = (th2^4*(1-cos(th1))-th1^4*(1-cos(th2)))/(th1^2*th2^2*(th2^2-th1^2));
        c3 = (th2*sin(th1)-th1*sin(th2))/(th1*th2*(th2^2-th1^2));
        c4 = (th2^2*(1-cos(th1))-th1^2*(1-cos(th2)))/(th1^2*th2^2*(th2^2-th1^2));
        eS = eye(5) + c1*S + c2*S^2 + c3*S^3 + c4*S^4;
    end
    
elseif n == 6
    a = coeff(3); b = coeff(5); c = coeff(7);
    % This code is complete for all subcases where there are distinct eigenvalues.
    % The remaining subcases for size 6 and above are to be included.
    p = 2*a^3-9*a*b+27*c; q = a^2-3*b; r = sqrt(q); phi = acos(0.5*p/q^1.5);
    th1 = sqrt((1/3)*(a-r*cos(phi/3)-sqrt(3)*r*sin(phi/3)));    
    th2 = sqrt((1/3)*(a-r*cos(phi/3)+sqrt(3)*r*sin(phi/3)));
    th3 = sqrt((1/3)*(a+2*r*cos(phi/3)));
    c0 = th2^2*th3^2*(th3^2-th2^2)*cos(th1)-th1^2*th3^2*(th3^2-th1^2)*cos(th2)+...
        th1^2*th2^2*(th2^2-th1^2)*cos(th3);
    c0 = c0/((th3^2-th1^2)*(th3^2-th2^2)*(th2^2-th1^2));
    c1 = th2^3*th3^3*(th3^2-th2^2)*sin(th1)-th1^3*th3^3*(th3^2-th1^2)*sin(th2)+...
        th1^3*th2^3*(th2^2-th1^2)*sin(th3);
    c1 = c1/(th1*th2*th3*(th3^2-th1^2)*(th3^2-th2^2)*(th2^2-th1^2));
    c2 = (th3^4-th2^4)*cos(th1)-(th3^4-th1^4)*cos(th2)+(th2^4-th1^4)*cos(th3);
    c2 = c2/((th3^2-th1^2)*(th3^2-th2^2)*(th2^2-th1^2));
    c3 = th2*th3*(th3^4-th2^4)*sin(th1)-th1*th3*(th3^4-th1^4)*sin(th2)+...
        th1*th2*(th2^4-th1^4)*sin(th3);
    c3 = c3/(th1*th2*th3*(th3^2-th1^2)*(th3^2-th2^2)*(th2^2-th1^2));
    c4 = (th3^2-th2^2)*cos(th1)-(th3^2-th1^2)*cos(th2)+(th2^2-th1^2)*cos(th3);
    c4 = c4/((th3^2-th1^2)*(th3^2-th2^2)*(th2^2-th1^2));
    c5 = th2*th3*(th3^2-th2^2)*sin(th1)-th1*th3*(th3^2-th1^2)*sin(th2)+...
        th1*th2*(th2^2-th1^2)*sin(th3);
    c5 = c5/(th1*th2*th3*(th3^2-th1^2)*(th3^2-th2^2)*(th2^2-th1^2));
    eS = c0*eye(6) + c1*S + c2*S^2 + c3*S^3 + c4*S^4 + c5*S^5;
    
elseif n == 7
    a = coeff(3); b = coeff(5); c = coeff(7);    
    p = 2*a^3-9*a*b+27*c; q = a^2-3*b; r = sqrt(q); phi = acos(0.5*p/q^1.5);
    th1 = sqrt((1/3)*(a-r*cos(phi/3)-sqrt(3)*r*sin(phi/3)));    
    th2 = sqrt((1/3)*(a-r*cos(phi/3)+sqrt(3)*r*sin(phi/3)));
    th3 = sqrt((1/3)*(a+2*r*cos(phi/3)));    
    c1 = th2^3*th3^3*(th3^2-th2^2)*sin(th1)-th1^3*th3^3*(th3^2-th1^2)*sin(th2)+...
        th1^3*th2^3*(th2^2-th1^2)*sin(th3);
    c1 = c1/(th1*th2*th3*(th3^2-th1^2)*(th3^2-th2^2)*(th2^2-th1^2));
    c2 = th2^4*th3^4*(th3^2-th2^2)*(1-cos(th1))-th1^4*th3^4*(th3^2-th1^2)*(1-cos(th2))+...
        th1^4*th2^4*(th2^2-th1^2)*(1-cos(th3));
    c2 = c2/(th1^2*th2^2*th3^2*(th3^2-th1^2)*(th3^2-th2^2)*(th2^2-th1^2));
    c3 = th2*th3*(th3^4-th2^4)*sin(th1)-th1*th3*(th3^4-th1^4)*sin(th2)+...
        th1*th2*(th2^4-th1^4)*sin(th3);
    c3 = c3/(th1*th2*th3*(th3^2-th1^2)*(th3^2-th2^2)*(th2^2-th1^2));
    c4 = th2^2*th3^2*(th3^4-th2^4)*(1-cos(th1))-th1^2*th3^2*(th3^4-th1^4)*(1-cos(th2))+...
        th1^2*th2^2*(th2^4-th1^4)*(1-cos(th3));
    c4 = c4/(th1^2*th2^2*th3^2*(th3^2-th1^2)*(th3^2-th2^2)*(th2^2-th1^2));
    c5 = th2*th3*(th3^2-th2^2)*sin(th1)-th1*th3*(th3^2-th1^2)*sin(th2)+...
        th1*th2*(th2^2-th1^2)*sin(th3);
    c5 = c5/(th1*th2*th3*(th3^2-th1^2)*(th3^2-th2^2)*(th2^2-th1^2));
    c6 = th2^2*th3^2*(th3^2-th2^2)*(1-cos(th1))-th1^2*th3^2*(th3^2-th1^2)*(1-cos(th2))+...
        th1^2*th2^2*(th2^2-th1^2)*(1-cos(th3));
    c6 = c6/(th1^2*th2^2*th3^2*(th3^2-th1^2)*(th3^2-th2^2)*(th2^2-th1^2));
    eS = eye(7) + c1*S + c2*S^2 + c3*S^3 + c4*S^4 + c5*S^5 + c6*S^6;
    
elseif n == 8
    a = coeff(3); b = coeff(5); c = coeff(7); d = coeff(9);
    p = -2*b^3+9*a*b*c+72*b*d-27*c^2-27*a^2*d; q = b^2-3*a*c+12*d; 
    phi = acos(-0.5*p/q^1.5); y = (1/3)*(b+2*sqrt(q)*cos(phi/3));
    R=0.5*sqrt(a^2-4*b+4*y);
    D=0.5*sqrt(3*a^2-4*R^2-8*b+(4*a*b-8*c-a^3)/R);
    E=0.5*sqrt(3*a^2-4*R^2-8*b-(4*a*b-8*c-a^3)/R);
    th1 = 0.5*sqrt(a-2*R-2*D);
    th2 = 0.5*sqrt(a-2*R+2*D);
    th3 = 0.5*sqrt(a+2*R-2*E);
    th4 = 0.5*sqrt(a+2*R+2*E);
    c0 = P(1)^2*T(3,2)*T(4,2)*T(4,3)*cos(th1)...
        -P(2)^2*T(3,1)*T(4,1)*T(4,3)*cos(th2)...
        +P(3)^2*T(2,1)*T(4,1)*T(4,2)*cos(th3)...
        -P(4)^2*T(2,1)*T(3,1)*T(3,2)*cos(th4);
    c0 = c0/(T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c1 = P(1)^3*T(3,2)*T(4,2)*T(4,3)*sin(th1)...
        -P(2)^3*T(3,1)*T(4,1)*T(4,3)*sin(th2)...
        +P(3)^3*T(2,1)*T(4,1)*T(4,2)*sin(th3)...
        -P(4)^3*T(2,1)*T(3,1)*T(3,2)*sin(th4);
    c1 = c1/(th1*th2*th3*th4*T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c2 = U(1)*T(3,2)*T(4,2)*T(4,3)*cos(th1)...
        -U(2)*T(3,1)*T(4,1)*T(4,3)*cos(th2)...
        +U(3)*T(2,1)*T(4,1)*T(4,2)*cos(th3)...
        -U(4)*T(2,1)*T(3,1)*T(3,2)*cos(th4);
    c2 = c2/(T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c3 = P(1)*U(1)*T(3,2)*T(4,2)*T(4,3)*sin(th1)...
        -P(2)*U(2)*T(3,1)*T(4,1)*T(4,3)*sin(th2)...
        +P(3)*U(3)*T(2,1)*T(4,1)*T(4,2)*sin(th3)...
        -P(4)*U(4)*T(2,1)*T(3,1)*T(3,2)*sin(th4);
    c3 = c3/(th1*th2*th3*th4*T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c4 = V(2,3,4)*T(3,2)*T(4,2)*T(4,3)*cos(th1)...
        -V(1,3,4)*T(3,1)*T(4,1)*T(4,3)*cos(th2)...
        +V(1,2,4)*T(2,1)*T(4,1)*T(4,2)*cos(th3)...
        -V(1,2,3)*T(2,1)*T(3,1)*T(3,2)*cos(th4);
    c4 = c4/(T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c5 = P(1)*V(2,3,4)*T(3,2)*T(4,2)*T(4,3)*sin(th1)...
        -P(2)*V(1,3,4)*T(3,1)*T(4,1)*T(4,3)*sin(th2)...
        +P(3)*V(1,2,4)*T(2,1)*T(4,1)*T(4,2)*sin(th3)...
        -P(4)*V(1,2,3)*T(2,1)*T(3,1)*T(3,2)*sin(th4);
    c5 = c5/(th1*th2*th3*th4*T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c6 = T(3,2)*T(4,2)*T(4,3)*cos(th1)...
        -T(3,1)*T(4,1)*T(4,3)*cos(th2)...
        +T(2,1)*T(4,1)*T(4,2)*cos(th3)...
        -T(2,1)*T(3,1)*T(3,2)*cos(th4);
    c6 = c6/(T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c7 = P(1)*T(3,2)*T(4,2)*T(4,3)*sin(th1)...
        -P(2)*T(3,1)*T(4,1)*T(4,3)*sin(th2)...
        +P(3)*T(2,1)*T(4,1)*T(4,2)*sin(th3)...
        -P(4)*T(2,1)*T(3,1)*T(3,2)*sin(th4);
    c7 = c7/(th1*th2*th3*th4*T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));    
    eS = c0*eye(8)+c1*S+c2*S^2+c3*S^3+c4*S^4+c5*S^5+c6*S^6+c7*S^7;
    
elseif n ==9
    a = coeff(3); b = coeff(5); c = coeff(7); d = coeff(9);
    p = -2*b^3+9*a*b*c+72*b*d-27*c^2-27*a^2*d; q = b^2-3*a*c+12*d; 
    phi = acos(-0.5*p/q^1.5); y = (1/3)*(b+2*sqrt(q)*cos(phi/3));
    R=0.5*sqrt(a^2-4*b+4*y);
    D=0.5*sqrt(3*a^2-4*R^2-8*b+(4*a*b-8*c-a^3)/R);
    E=0.5*sqrt(3*a^2-4*R^2-8*b-(4*a*b-8*c-a^3)/R);
    th1 = 0.5*sqrt(a-2*R-2*D);
    th2 = 0.5*sqrt(a-2*R+2*D);
    th3 = 0.5*sqrt(a+2*R-2*E);
    th4 = 0.5*sqrt(a+2*R+2*E);    
    c1 = P(1)^3*T(3,2)*T(4,2)*T(4,3)*sin(th1)...
        -P(2)^3*T(3,1)*T(4,1)*T(4,3)*sin(th2)...
        +P(3)^3*T(2,1)*T(4,1)*T(4,2)*sin(th3)...
        -P(4)^3*T(2,1)*T(3,1)*T(3,2)*sin(th4);
    c1 = c1/(th1*th2*th3*th4*T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c2 = P(1)^4*T(3,2)*T(4,2)*T(4,3)*(1-cos(th1))...
        -P(2)^4*T(3,1)*T(4,1)*T(4,3)*(1-cos(th2))...
        +P(3)^4*T(2,1)*T(4,1)*T(4,2)*(1-cos(th3))...
        -P(4)^4*T(2,1)*T(3,1)*T(3,2)*(1-cos(th4));
    c2 = c2/(th1^2*th2^2*th3^2*th4^2*T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c3 = P(1)*U(1)*T(3,2)*T(4,2)*T(4,3)*sin(th1)...
        -P(2)*U(2)*T(3,1)*T(4,1)*T(4,3)*sin(th2)...
        +P(3)*U(3)*T(2,1)*T(4,1)*T(4,2)*sin(th3)...
        -P(4)*U(4)*T(2,1)*T(3,1)*T(3,2)*sin(th4);
    c3 = c3/(th1*th2*th3*th4*T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c4 = P(1)^2*U(1)*T(3,2)*T(4,2)*T(4,3)*(1-cos(th1))...
        -P(2)^2*U(2)*T(3,1)*T(4,1)*T(4,3)*(1-cos(th2))...
        +P(3)^2*U(3)*T(2,1)*T(4,1)*T(4,2)*(1-cos(th3))...
        -P(4)^2*U(4)*T(2,1)*T(3,1)*T(3,2)*(1-cos(th4));
    c4 = c4/(th1^2*th2^2*th3^2*th4^2*T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c5 = P(1)*V(2,3,4)*T(3,2)*T(4,2)*T(4,3)*sin(th1)...
        -P(2)*V(1,3,4)*T(3,1)*T(4,1)*T(4,3)*sin(th2)...
        +P(3)*V(1,2,4)*T(2,1)*T(4,1)*T(4,2)*sin(th3)...
        -P(4)*V(1,2,3)*T(2,1)*T(3,1)*T(3,2)*sin(th4);
    c5 = c5/(th1*th2*th3*th4*T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c6 = P(1)^2*V(2,3,4)*T(3,2)*T(4,2)*T(4,3)*(1-cos(th1))...
        -P(2)^2*V(1,3,4)*T(3,1)*T(4,1)*T(4,3)*(1-cos(th2))...
        +P(3)^2*V(1,2,4)*T(2,1)*T(4,1)*T(4,2)*(1-cos(th3))...
        -P(4)^2*V(1,2,3)*T(2,1)*T(3,1)*T(3,2)*(1-cos(th4));
    c6 = c6/(th1^2*th2^2*th3^2*th4^2*T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    c7 = P(1)*T(3,2)*T(4,2)*T(4,3)*sin(th1)...
        -P(2)*T(3,1)*T(4,1)*T(4,3)*sin(th2)...
        +P(3)*T(2,1)*T(4,1)*T(4,2)*sin(th3)...
        -P(4)*T(2,1)*T(3,1)*T(3,2)*sin(th4);
    c7 = c7/(th1*th2*th3*th4*T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1)); 
    c8 = P(1)^2*T(3,2)*T(4,2)*T(4,3)*(1-cos(th1))...
        -P(2)^2*T(3,1)*T(4,1)*T(4,3)*(1-cos(th2))...
        +P(3)^2*T(2,1)*T(4,1)*T(4,2)*(1-cos(th3))...
        -P(4)^2*T(2,1)*T(3,1)*T(3,2)*(1-cos(th4));
    c8 = c8/(th1^2*th2^2*th3^2*th4^2*T(4,1)*T(4,2)*T(4,3)*T(3,1)*T(3,2)*T(2,1));
    eS = eye(9) + c1*S + c2*S^2 + c3*S^3 + c4*S^4 + c5*S^5 + c6*S^6 + c7*S^7 + c8*S^8; 
end   
    function t = T(j,k)
        th = [th1,th2,th3,th4];
        t = th(j)^2-th(k)^2;
    end
    function u = U(j)
        th = [th1,th2,th3,th4];
        nth = th(th~=th(j));
        u = nth(1)^2*nth(2)^2+nth(1)^2*nth(3)^2+nth(2)^2*nth(3)^2;
    end
    function v = V(a,b,c)
        th = [th1,th2,th3,th4];
        v = th(a)^2+th(b)^2+th(c)^2;
    end
    function p = P(j)
        th = [th1,th2,th3,th4];
        nth = th(th~=th(j));
        p = nth(1)*nth(2)*nth(3);
    end
end