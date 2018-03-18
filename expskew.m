function E = expskew(S)
[~, n] = size(S);
thetas = unique(abs(imag(eig(S))));
thetas = thetas(thetas>0);
[p, ~] = size(thetas);
THETA = zeros(p,p);
for j = 1:p
    sig = mod(j,4)+1;
    sig = (-1)^sig;
    THETA(j,:) = sig*thetas.^(2*j-1);
end
coefficients = inv(THETA); %inverting this probably causes a lot of error
Sp = zeros(n,n*p);
for j = 1:p
    for k = 1:p
        Sp(:,(n*j-n+1:n*j))= coefficients(k,j)*S^(2*j-1);
    end
end
E = eye(n);
for j = 1:p
    E = E + sin(thetas(j))*Sp(:,(n*j-n+1):n*j) + (1-cos(thetas(j)))*Sp(:,(n*j-n+1):n*j)^2;
end
