function [M, I, m, i, successratio, eigratio] = thesisresults(matrix_size, entries_size, number_of_matrices) 
close all
n = matrix_size; s = entries_size; N = number_of_matrices;
Myabserrors = zeros(1,N); 

Matabserrors = zeros(1,N); 

relerror = zeros(1,N);
Matrices = zeros(n,n*N); eigcomparison = zeros(1,N);
count = 0;
for j = 1:N
    S = skewsymgenerator(n,s); 
    eS = skewsymexpm(S); MateS = expm(S);
    Myabserrors(j) = norm((eye(n)-eS'*eS),'fro');
    Matabserrors(j) = norm((eye(n)-MateS'*MateS),'fro');    
    relerror(j) = norm(MateS-eS,'fro');
    Matrices(:,n*j-n+1:n*j) = S;
    if Myabserrors(j) <= Matabserrors(j)
        count = count + 1;
    else
        th = unique(abs(imag(eig(S)))); th = th(th~=0);
        eigcomparison(j) = min(th)/max(th);        
    end
    
end
eigcomparison = eigcomparison(eigcomparison~=0); h=length(eigcomparison);
eigratio = mean(eigcomparison);
successratio = count/N;
[~, I] = max(Matabserrors-Myabserrors); I=I(1);
[~, i] = min(Matabserrors-Myabserrors); i=i(1);
M = Matrices(:,n*I-n+1:n*I);
m = Matrices(:,n*i-n+1:n*i);

figure(1)
semilogy(1:N, Matabserrors, 'ob',1:N, Myabserrors,'*r');
title('Measure of Absolute Error Against the Identity')
xlabel(['Generated Matrices S of size n = ', num2str(n)])
ylabel('$||I-(e^S)^Te^S||$','Interpreter','latex')
legend('Matlab''s expm','Closed Formula')
figure(2)
plot(1:h,eigcomparison, 'ok');
title('Eigenvalue ratio')
