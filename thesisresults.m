function [] = thesisresults(matrix_size, entries_size, number_of_matrices) 
close all
n = matrix_size; s = entries_size; N = number_of_matrices;
Myabserrors = zeros(1,N); Myorthoerrors = zeros(1,N);
Matabserrors = zeros(1,N); Matorthoerrors = zeros(1,N);
for j = 1:N
    S = skewsymgenerator(n,s); 
    eS = skewsymexpm(S); MateS = expm(S);
    Myabserrors(j) = norm((eye(n)-eS*eS'),'fro');
    Matabserrors(j) = norm((eye(n)-MateS*MateS'),'fro');
    
    Myorthoerrors(j) = abs(1-det(eS));
    Matorthoerrors(j) = abs(1-det(MateS));
end

figure(1)
plot(1:N, Matabserrors, '--ob',1:N, Myabserrors,'-*r');
title('Measure of Absolute Error Against the Identity')
xlabel(['Generated Matrices S of size n = ', num2str(n)])
ylabel('$||I-(e^S)^Te^S||$','Interpreter','latex')
legend('Matlab''s expm','Closed Formula')
figure(2)
plot(1:N, Matorthoerrors, '--ob',1:N, Myorthoerrors,'-*r');
title('Measure of Orthogonality Error')
xlabel(['Generated Matrices S of size n = ', num2str(n)])
ylabel('$|1-\det(e^S)|$','Interpreter','latex')
legend('Matlab''s expm','Closed Formula')