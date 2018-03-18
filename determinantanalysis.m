function [interestingmatrix, I, notsointerestingmatrix, i] = determinantanalysis(N,s) 
Matrices = zeros(5,5*N);
close all
dets = zeros(1,N); err = zeros(1,N);
Dets = zeros(1,N); Err = zeros(1,N);

for j = 1:N
    r=randi([-s s],1,10);
    [M, m, ~, S] = expskew5(r);
    dets(j)=abs(1-det(m)); err(j)=norm((eye(5)-m'*m),'fro');
    Dets(j)=abs(1-det(M)); Err(j)=norm((eye(5)-M'*M),'fro');
    Matrices(:,5*j-4:5*j) = S;
end
[~, I] = max(err-Err); I=I(1);
[~, i] = min(err-Err); i=i(1);
interestingmatrix = Matrices(:,5*I-4:5*I);
notsointerestingmatrix = Matrices(:,5*i-4:5*i);
figure(1)
plot(1:N, dets, '--ob',1:N, Dets,'-*r');
title('Measure of Orthogonality Error')
xlabel('Generated Matrices S')
ylabel('$|1-\det(S)|$','Interpreter','latex')
legend('Matlab''s expm','Closed Formula')
figure(2)
plot(1:N, err, '--ob',1:N, Err,'-*r');
title('Measure of Error Against The Identity')
xlabel('Generated Matrices S')
ylabel('$||I-(e^S)^Te^S||_F$','Interpreter','latex')
legend('Matlab''s expm','Closed Formula')
