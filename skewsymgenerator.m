function S = skewsymgenerator(n,size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Skew-Symmetric Matrix Generator   
%                            by Diego Avalos
%                      http://www.diegoavalos.net/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The skwesymgenerator function randomly generates a skew-symmetric matrix
% S by taking two inputs:
% n    = dimension of matrix, and
% size = the integer bound for the entries S_{ij}.
% The ouptput is an n-by-n skew-symmetric matrix S, such that
% |S_{ij}| \leq size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Commented Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = n(n-1)/2 is the number of entries of S which we need to generate. We
% store the generated integer numbers in a row vector r.
N = n*(n-1)/2;
rng('shuffle');
r = randi([-size, size], 1, N);
% We preallocate our matrix S.
S = zeros(n,n);
% We now assign the entries in the vector r to the upper-triangular entries
% of S. The assingment of S is rowwise (from first to nth row of S).
for  rnum = 1:n
% We are in the rnum-th row assingment of S. We only need to assing the
% nonzero values of this row: [0 ... 0 S_{rnum,rnum+1} ... S_{rnum, n}], 
% which is a total of rsub = (\sigma_{jj=n-1}^{n-rnum}jj) assingments. This
% assings the following r values: r(rsub-jj+1), ..., r(rsub) to S.
    rsub = 0;
    for jj = (n-1):-1:(n-rnum)
        rsub = rsub +jj;
    end
    S(rnum,(rnum+1):n) = r((rsub-jj+1):rsub);
end
% We now have the strictly upper-triangular matrix S. To get the desired
% randomly-generated skew-symmetric matrix S, perform the operation:
S = S  - S';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 FIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
