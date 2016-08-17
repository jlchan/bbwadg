i = sym('i');
j = sym('j');
k = sym('k');

N = sym('N');


% id = \suj_{a=0}^{a=i-1} (N+1-a)*(N+2-a)/2 + ...
%      + \suj_{a=0}^{a=j-1} (N+1-i-a)
%      + k  

id1 = (N+1)*(N+2)/2*i - (2*N+3)*(i-1)*(i)/4 + (i-1)*(i)*(2*i-1)/12
id2 = (N+1-i)*j - (j-1)*j/2
id3 = k

id = id1+id2+id3

% check that the indexing is contiguous and matches this choice of ordering
N = 3
for i=0:N
    for j=0:N-i
        for k=0:N-i-j
            [i,j,k,eval(id)]
        end
    end
end

