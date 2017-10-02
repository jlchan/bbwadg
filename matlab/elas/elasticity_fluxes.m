syms n1 n2 n3 u1 u2 u3
B = [n1 0 0 ; 0 n2 0; 0 0 n3;0 n3 n2; n3 0 n1; n2 n1 0];

v = B*[u1;u2;u3];
subs(simplify(conj(v')*v),n1^2+n2^2+n3^2,1)

B = [ n1 0; 0 n2; n2 n1];
conj(B')*B*[u1;u2]

% subs(simplify(conj(B')*B*[u1;u2;u3]),n1^2+n2^2+n3^2,1)
% simplify(subs(expand(cross(-[n1;n2;n3],cross([n1;n2;n3],[u1;u2;u3]))),n1^2+n2^2+n3^2,1))