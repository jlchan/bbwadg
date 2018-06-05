%inputs are two matrices: the initial conditions at time = 0 and at time = 1
%C is another input, that indicates certain properties of the intensity of
%the initial disruption 

% Dirichlet boundary conditions

function explicit_time(A0, A1, C, x_length, y_length)

%A0 and A1 area mxn matrices 
%initialize time step 
t = 2; 

%initialize length of area mapped
a = x_length;       %input actual length here 
b = y_length;       %input actual length here

%initialize height of area mapped
(n,n) = size(A0);   %square matrices 

%discretize both direction 
%m is the discretization in the x-direction, n is the discretization in the
%y-direction
maxt = 100;         %max number of steps
hx = a/n;
hy = b/n;
ht = max/n
%using dirichlet boundary conditions

%initialize U

%initialize matrices
U1 = A1; %this will be updated as we go on 
U0 = A0; %as will thi
%initialize U
U_next = zeros(n,n);
for t = 2:maxt;
    for i = 1:n
        for j = 1:n
            if i == 1 
                if j == 1
                    U_next(i,j) = (C(i,j)^2*((U1(i+1,j) - 2*U1(i,j) + 0)/(hx^2) ...
                        + (U1(i,j+1) - 2*U1(i,j) + 0)/(hy^2)))*(ht^2) + ...
                        2*U1(i,j) - U0(i,j);
                end
                if j == n
                    U_next(i,j) = (C(i,j)^2*((U1(i+1,j) - 2*U1(i,j) + 0)/(hx^2) ...
                        + 0 - 2*U1(i,j) + U1(i,j-1))/(hy^2)))*(ht^2) + ...
                        2*U1(i,j) - U0(i,j);
                end
                U_next(i,j) = (C(i,j)^2*((U1(i+1,j) - 2*U1(i,j) + 0)/(hx^2) ...
                    + (U1(i,j+1) - 2*U1(i,j) + U1(i,j-1))/(hy^2)))*(ht^2) + ...
                    2*U1(i,j) - U0(i,j);
                
            end
            if i == n
                if j == 1
                    U_next(i,j) = (C(i,j)^2*(0 - 2*U1(i,j) + U1(i-1,j))/(hx^2) ...
                        + (U1(i,j+1) - 2*U1(i,j) + 0)/(hy^2)))*(ht^2) + ...
                        2*U1(i,j) - U0(i,j);
                end
                if j == n
                    U_next(i,j) = (C(i,j)^2*((0 - 2*U1(i,j) + U1(i-1,j))/(hx^2) ...
                        + 0 - 2*U1(i,j) + U1(i,j-1))/(hy^2)))*(ht^2) + ...
                        2*U1(i,j) - U0(i,j);
                end
                U_next(i,j) = (C(i,j)^2*((0 - 2*U1(i,j) + U1(i-1,j))/(hx^2) ...
                    + (U1(i,j+1) - 2*U1(i,j) + U1(i,j-1))/(hy^2)))*(ht^2) + ...
                    2*U1(i,j) - U0(i,j);    
            end
            
            U_next(i,j) = (C(i,j)^2*((U1(i+1,j) - 2*U1(i,j) + U1(i-1,j))/(hx^2) ...
                + (U1(i,j+1) - 2*U1(i,j) + U1(i,j-1))/(hy^2)))*(ht^2) + ...
                2*U1(i,j) - U0(i,j);
        end
    end
    U0 = U1;
    U1 = U_next;
end
    


end