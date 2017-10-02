clear
N = 50;
x = randn(N,1);
y = randn(N,1);
x = x/max(abs(x));
y = y/max(abs(y));
ids = 1:length(x);
x = [x; linspace(-1,1,10)'];
y = [y; -1*ones(10,1)];
x = [x; linspace(-1,1,10)'];
y = [y; 1*ones(10,1)];

y = [y; linspace(-1,1,10)'];
x = [x; -1*ones(10,1)];
y = [y; linspace(-1,1,10)'];
x = [x; 1*ones(10,1)];

tri = delaunayFixArea(x',y');

triplot(tri,x,y)

% adjacency matrix
A = sparse(N,N);
for e = 1:size(tri,1)
    A(tri(e,:),tri(e,:)) = 1; % adjacency matrix
end
nbrs = {};
for i = 1:N
    nbrs{i} = setdiff(find(A(i,:)),i);
end

maxit = 100;
h = 1/N;
for iter = 1:maxit
    
    alpha = iter/maxit;
    %y(Fmask(:,1),:) = y0(Fmask(:,1))*(1-alpha) + alpha*yf;
    
    clf;
    hold on;
    triplot(tri,x,y);plot(x,y,'s')
    title(sprintf('iteration %d',iter))
    drawnow
    
    f = zeros(N,2);
    for i = 1:N
        
        for j = 1:length(nbrs{i})
            nbr = nbrs{i}(j);
            v = [x(i)-x(nbr); y(i)-y(nbr)];
            d = norm(v);            
            f(i,:) = f(i,:) + v(:)'*(h > d)*(h-d) + v(:)'*(h < d)*(h-d); % makes compression forces stronger
            %                 f(i,:) = f(i,:) + v(:)'*(h-d); % compression terms stronger
        end
    end
    x(ids) = x(ids) + h*f((ids),1);
    y(ids) = y(ids) + h*f((ids),2);
    
end