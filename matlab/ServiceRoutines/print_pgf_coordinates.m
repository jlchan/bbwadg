function print_pgf_coordinates(x,y)

if nargin<2
    y = x;
    x = 1:length(y);
end

fprintf('coordinates{');
for i = 1:length(x)
    fprintf('(%g,%g)',x(i),y(i));
end
fprintf('};\n');