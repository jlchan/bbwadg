function [r s t] = pyr_nodes(N)

[r s t] = pyramidWBNodes3D(N); 
[r s t] = sym_to_biunit(r,s,t);
