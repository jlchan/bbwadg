% first the integration tables for the local mass and stiffness matrices
%   must be generated. This will take a while since it calculates the 
%   entries for p=1-129
createIntegrationTables


%% a few definitions
ghostBasis    = 1;
compatibility = 2;

%% plot eigenvalues
plotEigsByOrder( 201,101,compatibility );
plotEigsByOrder( 201,21,ghostBasis );

%% do a number of convergence studies
plotConvCompatibility
plotConvGhostBasis
plotConvCompatibilityExtreme

