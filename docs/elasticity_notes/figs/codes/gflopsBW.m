N = (1:7)';
Np = (N+1).*(N+2).*(N+3)/6; 
Nfields = 9;
K = 9918;

flopsU = [9346260 45804528  188283942 587290704 1542030798 3524953740 7343819868]';
flopsV = [9896040 50187060  185197320 546756210 1373570352 3057876360 6200575920]';
flopsS = [14796936 50092812 140272440 337996890 723070656  1408096536 2544130512]';

bwU = [70.251+39.375, 76.168+42.614,55.421+30.569, 37.920+20.755, 33.752+18.348, 28.664+15.520, 16.265+8.7676]';
bwV = [40.408+56.772, 57.966+52.889,49.888+47.655, 33.635+32.756, 23.767+23.353, 16.553+16.343, 11.601+11.483]';
bwS = [55.879+26.452, 75.624+22.578,95.730+23.234, 101.88+22.483, 98.656+21.492, 74.419+16.396, 66.842+14.767]';

dofTimeU = min(KblkU,[],2)/2.5;
dofTimeV = min(KblkV,[],2)/2.5;
dofTimeS = min(KblkS,[],2)/2.5;


timeU = dofTimeU.*K.*Np*Nfields;
timeV = dofTimeV.*K.*Np*Nfields;
timeS = dofTimeS.*K.*Np*Nfields;
gflops = [flopsV flopsS flopsU]./[timeV timeS timeU]/1e9
gflops = [N gflops]
bw = [N bwV bwS bwU]
time = [timeV timeS]
timeTotal = dofTimeV + dofTimeS + dofTimeU
% print_pgf_coordinates(N,timeTotal)
print_pgf_coordinates(N,dofTimeV)
print_pgf_coordinates(N,dofTimeS)
print_pgf_coordinates(N,dofTimeU)

