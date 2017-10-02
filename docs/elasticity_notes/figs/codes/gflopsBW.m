
N = 2:5;
Np = (N+1).*(N+2).*(N+3)/6; 
Nfields = 4;
K = 49748; %10087; 

gflopsU = [ 124867480, 521110300,1636957940,4313897820];
bwU = [89.283+48.331,47.350+25.108, 55.729+29.223,40.666+21.136];
dofTimeU = [1.709e-10                 2.856e-10                 2.697e-10                 3.652e-10];
timeU = dofTimeU.*K.*Np*Nfields;

% gflopsU = gflopsU./timeU/1e9
bwU

%% strong

bwV = [121.86+28.377,64.829+8.2766,34.887+3.1257,12.326+275.17*1e-3];
gflopsV = [249933952,1920869776,8403233168,9.0950e+10];
bwS = [133.40+13.911,110.80+10.833,60.906+5.7574, 20.539+1.9310];
gflopsS = [339082368,1261609280,3783235904,9139304576];
dofTimeV = [1.588e-10                  5.32e-10                 1.313e-09                 1.483e-08];
dofTimeS = [3.223e-10                 4.003e-10                 7.079e-10                 2.105e-09];
timeV = dofTimeV.*K.*Np*Nfields;
timeS = dofTimeS.*K.*Np*Nfields;
gflopsStrong = [gflopsV;gflopsS; gflopsU]./[timeV;timeS;timeU]/1e9
gflopsStrong = [2:5; gflopsStrong]'
bwStrong = [2:5;bwV;bwS;bwU]'
timeStrong = [timeV;timeS]
timeStrongTotal = dofTimeV + dofTimeS + dofTimeU
print_pgf_coordinates(2:5,timeStrongTotal)
print_pgf_coordinates(2:5,dofTimeV)
print_pgf_coordinates(2:5,dofTimeS)


%% skew

bwS = [146.12+13.499,141.33+12.802,118.56+14.328,119.02+16.076];
gflopsS = [149086952,558262440,1160834184,2646761376];
bwQ = [26.534+72.564,34.780+103.08, 43.400+93.589,42.090+83.818];
gflopsQ = [66861312,238790400,453701760,936058368];
bwV = [109.60+25.572, 61.286+12.856, 43.482+8.6943,27.571+5.3064];
gflopsV = [196305608,804375412,2499488764,6547284532];

dofTimeV = [ 1.754e-10                 3.448e-10                 4.656e-10                 7.611e-10];
dofTimeS = [ 3.298e-10                  3.34e-10                  2.83e-10                 2.518e-10];
dofTimeQ = [ 1.652e-10                 1.285e-10                 9.354e-11                 9.656e-11];
timeV = dofTimeV.*K.*Np*Nfields;
timeS = dofTimeS.*K.*Np*Nfields;
timeQ = dofTimeQ.*K.*Np*Nfields;

gflopsSkew = [gflopsV;gflopsS;gflopsQ]./[timeV;timeS;timeQ]/1e9
gflopsSkew = [2:5;gflopsSkew]'
bwSkew = [2:5;bwV;bwS;bwQ]'

timeSkew = [timeV;timeS+timeQ]

timeSkewTotal = dofTimeV + dofTimeS + dofTimeQ + dofTimeU
print_pgf_coordinates(2:5,timeSkewTotal)

print_pgf_coordinates(2:5,dofTimeV)
print_pgf_coordinates(2:5,dofTimeS)
print_pgf_coordinates(2:5,dofTimeQ)

% print_pgf_coordinates(2:5,timeSkewTotal)