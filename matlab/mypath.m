
pw = pwd;
addpath([pw,'/Codes1D']);
if useQuads
    addpath([pw,'/Codes2DQuad']);
else
    addpath([pw,'/Codes2D']);
end
addpath([pw,'/Codes3D']);
addpath([pw,'/Grid']);
addpath([pw,'/ServiceRoutines']);
% addpath([pw,'/wedge']);
% addpath([pw,'/pyr']);
% addpath([pw,'/pyr/sym_pyr']);
% addpath([pw,'/hex']);
% addpath([pw,'/tet']);
% addpath([pw,'/hybrid']);
addpath([pw,'/scratch']);
addpath([pw,'/bern']);
addpath([pw]);