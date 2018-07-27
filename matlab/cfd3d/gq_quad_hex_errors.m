%% 2D Gpu runs N = 4

K1D = [8 16 32 64];
h = 1./K1D;

errGLLa{2} = [1.22861 0.236646 0.0258883 0.00406122];
errGLLa{3} = [0.414449 0.0412224 0.00293039 0.000145701];
errGLLa{4} = [0.123516 0.00439285 0.000126852 5.61786e-06];
errGLLa{5} = [0.0332509 0.000970148 2.10334e-05 2.60008e-07];
errGLLa{6} = [0.0106752 0.000133284 7.97032e-07 nan];
errGLLa{7} = [0.00330917 1.67825e-05 1.32696e-07 nan];

errGQa{2}  = [0.553435 0.0625441 0.0117156 .0022841];
errGQa{3}  = [0.25082 0.0198863 0.000673152 2.98263e-05];
errGQa{4}  = [0.0768018 0.0020579 5.55322e-05  2.36108e-06];
errGQa{5}  = [0.0151231 0.000454595 1.00203e-05 6.57796e-08];
errGQa{6}  = [0.00657096 9.22845e-05 3.60025e-07 nan];
errGQa{7}  = [0.00185273 1.21465e-05 1.14289e-07 nan];

% errGLLt{2} = [1.54813 0.434129 0.075089]; % bilinear
% errGLLt{3} = [0.708886 0.134981 0.024557];
% errGLLt{4} = [0.288314 0.0507278 0.0022585];
% errGLLt{5} = [0.159946 0.00910103 0.000373439];
% errGLLt{6} = [];
% errGLLt{7} = [];
% errGQt{2} = [0.779976 0.1729 0.0351796]; % bilinear
% errGQt{3} = [0.346929 0.0734158 0.00386658];
% errGQt{4} = [0.168446 0.0119508 0.000311064];
% errGQt{5} = [0.0744237 0.0030628 6.51708e-05];
% errGQt{6} = [];
% errGQt{7} = [];

% lightly curved (a = .125/4 = .0312)
errGLLc1{2} = [1.2082 0.265357 0.0301118 0.00467693]; 
errGLLc1{3} = [0.448878 0.0485323 0.00385418 0.000225396]; 
errGLLc1{4} = [0.122081 0.00694351 0.000224524 1.14555e-05]; 
errGLLc1{5} = [0.0424351 0.00114569 3.59697e-05 5.09643e-07]; 
errGLLc1{6} = [0.0112838 0.000219049 1.84304e-06 nan]; 
errGLLc1{7} = [0.00388316 2.55074e-05 2.68102e-07 nan]; 
errGQc1{2}  = [0.555586 0.0773534 0.0129171 0.0024363];
errGQc1{3}  = [0.287685 0.0190431 0.000721055 3.19095e-05];
errGQc1{4}  = [0.0731583 0.00213914 6.32128e-05 2.73037e-06];
errGQc1{5}  = [0.0215298 0.000523714 9.27521e-06 7.14031e-08];
errGQc1{6}  = [0.00809138 8.23415e-05 3.95565e-07 nan];
errGQc1{7}  = [0.00222029 1.33111e-05 9.72405e-08 nan];

% % lightly curved (a = .125/2 = .0625)
% errGLLc1{2}  = [1.35014 0.330866 0.0355031 0.00582492];
% errGLLc1{3}  = [0.488452 0.0525712 0.00621317 0.0004342];
% errGLLc1{4}  = [0.143136 0.0105404 0.000425336 2.27161e-05];
% errGLLc1{5}  = [0.0541884 0.00166351 7.04547e-05 1.23574e-06];
% errGLLc1{6}  = [0.0148569 0.000438661 3.99603e-06 nan];
% errGLLc1{7}  = [0.00595663 5.74999e-05 5.8022e-07 nan];
% 
% errGQc1{2}  = [0.631228 0.0844294 0.0153232 0.00270066];
% errGQc1{3}  = [0.259987 0.0206532 0.000899763 3.92541e-05];
% errGQc1{4}  = [0.0901608 0.00269128 7.66298e-05 3.14632e-06];
% errGQc1{5}  = [0.0289158 0.000579468 1.17134e-05 9.30656e-08];
% errGQc1{6}  = [0.0106526 0.000122192 5.39937e-07 nan];
% errGQc1{7}  = [0.0035827 1.67731e-05 1.23304e-07 nan];

errGLLc2{2} = [1.64274 0.445945 0.0720148 0.0091066]; % curved (a = .125)
errGLLc2{3} = [0.662165 0.146131 0.0234934 0.00105368]; % curved (a = .125)
errGLLc2{4} = [0.349184 0.0508268 0.00234785 5.36192e-05]; % curved (a = .125)
errGLLc2{5} = [0.171986 0.0103578 0.000386834 5.76284e-06]; % curved (a = .125)
errGLLc2{6} = [0.0729963 0.00276104 4.11305e-05 2.36125e-07]; % curved (a = .125)
errGLLc2{7} = [0.026577 0.000589552 5.3039e-06 nan]; % curved (a = .125)

% curved (a = .125)
errGQc2{2}  = [0.882247 0.138116 0.0268145 0.0039945];
errGQc2{3}  = [0.368084 0.0734059 0.00420438 9.19378e-05];
errGQc2{4}  = [0.183472 0.0127617 0.000251929 7.6013e-06]; % curved (a = .125)
errGQc2{5}  = [0.0924215 0.00276519 6.80528e-05 3.27246e-07];
errGQc2{6}  = [0.0378013 0.0007649 4.15848e-06 nan];
errGQc2{7}  = [0.0154563 0.000146207 8.83156e-07 nan];

% curved (a = .25)
errGLLc3{2} = [2.08578 1.33668 0.433653 0.0721002]; % curved (a = .25)
errGLLc3{3} = [1.74649 0.546097 0.0879841 0.0090851]; % curved (a = .25)
errGLLc3{4} = [1.06289 0.218361 0.0209471 0.000665425]; % curved (a = .25)
errGLLc3{5} = [0.695283 0.0996971 0.00390099 0.000118992]; % curved (a = .25)
errGLLc3{6} = [0.449493 0.0374955 0.000909834 7.70064e-06]; % curved (a = .25)
errGLLc3{7} = [0.294238 0.0162153 0.000166829 1.19661e-06]; % curved (a = .25)

errGQc3{2}  = [1.82833 0.746283 0.11006 0.0110115]; % curved (a = .25)
errGQc3{3}  = [1.39268 0.245852 0.0255562 0.00165226]; % curved (a = .25)
errGQc3{4}  = [0.740071 0.105075 0.003505 8.26421e-05]; % curved (a = .25)
errGQc3{5}  = [0.457324 .0354479 0.000767746 1.95823e-05]; % curved (a = .25)
errGQc3{6}  = [0.274721 0.0120499 0.000157175 9.61243e-07]; % curved (a = .25)
errGQc3{7}  = [0.150916 0.00538908 2.79182e-05 nan]; % curved (a = .25)

for N = 2:7 
%     loglog(h,errGLLa{N},'o-');hold on    
%     loglog(h,errGQa{N},'x--')
    
    loglog(h,errGLLc1{N},'o-');hold on    
    loglog(h,errGQc1{N},'x--')
    
%     loglog(h,errGLLc2{N},'o-');hold on    
%     loglog(h,errGQc2{N},'x--')

% loglog(h,errGLLc3{N},'o-');hold on    
%     loglog(h,errGQc3{N},'x--')
    
%     print_pgf_coordinates(h,errGLLa{N})
%     print_pgf_coordinates(h,errGQa{N})
    print_pgf_coordinates(h,errGLLc1{N})
    print_pgf_coordinates(h,errGQc1{N})
%     print_pgf_coordinates(hcc,errGLLc2{N})
%     print_pgf_coordinates(hcc,errGQc2{N})
%     print_pgf_coordinates(hcc,errGLLc3{N})
%     print_pgf_coordinates(hcc,errGQc3{N})
end
    
return

%% 3d

h = 1./[4 8 16 32];
errGLLa{2} = [2.62759 .728523 0.112225 0.018711];
errGLLa{3} = [0.934758 0.191949 0.0119788 .000558681];
errGLLa{4} = [0.460425 0.02445 0.000821249 4.67491e-05];
errGLLa{5} = [0.138105 0.00744055 0.000126057 nan]
errGLLa{6} = [0.0489719 0.000897467 2.55644e-05 nan]
errGLLa{7} = [0.023319 0.00021951 nan nan]
errGQa{2} = [1.53203 0.23283 0.0417033 0.00806707];
errGQa{3} = [0.562488 0.0864096 0.00268061 0.000136725];
errGQa{4} = [0.239027 0.0077433 0.000219705 2.67589e-05];
errGQa{5} = [0.0547127 0.00242675 4.2962e-05 nan];
errGQa{6} = [0.0246519 0.000324802 2.44009e-05 nan];
errGQa{7} = [0.00801774 7.39203e-05 nan nan];

% % a = .5 
% errGLLc{2} = [2.70281 0.808467 0.13359 0.0224913];
% errGLLc{3} = [1.02978 0.214321 0.025052 0.00153341];
% errGLLc{4} = [0.506958 0.0514311 0.00144395 7.83588e-05];
% 
% errGQc{2} = [1.65196 0.288063 0.0573897 0.0125185];
% errGQc{3} = [0.592762 0.110257 0.00576595 nan];
% errGQc{4} = [0.275171 0.0127005 0.000331244 nan];

% a = 1.0
errGLLc{2} = [2.93541 1.00445 0.18488 0.0323814];
errGLLc{3} = [1.23231 0.266638 0.0513328 0.00390276];
errGLLc{4} = [0.599851 0.101098 0.00511836 0.00018702];
errGLLc{5} = [0.310649 0.0199095 0.00157312 nan];
errGLLc{6} = [0.130579 0.0108408 0.000188745 nan];
errGLLc{7} = [0.0558635 0.00251764 nan nan];
errGQc{2} = [1.99429 0.423426 .0787544 .0173475];
errGQc{3} = [0.727022 0.148442 0.0158392 0.000446109];
errGQc{4} = [0.351163 0.0353206 0.000705351 4.41811e-05];
errGQc{5} = [0.144996 0.00619848 0.000311659 nan];
errGQc{6} = [0.0517392 0.00290387 3.71606e-05 nan];
errGQc{7} = [0.020886 0.000421991 nan nan];

% a = 1.0, Ngeo = floor(N/2) + 1 approximation
errGLLc2{2} = [2.80762 0.951269 0.173495 0.0295367];
errGLLc2{3} = [1.21946 0.257914 0.0494569 0.00380144];
errGLLc2{4} = [0.594597 0.100782 0.0051111 0.000186621];
errGLLc2{5} = [0.306622 0.0198205 0.00156957 nan];
errGLLc2{6} = [0.130865 0.0108333 0.000188656 nan ];
errGLLc2{7} = [0.0558827 0.00251998 nan nan ];
errGQc2{2} = [1.89337 0.366475 0.0599646 0.0133605];
errGQc2{3} = [0.708125 0.139733 0.015052 0.000419943];
errGQc2{4} = [0.347316 0.0351013 0.000702893 4.39052e-05];
errGQc2{5} = [0.142279 0.00620177  0.000311399 nan ];
errGQc2{6} = [0.0518407 0.00289859 3.71549e-05 nan ];
errGQc2{7} = [0.0208791 0.000423102 nan nan ];
for N = 2:7
%     loglog(h,errGLLa{N},'o-');hold on    
%     loglog(h,errGQa{N},'x--')
    
    loglog(h,errGLLc{N},'o-');hold on    
    loglog(h,errGQc{N},'x--')        
    
%     loglog(h,errGLLc2{N},'o-');hold on    
%     loglog(h,errGQc2{N},'x--')        
    
% disp('\addplot+[semithick, mark options={solid, fill=markercolor}]')
%     print_pgf_coordinates(h,errGLLa{N})
%     disp('\addplot+[semithick, mark options={solid, fill=markercolor}]')
%     print_pgf_coordinates(h,errGQa{N})
    
% disp('\addplot+[semithick, mark options={solid, fill=markercolor}]')    
%     print_pgf_coordinates(h,errGLLc2{N})
% disp('\addplot+[semithick, mark options={solid, fill=markercolor}]')        
%     print_pgf_coordinates(h,errGQc2{N})


%     disp('\addplot+[semithick, mark options={solid, fill=markercolor}]')
%     print_pgf_coordinates(h,errGLLc{N})
%     disp('\addplot+[semithick, mark options={solid, fill=markercolor}]')
%     print_pgf_coordinates(h,errGQc{N})

end


















%% old results

errGLL = [0.190411 0.0168452 0.00104654 2.71597e-05 1.06648e-06]; % affine
errGQ = [0.142168 0.0118826 0.000397138 9.60209e-06 4.70206e-07]; % affine
ratio_affine = errGLL./errGQ

loglog(h,errGLL,'o-')
hold on
loglog(h,errGQ,'x-')
loglog(h,1e3*h.^5,'k--')

errGLL = [0.518557 0.241366 0.0485364 0.00413028 0.000262928]; % curved (a = .25)
errGQ = [0.482721 0.158578 0.0222502 0.000690727 4.09657e-05]; % curved (a = .25)
ratio_curved = errGLL./errGQ

% really curved: a = .4
% K1D = 16: errGLL = 0.2601, errGQ = 0.15862
% K1D = 32: errGLL = 0.0666038, errGQ = 0.0261791

loglog(h,errGLL,'o-')
hold on
loglog(h,errGQ,'x-')
loglog(h,5e3*h.^4,'k--')

%% 3D gpu runs

h = 1./[4 8 16 32];

% N = 4
errGLL  = [0.0765747 0.00404953 0.000134603]; % affine
errGQ = [0.0280623 0.00130065 3.7327e-05]; % affine

errGLL = [0.099558 0.010684 0.000799367];
errGQ = [0.0554412 0.00943044 0.000201176];

ratio_affine = errGLL./errGQ

loglog(h,errGLL,'o-')
hold on
loglog(h,errGQ,'x-')
loglog(h,1e2*h.^5,'k--')
