K = [4 8 16 32 64];

% knots
print_pgf_coordinates(K,[0.012097 0.007485 0.004234 0.002268 0.001176]./[0.057496 0.045147 0.027045 0.014692 0.007646])
print_pgf_coordinates(K,[0.008967 0.004076 0.002625 0.001981 0.001197]./[0.101090 0.083937 0.051790 0.028574 0.014989])
print_pgf_coordinates(K,[0.003910 0.009638 0.007409 0.004827 0.002777]./[0.134560 0.117629 0.074592 0.041765 0.022080]);

% greville
print_pgf_coordinates(K,[0.006048 0.003743 0.002117 0.001134 0.000588]./[0.028748 0.037257 0.024875 0.014071 0.007453])
print_pgf_coordinates(K,[0.002989 0.002040 0.001727 0.001166 0.000677]./[0.033697 0.055806 0.044250 0.026584 0.014460])
print_pgf_coordinates(K,[0.000978 0.002409 0.001611 0.000895 0.000637]./[0.033640 0.059060 0.058858 0.037744 0.021095])

%%

% N = 2, K = 4: knot error = 0.012097, , original knot error = 0.057496, Greville error = 0.006048, original Greville error = 0.028748
% N = 2, K = 8: knot error = 0.007485, , original knot error = 0.045147, Greville error = 0.003743, original Greville error = 0.037257
% N = 2, K = 16: knot error = 0.004234, , original knot error = 0.027045, Greville error = 0.002117, original Greville error = 0.024875
% N = 2, K = 32: knot error = 0.002268, , original knot error = 0.014692, Greville error = 0.001134, original Greville error = 0.014071
% N = 2, K = 64: knot error = 0.001176, , original knot error = 0.007646, Greville error = 0.000588, original Greville error = 0.007453
% 
% N = 3, K = 4: knot error = 0.008967, , original knot error = 0.101090, Greville error = 0.002989, original Greville error = 0.033697
% N = 3, K = 8: knot error = 0.004076, , original knot error = 0.083937, Greville error = 0.002040, original Greville error = 0.055806
% N = 3, K = 16: knot error = 0.002625, , original knot error = 0.051790, Greville error = 0.001727, original Greville error = 0.044250
% N = 3, K = 32: knot error = 0.001981, , original knot error = 0.028574, Greville error = 0.001166, original Greville error = 0.026584
% N = 3, K = 64: knot error = 0.001197, , original knot error = 0.014989, Greville error = 0.000677, original Greville error = 0.014460
% 
% N = 4, K = 4: knot error = 0.003910, , original knot error = 0.134560, Greville error = 0.000978, original Greville error = 0.033640
% N = 4, K = 8: knot error = 0.009638, , original knot error = 0.117629, Greville error = 0.002409, original Greville error = 0.059060
% N = 4, K = 16: knot error = 0.007409, , original knot error = 0.074592, Greville error = 0.001611, original Greville error = 0.058858
% N = 4, K = 32: knot error = 0.004827, , original knot error = 0.041765, Greville error = 0.000895, original Greville error = 0.037744
% N = 4, K = 64: knot error = 0.002777, , original knot error = 0.022080, Greville error = 0.000637, original Greville error = 0.021095
