clc
% GLL EC
h = [.5.^(1:3)]';
err = [0.203599,0.0947163,0.0463705]';
CC = [ones(size(h)) log(h)]\log(err); C(1) = CC(2);
err = [0.0137392,0.000701926,8.64531e-05]';
CC = [ones(size(h)) log(h)]\log(err); C(2) = CC(2);
err = [0.00073573,9.05975e-05,1.13596e-05]';
CC = [ones(size(h)) log(h)]\log(err); C(3) = CC(2);
err = [0.000176663,1.64135e-06,2.66024e-08]';
CC = [ones(size(h)) log(h)]\log(err); C(4) = CC(2);
err = [2.1282e-06,4.49046e-08,9.99912e-10]';
CC = [ones(size(h)) log(h)]\log(err); C(5) = CC(2);
C

% GQ (N+2) EC
h = [.5.^(1:3)]';
err = [0.106574,0.058359,0.0298728]';
CC = [ones(size(h)) log(h)]\log(err); C(1) = CC(2);
err = [0.00180028,0.000194939,2.4045e-05]';
CC = [ones(size(h)) log(h)]\log(err); C(2) = CC(2);
err = [0.000260859,3.76182e-05,3.86238e-06]';
CC = [ones(size(h)) log(h)]\log(err); C(3) = CC(2);
err = [3.57758e-06,9.58294e-08,2.94985e-09]';
CC = [ones(size(h)) log(h)]\log(err); C(4) = CC(2);
err = [3.6544e-07,4.95271e-09,2.42763e-10]';
CC = [ones(size(h)) log(h)]\log(err); C(5) = CC(2);
C


% GLL LF
err = [0.0562398,0.0151873];
C(1) = -(log(err(end))-log(err(end-1)))/(log(2));
err = [0.000748842,9.25456e-05];
C(2) = -(log(err(end))-log(err(end-1)))/(log(2));
err = [1.99064e-05,1.22357e-06];
C(3) = -(log(err(end))-log(err(end-1)))/(log(2));
err = [7.03865e-07,2.10265e-08];
C(4) = -(log(err(end))-log(err(end-1)))/(log(2));
err = [3.20127e-08,4.63639e-10];
C(5) = -(log(err(end))-log(err(end-1)))/(log(2));
C

% GQ N+2 LF
err = [0.00974763,0.00244539]; % N = 1
C(1) = -(log(err(end))-log(err(end-1)))/(log(2));
err = [0.00028425,3.54865e-05]; % N = 2
C(2) = -(log(err(end))-log(err(end-1)))/(log(2));
err = [7.39839e-06,4.6305e-07];
C(3) = -(log(err(end))-log(err(end-1)))/(log(2));
err = [1.37151e-07,4.16335e-09];
C(4) = -(log(err(end))-log(err(end-1)))/(log(2));
err = [2.34345e-09,3.65306e-11];
C(5) = -(log(err(end))-log(err(end-1)))/(log(2));
C




