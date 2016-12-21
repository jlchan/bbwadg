% compute quadrature nodes of Williams, et al

function [x,y,w,Ncub] = QNodes2D(N,cubDeg)

degrees = [2 3 5 6 8 9 11 13]; % degrees of integration
if nargin>1
    N = find(degrees >= cubDeg,1) - 1;   % -1 to index from 0    
    Ncub = (N+1)*(N+2)/2;
end

% total number of nodes
Np = (N+1)*(N+2)/2;

v(1,:) = 2*[-.5 -sqrt(3)/6];
v(2,:) = 2*[.5 -sqrt(3)/6];
v(3,:) = 2*[0 sqrt(3)/3];

switch N
    case 0
        a = [1/3 1/3 1/3];
        w = 1;
    case 1
        a = [2/3 1/6 1/6;
            1/6 2/3 1/6;
            1/6 1/6 2/3];
        w = [1/3 1/3 1/3];
    case 2
        a = [0.816847572980440 0.091576213509780 0.091576213509780;
            0.091576213509780 0.816847572980440 0.091576213509780;
            0.091576213509780 0.091576213509780 0.816847572980440;
            0.445948490915964 0.445948490915964 0.108103018168071;
            0.445948490915964 0.108103018168071 0.445948490915964;
            0.108103018168071 0.445948490915964 0.445948490915964];                

        w = [.109951743655333 0.109951743655333 0.109951743655333
            0.223381589678000 0.223381589678000 0.223381589678000];
    case 3        
        a = [0.888871894660413 0.055564052669793 0.055564052669793;
            0.055564052669793 0.888871894660413 0.055564052669793;
            0.055564052669793 0.055564052669793 0.888871894660413;
            0.295533711735893 0.634210747745723 0.070255540518384;
            0.295533711735893 0.070255540518384 0.634210747745723;
            0.070255540518384 0.295533711735893 0.634210747745723;
            0.634210747745723 0.295533711735893 0.070255540518384;
            0.634210747745723 0.070255540518384 0.295533711735893;
            0.070255540518384 0.634210747745723 0.295533711735893;
            0.333333333333333 0.333333333333333 0.333333333333333];            
        w = [0.041955512996649 
            0.041955512996649 
            0.041955512996649
            0.112098412070887 
            0.112098412070887 
            0.112098412070887 
            0.112098412070887 
            0.112098412070887 
            0.112098412070887
            0.201542988584730];

    case 4
        a = [0.928258244608533 0.035870877695734 0.035870877695734;
            0.035870877695734 0.928258244608533 0.035870877695734;
            0.035870877695734 0.035870877695734 0.928258244608533;
            0.516541208464066 0.241729395767967 0.241729395767967;
            0.241729395767967 0.516541208464066 0.241729395767967;
            0.241729395767967 0.241729395767967 0.516541208464066;
            0.474308787777079 0.474308787777079 0.051382424445843;
            0.474308787777079 0.051382424445843 0.474308787777079;
            0.051382424445843 0.474308787777079 0.474308787777079;
            0.201503881881800 0.751183631106484 0.047312487011716;
            0.201503881881800 0.047312487011716 0.751183631106484;
            0.047312487011716 0.201503881881800 0.751183631106484;
            0.751183631106484 0.201503881881800 0.047312487011716;
            0.751183631106484 0.047312487011716 0.201503881881800;
            0.047312487011716 0.751183631106484 0.201503881881800];
        w = [0.017915455012303 0.017915455012303 0.017915455012303
            0.127712195881265 0.127712195881265 0.127712195881265
            0.076206062385535 0.076206062385535 0.076206062385535
            0.055749810027115 0.055749810027115 0.055749810027115
            0.055749810027115 0.055749810027115 0.055749810027115];

    case 5
        a = [0.943774095634672 0.028112952182664 0.028112952182664;
            0.028112952182664 0.943774095634672 0.028112952182664;
            0.028112952182664 0.028112952182664 0.943774095634672;
            0.645721803061365 0.177139098469317 0.177139098469317;
            0.177139098469317 0.645721803061365 0.177139098469317;
            0.177139098469317 0.177139098469317 0.645721803061365;
            0.405508595867433 0.405508595867433 0.188982808265134;
            0.405508595867433 0.188982808265134 0.405508595867433;
            0.188982808265134 0.405508595867433 0.405508595867433;
            0.148565812270887 0.817900980028499 0.033533207700614;
            0.148565812270887 0.033533207700614 0.817900980028499;
            0.033533207700614 0.148565812270887 0.817900980028499;
            0.817900980028499 0.148565812270887 0.033533207700614;
            0.817900980028499 0.033533207700614 0.148565812270887;
            0.033533207700614 0.817900980028499 0.148565812270887;
            0.357196298615681 0.604978911775132 0.037824789609186;
            0.357196298615681 0.037824789609186 0.604978911775132;
            0.037824789609186 0.357196298615681 0.604978911775132;
            0.604978911775132 0.357196298615681 0.037824789609186;
            0.604978911775132 0.037824789609186 0.357196298615681;
            0.037824789609186 0.604978911775132 0.357196298615681];
        w = [0.010359374696538 0.010359374696538 0.010359374696538
            0.075394884326738 0.075394884326738 0.075394884326738
            0.097547802373242 0.097547802373242 0.097547802373242
            0.028969269372473 0.028969269372473 0.028969269372473
            0.028969269372473 0.028969269372473 0.028969269372473
            0.046046366595935 0.046046366595935 0.046046366595935
            0.046046366595935 0.046046366595935 0.046046366595935];
    case 6
        a = [0.960045625755613 0.019977187122193 0.019977187122193;
            0.019977187122193 0.960045625755613 0.019977187122193;
            0.019977187122193 0.019977187122193 0.960045625755613;
            0.736556464940005 0.131721767529998 0.131721767529998;
            0.131721767529998 0.736556464940005 0.131721767529998;
            0.131721767529998 0.131721767529998 0.736556464940005;
            0.333333333333333 0.333333333333333 0.333333333333333;
            0.485135346793461 0.485135346793461 0.029729306413079;
            0.485135346793461 0.029729306413079 0.485135346793461;
            0.029729306413079 0.485135346793461 0.485135346793461;
            0.107951981846011 0.867911210117951 0.024136808036039;
            0.107951981846011 0.024136808036039 0.867911210117951;
            0.024136808036039 0.107951981846011 0.867911210117951;
            0.867911210117951 0.107951981846011 0.024136808036039;
            0.867911210117951 0.024136808036039 0.107951981846011;
            0.024136808036039 0.867911210117951 0.107951981846011;
            0.270840772921567 0.700872570380723 0.028286656697710;
            0.270840772921567 0.028286656697710 0.700872570380723;
            0.028286656697710 0.270840772921567 0.700872570380723;
            0.700872570380723 0.270840772921567 0.028286656697710;
            0.700872570380723 0.028286656697710 0.270840772921567;
            0.028286656697710 0.700872570380723 0.270840772921567;
            0.316549598844617 0.536654684206138 0.146795716949245;
            0.316549598844617 0.146795716949245 0.536654684206138;
            0.146795716949245 0.316549598844617 0.536654684206138;
            0.536654684206138 0.316549598844617 0.146795716949245;
            0.536654684206138 0.146795716949245 0.316549598844617;
            0.146795716949245 0.536654684206138 0.316549598844617];
        w = [0.005272170280495
            0.005272170280495
            0.005272170280495
            0.044552936679504
            0.044552936679504
            0.044552936679504
            0.083608212215637
            0.033815712804198
            0.033815712804198
            0.033815712804198
            0.015710461340183
            0.015710461340183
            0.015710461340183
            0.015710461340183
            0.015710461340183
            0.015710461340183
            0.028205136280616
            0.028205136280616
            0.028205136280616
            0.028205136280616
            0.028205136280616
            0.028205136280616
            0.066995957127830
            0.066995957127830
            0.066995957127830
            0.066995957127830
            0.066995957127830
            0.066995957127830];
    case 7
        a = [0.957657154441070 0.021171422779465 0.021171422779465;
            0.021171422779465 0.957657154441070 0.021171422779465;
            0.021171422779465 0.021171422779465 0.957657154441070;
            0.798831205208225 0.100584397395888 0.100584397395888;
            0.100584397395888 0.798831205208225 0.100584397395888;
            0.100584397395888 0.100584397395888 0.798831205208225;
            0.457923384576135 0.271038307711932 0.271038307711932;
            0.271038307711932 0.457923384576135 0.271038307711932;
            0.271038307711932 0.271038307711932 0.457923384576135;
            0.440191258403832 0.440191258403832 0.119617483192335;
            0.440191258403832 0.119617483192335 0.440191258403832;
            0.119617483192335 0.440191258403832 0.440191258403832;
            0.101763679498021 0.879979641427232 0.018256679074748;
            0.101763679498021 0.018256679074748 0.879979641427232;
            0.018256679074748 0.101763679498021 0.879979641427232;
            0.879979641427232 0.101763679498021 0.018256679074748;
            0.879979641427232 0.018256679074748 0.101763679498021;
            0.018256679074748 0.879979641427232 0.101763679498021;
            0.394033271669987 0.582562022863673 0.023404705466341;
            0.394033271669987 0.023404705466341 0.582562022863673;
            0.023404705466341 0.394033271669987 0.582562022863673;
            0.582562022863673 0.394033271669987 0.023404705466341;
            0.582562022863673 0.023404705466341 0.394033271669987;
            0.023404705466341 0.582562022863673 0.394033271669987;
            0.226245530909229 0.751530614542782 0.022223854547989;
            0.226245530909229 0.022223854547989 0.751530614542782;
            0.022223854547989 0.226245530909229 0.751530614542782;
            0.751530614542782 0.226245530909229 0.022223854547989;
            0.751530614542782 0.022223854547989 0.226245530909229;
            0.022223854547989 0.751530614542782 0.226245530909229;
            0.635737183263105 0.249079227621332 0.115183589115563;
            0.635737183263105 0.115183589115563 0.249079227621332;
            0.115183589115563 0.635737183263105 0.249079227621332;
            0.249079227621332 0.635737183263105 0.115183589115563;
            0.249079227621332 0.115183589115563 0.635737183263105;
            0.115183589115563 0.249079227621332 0.635737183263105];
        w = [0.005639123786910
            0.005639123786910
            0.005639123786910
            0.027148968192278
            0.027148968192278
            0.027148968192278
            0.063100912533359
            0.063100912533359
            0.063100912533359
            0.051752795679899
            0.051752795679899
            0.051752795679899
            0.009866753574646
            0.009866753574646
            0.009866753574646
            0.009866753574646
            0.009866753574646
            0.009866753574646
            0.022008204800147
            0.022008204800147
            0.022008204800147
            0.022008204800147
            0.022008204800147
            0.022008204800147
            0.016644570076736
            0.016644570076736
            0.016644570076736
            0.016644570076736
            0.016644570076736
            0.016644570076736
            0.044326238118914
            0.044326238118914
            0.044326238118914
            0.044326238118914
            0.044326238118914
            0.044326238118914
            ];            
    case default
        error('No quad pts coded in')
end
XY = a*v;
x = XY(:,1);
y = XY(:,2);
w = 2.0*w/sum(w(:));
w = w(:);
return;
