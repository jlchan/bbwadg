function [rfq sfq tfq wfq] = tet_surface_cubature(N)

[rqt sqt wqtri] = Cubature2D(N);

rfq = []; sfq = []; tfq = [];
e = ones(size(rqt,1),1);

rfqf = rqt; sfqf = sqt; tfqf = -e;
rfq = [rfq rfqf]; sfq = [sfq sfqf]; tfq = [tfq tfqf];

rfqf = rqt; sfqf = -e; tfqf = sqt;
rfq = [rfq rfqf]; sfq = [sfq sfqf]; tfq = [tfq tfqf];

rfqf = -(1 + rqt + sqt); sfqf = rqt; tfqf = sqt; 
rfq = [rfq rfqf]; sfq = [sfq sfqf]; tfq = [tfq tfqf];

rfqf = -e; sfqf = rqt; tfqf = sqt;
rfq = [rfq rfqf]; sfq = [sfq sfqf]; tfq = [tfq tfqf];

rfq = rfq(:); sfq = sfq(:); tfq = tfq(:);
wfq = [wqtri; wqtri; wqtri*sqrt(2); wqtri];

rfq = rfq(:);
sfq = sfq(:);
tfq = tfq(:);
wfq = wfq(:);