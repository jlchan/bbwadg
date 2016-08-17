function [rq sq tq wq Fmaskq wfq] = tet_cubature_TW(N)

a = JacobiGL(0,0,N+1); wa1D = sum(inv(Vandermonde1D(N+1,a)*Vandermonde1D(N+1,a)'),2);
%     [a wa1D] = JacobiGQ(0,0,N);
[b wb1D] = JacobiGR(1,0,N);
[c wc1D] = JacobiGR(2,0,N);
[aq bq cq] = meshgrid(a,b,c);
aq = aq(:); bq = bq(:); cq = cq(:);
[rq sq tq] = abctorst(aq,bq,cq);
[wa wb wc] = meshgrid(wa1D,wb1D,wc1D);
wq = wa(:).*wb(:).*wc(:);
wq = (4/3)*wq/sum(wq);

% find all the nodes that lie on each edge
NODETOL = 1e-8;
fmask1   = find( abs(1+tq) < NODETOL)'; 
fmask2   = find( abs(1+sq) < NODETOL)';
fmask3   = find( abs(1+rq+sq+tq) < NODETOL)';
fmask4   = find( abs(1+rq) < NODETOL)';
Fmaskq  = [fmask1(:);fmask2(:);fmask3(:);fmask4(:)];
% Fmaskq  = [fmask1(:) fmask2(:) fmask3(:) fmask4(:)];

% f = 1: a,b, c = -1
[wfa wfb] = meshgrid(wa1D,wb1D);
wf1 = wfa(:).*wfb(:)*wc1D(1);
wf1 = 2*wf1/sum(wf1);

% f = 2: a,c, b = -1
[wfa wfc] = meshgrid(wa1D,wc1D);
wf2 = wfa(:).*wfc(:)*wb1D(1);
wf2 = 2*wf2/sum(wf2);

% f = 3: b,c, a = 1
[wfb wfc] = meshgrid(wb1D,wc1D);
wf3 = wfb(:).*wfc(:)*wa1D(end);
wf3 = 2*sqrt(2)*wf3/sum(wf3);

% f = 4: b,c, a = -1
[wfb wfc] = meshgrid(wb1D,wc1D);
wf4 = wfb(:).*wfc(:)*wa1D(1);
wf4 = 2*wf4/sum(wf4);

wfq = [wf1(:);wf2(:);wf3(:);wf4(:)];

% plot3(aq,bq,cq,'o')
% rfq = rq(fmask1); sfq = sq(fmask1); tfq = tq(fmask1);
% keyboard
