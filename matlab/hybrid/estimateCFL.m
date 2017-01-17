% estimate CFL for a given mesh
function [C CK] = estimateCFL()

hybridgGlobals3D
hybridgGlobalFlags

% precomputed trace constants
CH = [ 9.0000   18.0000   30.0000   45.0000   63.0000   84.0000  108.0000];
CW = [ 9.926135933275306  18.563576705381948  29.030325215439685  42.988345972840094 58.802145509223905  78.006158337859318  99.271490513773074];
CP = 1.0e+02 *[ 0.116837586022311   0.208869058826159   0.328373010378272   0.475870076602500    0.651676634297143   0.855997828354466   1.088964369062309];
CT = [ 12.222306675423649  20.461046774606864  29.175757913298817  41.653800010846517 54.447306071124011  71.104042033768295  88.317462093087173];

CK = zeros(K,1);
for ee = 1:length(hexK)    
    e = hexK(ee);
    JsH = JsB(1:NfcH,e); JH = J(1:NcH,e);
    bCH = CH(N) * max(JsH(:)) * max(1./JH(:));
    CK(e) = bCH;
end

for ee = 1:length(wedgK)
    e = wedgK(ee);    
    JsW = JsB(1:NfcW,e);  JW = J(1:NcW,e);
    if useLSC
        % figure out trace inequality here
        bCW = CW(N) * max(JsW(:))*max(1./JW(:)); % wrong one!
    else
        bCW = CW(N) * max(JsW(:))*max(1./JW(:));
    end
    CK(e) = bCW;
end

for ee = 1:length(pyrK)
    e = pyrK(ee);        
    JsP = JsB(1:NfcP,e); JP = J(1:NcP,e);
    bCP = CP(N) * max(JsP(:)) * max(1./JP(:));
    CK(e) = bCP;
end

for ee = 1:length(tetK)
    e = tetK(ee);        
    JsT = JsB(1:NfcT,e); JT = J(1:NcT,e);
    bCT = CT(N) * max(JsT(:)) * max(1./JT(:));
    CK(e) = bCT;
end

C = max(CK);