function plotlabel(r,s,t)

dr = abs(max(r(:))-min(r(:)));
if nargin==2
    plot(r,s,'o');
    text(r+.05*dr,s,num2str((1:length(r(:)))'))
else
    plot3(r,s,t,'o');
    text(r+.05*dr,s,t,num2str((1:length(r(:)))'))
end