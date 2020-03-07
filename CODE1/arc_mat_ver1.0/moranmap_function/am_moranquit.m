function am_moranquit(hObject, eventdata)
% PURPOSE: closes figure windows
%

hfig = gcf;
poly = get(hfig,'UserData');

results = guidata(gcbo);
qpop = results.qpop;
lfig = results.moran_fig;

val = get(qpop,'Value');
if val == 2
close(hfig);
close(lfig);
end;
%clf reset;
close all;