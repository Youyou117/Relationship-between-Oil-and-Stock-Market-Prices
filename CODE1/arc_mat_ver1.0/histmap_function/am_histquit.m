function am_histquit(hObject, eventdata)
% PURPOSE: closes histmap figure windows
%

hfig = gcf;
poly = get(hfig,'UserData');

results = guidata(gcbo);
qpop = results.qpop;
cpop = results.cpop;
hpop = results.hpop;
lfig = results.legend_fig;

val = get(qpop,'Value');
if val == 2
close(hfig);
close(lfig);
end;
%clf reset;
close all;