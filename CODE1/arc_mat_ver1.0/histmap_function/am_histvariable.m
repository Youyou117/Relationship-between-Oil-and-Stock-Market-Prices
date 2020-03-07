function am_histvariable(hObject, eventdata)
% PURPOSE: sets the variable for use by histmap
%
%       'String', 'variables'

hfig = gcf;
poly = get(hfig,'UserData');

results = guidata(gcbo);

vpop = results.vpop;
vnames = results.vnames;
variable = results.svariable; % all variables possibly zoomed

vindex = get(vpop,'Value');
vindex = vindex - 1;
vflag = results.vflag;

results.cvariable = variable(:,vindex); % current variable selection
results.vindex = vindex;
% uses results.cvariable
results = am_histlegend(results);

    cnt = 1;
	for i=1:results.npoly;
      for k = 1:results.nparts(i);
      tst = get(poly(i).handles(k),'Visible');
          if strcmp(tst,'on');
            set(poly(i).handles(k),'FaceColor',results.map_colors(cnt,:));
            hi = get(poly(i).handles(k),'UserData');
            set(hi,'Label',[num2str(i) ') ' num2str(results.cvariable(cnt))]);
          end;
      end;
     if strcmp(tst,'on');
     cnt = cnt + 1;
     end;
	end;

if vflag == 1
set(poly(1).fig_handle,'Name',['Map of ' results.vnames(vindex,:)]);
elseif vflag == 0
set(poly(1).fig_handle,'Name',['Map of ' results.vnames(vindex+1,:)]);
end;

set(poly(1).fig_handle,'UserData',poly);

guidata(results.spop,results);
guidata(results.cpop,results);
guidata(results.hpop,results);
guidata(results.vpop,results);

% bring legend figure back to the front
figure(results.legend_fig);

