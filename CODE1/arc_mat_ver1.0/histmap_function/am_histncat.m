function am_histncat(hObject, eventdata)
% PURPOSE: sets the # of histogram bars in histmap
%
%       'String', '# categories|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20',...

hfig = gcf;
poly = get(hfig,'UserData');

results = guidata(gcbo);
cpop = results.cpop;

val = get(cpop,'Value');
if val == 2
  nbc = 2;
elseif val == 3
  nbc = 3;
elseif val == 4
  nbc = 4;
elseif val == 5
  nbc = 5;
elseif val == 6
  nbc = 6;
elseif val == 7
  nbc = 7;
elseif val == 8
  nbc = 8;
elseif val == 9
  nbc = 9;
elseif val == 10
  nbc = 10;
elseif val == 11
  nbc = 11;
elseif val == 12
  nbc = 12;
elseif val == 13
  nbc = 13;
elseif val == 14
  nbc = 14;
elseif val == 15
  nbc = 15;
elseif val == 16
  nbc = 16;
elseif val == 17
  nbc = 17;
elseif val == 18
  nbc = 18;
elseif val == 19
  nbc = 19;
elseif val == 20
  nbc = 20;
else 
  nbc = 5;
end;

results.nbc = nbc;
results = am_histlegend(results);

    cnt = 1;
	for i=1:results.npoly;
      for k = 1:results.nparts(i);
      tst = get(poly(i).handles(k),'Visible');
          if strcmp(tst,'on');
            set(poly(i).handles(k),'FaceColor',results.map_colors(cnt,:));
          end;
      end;
     if strcmp(tst,'on');
     cnt = cnt + 1;
     end;

	end;


set(poly(1).fig_handle,'UserData',poly);

guidata(results.spop,results);
guidata(results.cpop,results);
guidata(results.hpop,results);
guidata(results.vpop,results);

% bring legend figure back to the front
figure(results.legend_fig);

