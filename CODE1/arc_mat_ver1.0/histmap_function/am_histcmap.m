function am_histcmap(hObject, eventdata)
% PURPOSE: sets the colormap for the histmap figure
%
%   Color maps.
%     hsv        - Hue-saturation-value color map.
%     hot        - Black-red-yellow-white color map.
%     cool       - Shades of cyan and magenta color map.
%     gray       - Linear gray-scale color map.
%     bone       - Gray-scale with tinge of blue color map.
%     copper     - Linear copper-tone color map.
%     jet        - Variant of HSV.
%     pink       - Pastel shades of pink color map.
%     white      - All white color map.
%     autumn     - Shades of red and yellow color map.
%     spring     - Shades of magenta and yellow color map.
%     winter     - Shades of blue and green color map.
%     summer     - Shades of green and yellow color map.
 
%       'String', 'hsv|hot|cool|gray|bone|copper|jet|pink|white|autumn|spring|winter|summer',...

hfig = gcf;
poly = get(hfig,'UserData');

results = guidata(gcbo);
hpop = results.hpop;

val = get(hpop,'Value');
if val == 2
  cmap = 'hsv';
elseif val == 3
  cmap = 'hot';
elseif val == 4
  cmap = 'cool';
elseif val == 5
  cmap = 'gray';
elseif val == 6
  cmap = 'bone';
elseif val == 7
  cmap = 'copper';
elseif val == 8
  cmap = 'jet';
elseif val == 9
  cmap = 'pink';
elseif val == 10
  cmap = 'white';
elseif val == 11
  cmap = 'autumn';
elseif val == 12
  cmap = 'spring';
elseif val == 13
  cmap = 'winter';
elseif val == 14
  cmap = 'summer';
else
  cmap = results.cmap;
end;

results.cmap = cmap;
results = am_histlegend(results);

    cnt = 1;
	for i=1:results.npoly;
      for k = 1:results.nparts(i);
      tst = get(poly(i).handles(k),'Visible');
          if strcmp(tst,'on');
            set(poly(i).handles(k),'FaceColor',results.map_colors(cnt,:));
          end;
      end;
     tst = get(poly(i).handles(1),'Visible');
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

