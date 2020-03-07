function am_selectmap(hObject, eventdata)
% PURPOSE: uses rubber band rectangle to select a sub-region of the map

hfig = gcf;
poly = get(hfig,'UserData');

results = guidata(gcbo);
spop = results.spop;
val = get(spop,'Value');
thandles = results.thandles;
%missing = results.cmissing;

% ============================================================
if val == 2 % we zoom in
% ============================================================
	
figure(poly(1).fig_handle);

    k = waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
    y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];

	xlim1 = min(x);
	xlim2 = max(x);
	ylim1 = min(y);
	ylim2 = max(y);

% =======================================
maxx = xlim2;
maxy = ylim2;
minx = xlim1;
miny = ylim1;

figure(poly(1).fig_handle);

% =======================================

obs_selected = [];

	for i=1:results.npoly;
	ind = find(results.xmin(i) > xlim1 & results.xmax(i) < xlim2 & results.ymin(i) > ylim1 & results.ymax(i) < ylim2);
     if length(ind) == 0
      for k = 1:results.nparts(i);
      set(poly(i).handles(k),'Visible','off');
      end;
     else
     obs_selected = [obs_selected
                     i];
     end;
	end;


if isempty(obs_selected); % the user botched a selection
tmp = 1:results.npoly;    % so we save them by resetting to all of the data
obs_selected = tmp';
	for i=1:results.npoly;
      for k = 1:results.nparts(i);
      set(poly(i).handles(k),'Visible','on');
      end;
	end;
end;

results.obs_selected = obs_selected; % we send down an index to the selected obs
results = am_moran(results);

	for i=1:results.npoly;
      for k = 1:results.nparts(i);
      tst = get(poly(i).handles(k),'Visible');
          if strcmp(tst,'on');
          set(poly(i).handles(k),'FaceColor',results.map_colors(i,:));
          end;
      end;
      if results.labels == 1
          tst = get(poly(i).handles(k),'Visible');
          if strcmp(tst,'on');
          % do nothing
          else
          set(thandles(i),'Visible','off');
          end;
      end;
	end;

% update data structure
set(poly(1).fig_handle,'UserData',poly);
guidata(results.spop,results);

% bring legend figure back to the front
figure(results.moran_fig);


% ============================================================
elseif val == 3 % we zoom back out
% ============================================================

tmp = 1:results.npoly;
obs_selected = tmp';
results.obs_selected = obs_selected;
results.cmissing = results.missing;
results = am_moran(results);

figure(poly(1).fig_handle);

	for i=1:results.npoly;
      for k = 1:results.nparts(i);
          set(poly(i).handles(k),'FaceColor',results.map_colors(i,:),'Visible','on');
          set(poly(i).handles(k),'EdgeColor',[0 0 0]);
      end;
      if results.labels == 1
      set(thandles(i),'Visible','on');
      end;
	end;

% update the results structure passed along in UserData
set(poly(1).fig_handle,'UserData',poly);
guidata(results.spop,results);

% bring legend figure back to the front
figure(results.moran_fig);

end; 