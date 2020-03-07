function am_histmap(hObject, eventdata)
% PURPOSE: uses rubber band rectangle to select a sub-region of the histmap

hfig = gcf;
poly = get(hfig,'UserData');

results = guidata(gcbo);
spop = results.spop;
val = get(spop,'Value');
vindex = results.vindex;
missing = results.missing;

% ============================================================
if val == 2 % we zoom in
% ============================================================

% recover stuff from the userdata structure
	
variable = results.variable; % we need the original variable vector here

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

set(gca,'Units','Normalized');

new_variable = [];
new_missing = [];

	for i=1:results.npoly;
	ind = find(results.xmin(i) > xlim1 & results.xmax(i) < xlim2 & results.ymin(i) > ylim1 & results.ymax(i) < ylim2);
     if length(ind) == 0
      for k = 1:results.nparts(i);
      set(poly(i).handles(k),'Visible','off');
      end;
      
     else
     new_variable = [new_variable
                     variable(i,:)]; % we zoom all variables so the user can switch variables for the zoomed map
     new_missing = [new_missing
                    missing(i,1)];   % we update the missing varables index to reflect zooming
     end;
	end;
    
% the user botched the selection, so we reset to the full map
if isempty(new_variable);
new_variable = variable;
new_missing = missing;
thandles = zeros(results.npoly,1);
hold on;
	for i=1:results.npoly;
      for k = 1:results.nparts(i);
              set(poly(i).handles(k),'Visible','on');
              hi = get(poly(i).handles(k),'UserData');
              set(hi,'Label',[num2str(i) ') ' num2str(new_variable(i,vindex))]);
      end;
	end;
end;
hold off
results.cvariable = new_variable(:,vindex); % we send down the current variable to create the legend
results.svariable = new_variable; % contains all variables zoomed
results.cmissing = new_missing;

results = am_histlegend(results);


    cnt = 1;
	for i=1:results.npoly;
      for k = 1:results.nparts(i);
      tst = get(poly(i).handles(k),'Visible');
          if strcmp(tst,'on');
              set(poly(i).handles(k),'FaceColor',results.map_colors(cnt,:));
              hi = get(poly(i).handles(k),'UserData');
              set(hi,'Label',[num2str(cnt) ') ' num2str(results.svariable(cnt,vindex))]);
              if k == results.nparts(i)
              cnt = cnt+1;
              end;
          end;
      end;
	end;

% update data structure
set(poly(1).fig_handle,'UserData',poly);

guidata(results.spop,results);
guidata(results.cpop,results);
guidata(results.hpop,results);
guidata(results.vpop,results);

% bring legend figure back to the front
figure(results.legend_fig);

% ============================================================
elseif val == 3 % we zoom back out
% ============================================================

results.cvariable = results.variable(:,vindex);
results.svariable = results.variable;
results.cmissing = results.missing;

results = am_histlegend(results);

figure(poly(1).fig_handle);

set(gca,'Units','Normalized');

	for i=1:results.npoly;
      for k = 1:results.nparts(i);
          set(poly(i).handles(k),'FaceColor',results.map_colors(i,:),'Visible','on');
          hi = get(poly(i).handles(k),'UserData');
          set(hi,'Label',[num2str(i) ') ' num2str(results.cvariable(i))]);
      end;
	end;

% update the results structure passed along in UserData
set(poly(1).fig_handle,'UserData',poly);
guidata(results.spop,results);
guidata(results.cpop,results);
guidata(results.hpop,results);
guidata(results.vpop,results);

% bring legend figure back to the front
figure(results.legend_fig);

end; 