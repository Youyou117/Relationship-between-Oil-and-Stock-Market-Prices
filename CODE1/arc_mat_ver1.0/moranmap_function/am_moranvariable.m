function am_moranvariable(hObject, eventdata)
% PURPOSE: sets the variable for use by moranmap
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
cvariable = results.cvariable;
WX = results.WX;
% define the quadrants
Q0=find(cvariable == 0 & WX(:,vindex)== 0);
Q1=find(cvariable>0 & WX(:,vindex)>0);
Q2=find(WX(:,vindex)>0 & cvariable<= 0);
Q3=find(WX(:,vindex)<= 0 & cvariable<= 0);
Q4=find(cvariable>0 & WX(:,vindex)<= 0);

%%%%%%%%%%%%%%%%%%%%%%%
results.Q1 = Q1;
results.Q2 = Q2;
results.Q3 = Q3;
results.Q4 = Q4;


results.vindex = vindex;
results.vflag = vflag;
% uses results.cvariable
results = am_moran(results);
thandles = results.thandles;


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
%     % cnt = 1;
% 	for i=1:results.npoly;
%       for k = 1:results.nparts(i);
%       tst = get(poly(i).handles(k),'Visible');
%           if strcmp(tst,'on');
%             set(poly(i).handles(k),'FaceColor',results.map_colors(i,:));
%             hi = get(poly(i).handles(k),'UserData');
%             set(hi,'Label',[num2str(i) ') ' num2str(cvariable(i,1))]);
%           end;
%       end;
%      if results.labels == 1
%        if strcmp(tst,'on');
%      %  thandles(i,1) = text(results.xc(i),results.yc(i),num2str(i));
% 	%   set(thandles(i,1),'fontsize',6);
%        end;
%      end;
% 	end;

if vflag == 1
set(poly(1).fig_handle,'Name',['Map of ' results.vnames(vindex,:)]);
elseif vflag == 0
set(poly(1).fig_handle,'Name',['Map of ' results.vnames(vindex+1,:)]);
end;

set(poly(1).fig_handle,'UserData',poly);

guidata(results.spop,results);
guidata(results.qpop,results);
guidata(results.vpop,results);

% bring moran figure back to the front
figure(results.moran_fig);

