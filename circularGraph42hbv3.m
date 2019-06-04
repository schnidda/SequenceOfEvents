classdef circularGraph42hbv3 < handle
% circularGraph42hbv3 Plot an interactive circular graph to illustrate connections in a network.
%adopted from Matlab
%% Syntax
% circularGraph42hbv3(X)
% circularGraph42hbv3(X,'PropertyName',propertyvalue,...)
% h = circularGraph42hbv3(...)
%
%% Description
% A 'circular graph' is a visualization of a network of nodes and their
% connections. The nodes are laid out along a circle, and the connections
% are drawn within the circle. Click on a node to make the connections that
% emanate from it more visible or less visible. Click on the 'Show All'
% button to make all nodes and their connections visible. Click on the
% 'Hide All' button to make all nodes and their connections less visible.
%
% Required input arguments.
% X : A symmetric matrix of numeric or logical values.
%
% Optional properties.
% Colormap : A N by 3 matrix of [r g b] triples, where N is the 
%            length(adjacenyMatrix).
% Label    : A cell array of N strings.
%%
% Copyright 2016 The MathWorks, Inc.
  properties
    Node = node(0,0); % Array of nodes
    ColorMap;         % Colormap
    Label;            % Cell array of strings
    ShowButton;       % Turn all nodes on
    HideButton;       % Turn all nodes off
    
    Incorp;           %incorporation times for modificatin
  end
  
  methods
    function this = circularGraph42hbv3(adjacencyMatrix,varargin)
      % Constructor
      p = inputParser;
      %this is modified to fit  the 42hbv3 part
      defaultColorMap = jet(140);%length(adjacencyMatrix));
      maxT=max(adjacencyMatrix);
      minT=min(adjacencyMatrix);
      
      
      defaultLabel = cell(length(adjacencyMatrix));
      for i = 1:length(defaultLabel)
        if mod(i,500)==0
        defaultLabel{i} = num2str(i);
        elseif sum(adjacencyMatrix(i,:))~=0
        defaultLabel{i} = num2str(i);            
        else
        defaultLabel{i} = ''; %deleted space as empty label
        end
      end
      
      addRequired(p,'adjacencyMatrix',@(x)(isnumeric(x) || islogical(x)));
      addParameter(p,'ColorMap',defaultColorMap,@(colormap)length(colormap) == length(adjacencyMatrix));
      addParameter(p,'Label'   ,defaultLabel   ,@iscell);
      addParameter(p,'Incorp'   ,@iscell);
    
      parse(p,adjacencyMatrix,varargin{:});
      this.ColorMap = p.Results.ColorMap;
      this.Label    = p.Results.Label;
      
      %parse(p, Incorp,varargin{:});
      this.Incorp =p.Results.Incorp;
      
      
      
      this.ColorMap;
      
      this.ShowButton = uicontrol(...
        'Style','pushbutton',...
        'Position',[0 40 80 40],...
        'String','Show All',...
        'Callback',@circularGraph42hbv3.showNodes,...
        'UserData',this);
      
      this.HideButton = uicontrol(...
        'Style','pushbutton',...
        'Position',[0 0 80 40],...
        'String','Hide All',...
        'Callback',@circularGraph42hbv3.hideNodes,...
        'UserData',this);
      
      fig = gcf;
      set(fig,...
        'UserData',this,...
        'CloseRequestFcn',@circularGraph42hbv3.CloseRequestFcn);
      
      % Draw the nodes
      delete(this.Node);
      t = linspace(-pi,pi,length(adjacencyMatrix) + 1).'; % theta for each node
      extent = zeros(length(adjacencyMatrix),1);
      for i = 1:length(adjacencyMatrix)
            this.Node(i) = node(cos(t(i)),sin(t(i)));
            %auskommentiert um die nodes schwarz zu machen this.Node(i).Color = this.ColorMap(i,:);
          %if mod(i,10)==0
            this.Node(i).Label = this.Label{i};
         %end
      end
      
      % Find non-zero values of s and their indices
      [row,col,v] = find(adjacencyMatrix);
      % Calculate line widths based on values of s (stored in v).
      minLineWidth  = 0.5;
      lineWidthCoef = 1; %5;
      lineWidth = v./max(v);
      if sum(lineWidth) == numel(lineWidth) % all lines are the same width.
        lineWidth = repmat(minLineWidth,numel(lineWidth),1);
      else % lines of variable width.
        lineWidth = lineWidthCoef*lineWidth + minLineWidth;
      end
      
      % Draw connections on the Poincare hyperbolic disk.
      %
      % Equation of the circles on the disk:
      % x^2 + y^2 
      % + 2*(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1))*x 
      % - 2*(u(1)-v(1))/(u(1)*v(2)-u(2)*v(1))*y + 1 = 0,
      % where u and v are points on the boundary.
      %
      % Standard form of equation of a circle
      % (x - x0)^2 + (y - y0)^2 = r^2
      %
      % Therefore we can identify
      % x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
      % y0 = (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
      % r^2 = x0^2 + y0^2 - 1
     
      
      
      
      
      %this.ColorMap
      %round(v/max(v)*140)
      global actualcolors;
      %actualcolors=this.ColorMap(141-round(v/max(v)*140),:);
      
      %modification here: overwrite min and max with global min/max from
      %both 5p and 3p
      
      % case A: global values
      %min(min(this.Incorp)), max(max(this.Incorp))
      %140-round(((v-min(min(this.Incorp)))/(max(max(this.Incorp))-min(min(this.Incorp))))*139)
      %actualcolors=this.ColorMap(140-round(((v-min(min(this.Incorp)))/(max(max(this.Incorp))-min(min(this.Incorp))))*139),:);
     
      %listforcolorbar=140-round(((v-min(min(this.Incorp)))/(max(max(this.Incorp))-min(min(this.Incorp))))*139)
                    %global incorporation;
 
      
      % case B: each on its own
      v
      min(v), max(v)
%%this part is changed to have the same coloring as in the python script,
%%slow incorporation should be red (jet colormap starts blue to red)

      actualcolors=this.ColorMap(1+round(((v-min(v))/(max(v)-min(v)))*139),:)
      listforcolorbar=1+round(((v-min(v))/(max(v)-min(v)))*139);
                %listforcolorbar=141-round(((v)/(max(v)))*140);

      %actualcolors=this.ColorMap(140-round(((v-min(v))/(max(v)-min(v)))*139),:)
      %listforcolorbar=140-round(((v-min(v))/(max(v)-min(v)))*139);
                %listforcolorbar=141-round(((v)/(max(v)))*140);
            
         
      colors_sorted=actualcolors;     
      global listforcolorbar;
      listforcolorbar
      
      
      
      
      for i = 1:length(v)
        if row(i) ~= col(i)
          if abs(row(i) - col(i)) - length(adjacencyMatrix)/2 == 0 
            % points are diametric, so draw a straight line
            u = [cos(t(row(i)));sin(t(row(i)))];
            v = [cos(t(col(i)));sin(t(col(i)))];
            this.Node(row(i)).Color = actualcolors(i,:);
            this.Node(row(i)).MarkerSize=20;
            %this.Node(i).Label = num2str(i);
            this.Node(row(i)).Connection(end+1) = line(...
              [u(1);v(1)],...
              [u(2);v(2)],...
              'LineWidth', lineWidth(i),... 
               'Color', actualcolors(i,:),...
              'PickableParts','none');
          else % points are not diametric, so draw an arc
            u  = [cos(t(row(i)));sin(t(row(i)))];
            v  = [cos(t(col(i)));sin(t(col(i)))];
            this.Node(row(i)).Color = actualcolors(i,:);
            this.Node(row(i)).MarkerSize=20;
            %this.Node(i).Label = num2str(i);
            x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            r  = sqrt(x0^2 + y0^2 - 1);
            thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
            thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
            
            if u(1) >= 0 && v(1) >= 0 
              % ensure the arc is within the unit disk
              theta = [linspace(max(thetaLim),pi,50),...
                       linspace(-pi,min(thetaLim),50)].';
            else
              theta = linspace(thetaLim(1),thetaLim(2)).';
            end
            this.Node(row(i)).Connection(end+1) = line(...
              r*cos(theta)+x0,...
              r*sin(theta)+y0,...
              'LineWidth', lineWidth(i),...
               'Color', actualcolors(i,:),...
              'PickableParts','none');
          end
        end
      end
      
      axis image;
      ax = gca;
      for i = 1:length(adjacencyMatrix)
        extent(i) = this.Node(i).Extent;
      end
      extent = max(extent(:));
      ax.XLim = ax.XLim + extent*[-1 1];
      fudgeFactor = 1.75; % Not sure why this is necessary. Eyeballed it.
      ax.YLim = ax.YLim + fudgeFactor*extent*[-1 1];
      ax.Visible = 'off';
      ax.SortMethod = 'depth';
 
      
      
      [row2,col2,v2] = find(adjacencyMatrix);

      fig = gcf;
      fig.Color = [1 1 1];
      %% flipud changes the direction of the colormap and therefore the matching betweening incorporation and color
      %colormap (flipud(jet))
      colormap(jet)
      cb = colorbar('Direction','reverse')%(listforcolorbar);%(colors_sorted);%(actualcolors)
      set(cb,'position',[.8 .3 .01 .4])
      %A)this is changed to compare 5' and 3'
      %caxis([min(min(this.Incorp)),max(max(this.Incorp))]);

      %B) 5' and 3' individual
      caxis([min(v2),max(v2)]);

    end
    
  end
  
  methods (Static = true)
    function showNodes(this,~)
      % Callback for 'Show All' button
      n = this.UserData.Node;
      for i = 1:length(n)
        n(i).Visible = true;
      end
    end
    
    function hideNodes(this,~)
      % Callback for 'Hide All' button
      n = this.UserData.Node;
      for i = 1:length(n)
        n(i).Visible = false;
      end
    end
    
    function CloseRequestFcn(this,~)
      % Callback for figure CloseRequestFcn
      c = this.UserData;
      for i = 1:length(c.Node)
        delete(c.Node(i));
      end
      delete(gcf);
    end
    
  end
  
end