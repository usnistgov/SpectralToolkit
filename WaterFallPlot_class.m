classdef WaterFallPlot_class
    % A waterfall plot is a figure object.  It is a 2D heatmap with time on the T axis and
    % the independent variable(such as frequency) on the X axis.  Data is added in strips
    
    properties
        figHandle      % Hanfle to the figure object
        axHandle
        timeWidth      % (seconds) Time will be on the Y axis, timeWidth duration to be displayed in the plot
        stripWidth     % (seconds) The width of each strip of heatmap data
        currentTime = 0;    %  time of the last written strip
        stripLength = -1;    % number of elements in a strip
        tData          % the time vector
        xyData          % 2D array of the xData over time
        colorMap
        flipColorMap
        colorScale
        yLabel
        xLabel
    end
    
    methods
        % Class constructor and destructor
        function obj = WaterFallPlot_class(varargin)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            
            defaultFigNum = [];
            defaultTimeWidth = 1;
            defaultStripWidth = 1/60;  
            defaultColorMap = gray;
            defaultFlipColorMap = true;
            defaultYLabel = 'Time (s)';
            defaultXLabel = 'Frequency (Hz)';
            defaultColorScale = 'linear';
            
            p = inputParser;
            
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            
            
            addParameter(p,'FigNum',defaultFigNum,validScalarPosNum)
            addParameter(p,'TimeWidth',defaultTimeWidth,validScalarPosNum)
            addParameter(p,'StripWidth',defaultStripWidth,validScalarPosNum)
            addParameter(p,'ColorMap',defaultColorMap,@ischar)
            addParameter(p,'FlipColorMap',defaultFlipColorMap,@islogical)
            addParameter(p,'YLabel',defaultYLabel,@ischar)
            addParameter(p,'XLabel',defaultXLabel,@ischar)
            addParameter(p,'ColorScale',defaultColorScale,@ischar)
            

            
            
            parse(p,varargin{:})
            
            % create the figure object
            if isempty(p.Results.FigNum)
                obj.figHandle = figure();
            else
                obj.figHandle = figure(p.Results.FigNum);
            end
            obj.axHandle = gca;
            obj.yLabel = p.Results.YLabel;
            obj.xLabel = p.Results.XLabel;
            ytickformat(obj.axHandle,'%3.2f')
            obj.timeWidth = p.Results.TimeWidth;
            obj.stripWidth = p.Results.StripWidth;
            obj.colorMap = p.Results.ColorMap;
            obj.flipColorMap = p.Results.FlipColorMap;
            if obj.flipColorMap
                colormap(flipud(obj.colorMap));
            else
                colormap(obj.colorMap);
            end
            obj.colorScale = p.Results.ColorScale;

            
            % initialize the table that will be shown as a heatmap
            obj.tData = -obj.timeWidth:obj.stripWidth:0-obj.stripWidth;
            
        end
        
        function delete(obj)
           close(obj.figHandle) 
        end
                
    end
    
    methods (Access = public)
        % public methods
        function obj = addStrip(obj,stripData,xData)
            % add a strip of data to the cuurentTime + stripWidth
            
            % check the length of xData equals obj.stripLength
            %   if obj.stripLength <=  0, set stripLength to the length of xData
            if obj.stripLength <= 0
                obj.stripLength = numel(stripData); 
                obj.xyData = ones(length(obj.tData),length(stripData))*min(stripData);                
            end
            if ~(numel(stripData)==obj.stripLength)
                error('xData length %d must equal stripLengh property %d',numel(xData),obj.stripLength)
            end
            
            % append the stripData to the bottom of xData
            obj.xyData(1:end-1,:) = obj.xyData(2:end,:);
            obj.xyData(end,:) = stripData;
            imagesc(obj.xyData)
            
            % Y axis tick marks
            obj.tData = obj.tData+obj.stripWidth;
            ticks = obj.tData(obj.axHandle.YTick);
            obj.axHandle.YTickLabel = string(round(ticks*100)/100);
            
            % varargin will hold the X axis tickmarks if it is there
            ticks = xData(obj.axHandle.XTick);
            obj.axHandle.XTickLabel = string(ticks);
               
            % axis labels
            ylabel(obj.yLabel);
            xlabel(obj.xLabel);
            
            set(obj.axHandle,'ColorScale',obj.colorScale)
            
            
        end
        
    end
    
end

