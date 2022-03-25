function plot(obj,varargin)
%PLOTS plot override for the HilbertHuang_class

i = 1;
while i<= numel(varargin)
    switch varargin{i}
        case 'IMFs'
            fig=figure(obj.fig); obj.fig = obj.fig+1;
            set(fig,'NumberTitle','off','Name',sprintf('%s All IMFs #%d',obj.IMFs.Name))            
            nImfs = size(obj.IMFs.Data,2);
            for j = 1:nImfs
                subplot(nImfs+1,1,j)
                plot(obj.IMFs.Time,obj.IMFs.Data(:,j))
                title(sprintf('%s %d',obj.IMFs.Name,j))
                xlabel('time (s)')
                ylabel('ampl')
            end
            subplot(nImfs+1,1,j+1)
            plot(obj.IMFs.Time,obj.Residual)
            title('Residual')
            xlabel('time (s)')
            ylabel('ampl')
            i = i+1;
            
        case 'Hilbert'
            
            for j=1:size(obj.Hilbert,2)
                fig=figure(obj.fig); obj.fig = obj.fig+1;
                set(fig,'NumberTitle','off','Name',sprintf('%s IMF #%d',obj.IMFs.Name,j))
                
                
                subplot(4,1,1)                
                plot(obj.IMFs.Time,real(obj.Hilbert{1,j}),'-b',...
                     obj.IMFs.Time,imag(obj.Hilbert{1,j}),'-r')
                ylabel('analytic')
                xlabel('time (s)')
                title(sprintf('IMF %d', j))
                
                subplot(4,1,2)
                plot(obj.IMFs.Time,obj.Hilbert{2,j},'-k')
                ylabel('energy')
                xlabel('time (s)')
                
                subplot(4,1,3)
                plot(obj.IMFs.Time,obj.Hilbert{3,j},'-k')                
                ylabel('angle')
                xlabel('time (s)')                
                
                subplot(4,1,4)
                plot(obj.IMFs.Time,obj.Hilbert{4,j},'-k')                
                ylabel('freq (Hz)')
                xlabel('time (s)') 
            end
            i=i+1;
            
        case 'Spectrum'
            % creates a spectrum plot of all of the IMFs
            fig=figure(obj.fig); obj.fig = obj.fig+1;     
            set(fig,'NumberTitle','off','Name',sprintf('%s Hilbert Spectrum',obj.IMFs.Name))
                        
            tVec = obj.IMFs.Time; %time vector
            fRange = obj.hilOpts.FreqLimits;
            xyrange = [tVec(1),tVec(end),fRange(1),fRange(2)];
            axis(xyrange);
            cla
            
            
            for j=1:size(obj.Hilbert,2)            
                patch([tVec(1);tVec;tVec(end)],[0;obj.Hilbert{4,j};0],[nan;obj.Hilbert{2,j};nan],...
                    'EdgeColor','interp','EdgeAlpha','interp',...
                    'FaceColor', 'none', 'FaceVertexAlphaData',[nan;obj.Hilbert{2,j};nan],...
                    'LineWidth', 2, 'FaceAlpha', 'interp');
                    
            end
                        
            xlabel('time(s)')
            ylabel('frequency (Hz)')
            title('Hilbert Spectrum')
            colormap parula
            colorbar
            i = i + 1;
            
         otherwise
            warning('unrecognized plot argument: %s',varargin{i})
            i = i+1;
    end
                



end

