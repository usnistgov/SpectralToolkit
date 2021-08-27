function plot(obj,varargin)
%PLOTS plot override for the HilbertHuang_class

i = 1;
while i<= numel(varargin)
    switch varargin{i}
        case 'IMFs'
            figure(obj.fig); obj.fig = obj.fig+1;
            nImfs = size(obj.IMFs.Data,2);
            for i = 1:nImfs
                subplot(nImfs+1,1,i)
                plot(obj.IMFs.Time,obj.IMFs.Data(:,i))
                title(sprintf('%s %d',obj.IMFs.Name,i))
                xlabel('time (s)')
                ylabel('ampl')
            end
            subplot(nImfs+1,1,i+1)
            plot(obj.IMFs.Time,obj.Residual)
            title('Residual')
            xlabel('time (s)')
            ylabel('ampl')
            i = i+1;
            
        case 'Hilbert'
            fig=figure(obj.fig); obj.fig = obj.fig+1;
            set(fig,'NumberTitle','off','Name',sprintf('%s IMF #%d',obj.IMFs.Name,i))
            
            for i=1:size(obj.Hilbert,2)
                subplot(4,1,1)                
                plot(obj.IMFs.Time,real(obj.Hilbert{1,i}),'-b')
                hold on
                plot(obj.IMFs.Time,imag(obj.Hilbert{1,i}),'-r')
                hold off
                ylabel('analytic')
                xlabel('time (s)')
                title(sprintf('IMF %d', i))
                
                subplot(4,1,2)
                plot(obj.IMFs.Time,obj.Hilbert{2,i},'-k')
                ylabel('energy')
                xlabel('time (s)')
                
                subplot(4,1,3)
                plot(obj.IMFs.Time,obj.Hilbert{3,i},'-k')                
                ylabel('angle')
                xlabel('time (s)')                
                
                subplot(4,1,4)
                plot(obj.IMFs.Time,obj.Hilbert{4,i}/(2*pi),'-k')                
                ylabel('freq (Hz)')
                xlabel('time (s)') 
            end
            
            
            i=i+1;
            
            
            
        otherwise
            warning('unrecognized plot argument: %s',varargin{i})
            i = i+1;
    end
                



end

