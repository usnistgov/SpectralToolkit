function plot(obj, varargin)
%PLOT functions for the Hook_Jeeves_Student3 class
%  Examples: 
%       HJS.plot('NLS') % plots modulating signal NLS analysis results
%       HJS.plot('BSA') % plots modulated signal BSA analysis results
%       HJS.plot('NLS','BSA') % plots both analysis results in separate figures

i = 1;
while i<= numel(varargin)
    switch varargin{i}
        case 'NLS'
            f = figure(obj.fig); obj.fig = obj.fig+1;
            f.Name = sprintf('%s_NLS',obj.Name);
            subplot(2,1,1)
            t=1:length(obj.data);
            plot(t,obj.data,'-b',t,obj.Result_NLS,'-r')
            title('NLS Best Fit')
            ylabel('Amplitude')
            xlabel('Samples')
            legend('data','Result\_NLS')
            
            subplot(2,1,2)
            plot(t,obj.Residue_NLS,'.r')
            msg = sprintf('Residue, THD = %f',obj.THD_NLS);
            title(msg)
            ylabel('Amplitude')
            xlabel('Samples')
            
        case('BSA')
             f = figure(obj.fig); obj.fig = obj.fig+1;
            f.Name = sprintf('%s_BSA',obj.Name);
            t=1:length(obj.data);            
            subplot(2,1,1)
            plot(t,obj.data1,'-b',t,obj.Result_BSA,'-r')
            title('BSA Result')
            xlabel('Samples')
            ylabel('Amplitude')
            legend('data','Result\_BSA')
            
            subplot(2,1,2)
            plot(t,obj.Residue_BSA,'.k')
            title('BSA Residue')
            xlabel('Samples')
            ylabel('Amplitude')                             
        
         case('hookeContour')             
             x = cell2mat(obj.hookeContour(:,1));
             fx = cell2mat(obj.hookeContour(:,2));
             f = figure(obj.fig); obj.fig = obj.fig+1;
             plot(x,fx,'.k')
             w0 = cell2mat(obj.hookeContour(:,3));
             ka = cell2mat(obj.hookeContour(:,4));
             phi = cell2mat(obj.hookeContour(:,5));
             
             f = figure(obj.fig); obj.fig = obj.fig+1;
             subplot(3,1,1)
             plot(x,w0,'.k')
             %ylim([.0204,.0206])
             title('omega(0)')
            
             subplot(3,1,2)
             plot(x,ka,'.k')
             %ylim([1.568,1.574])
             title('omega(1)')
              
             subplot(3,1,3)
             plot(x,phi,'.k')
             %ylim([0.008,0.012])
             title('omega(2)')
             
             %plot3(w0-mean(w0),ka-mean(ka),phi-mean(phi))
             
             
        otherwise
            warning('unrecognized plot argument: %s',varargin{i})
    end
    i=i+1;


end

