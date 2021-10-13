function plot(obj, varargin)
%PLOT functions for the Hook_Jeeves_Student3 class
%   Detailed explanation goes here

i = 1;
while i<= numel(varargin)
    switch varargin{i}
        case 'NLS'
            figure(obj.fig); obj.fig = obj.fig+1;
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
            figure(obj.fig); obj.fig = obj.fig+1;
            t=1:length(obj.data);            
            subplot(2,1,1)
            plot(t,obj.data1,'.b',t,obj.Result_BSA,'.r')
            title('BSA Result')
            xlabel('Samples')
            ylabel('Amplitude')
            legend('data','Result\_BSA')
            
            subplot(2,1,2)
            plot(t,obj.Residue_BSA,'.k')
            title('BSA Residue')
            xlabel('Samples')
            ylabel('Amplitude')                             
        
        otherwise
            warning('unrecognized plot argument: %s',varargin{i})
    end
    i=i+1;


end

