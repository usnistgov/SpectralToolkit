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
        
        case ('fcontour3')
            % plots a 3D contour map of the objective function obj.f(x)
            % obj.plot('fcontour3',omega,resolution) 
            %   omega:  
            % an N x 2 matrix with N dimensions and the columns being the begin and
            % end points for each dimension.  All dimensions except 2 of
            % them must be the same begin and end point.
            %   resolution:
            % The number f points to plot in the x and y axis.  The Z axis
            % will have resolution^2 points
            
            % Later, can add a function handle to plot other function contour maps.
            %   fun:
            % function handle
            
            % check for the correct number of input arguments
            if (nargin - i) < 2
                fprintf ('Not enough arguments: obj.plot(''fcontour3'',omega,resolution,fun)\n')
                break
            end
            i = i+1; omega = varargin{i};
            i = i+1; resolution = varargin{i};
            
            % Later, can add a function handle to plot other function contour maps.
            %i = i+1; fun = varargin{i};
            %OMEGA = fcontour3(omega,resolution);
            
                % check the input arguments
                if (size(omega,2)~= 2)
                    warning('fcontour3 argument omega must be a matrix with 2 columns')
                    return
                else
                    if ~(isscalar(resolution))
                        warning('fcontour3 argument resolution must be a scalar integer')
                        return
                        %        else
                        %             if ~isa(fun,'function_handle')
                        %                 warning('fcontour3 argument fun must be a function handle')
                        %                 return
                        %             end
                    end
                end
                % final check, only two of the rows in omega can have unequal columns
                count = 0;
                for i = 1 : size(omega,1)
                    if omega(i,1)~=omega(i,2)
                        count = count+1;
                    end
                end
                if count~=2
                    warning('omega must have exactly 2 rows of unequal ranges')
                    return
                end
                
                % name the existing figure
                f = gcf;
                set(gca,'FontSize',18)
                f.Name = sprintf('%s_BSA',obj.Name);

                % Create vectors of omega values to check.  Also x and y will be the 
                % values to be plotted
                colOMEGA = size(omega,1);
                x = []; y = [];
                OMEGA = [];
                for i = 1 : colOMEGA
                    OMEGA = horzcat(OMEGA,linspace(omega(i,1),omega(i,2),resolution)');
                    if omega(i,1)~=omega(i,2)
                        if isempty(x)
                            x = OMEGA(:,i);
                            if i==1,xLbl='Carrier Freq (Hz)';else,xLbl='Phase (rad)';end
                        else
                            y = OMEGA(:,i);
                            if i==2,yLbl='Phase (rad)';else,yLbl='Delta Freq (Hz/s)';end
                            
                        end
                    end
                end
                
            count = 0;
            wb=waitbar(count,'Drawing Contour');
            z = zeros(resolution,resolution);
            for i = 1:resolution
                for j = 1:resolution
                    z(i,j) = obj.f([OMEGA(i,1),OMEGA(j,2),OMEGA(i,3)]);
                    count = count+1;
                    waitbar(count/resolution^2)
                end
            end
            close(wb)
            contour3(x,y,z,30,'-k')
            xlabel(xLbl),ylabel(yLbl),zlabel('f eval')
        
        % This may be obsolete if we are no longer debugging the hooke-jeeves search 
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


end
