% contour plots for debugging
    function fcontour3(omega,resolution,fun)
        % plots a 3D contour map of the objective function f(x)
        % fcontour3(omega,fun, res)
        %
        %   omega:
        % an N x 2 matrix with N dimensions and the columns being the begin and
        % end points for each dimension.  All dimensions except 2 of
        % them must be the same begin and end point.
        % omega(:,1): carrier frequency range
        % omega(:,2): modulation phase range
        % omega(:,3): pean frequency deviation
        %
        %   res:
        % The number f points to plot in the x and y axis.  The Z axis
        % will have resolution^2 points
        %
        %   fun:
        % Handle to the objective function
        % Create vectors of omega values to check.  Also x and y will be the
        % values to be plotted
        zp = zeros(resolution,resolution);
        colOMEGA = size(omega,1);
        x = []; y = [];
        OMEGA = [];
        f_dF = false;
        for ip = 1 : colOMEGA
            OMEGA = horzcat(OMEGA,linspace(omega(ip,1),omega(ip,2),resolution)');
            if omega(ip,1)~=omega(ip,2)
                if isempty(x)
                    x = OMEGA(:,ip);
                    if ip==1
                        xLbl='Carrier Freq (Hz)';
                        f_dF = true;
                    else
                        xLbl='Phase (rad)';
                    end
                else
                    y = OMEGA(:,ip);
                    if ip==2
                        yLbl='Phase (rad)';
                    else
                        yLbl='Frequency Deviation (Hz)';
                    end
                end
            end
        end
        count = 0;
        wb=waitbar(count,'Drawing Contour');
        for ip = 1:resolution
            for jp = 1:resolution
                if f_dF
                    zp(ip,jp) = fun([OMEGA(jp,1),OMEGA(ip,2),OMEGA(ip,3)]);
                else
                    zp(ip,jp) = fun([OMEGA(ip,1),OMEGA(jp,2),OMEGA(ip,3)]);
                end
                count = count+1;
                waitbar(count/resolution^2)
            end
        end
        close(wb)
        contour3(x,y,zp,resolution,'-k')
        xlabel(xLbl),ylabel(yLbl),zlabel('objective function value')
    end
