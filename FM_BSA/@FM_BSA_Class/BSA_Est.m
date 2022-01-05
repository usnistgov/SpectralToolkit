function [endpt_BSA] = BSA_Est(obj,startPt)
% uses Nelder-Meade to optimize the objective BSA function 

% Using fminsearch(problem)where problem is a structure.  see MATLAB "doc fminsearch"
opts = optimset(@fminsearch);    % configures NM options with default values
opts.TolX = 1e-8;
opts.TolFun = 1e8;
PROBLEM = struct('objective',[],'X0',[],'options',opts,'solver','fminsearch');
X0 = startPt;
PROBLEM.X0 = X0;
PROBLEM.objective = @(X0) obj.objFun(X0);

if obj.debug
    PROBLEM.options.OutputFcn = @out;             % plots the contour map and points
    [endpt_BSA,~,~,~] = fminsearch(PROBLEM);
    hold off
else
    [endpt_BSA] = fminsearch(PROBLEM);
end

% verbose status ------------------------------------------------------
if obj.verbose
    fprintf('Total iterations: %d\n',funEvals(phase))
    nfun = (iFun-1)/2;
    fList = '';
    for m = 1:nfun
        fItem = sprintf('%d: %f \n',m,endpt_BSA(m));
        fList = sprintf('%s %s',fList,fItem);
    end
    fprintf('BSA Frequencies (rad/unit):\n%s',fList);
end
% end verbose ---------------------------------------------------------

if obj.debug,hold off,end

 function stop = out(x, optimValue, state)
        % fminsearch output function for ploting the trajectory during NM minimization
        stop = false;
        switch state
            case 'init'
                %hold off
                %figure(obj.fig);obj.fig=obj.fig+1;
                f = 2*pi*obj.Fcarr(obj.phase)*obj.dT;
                dF = 2*pi*obj.Fm(obj.phase)*obj.Km(obj.phase)*obj.dT;
                obj.fcontour3([f,f;-pi,pi;0,2*dF],30,@obj.objFun)
                hold on;
            case 'iter'
                plot3((x(2)),abs(x(3)),optimValue.fval,'.');                
                drawnow
        end
    end    

end
