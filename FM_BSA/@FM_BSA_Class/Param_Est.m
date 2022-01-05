function [estParams] = Param_Est(obj,endpt_BSA)
% determines the FM modulaton signal parameters from the endpoint of the
% BSA estimation and some of the intermediate values from the final
% objective function evaluation


% We are going to have to evaluate the objective functon one more time to
% get some of the intermediate values

nSamples = obj.nSamples;
mod_Freq = 2*pi*obj.Fm(obj.phase)*obj.dT;

t = double(0:nSamples-1)';   % time vector

% The hypothesis
GIJ = zeros(nSamples,3);
GIJ(:,1) = 1;       %DC
GIJ(:,2) = cos(endpt_BSA(1)*t+endpt_BSA(3)/mod_Freq*sin(mod_Freq*t+endpt_BSA(2)));
GIJ(:,3) = sin(endpt_BSA(1)*t+endpt_BSA(3)/mod_Freq*sin(mod_Freq*t+endpt_BSA(2)));


% pre allocate some shared variables to get the values we need
nFun = size(GIJ,2);
V_norm = [];
hi = zeros(nFun,1);
invD_norm = [];

% call the objective function one last time
y = -prob(GIJ); 

% Calculate the parameters
B = V_norm*invD_norm;
bi = B*hi;
a(1) = bi(1);
k = 2;
while k <= nFun
    if k < (nFun+1)/2+1
        a(k) = bi(k);
    else
        b(k-(nFun-1)/2) = bi(k);
    end
    k = k+1;
end

cBSA = complex(a,b);        % complex
Modulo_BSA = abs(cBSA);
Phi_BSA = -angle(cBSA);
Fcarr_BSA = endpt_BSA(1);
Phim_BSA = endpt_BSA(2);
dF_BSA = endpt_BSA(3);


% if in debug mode, show the best fit and residual
if obj.debug
    % determine the best-fit signal
    n = double(0:nSamples-1)'*obj.dT;   % time vector
    Result = Modulo_BSA(2).*cos(Fcarr_BSA/obj.dT.*n + dF_BSA/mod_Freq * sin(mod_Freq/obj.dT.*n+Phim_BSA)+Phi_BSA(2));
    %figure(obj.fig);obj.fig=obj.fig+1;
    subplot(2,1,1)
    plot(n,obj.Samples(:,obj.phase),n,Result)
    subplot(2,1,2)
    plot(n,obj.Samples(:,obj.phase)-Result)
end

estParams = struct('Modulo_BSA',Modulo_BSA,...
                   'Phi_BSA',Phi_BSA,...
                   'Fcarr_BSA',Fcarr_BSA,...
                   'Phim_BSA',Phim_BSA,...
                   'dF_BSA',dF_BSA);

%%=========================================================================
% nested functions
    function stloge = prob(GIJ)
        iNo = double(nSamples);
        HIJ = ortho(GIJ);    % orthogonolize the objective function
        %nFun = size(HIJ,2);
        %hi = zeros(nFun,1);
        for j=1:nFun
            h1 = sum(obj.Samples(:,obj.phase).*HIJ(:,j));
            hi(j)=h1;
        end
        h2 = sum(hi.^2);
        h2bar = h2/nFun;
        
        y2 = sum(obj.Samples(:,obj.phase).^2);
        y2bar = y2/iNo;
        
        qq = 1-h2/y2;
        if qq<=0
            qq = 1e-6;
        end
        
        stloge = log(qq)*(nFun-iNo)/2;
        dif = y2bar-nFun*h2bar/iNo;
        sigma = sqrt(abs(dif)*iNo/(iNo-nFun-2));
    end

    function HIJ = ortho(GIJ)
        % Orthogonalizes the hypothesis function
        M = GIJ'*GIJ;
        M = (M+M.')/2; % The matlab eig function must recognize the matrix as symetrical
        [V,D_vector] = eig(M,'vector');         % eiganvalues and eiganvectors
        SqrSumCol = sum(V.^2);
        norm = sqrt(SqrSumCol);
        %V_norm = V./norm;                       % matlab 2015 has issues with elementwise divide
        V_norm = bsxfun(@rdivide,V,norm);       % this is the workaround
        D_vec = sqrt(abs(D_vector));
        D_norm = diag(D_vec);
        invD_norm = inv(D_norm);
        A = GIJ*V_norm;
        HIJ = A*invD_norm;
    end


end