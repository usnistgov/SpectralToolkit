function [y] = objFun(obj, x)
obj.funEvals = obj.funEvals+1;
nSamples = obj.nSamples;
mod_Freq = 2*pi*obj.Fm(obj.phase)*obj.dT;

omega = x;  % [carrier freq, modulation phase, delta_Freq;
t = double(0:nSamples-1)';   % time vector

% The hypothesis
GIJ = zeros(nSamples,3);
GIJ(:,1) = 1;       %DC
GIJ(:,2) = cos(omega(1)*t+omega(3)/mod_Freq*sin(mod_Freq*t+omega(2)));
GIJ(:,3) = sin(omega(1)*t+omega(3)/mod_Freq*sin(mod_Freq*t+omega(2)));

y = -prob(GIJ); 

%%=========================================================================
% nested functions
    function stloge = prob(GIJ)
        iNo = double(nSamples);
        HIJ = ortho(GIJ);    % orthogonolize the objective function
        nFun = size(HIJ,2);
        hi = zeros(nFun,1);
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