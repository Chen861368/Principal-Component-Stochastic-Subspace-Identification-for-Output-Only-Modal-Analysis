function [fn, zeta, phi, varargout] = PCSSI(y, dt, rows_f, cols_f, rows_p, cols_p, varargin)
% PCSSI - Identifies the modal parameters of an M-DOF system from response data.
%
% Syntax:
%   [fn, zeta, phi, varargout] = PCSSI(y, dt, rows_f, cols_f, rows_p, cols_p, varargin)
%
% Description:
%   This function extracts the modal properties of a multi-degree-of-freedom (M-DOF) 
%   system from ambient vibration response data using a data-driven principal component subspace identification method.
%
% Inputs:
%   y       - Time-series data of ambient vibrations, given as an [M × N] matrix,
%             where M is the number of sensors and N is the number of time steps.
%   dt      - Scalar representing the time step between consecutive samples.
%   rows_f  - Number of rows for the matrix Yf, which is constructed and used in subsequent calculations.
%   cols_f  - Number of columns for the matrix Yf, determining its dimensions along with rows_f.
%   rows_p  - Number of rows for the matrix Yp, a key matrix in the identification process.
%   cols_p  - Number of columns for the matrix Yp, defining its size along with rows_p.
%   varargin - Optional input parameters:
%       'Nmin'       - Scalar, minimum model order.
%       'Nmax'       - Scalar, maximum model order.
%       'eps_freq'   - Scalar, frequency accuracy threshold.
%       'eps_zeta'   - Scalar, damping ratio accuracy threshold.
%       'eps_MAC'    - Scalar, MAC (Modal Assurance Criterion) accuracy threshold.
%       'eps_cluster'- Scalar, maximum distance within each modal cluster.
%
% Outputs:
%   fn        - Identified natural frequencies.
%   zeta      - Identified modal damping ratios.
%   phi       - Identified mode shapes.
%   varargout - Structured output containing data for the stabilization diagram.
%
% -------------------------------------------------------------------------
% Author: Biqi Chen
% Note: The clustering and stabilization diagram algorithm is adapted from the original implementation by E. Cheynet.
% -------------------------------------------------------------------------

%% Parse optional input arguments and set default values
p = inputParser();
p.CaseSensitive = false;
p.addOptional('Nmin', 2);
p.addOptional('Nmax', 30);
p.addOptional('eps_freq', 1e-2);
p.addOptional('eps_zeta', 4e-2);
p.addOptional('eps_MAC', 5e-3);
p.addOptional('eps_cluster', 0.2);
p.parse(varargin{:});

% Ensure the number of output arguments is between 3 and 4
nargoutchk(3, 4)

% Get the number of sensors (rows) from input matrix y
[Nyy, ~] = size(y);

% Shorten variable names for convenience
eps_freq = p.Results.eps_freq;
eps_zeta = p.Results.eps_zeta;
eps_MAC = p.Results.eps_MAC;
eps_cluster = p.Results.eps_cluster;
Nmin = p.Results.Nmin;
Nmax = p.Results.Nmax;

% Transpose data for Hankel matrix construction
data = y';

% Adjust rows_p to ensure it meets the minimum model order requirement
if rows_p * Nyy < Nmax
    rows_p = floor(Nmax / Nyy) + 1;
end

% Ensure cols_p matches cols_f
if cols_p ~= cols_f
    cols_p = cols_f;
end

%% Step 1: Extract the required columns from the data matrix
Yp_data = data(1:rows_p + cols_p - 1, :);  % Extract rows 1 to (i+j-1)
Yf_data = data(rows_p + 1:end, :);         % Extract rows (i+1) to the last row

%% Step 2: Construct Hankel matrices
Yp = create_hankel_matrix(Yp_data, rows_p, cols_p);
Yf = create_hankel_matrix(Yf_data, rows_f, cols_f);

%% Step 3: Perform economical SVD decomposition on Yp and denoise
[~, ~, V_s] = svd_decompose(Yp, Nmax);

%% Step 4: Select instrumental variables and apply subspace-based denoising
Yf_IV = instrumental_variable_denoise(Yf, V_s', cols_f);

%% Step 5: Perform economical SVD decomposition on the denoised Yf_IV
[U_k, Sigma_k, V_k] = svd_decompose(Yf_IV, Nmax);


if isnan(U_k)
    fn = nan;
    zeta = nan;
    phi = nan;
    if nargout==4
        varargout = {nan};
    end
    return
end
% Stability check
kk=1;
for ii=Nmax:-1:Nmin % decreasing order of poles
    if kk==1
        [fn0,zeta0,phi0] = modalID(U_k,Sigma_k,V_k,ii,dt,Nyy);
    else
        [fn1,zeta1,phi1] = modalID(U_k,Sigma_k,V_k,ii,dt,Nyy);
        [a,b,c,d,e] = stabilityCheck(fn0,zeta0,phi0,fn1,zeta1,phi1);
        fn2{kk-1}=a;
        zeta2{kk-1}=b;
        phi2{kk-1}=c;
        MAC{kk-1}=d;
        stablity_status{kk-1}=e;
        fn0=fn1;
        zeta0=zeta1;
        phi0=phi1;
    end
    kk=kk+1;
end
% sort for increasing order of poles
stablity_status=fliplr(stablity_status);
fn2=fliplr(fn2);
zeta2=fliplr(zeta2);
phi2=fliplr(phi2);
MAC=fliplr(MAC);
% get only stable poles
[fnS,zetaS,phiS,MACS] = getStablePoles(fn2,zeta2,phi2,MAC,stablity_status);
if isempty(fnS)
    warning('No stable poles found');
    fn = nan;
    zeta = nan;
    phi = nan;
    if nargout==4
        varargout = {nan};
    end
    return
end
% Hierarchical cluster
[fn3,zeta3,phi3] = myClusterFun(fnS,zetaS,phiS);
if isnumeric(fn3)
    warning('Hierarchical cluster failed to find any cluster');
    fn = nan;
    zeta = nan;
    phi = nan;
    if nargout==4
        varargout = {nan};
    end
    return
end
% average the clusters to get the frequency and mode shapes
% Up to Nmax parameters are identified
fn = zeros(1,Nmax);
zeta = zeros(1,Nmax);
phi = zeros(Nmax,Nyy);
for ii=1:numel(fn3)
    fn(ii)=nanmean(fn3{ii});
    zeta(ii)=nanmean(zeta3{ii});
    phi(ii,:)=nanmean(phi3{ii},2);
end
phi(fn==0,:)=[];
zeta(fn==0)=[];
fn(fn==0)=[];
% sort the eigen frequencies
[fn,indSort]=sort(fn);
zeta = zeta(indSort);
phi = phi(indSort,:);
% varargout for stabilization diagram
if nargout==4
    paraPlot.status=stablity_status;
    paraPlot.Nmin = Nmin;
    paraPlot.Nmax = Nmax;
    paraPlot.fn = fn2;
    varargout = {paraPlot};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fn, zeta, phi] = modalID(U_k, Sigma_k, V_k, t, dt, m)
    % modalID - Identifies the modal properties of a system.
    %
    % Syntax:
    %   [fn, zeta, phi] = modalID(U_k, Sigma_k, V_k, t, dt, m)
    %
    % Inputs:
    %   U_k     - Left singular matrix obtained from the Singular Value Decomposition (SVD) of the denoised data Yf_IV.
    %   Sigma_k - Singular value matrix from the SVD of Yf_IV.
    %   V_k     - Right singular matrix from the SVD of Yf_IV.
    %   t       - Model order, used to extract a subset of U_k, Sigma_k, and V_k.
    %   dt      - Time step, used to compute the frequencies of the system poles.
    %   m       - Number of sensors, used for computing the system and observation matrices.
    %
    % Outputs:
    %   fn  - Identified natural frequencies.
    %   zeta - Identified damping ratios.
    %   phi  - Identified mode shapes.
    
    % Extract the first t columns of U_k
    Ut = U_k(:, 1:t);

    % Extract the first t singular values and form a diagonal matrix
    St = diag(Sigma_k(1:t));

    % Extract the first t rows of V_k
    Vt = V_k(1:t, :);

    % Compute system and observation matrices
    [A, C] = compute_system_matrices(Ut, St, Vt, m);

    % Compute eigenvalues and eigenvectors of the system matrix A
    [Vi, Di] = eig(A);

    % Compute poles
    mu = log(diag(Di)) ./ dt;

    % Extract natural frequencies
    fn = abs(mu(2:2:end)) / (2 * pi);

    % Compute damping ratios
    zeta = -real(mu(2:2:end)) ./ abs(mu(2:2:end));

    % Compute mode shapes
    phi = real(C * Vi);
    phi = phi(:, 2:2:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, C_1] = compute_system_matrices(Ut, St, Vt, m)
    % compute_system_matrices - Computes the system and observation matrices.
    %
    % Syntax:
    %   [A, C_1] = compute_system_matrices(Ut, St, Vt, m)
    %
    % Inputs:
    %   Ut - First t columns of the left singular matrix U.
    %   St - Diagonal matrix containing the first t singular values.
    %   Vt - First t rows of the right singular matrix V.
    %   m  - Number of sensors.
    %
    % Outputs:
    %   A   - System matrix.
    %   C_1 - First m rows of the observation matrix.
    
    % Compute the observation matrix C = Ut * St
    C = Ut * St;
    C_1 = C(1:m, :); % Extract the first m rows

    % Compute the system matrix A
    X_2 = Vt(:, 1:end-1); % Equivalent to Vt[:, :-1] in Python
    X_1 = Vt(:, 2:end);   % Equivalent to Vt[:, 1:] in Python
    A = X_2 * X_1';       % Compute the system matrix A
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [fn,zeta,phi,MAC,stablity_status] = stabilityCheck(fn0,zeta0,phi0,fn1,zeta1,phi1)
        % [fn,zeta,phi,MAC,stablity_status] = stabilityCheck(fn0,zeta0,phi0,fn1,zeta1,phi1)
        % calculate the stability status of each mode obtained for
        % two adjacent poles (i,j).
        %
        % Input:
        % fn0: eigen frequencies calculated for pole i: vetor of N-modes [1 x N]
        % zeta0: modal damping ratio for pole i: vetor of N-modes [1 x N]
        % phi0: mode shape for pole i: vetor of N-modes [Nyy x N]
        % fn1: eigen frequencies calculated for pole j: vetor of N-modes [1 x N+1]
        % zeta1: modal damping ratio for pole j: vetor of N-modes [1 x N+1]
        % phi1: mode shape for pole j: vetor of N-modes [Nyy x N+1]
        %
        % Output:
        % fn: eigen frequencies calculated for pole j
        % zeta:  modal damping ratio for pole i
        % phi:mode shape for pole i
        % MAC: Mode Accuracy
        % stablity_status: stabilitystatus
        %%
        
        % Preallocation
        stablity_status = [];
        fn = [];
        zeta = [];
        phi = [];
        MAC=[];
        % frequency stability
        N0 = numel(fn0);
        N1 = numel(fn1);
        for rr=1:N0
            for jj=1:N1
                stab_fn = errCheck(fn0(rr),fn1(jj),eps_freq);
                stab_zeta = errCheck(zeta0(rr),zeta1(jj),eps_zeta);
                [stab_phi,dummyMAC] = getMAC(phi0(:,rr),phi1(:,jj),eps_MAC);
                % get stability status
                if stab_fn==0,
                    stabStatus = 0; % new pole
                elseif stab_fn == 1 & stab_phi == 1 & stab_zeta == 1,
                    stabStatus = 1; % stable pole
                elseif stab_fn == 1 & stab_zeta ==0 & stab_phi == 1,
                    stabStatus = 2; % pole with stable frequency and vector
                elseif stab_fn == 1 & stab_zeta == 1 & stab_phi ==0,
                    stabStatus = 3; % pole with stable frequency and damping
                elseif stab_fn == 1 & stab_zeta ==0 & stab_phi ==0,
                    stabStatus = 4; % pole with stable frequency
                else
                    error('Error: stablity_status is undefined')
                end
                fn = [fn,fn1(jj)];
                zeta = [zeta,zeta1(jj)];
                phi = [phi,phi1(:,jj)];
                MAC = [MAC,dummyMAC];
                stablity_status = [stablity_status,stabStatus];
            end
        end
        
        [fn,ind] = sort(fn);
        zeta = zeta(ind);
        phi = phi(:,ind);
        MAC = MAC(ind);
        stablity_status = stablity_status(ind);
        
        function y = errCheck(x0,x1,eps)
            if or(numel(x0)>1,numel(x1)>1),
                error('x0 and x1 must be a scalar');
            end
            if abs(1-x0./x1)<eps % if frequency for mode i+1 is almost unchanged
                y =1;
            else
                y = 0;
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [fnS,zetaS,phiS,MACS] = getStablePoles(fn,zeta,phi,MAC,stablity_status)
        fnS = [];zetaS = [];phiS=[];MACS = [];
        for oo=1:numel(fn)
            for jj=1:numel(stablity_status{oo})
                if stablity_status{oo}(jj)==1
                    fnS = [fnS,fn{oo}(jj)];
                    zetaS = [zetaS,zeta{oo}(jj)];
                    phiS = [phiS,phi{oo}(:,jj)];
                    MACS = [MACS,MAC{oo}(jj)];
                end
            end
        end
        
        % remove negative damping
        fnS(zetaS<=0)=[];
        phiS(:,zetaS<=0)=[];
        MACS(zetaS<=0)=[];
        zetaS(zetaS<=0)=[];
        
        % Normalized mode shape
        for oo=1:size(phiS,2)
            phiS(:,oo)= phiS(:,oo)./max(abs(phiS(:,oo)));
            if diff(phiS(1:2,oo))<0
                phiS(:,oo)=-phiS(:,oo);
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [fn,zeta,phi] = myClusterFun(fn0,zeta0,phi0)
        
        [~,Nsamples] = size(phi0);
        pos = zeros(Nsamples,Nsamples);
        for i1=1:Nsamples
            for i2=1:Nsamples
                [~,MAC0] = getMAC(phi0(:,i1),phi0(:,i2),eps_MAC); % here, eps_MAC is not important.
                pos(i1,i2) = abs((fn0(i1)-fn0(i2))./fn0(i2)) +1-MAC0; % compute MAC number between the selected mode shapes
            end
            
        end
        
        if numel(pos)==1
            warning('linkage failed: at least one distance (two observations) are required');
            fn = nan;
            zeta = nan;
            phi = nan;
            return
        else
            Z =  linkage(pos,'single','euclidean');
            myClus = cluster(Z,'Cutoff',eps_cluster,'Criterion','distance');
            Ncluster = max(myClus);
            
            ss=1;
            fn = {}; zeta = {}; phi = {};
            for rr=1:Ncluster
                if numel(myClus(myClus==rr))>5
                    dummyZeta = zeta0(myClus==rr);
                    dummyFn = fn0(myClus==rr);
                    dummyPhi = phi0(:,myClus==rr);
                    valMin = max(0,(quantile(dummyZeta,0.25) - abs(quantile(dummyZeta,0.75)-quantile(dummyZeta,0.25))*1.5));
                    valMax =quantile(dummyZeta,0.75) + abs(quantile(dummyZeta,0.75)-quantile(dummyZeta,0.25))*1.5;
                    dummyFn(or(dummyZeta>valMax,dummyZeta<valMin)) = [];
                    dummyPhi(:,or(dummyZeta>valMax,dummyZeta<valMin)) = [];
                    dummyZeta(or(dummyZeta>valMax,dummyZeta<valMin)) = [];
                    fn{ss} = dummyFn;
                    zeta{ss} = dummyZeta;
                    phi{ss} = dummyPhi;
                    ss=ss+1;
                end
            end
            if isempty(fn)
                fn = nan;
                zeta = nan;
                phi = nan;
                return
            end
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [y,dummyMAC] = getMAC(x0,x1,eps)
        Num = abs(x0(:)'*x1(:)).^2;
        D1= x0(:)'*x0(:);
        D2= x1(:)'*x1(:);
        dummyMAC = Num/(D1.*D2);
        if dummyMAC >(1-eps)
            y = 1;
        else
            y = 0;
        end
    end
end


function hankel_matrix = create_hankel_matrix(data, rows, cols)
    % Generate a Hankel matrix from the given data.
    % A Hankel matrix is a structured matrix in which each ascending diagonal from left to right has constant values.
    % It is commonly used in system identification and signal processing.

    % Get the dimensions of the input data.
    [n_samples, n_features] = size(data);
    
    % Initialize the Hankel matrix with zeros.
    hankel_matrix = zeros(cols, rows * n_features);

    % Populate the Hankel matrix.
    for i = 1:cols
        for j = 1:rows
            if i + j - 1 <= n_samples
                hankel_matrix(i, (j - 1) * n_features + 1 : j * n_features) = data(i + j - 1, :);
            end
        end
    end

    % Transpose the matrix to match the standard Hankel matrix format.
    hankel_matrix = hankel_matrix';
end


function [U, Sigma, Vt] = svd_decompose(A, k)
    % Perform an economical Singular Value Decomposition (SVD) on matrix A, retaining the top k singular values.
    % This function applies SVD-based dimensionality reduction by keeping only the leading k singular values.
    
    [U, Sigma, Vt] = svd(A, 'econ'); % Compute the economical SVD of A.
    
    % Extract the first k components of U, Sigma, and Vt.
    U = U(:, 1:k);
    Sigma = diag(Sigma(1:k, 1:k)); % Convert the kxk diagonal matrix into a vector of singular values.
    Vt = Vt'; % Transpose Vt to match the standard convention.
    Vt = Vt(1:k, :);
end

function Yf_IV = instrumental_variable_denoise(Yf, V_s, j)
    % Denoise the signal Yf using instrumental variables.
    % The denoised signal is obtained by computing (1/j) * Yf * V_s * V_s^T.
    
    Yf_IV = (1 / j) * Yf * V_s * V_s';
end

