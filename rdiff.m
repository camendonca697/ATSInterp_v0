function x1 = rdiff(x, t, method, sigma, strategy, t1, x1max)
% RDIFF   Regularised numerical differentiation.
%    X1 = RDIFF(X, T) returns an estimate of the first derivative of X;
%    if T is a scalar, it is used as spacing between abscissa; if T is
%    a vector, its elements are used as the abscissa. If X is a matrix,
%    RDIFF differentiates each of its columns independently.
%    
%    X1 = RDIFF(X, T, METHOD) uses method selected among following ones:
%       'totalvar'  based on total-variation regularisation
%       'tikhonov'  based on Tikhonov regularisation
%       'savitzky'  based on Savitzky-Golay filter
%       'kalman'    based on Kalman filter
%       'central'   based on central-difference formula with extended
%                   differentiation step
%       'sinstep'   based on modified central-difference formula
%       'expstep'   based on modified central-difference formula
%       'nbgauss'   based on smoothing approximation using Gaussian basis
%                   functions
%       'nbcos'     based on smoothing approximation using even powers
%                   of cosine functions
%       'tsvd'      based on truncated singular value decomposition
%       'landweber' based on Landweber iteration
%    By default, the 'totalvar' method is used.
%    
%    X1 = RDIFF(X, T, METHOD, SIGMA), where SIGMA is a scalar, uses SIGMA
%    as an estimate of the variance of the errors corrupting the elements
%    of X under the assumption that those errors are independent and
%    identically distributed. If SIGMA is an NxN matrix, it is used as the
%    covariance matrix of those errors.
%    
%    X1 = RDIFF(X, T, METHOD, SIGMA, STRATEGY) uses a selected strategy for
%    the optimisation of the regularisation parameter(s); see documentation
%    of M-files corresponding to particular methods ('rdiff_[METHOD].m')
%    for lists of available strategies.
%    
%    X1 = RDIFF(X, T, METHOD, SIGMA, STRATEGY, T1), if METHOD is 'nbgauss'
%    or 'nbcos', estimates the derivative at abscissa specified by T1.
%    For methods other than 'nbgauss' and 'nbcos', T1 is ignored.
%    
%    X1 = RDIFF(X, T, METHOD, SIGMA, STRATEGY, [], X1MAX), if METHOD is
%    'landweber', uses X1MAX as an a priori upper bound for absolute value
%    of the first derivative. For methods other than 'landweber', X1MAX is
%    ignored.
% 
%    Jakub Wagner, February 3, 2020
%    Institute of Radioelectronics and Multimedia Technology
%    Warsaw University of Technology

% Check input arguments
if ~exist('method', 'var') || isempty(method)
    if exist('t1', 'var') && ~isempty(t1)
        method = 'nbgauss';
    elseif exist('x1max', 'var') && ~isempty(x1max)
        method = 'landweber';
    else
        method = 'totalvar';
    end
end
if ~exist('sigma', 'var')
    sigma = [];
end
if ~exist('strategy', 'var')
    strategy = '';
end
if ~exist('t1', 'var')
    t1 = [];
end
if ~any(strcmpi(method, {'nbgauss', 'nbcos'})) && ~isempty(t1)
    warning('Cannot use parameter ''t1'' with method ''%s''.', method)
end
if ~exist('x1max', 'var')
    x1max = [];
elseif ~strcmpi(method, 'landweber') && ~isempty(x1max)
    warning('Parameter ''x1max'' can only be used with method ''landweber''.')
end

% Differentiate
switch lower(method)
        
    case 'savitzky'
        x1 = rdiff_savitzky(x, t, sigma, strategy);
        
    case 'central'
        x1 = rdiff_central(x, t, sigma, strategy);
        
    case 'sinstep'
        x1 = rdiff_sinstep(x, t, sigma, strategy);
        
    case 'expstep'
        x1 = rdiff_expstep(x, t, sigma, strategy);
        
    case 'nbgauss'
        x1 = rdiff_nbgauss(x, t, sigma, strategy, t1);
        
    case 'nbcos'
        x1 = rdiff_nbcos(x, t, sigma, strategy, t1);
        
    case 'tikhonov'
        x1 = rdiff_tikhonov(x, t, sigma, strategy);
        
    case 'totalvar'
        x1 = rdiff_totalvar(x, t, sigma, strategy);
        
    case 'tsvd'
        x1 = rdiff_tsvd(x, t, sigma, strategy);
        
    case 'landweber'
        x1 = rdiff_landweber(x, t, sigma, strategy, x1max);
        
    case 'kalman'
        x1 = rdiff_kalman(x, t, sigma, strategy);
        
    otherwise
        error('Unknown method: %s', method)
end
%
function x1 = rdiff_tikhonov(x, t, sigma, strategy)
% RDIFF_TIKHONOV   Regularised numerical differentiation using method based
%                  on constraining the 2-norm of the third derivative.
%    X1 = RDIFF_TIKHONOV(X, T) returns an estimate of the first derivative
%    of X; if T is a scalar, it is used as spacing between abscissa; if T
%    is a vector, its elements are used as the abscissa. If X is a matrix,
%    RDIFF_TIKHONOV differentiates each of its columns independently.
%    
%    X1 = RDIFF_TIKHONOV(X, T, SIGMA), where SIGMA is a scalar, uses SIGMA
%    as an estimate of the variance of the errors corrupting the elements
%    of X under the assumption that those errors are independent and
%    identically distributed. If SIGMA is an NxN matrix, it is used as the
%    covariance matrix of those errors.
%    
%    X1 = RDIFF_TIKHONOV(X, T, SIGMA, STRATEGY) or
%    X1 = RDIFF_TIKHONOV(X, T, [], STRATEGY) uses the strategy for the
%    optimisation of the constraint on the 2-norm of the third derivative
%    selected among the following ones:
%       'dp'      based on discrepancy principle (requires SIGMA)
%       'gcv'     based on generalised cross validation
%       'ncp'     based on normalised cumulative periodogram
%       'sure'    based on Stein's unbiased risk estimator (requires SIGMA)
%    By default, the 'dp' strategy is used if SIGMA is provided, and the
%    'ncp' strategy otherwise.
% 
%    Jakub Wagner, February 2, 2020
%    Institute of Radioelectronics and Multimedia Technology
%    Warsaw University of Technology

% Check size of x
if isrow(x)
    x = x';
    xrow = true;
end
if isscalar(x)
    error('x must be a vector or a 2-D matrix.')
end
[N, R] = size(x);

% Check size of t
if isscalar(t), t = t*(0:N-1)'; end
if ~isvector(t) || length(t) ~= N
    error('for size(x) = [%d, %d], t should be a scalar or a vector of length %d; size(t) = [%d, %d]', size(x,1), size(x,2), N, size(t,1), size(t,2))
end
if isrow(t); t = t'; end

% Select default strategy
if ~exist('strategy', 'var') || isempty(strategy)
    if exist('sigma', 'var') && ~isempty(sigma)
        strategy = 'dp';
    else
        strategy = 'ncp';
    end
end

% Check availability and size of sigma
if any(strcmpi(strategy, {'dp', 'sure'}))
    if ~exist('sigma', 'var') || isempty(sigma)
        warning('Cannot use ''%s'' strategy when sigma is not specified; using ''ncp'' strategy instead.', strategy)
        strategy = 'ncp';
    else
        if numel(sigma) == 1
            sigma = sigma * eye(N);
        elseif any(size(sigma) ~= N)
            error('Wrong size of sigma: [%d, %d] (should be scalar or %dx%d matrix)', size(sigma, 1), size(sigma, 2), N, N)
        end
    end
end

% Get quadrature matrix
Q = qmatrix(t);

% Modify matrix Q to add equation x'(1) = x'(2)
Q(1, 1:2) = [1, -1];

% Get progressive-difference matrix
d = 3;  % degree of derivative whose norm is to be constrained
D = pdmatrix(t, d);
% Other d might be worth giving a try.

% Coefficient for scaling regularisation parameter
di = N^(d+2) / trace(D'*D);

% Considered values of the regularisation parameter (to be scaled later on)
ls = logspace(-10, 10, 21)';
Nl = numel(ls);

% Constants used during optimisation of regularisation parameters
switch lower(strategy)
    case 'dp'
        uDP = 1.3;
        trS = trace(sigma);
    case 'gcv'
        
    case 'ncp'
        L = floor(N/2);
        l = repmat((1:L)' / L, [1 R]);
    case 'sure'
        sigma12 = sqrt(sigma);
    otherwise
        error('Unknown strategy: %s', strategy)
end

% Criterion values for optimisation of regularisation parameter
crit = inf(Nl, R);

% Estimates of the derivative
X1 = zeros(N, R, Nl);

% Make sure that first data point = 0.
xp = x - repmat(x(1,:), [N 1]);

% Optimise regularisation parameter
for nl = 1:Nl   % for each considered value of regularisation parameter
    
    % Scale regularisation parameter
    a = ls(nl) * di;
    
    % Generate "pre-hat" matrix
    H10 = (Q'*Q + a*(D'*D));
    
    if cond(H10) < 1e10
        
        % Differentiate and approximate data
        H1 = H10 \ Q';
        X1(:, :, nl) = H1 * xp;
        xa = Q * X1(:, :, nl);
        
        % Compute criteria for optimisation of regularisation parameter
        switch lower(strategy)
            case 'dp'
                crit(nl,:) = abs(sum((xa - xp).^2, 1) - uDP*trS);
            case 'gcv'
                crit(nl,:) = sum((xa - xp).^2, 1) / trace(eye(N) - Q*H1)^2;
            case 'ncp'
                residual = xa - xp;
                rf = fft(residual,[],1);
                rp = abs(rf(1:L,:)).^2;
                rncp = cumsum(rp,1) ./ repmat(sum(rp,1), [L 1]);
                crit(nl,:) = sqrt(sum((rncp - l).^2, 1));
            case 'sure'
                if R > 1
                    crit(nl,:) = sum(xa.^2, 1) - 2*diag(xp'*xa)' + 2*trace(sigma12*(Q*H1)*sigma12);
                else
                    crit(nl,:) = sum(xa.^2, 1) - 2*xp'*xa + 2*trace(sigma12*(Q*H1)*sigma12);
                end
        end
    end
end

% Select optimum values of regularisation parameter
x1 = zeros(N,R);
for r = 1:R
    nlSel = find(crit(:,r) == min(crit(:,r)), 1, 'last');
    x1(:,r) = X1(:, r, nlSel);
end

% Return row vector if x was one
if exist('xrow', 'var') && xrow
    x1 = x1';
end
%
%
function Q = qmatrix(t)
% QMATRIX   Matrix of coefficients of 2nd-order Newton-Cotes quadrature.
%    Q = QMATRIX(N), where N is a scalar, returns an NxN matrix whose
%    elements are the coefficients of the second-order Newton-Cotes
%    quadrature corresponding to equidistant abscissa.
% 
%    Q = QMATRIX(T), where T is an N-dimensional vector, returns an NxN
%    matrix whose elements are the coefficients of the second-order
%    Newton-Cotes quadrature corresponding to the abscissa specified by T.
%    The elements of T do not need to be equidistant.
%    
%    The integral of an N-dimensional vector X can be estimated as:
%    I = DT * QMATRIX(N) * X
%    where DT is the spacing between the abscissa or:
%    I = QMATRIX(T) * X
%    where it is assumed that X represents the values of the integrated
%    function evaluated at the abscissa T.
%    
%    Jakub Wagner, February 2, 2020
%    Institute of Radioelectronics and Multimedia Technology
%    Warsaw University of Technology

if isscalar(t), t = (1:t)'; end     % number of points specified instead
                                    % of vector of abscissa
N = length(t);
Q = zeros(N,N);
for n = 2:N
    Q(n:N, [n-1, n]) = Q(n:N, [n-1, n]) + (t(n)-t(n-1))/2;
end
%
%
function D = pdmatrix(t, k)
% PDMATRIX   Matrix of coefficients of progressive-difference formula.
%    D = PDMATRIX(N) generates an NxN matrix whose elements are the
%    coefficients of the progressive-difference formula.
%    The first derivative of an N-dimensional vector X can be estimated as:
%    X1 = PDMATRIX(N) * X / DT
%    where DT is the spacing of the data contained in X.
%    
%    D = PDMATRIX(T), where T is an N-dimensional vector, generates such a
%    matrix for differentiation of vectors which correspond to abscissa
%    specified by T, i.e.:
%    X1 = PDMATRIX(T) * X
%    The elements of T do not need to be equidistant.
%    
%    D = PDMATRIX(N, K) or D = PDMATRIX(T, K) generates a matrix
%    representative of the progressive-difference formula applied K times.
%    The K-th derivative of an N-dimensional vector X can be estimated as:
%    XK = PDMATRIX(N) * X / (DT^K)
%    or:
%    XK = PDMATRIX(T) * X
%    
%    Based on: J. J. Stickel, "Data smoothing and numerical differentiation
%    by a regularization method", Computers & Chemical Engineering,
%    vol. 34, no. 4, pp. 467-475, 2010.
%    
%    Jakub Wagner, February 2, 2020
%    Institute of Radioelectronics and Multimedia Technology
%    Warsaw University of Technology

if ~exist('k', 'var'), k = 1; end   % default order of differentiation = 1
if isscalar(t), t = (1:t)'; end     % number of points specified instead
                                    % of vector of abscissa
N = length(t);

if k == 0, D = eye(N); return; end % "order = 0" -> identity matrix

% Matrix of inverses of spacings between abscissa
Vd = zeros(N-k,N-k);
for n = (ceil(k/2)+1) : (N-floor(k/2))
    i = n - ceil(k/2);
    Vd(i,i) = 1/(t(n+floor(k/2)) - t(n-ceil(k/2)));
end

% Matrix of "unscaled" progressive-difference formula coefficients
Dh1 = zeros(N-k,N-k+1);
for n = 1:N-k
    Dh1(n,n) = -1;
    Dh1(n,n+1) = 1;
end

% Recursive implementation for order higher than 1
if k == 1
    D = Vd*Dh1;
else
    D = k*Vd*(Dh1*pdmatrix(t,k-1));
end
