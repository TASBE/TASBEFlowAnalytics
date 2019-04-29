% Copyright (C) 2015 Lachlan Andrew <lachlanbis@gmail.com>
% Converted to old-style classes by Jacob Beal, 2017
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.
%
% -*- texinfo -*-
% @deftypefn {Function File} {@var{GMdist} =} gmdistribution (@var{mu}, @var{Sigma})
% @deftypefnx {Function File} {@var{GMdist} =} gmdistribution (@var{mu}, @var{Sigma}, @var{p})
% @deftypefnx {Function File} {@var{GMdist} =} gmdistribution (@var{mu}, @var{Sigma}, @var{p}, @var{extra})
% Create an object of the  gmdistribution  class which represents a Gaussian
% mixture model with k components of n-dimensional Gaussians.
%
% Input @var{mu} is a k-by-n matrix specifying the n-dimensional mean of each
% of the k components of the distribution.
%
% Input @var{Sigma} is an array that specifies the variances of the
% distributions, in one of four forms depending on its dimension.
% @itemize
%   @item n-by-n-by-k: Slice @var{Sigma}(:,:,i) is the variance of the
%         i'th component
%   @item 1-by-n-by-k: Slice diag(@var{Sigma}(1,:,i)) is the variance of the
%         i'th component
%   @item n-by-n: @var{Sigma} is the variance of every component
%   @item 1-by-n-by-k: Slice diag(@var{Sigma}) is the variance of every
%         component
% @end itemize
%
% If @var{p} is specified, it is a vector of length k specifying the proportion
% of each component.  If it is omitted or empty, each component has an equal
% proportion.
%
% Input @var{extra} is used by fitgmdist to indicate the parameters of the
% fitting process.
% @seealso{fitgmdist}
% @end deftypefn
function obj = gmdistribution(mu,sigma,p,extra)

obj.mu = [];                        % means
obj.Sigma = [];                     % covariances
obj.PComponents = [];               % mixing proportions
obj.DistributionName = 'gaussian mixture distribution';          % 'gaussian mixture distribution'
obj.NumComponents = [];             % Number of mixture components
obj.NumVariables = [];              % Dimension d of each Gaussian component

obj.CovarianceType = [];            % 'diagonal' if DiagonalCovariance, 'full' othw
obj.SharedCovariance = [];          % true if all components have equal covariance

% Set by a call to fitgmdist
obj.AIC = [];                       % Akaike Information Criterion
obj.BIC = [];                       % Bayes Information Criterion
obj.Converged = [];                 % true  if algorithm converged by MaxIter
obj.NegativeLogLikelihood = [];     % Negative of log-likelihood
obj.NlogL = [];                     % Negative of log-likelihood
obj.NumIterations = [];             % Number of iterations
obj.RegularizationValue = [];       % const added to diag of cov to make +ve def

% Private property, which is suggested not to be accessed
obj.DiagonalCovariance = [];        % bool summary of 'CovarianceType'

if nargin>0 
    % Construction
    if nargin<3, p = []; end;
    if nargin<4, extra = []; end;
    obj.DistributionName = 'gaussian mixture distribution';
    obj.mu = mu;
    obj.Sigma = sigma;
    [row_n, col_n] = size(mu); 
    obj.NumComponents = row_n; %rows (mu);
    obj.NumVariables = columns (mu);
    if (isempty (p))
      obj.PComponents = ones (1,obj.NumComponents) / obj.NumComponents;
    else
      if any (p < 0)
        error ('gmmdistribution: component weights must be non-negative');
      end
      s = sum(p);
      if (s == 0)
        error ('gmmdistribution: component weights must not be all zero');
      elseif (s ~= 1)
        p = p / s;
      end
      obj.PComponents = p(:)';
    end
    if (length (size (sigma)) == 3)
      obj.SharedCovariance = false;
    else
      obj.SharedCovariance = true;
    end
    [r_sig, c_sig] = size(sigma); 
    [r_mu, c_mu] = size(mu);
    if (r_sig == 1 && c_mu > 1) %(rows (sigma) == 1 && columns (mu) > 1)
      obj.DiagonalCovariance = true;
      obj.CovarianceType = 'diagonal';
    else
      obj.DiagonalCovariance = false;       %% full
      obj.CovarianceType = 'full';
    end

    if (~isempty (extra))
      obj.AIC                   = extra.AIC;
      obj.BIC                   = extra.BIC;
      obj.Converged             = extra.Converged;
      obj.NegativeLogLikelihood = extra.NegativeLogLikelihood;
      obj.NlogL                 = extra.NegativeLogLikelihood;
      obj.NumIterations         = extra.NumIterations;
      obj.RegularizationValue   = extra.RegularizationValue;
    end
end

obj = class(obj,'gmdistribution');
