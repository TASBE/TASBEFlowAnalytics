% PlasmidExpressionModel is the class that maps ERFs to constitutive
% fluorescence to plasmid count with the following properties:

        %plot_pem = false
        %CFP_counts              % source for model
        % piecewise-gaussian model of plasmid distribution
        %plasmid_mean            % max point of distribution
        %lower_gaussian_sigma    % gaussian std.dev. for CFP<mean 
        %upper_gaussian_sigma    % gaussian std.dev. for CFP>mean
        % plasmid -> ERF model
        %noise                   % std noise
        %ERF_per_plasmid
        % ERF -> plasmid function
        %bins                    % ERF bins
        %estimated_plasmids      % corresponding estimates per bin
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

 function PEM = PlasmidExpressionModel(cfp_data,CFP_af,CFP_noise,ERF_per_plasmid,drop_threshold,n_components)
        
        % init
        PEM.plot_pem = false;
        PEM.CFP_counts =[];             % source for model
        %%% OLD:
        %%% piecewise-gaussian model of plasmid distribution
        %%%PEM.plasmid_mean=[];            % max point of distribution
        %%%PEM.lower_gaussian_sigma=[];    % gaussian std.dev. for CFP<mean 
        %%%PEM.upper_gaussian_sigma=[];    % gaussian std.dev. for CFP>mean
        % model of plasmid distribution
        PEM.fp_dist = [];                % 2-component gaussian mixture model of log-scale fluorescence (active vs. inactive cells)
        PEM.active_component = [];       % which component is the active one?
        PEM.inactive_component = [];     % which component is the inactive one?
        PEM.plasmid_mean = [];           % mean of plasmids for active cells
        PEM.plasmid_std = [];            % std.dev of plasmids for active cells
        PEM.fraction_active = [];        % corresponding portion active cells per bin
        % plasmid -> ERF model
        PEM.noise=[];                   % std noise
        PEM.ERF_per_plasmid=[];
        % ERF -> plasmid function
        PEM.bins =BinSequence();                   % ERF bins
        PEM.estimated_plasmids=[];      % corresponding estimates per bin

        if nargin > 0
            
            % default to two gaussian components of the model
            if nargin < 6, n_components = 2; end;
 
            PEM.ERF_per_plasmid = ERF_per_plasmid;
            PEM.noise = CFP_noise;
            PEM.bins = BinSequence(0,0.1,12,'log_bins');
            
            %%%% GAUSSIAN MIXTURE MODEL METHOD
            which = cfp_data>drop_threshold;
            if(sum(which)<10)
                TASBESession.warn('Model:CFPDistributionSeparation','NotEnoughData','Not enough data to model');
                PEM=class(PEM,'PlasmidExpressionModel');
                return;
            end
            
            % Seeded centers: [8], [5 8], [5 6.5 8], [5 6 7 8], ...
            if(~isempty(CFP_af)),
                min_center = log10(getMeanERF(CFP_af));
            else
                min_center = 1;
            end
            if n_components == 1, 
                initial_centers = 1;
            else
                initial_centers = min_center + (0:(3/(n_components-1)):3);  % spread on a 1000x scale from AF
            end;
            [label PEM.fp_dist] = emgm(log10(cfp_data(which))',initial_centers); % compute on log scale

            if (numel(PEM.fp_dist.mu) == 1 || max(PEM.fp_dist.weight)>0.999) && n_components>1, % failure to separate -> duplicate
                TASBESession.warn('Model:CFPDistributionSeparation','CannotSeparateDistribution','Distribution could not be separated');
                PEM.fp_dist.mu(2) = PEM.fp_dist.mu(1);
                PEM.fp_dist.Sigma(:,:,2) = PEM.fp_dist.Sigma(:,:,1); % This Sigma is a *variance*, and not a std. dev.
                PEM.fp_dist.weight = [0.5 0.5];
            else if numel(PEM.fp_dist.mu) < n_components
                TASBESession.warn('Model:CFPDistributionSeparation','NotEnoughComponents','Less components than expected: %i vs. %i',numel(PEM.fp_dist.mu),n_components);
                % add zero-weight components
                for i = (numel(PEM.fp_dist.mu)+1):n_components
                    PEM.fp_dist.mu(i) = PEM.fp_dist.mu(1);
                    PEM.fp_dist.Sigma(:,:,i) = PEM.fp_dist.Sigma(:,:,1); % This Sigma is a *variance*, and not a std. dev.
                    PEM.fp_dist.weight(i) = 0;
                end
                end
            end
            
            % Transform from CFP to plasmids
            high_dist = find(PEM.fp_dist.mu == max(PEM.fp_dist.mu),1);
            PEM.plasmid_mean = 10.^(PEM.fp_dist.mu(high_dist))/ERF_per_plasmid;
            PEM.plasmid_std = 10.^sqrt(max(0,PEM.fp_dist.Sigma(high_dist)-log10(PEM.noise)^2)); % No square on fp_dist.Sigma since it's variance
            PEM.active_component = high_dist;
            PEM.inactive_component = find(PEM.fp_dist.mu ~= max(PEM.fp_dist.mu));

            
            PEM=class(PEM,'PlasmidExpressionModel');
 
            [est_plasmids PEM] = compute_plasmid_estimates(PEM);
            PEM.estimated_plasmids = est_plasmids;
            
        else
            PEM=class(PEM,'PlasmidExpressionModel');
        end
