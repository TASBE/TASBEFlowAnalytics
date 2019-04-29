% COMPUTE_PLASMID_ESTIMATES computes the active and inactive estimates for
% a PlasmidExpressionModel object.
%%%%%%
% WARNING: THIS FUNCTIONALITY IS EXPERIMENTAL AND SHOULD NOT BE TRUSTED
%%%%%%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [estimates PEM] = compute_plasmid_estimates(PEM)
    ERF_off = round(log10(PEM.ERF_per_plasmid)/log10(get_bin_widths(PEM.bins)));
%    pdist = plasmid_distribution_model(PEM); % get underlying plasmid model

    logbin = log10(get_bin_widths(PEM.bins));
    offsets = [0 logbin/2:logbin:6*log10(PEM.noise)];
    offset_gaussian = zeros(numel(offsets)-1,1);
        
    for j=1:(numel(offsets)-1) % ... adding probability mass to distribution
        offset_gaussian(j) = p_gaussian(offsets(j),offsets(j+1),0,log10(PEM.noise));
    end

    bin_edges = get_bin_edges(PEM.bins);
    num_bins = get_n_bins(PEM.bins);
    pdist = p_gaussian(log10(bin_edges(1:num_bins)),log10(bin_edges(2:(num_bins+1))), ...
        log10(PEM.plasmid_mean),log10(PEM.plasmid_std));


    % E(p|c) = sum_p p * P(c and p)/P(c) = 1/P(c) * E(p and c)
    E_p_and_c = zeros(get_n_bins(PEM.bins),1); % E(p and c) = sum_p p P(c and p)
    p_cfp = zeros(get_n_bins(PEM.bins),1);     % P(c) = sum_p P(c and p)
    
    bin_centers = get_bin_centers(PEM.bins);
    for i=1:get_n_bins(PEM.bins) % walk through the plasmid model...
        for j=1:(numel(offsets)-1) % ... adding probability mass to distribution
            p_p_and_c = pdist(i)*offset_gaussian(j);
            E_frag = p_p_and_c*log10(bin_centers(i));

            idx = (i+ERF_off)+(j-1); % add with upward offset
            if(idx >= 1 && idx <= numel(E_p_and_c)), p_cfp(idx) = p_cfp(idx) + p_p_and_c; E_p_and_c(idx) = E_p_and_c(idx) + E_frag; end;
            idx = (i+ERF_off)-(j-1); % add with downward offset
            if(idx >= 1 && idx <= numel(E_p_and_c)), p_cfp(idx) = p_cfp(idx) + p_p_and_c; E_p_and_c(idx) = E_p_and_c(idx) + E_frag; end;
        end
    end
    active_estimates = 10.^(E_p_and_c./p_cfp);
  
    % finally, blend with inactive data points
    p_active = p_gaussian(log10(bin_edges(1:num_bins)),log10(bin_edges(2:(num_bins+1))), ...
        PEM.fp_dist.mu(PEM.active_component),sqrt(PEM.fp_dist.Sigma(:,:,PEM.active_component)));
    w_a = PEM.fp_dist.weight(PEM.active_component);
    p_inactive = 0;
    for i=1:numel(PEM.inactive_component),
        p_inactive = p_inactive + p_gaussian(log10(bin_edges(1:num_bins)),log10(bin_edges(2:(num_bins+1))), ...
            PEM.fp_dist.mu(PEM.inactive_component(i)),sqrt(PEM.fp_dist.Sigma(:,:,PEM.inactive_component(i))));
    end
    w_i = sum(PEM.fp_dist.weight(PEM.inactive_component));
    PEM.fraction_active = (w_a*p_active./(w_a*p_active+w_i*p_inactive))';
    if numel(PEM.inactive_component)
        PEM.fraction_active(bin_centers<10.^(min(PEM.fp_dist.mu(PEM.inactive_component)))) = 0; % Assert that there's nothing usable below the inactive component
    end    
    estimates = active_estimates.*PEM.fraction_active;
