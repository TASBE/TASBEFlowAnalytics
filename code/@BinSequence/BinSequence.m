% BinSequence holds geometric or arithmetic sequences of bins
% BinSequence(min_edge, step, max_edge, mode)
%    Mode = 'geometric': centers at geometric mean of edges
%    Mode = 'log_bins': geometric, except input values are assumed to be log10 scale
%    Mode = 'arithmetic': centers at arithmetic mean of edges
%
% mode        % 'geometric' or 'arithmetic'
% n_bins      % number of bins
% bin_edges   % array n_bins+1
% bin_centers % array
% bin_widths  % width per bin
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function B = BinSequence(min_edge, step, max_edge, mode)
    if nargin == 0
        B.mode = '';
        B.bin_edges = 0:0.25:1;
        B.bin_centers = mean([B.bin_edges(1:end-1);B.bin_edges(2:end)]);
        B.bin_widths = 0.25;
    elseif nargin == 4
        B.mode = mode;
        
        switch mode
            case 'log_bins'
                B.mode = 'geometric';
                B.bin_edges = 10.^(min_edge:step:max_edge);
                B.bin_centers = geomean([B.bin_edges(1:end-1);B.bin_edges(2:end)]);
                B.bin_widths = 10^step;
            case 'geometric'
                B.mode = 'geometric';
                B.bin_edges = 10.^(log10(min_edge):log10(step):log10(max_edge));
                B.bin_centers = geomean([B.bin_edges(1:end-1);B.bin_edges(2:end)]);
                B.bin_widths = step;
            case 'arithmetic'
                B.mode = 'arithmetic';
                B.bin_edges = min_edge:step:max_edge;
                B.bin_centers = mean([B.bin_edges(1:end-1);B.bin_edges(2:end)]);
                B.bin_widths = step;
            otherwise
                TASBESession.error('TASBE:BinSequence','UnknownBinMode','Unknown bin mode %s',mode);
        end
    end
            
    B.n_bins = numel(B.bin_centers);
    
    B=class(B,'BinSequence');
