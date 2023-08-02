function [k mu err] = estimator_RSK(ROI, display_plot)
% ESTIMATOR_RSK_V4: homodyned k distribution parameter estimation using
% level curves based on:
%                     R = signal-to-noise ratio (SNR),
%                     S = skewness, and
%                     K = kurtosis
% of arbitrary moments of the envelope echo samples.
%
% Usage: [k mu err] = estimator_RSK(ROI)
%        [k mu err] = estimator_RSK(ROI, display_plot)
% Inputs:  -- ROI = envelope samples in region of interest
%          -- display_plot = Boolean option: 1 = show level curves plot
%                                            0 = don't show plot (default)
% Outputs: -- k, mu = estimated model parameters: k = structure/periodicity
%                     parameter; mu = effective scatterer number density
%                     (effective number of scatterers per resolution cell)
%          -- err = sum of squared distances from estimated point
%                   to level curves (i.e., a measure of goodness of fit)
% Estimator will return an error if an estimate is not possible using the
% provided data.

% Code makes use of precomputed RSK values in file RSK_data_501.mat. If
% this file name be changed, this code must be modified in FOUR locations.

% Important Variables:
% RSK_curves(l,m,n,p) = precomputed RSK data matrix. Indexing convention:
%  l <--> k, m <--> mu, n <--> nu, p <--> 1,2,3 for R,S,K classifiers.
%  k, mu = model parameters; nu = moment order; p = which summary statistic
% k_values(l), mu_values(m), nu_values(n) = values where RSK_curves are
%  calculated (i.e., axis labels of 3D matrix volume)
% RSK_nu_values(n,p) = Estimated summary statistic values as a function of
%  moment order index n and summary statistic index p
% curve_pts(:,n,p) = cell array that holds a representation of the level
%  curve for each n and p.

% Some potentially useful references
%
% M. Martin-Fernandez, R. Cardenes, and C. Alberola-Lopez, "Parameter
% estimation of the homodyned K distribution based on signal to noise
% ratio," in Proceedings of the IEEE Ultrasonics Symposium, 2007. pp.
% 158-161.
%
% R. W. Prager, H. Gee, G. M. Treece, and L. H. Berman, "Analysis of
% speckle in ultrasound images using fractional order statistics and the
% homodyned k distribution," Ultrasonics, vol. 40, pp. 133-137, 2002.
%
% V. Dutt and J. F. Greenleaf, "Speckle analysis using signal to noise
% ratios based on fractional order moments," Ultrasonic Imaging, vol. 17
% pp. 251-268, 1995.
%
% V. Dutt and J. F. Greenleaf, "Ultrasound echo envelope analysis using a
% homodyned K distribution signal model," Ultrasonic Imaging, vol. 16, pp.
% 265-287, 1994.

% Estimator developed July 2008 - March 2009 at the Bioacoustics Research
% Laboratory (BRL) at the University of Illinois at Urbana-Champaign by
% David P. Hruska (email: david.hruska@alumni.illinois.edu)

% Based on optimization of the angle of intersection of level curves,
% Version 4 (March 2009) uses:
% -- moment orders nu = 0.72 and 0.88
% -- k on interval [0 5]
% -- mu on interval [0.001 100]
% Version 4 uses a grid search method to locate the global L2-norm minimum
% of the distance from level curves.

% --- Handle inputs ---

if (nargin < 2)   % by default (i.e., if not specified), don't display plot
    display_plot = 1;
end

% --- Load RSK data file ---
persistent RSK_curves;  %       use persistent variables to prevent
persistent k_values;    %#ok    repeatedly loading large RSK data
persistent mu_values;   %#ok    file when analyzing many ROIs in
persistent nu_values;   %#ok    the same image
if isempty(RSK_curves)
    load RSK_data_501.mat
end
grid_size = length(k_values);  % assume square grid
num_nu = length(nu_values);    % number of moment orders

% --- Estimate and check validity of RSK classifier functions ---
RSK_nu_values = zeros(num_nu,3);
for n = 1:num_nu
    nu = nu_values(n);
    ROI_nu = ROI.^nu;  % for convenience
    
    RSK_nu_values(n,1) = mean(ROI_nu)/std(ROI_nu);  % estimate of SNR
    % check if estimated SNR is within range
    wee_SNR = check_params(RSK_nu_values(n,1), nu_values(n), 1, ...
        min(min(RSK_curves(:,:,n,1))), max(max(RSK_curves(:,:,n,1))));

    RSK_nu_values(n,2) = skewness(ROI_nu, 0);  % unbiased estimate skewness
    % check if estimated skewness is within range
    wee_skewness = check_params(RSK_nu_values(n,2), nu_values(n), 2, ...
        min(min(RSK_curves(:,:,n,2))), max(max(RSK_curves(:,:,n,2))));
 
    RSK_nu_values(n,3) = kurtosis(ROI_nu, 0);  % unbiased estimate kurtosis
    % check if estimated kurtosis is within range
    %%%PL added for error problem
    wee_kurtosis = check_params(RSK_nu_values(n,3), nu_values(n), 3, ...
        min(min(RSK_curves(:,:,n,3))), max(max(RSK_curves(:,:,n,3))));
    if wee_kurtosis == -1111
        k = -1111;
        mu = -1111;
        err = -1111;
        return
    end
    %%%
    %check_params(RSK_nu_values(n,3), nu_values(n), 3, ...
    %    min(min(RSK_curves(:,:,n,3))), max(max(RSK_curves(:,:,n,3))));
end

% --- Find RSK level curves ---
curve_pts = cell(num_nu, 3);
for n = 1:num_nu
    for p = 1:3   % iterate over R, S, and K
        RSK_nu = RSK_nu_values(n,p);
        % need to search both horiz & vert directions to get complete curve
        h_ind = search_horiz(n, p, RSK_nu, grid_size);
        
        %%%PL added
        if isempty(h_ind)
            k = -3333;
            mu = -3333;
            err = -3333;
            return
        end
        %%%
        
        v_ind = search_vert(n, p, RSK_nu, grid_size);
        % keep only unique points
        curve_pts{n, p} = unique([h_ind; v_ind], 'rows');
    end
end

% --- Display level curve plot ---
if display_plot
    accumulative = zeros(grid_size);  % pixelwise sum of all level curves
    for n = 1:num_nu
        for p = 1:3
            % convert coordinates to linear array indices
            curve_indices = sub2ind([grid_size grid_size], ...
                curve_pts{n,p}(:,1), curve_pts{n,p}(:,2));
            accumulative(curve_indices) = accumulative(curve_indices) + 1;
        end
    end
    
    figure; imagesc(10*log10(mu_values), k_values, accumulative);
    xlabel('10log_1_0\mu'); ylabel('k'); title('Sum of all level curves');
    set(gca,'YDir','normal'); axis square; colorbar;
    
end

% --- Find point that minimizes total squared distance to level curves ---
% initial grid search parameters
num_grid_pts = 11;
l_min = 1;
l_max = grid_size;
m_min = 1;
m_max = grid_size;

go = 1;
while go
    % determine L2-norm distance metric at a small number of sample points
    l_vals = round(linspace(l_min, l_max, num_grid_pts));
    m_vals = round(linspace(m_min, m_max, num_grid_pts));
    [l_grid m_grid] = meshgrid(l_vals, m_vals);
    cell_grid = mat2cell([l_grid(:) m_grid(:)], ones(1, numel(l_grid)), 2);
    dist2 = cellfun(@(x) level_curve_dist(x,curve_pts,num_nu), cell_grid);
    
    % find coordinates of the point with least L2-norm
    [err index] = min(dist2);
    [m_ind l_ind] = ind2sub(size(l_grid), index);
    l = l_vals(l_ind);
    m = m_vals(m_ind);
    
    % stop search once we have sampled on a dense grid
    if ((range(l_vals) <= num_grid_pts) && (range(m_vals) <= num_grid_pts))
        go = 0;
    else
        % update search space for next iteration (make sure to keep the
        % search space stays a square when we bump up against the edges)
        new_box_size = round(range(l_vals)/4);   % half-size of new box
        
        if (l <= new_box_size)                   % hit bottom edge
            l_min = 1;
            l_max = 2*new_box_size + 1;
        elseif (l + new_box_size >= grid_size)   % hit top edge
            l_min = grid_size - 2*new_box_size;
            l_max = grid_size;
        else                  % somewhere in the middle of the search space
            l_min = l - new_box_size;
            l_max = l + new_box_size;
        end
        
        if (m <= new_box_size)                   % hit left edge
            m_min = 1;
            m_max = 2*new_box_size + 1;
        elseif (m + new_box_size >= grid_size)   % hit right edge
            m_min = grid_size - 2*new_box_size;
            m_max = grid_size;
        else                  % somewhere in the middle of the search space
            m_min = m - new_box_size;
            m_max = m + new_box_size;
        end
    end
end

% --- Return parameter estimates ---
k = k_values(l);
mu = mu_values(m);
err = sqrt(sum(level_curve_dist([l m], curve_pts, num_nu).^2));

% %%%
% SNR072=RSK_curves(:,:,1,1);
% SNR088=RSK_curves(:,:,2,1);
% kurtosis072=RSK_curves(:,:,1,3);
% kurtosis088=RSK_curves(:,:,2,3);
% skewness072=RSK_curves(:,:,1,2);
% skewness088=RSK_curves(:,:,2,2);
% 
% K=0:.01:5;
% LogMU=-3:.01:2;
% 
% figure
% hold on
% contour(LogMU,K,SNR072(:,:),'r');
% contour(LogMU,K,kurtosis072(:,:),'k');
% contour(LogMU,K,skewness072(:,:),'g');
% contour(LogMU,K,SNR088,'r-.');
% contour(LogMU,K,kurtosis088,'k-.');
% contour(LogMU,K,skewness088,'g-.');
% plot(log10(mu),k,'m+')
% xlabel('log10(mu)', 'FontSize', 15, 'FontWeight', 'bold')
% ylabel('k', 'FontSize', 15, 'FontWeight', 'bold'); 
% set(gca,'xlim',[-3 2],'ylim',[0 5],'FontSize', 15)
% set(gca, 'XTick', [-3 -2 -1 0 1 2], 'YTick', [01 2 3 4 5],'FontSize', 15)
% axis square; colorbar;
% hold off
% %%%

end

function [temp_err] = check_params(param, nu, p, param_min, param_max) %changed for error problem should be: [] = check_params(param, nu, p, param_min, param_max)
% Check validity of estimated classifier function values versus the minimum
% and maximum pre-computed values.
% Inputs: -- param = classifier function (R, S, or K)
%         -- nu = moment order
%         -- p = identifies classifier: R (p=1), S (p=2), or K (p=3)
%         -- param_min = minimum allowable parameter value
%         -- param_max = maximum allowable parameter value
% No Output (estimator_RSK function will terminate on an error)

% determine which classifier function (cf) we are dealing with
switch p
    case 1; cf = 'R';
    case 2; cf = 'S';
    case 3; cf = 'K';
end
%%%PL added for error problem
temp_err = 0;
%%%
% check for and display errors
if isnan(param)
    error(['estimator_RSK:' cf 'isNaN'], [cf ' parameter estimate ' ...
        'is NaN for moment order ' num2str(nu)]);
elseif (param < param_min)
    %%%PL added for error problem
    temp_err = -1111;
    %%%
    %     error(['estimator_RSK:' cf 'outOfRange'], [cf ' parameter ' ...
    %         'estimate is out of range for moment order ' num2str(nu) ...
    %         '\nEstimate is ' num2str(param) '; minimum allowable is ' ...
    %         num2str(param_min)]);
elseif (param > param_max)
    %%%PL added for error problem
    temp_err = -1111;
    %%%
    %    error(['estimator_RSK:' cf 'outOfRange'], [cf ' parameter ' ...
    %        'estimate is out of range for moment order ' num2str(nu) ...
    %        '\nEstimate is ' num2str(param) '; maximum allowable is ' ...
    %        num2str(param_max)]);
end

end

function [indices] = search_vert(n, p, RSK_nu, grid_size)
% Find the particular R, S, or K level curve by searching vertically (to
% find an l value) for each m value.
% Inputs: -- n = moment order index
%         -- p = specifies R (p=1), S (p=2), or K (p=3)
%         -- RSK_nu = R, S, or K value calculated from data from ROI
%         -- grid_size = size of grid of precomputed RSK values
% Output: -- indices = 2D array of matrix indices [l,m] that define the
%                      R, S, or K level curve

persistent RSK_curves;
if isempty(RSK_curves)
    load RSK_data_501.mat
end

indices = [];

for m = 1:grid_size
    sign_changes = find(diff(sign(RSK_curves(:,m,n,p) - RSK_nu)));
    
    % test each side of the sign change to find which is better
    for i = 1:length(sign_changes)
        l = closer_RSK_pt([sign_changes(i) m n p], ...
            [sign_changes(i)+1 m n p], RSK_nu, 1);
        indices = [indices; [l m]]; %#ok
    end
end

% earlier error checking should have identified any problems but just in
% case, perform another simple error check

if isempty(indices)
    error('estimator_RSK:GenOutOfRange', 'General parameter out of range');
end

end

function [indices] = search_horiz(n, p, RSK_nu, grid_size)
% Find the particular R, S, or K level curve by searching horizontally (to
% find an m value) for each l value.
% Inputs: -- n = moment order index
%         -- p = specifies R (p=1), S (p=2), or K (p=3)
%         -- RSK_nu = R, S, or K value calculated from data from ROI
%         -- grid_size = size of grid of precomputed RS values
% Output: -- indices = 2D array of matrix indices [l,m] that define the
%                      R, S, or K level curve

persistent RSK_curves;
if isempty(RSK_curves)
    load RSK_data_501.mat
end

indices = [];

for l = 1:grid_size
    sign_changes = find(diff(sign(RSK_curves(l,:,n,p) - RSK_nu)));
    
    % test each side of the sign change to find which is better
    for i = 1:length(sign_changes)
        m = closer_RSK_pt([l sign_changes(i) n p], ...
            [l sign_changes(i)+1 n p], RSK_nu, 2);
        indices = [indices; [l m]]; %#ok
    end
end

% earlier error checking should have identified any problems but just in
% case, perform another simple error check

%%% PL commented out
%if isempty(indices)
%    error('estimator_RSK:GenOutOfRange', 'General parameter out of range');
%end
%%%
end

function [closer] = closer_RSK_pt(pt1, pt2, RSK_nu, index)
% Given 2 points in the form [l1,m1,n1,p1], [l2,m2,n2,p2] and an R, S, or
% K value RSK_nu, determine which point on RSK_curves is closer to RSK_nu
% Inputs: -- pt1, pt2 = points of RSK_curves
%         -- RSK_nu = R, S, or K value calculated from data from ROI
%         -- index = index of the point to return, for example index = 1
%                    would return either l1 or l2.
% Output: -- closer = one of l1, l2, m1, m2, n1, n2, p1, p2 corresponding
%                     to point with R, S, or K closer to RSK_nu

persistent RSK_curves;  %#ok
if isempty(RSK_curves)
    load RSK_data_501.mat
end

% find closer point in terms of absolute error
if (abs(RSK_curves(pt2(1), pt2(2), pt2(3), pt2(4)) - RSK_nu) > ...
        abs(RSK_curves(pt1(1), pt1(2), pt1(3), pt1(4)) - RSK_nu))
    closer = pt1(index);
else
    closer = pt2(index);
end

end

function [dist2] = level_curve_dist(lm, curve_pts, num_nu)
% Total squared distance between the point defined by indices lm = [l m]
% and all level curves
% Inputs: -- lm = [l m] = indices of point in question
%         -- curve_pts = points that make up all level curves to consider
%         -- num_nu = number of moment orders
% Output: -- dist2 = total squared distance from all level curves

err_vect = zeros(num_nu, 3);

for n = 1:num_nu
    for p = 1:3
        % find minimum distance to curve by checking all points on curve
        err_vect(n,p) = min(sum((repmat(lm, length(curve_pts{n,p}),1) - ...
            curve_pts{n,p}).^2, 2));
    end
end
dist2 = sum(err_vect(:));

end






