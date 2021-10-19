% Density estimation for a beta(1,4) distribution

% See June 5, 2013 entry in lab book

clear all
close all

% The number of trials
N_trials = 2;

% The number of points between 0 and 1 where we will estimate the density
N_points = 256;

% The kth nearest neighbor to use
knn = 1;

% Exponent to decide how many subsets
a = 1/3;

% Vector to hold the density estimates
density_estimate_s_nearest = zeros(N_trials, N_points);
density_estimate_s_nearest_normalized = zeros(N_trials, N_points);

% Vector to hold the density estimates from Scott's rule
density_estimate_scott = zeros(N_trials, N_points);

% Vector to hold the density estimate from Botev-Grotowski-Kroese
density_estimate_kernel = zeros(N_trials, N_points);

% Vector to hold the errors for each trial
error_s_nearest = zeros(N_trials, N_points);
error_s_nearest_normalized = zeros(N_trials, N_points);
error_L1_s_nearest = zeros(N_trials,1);
error_L2_s_nearest = zeros(N_trials,1);

% Vector to hold the error from Scott's rule
error_scott = zeros(N_trials, N_points);

% Vector to hold the error from B-G-K KDE
error_kde = zeros(N_trials, N_points);
error_L1_kde = zeros(N_trials,1);

% x locations where the kde estimator returns densities
xmesh = zeros(N_trials, N_points);

% area vector
area = zeros(N_trials,1);

% Vector of points where we will estimate the density
x_vec = linspace(-3,3,N_points);

% The total number of samples
N_total_samples = 1000;
    
% Factor by which to reduce the number of representations
representation_factor = 1;

% The number of different representations
N_representations = round(representation_factor*N_total_samples^(1-a));
    
% The number of samples in each representation
N_samples = round(N_total_samples/N_representations);

N_total_samples = N_representations*N_samples;

% For each trial
for t = 1:N_trials
    
    % Print out which trial w e are at
    sprintf('Trial: %d', t)
                
    % First generate N_samples raw samples
    raw_y_1 = rand(N_total_samples, 1);
    raw_y_2 = -pi/2 + (pi*raw_y_1);
    y = tan(raw_y_2);  
    
    y_min = min(y);
    y_max = max(y);
    
    % Reshape y into N_R representations with N_samples in each
    % representation
    X_N_s_N_R = reshape(y, N_samples, N_representations);                
    
    % Loop over points at which we will estimate density
    for h = 1:length(x_vec)
    
    % The current point where the density is being estimated
    x_est = x_vec(h);        
    
    % Vector to hold the estimates of the density^-1
    density_inverse_estimate = zeros(N_representations, 1);
        
    % If outside the range of samples, the density estimate is just zero
    if(x_est < y_min || x_est > y_max)
        density_estimate_s_nearest(t,h) = 0;
    % Otherwise do the MLD-DE density estimation   
    else         
        % Loop over every representation
        for i = 1:N_representations
            
            % The current representation
            y_current = zeros(N_samples, 1);
            y_current = X_N_s_N_R(:,i);
                        
            % Get the distances of the current samples from the point we are
            % estimating the density at
            dist = y_current - x_est;
                    	                                                                              
            % Sort to get the nearest neighbor distance
            y_current_dist_sorted = sort(abs(dist));
                           
            % The first element of this sorted vector gives us our statistic
            density_inverse_estimate(i) = (N_samples+1)*y_current_dist_sorted(knn)*2;
                                                             
        end
        % End loop over representations
    
        % Now average over the representations to get an esimate of the
        % expectation
        averaged_density_inverse = mean(density_inverse_estimate);
    
        % The density estimate at trial t
        density_estimate_s_nearest(t,h) = 1/averaged_density_inverse;
        
    end
    % End if for x_est outside sample range
    
        % Error for trial t
        error_s_nearest(t,h) = abs(cauchypdf(x_est) - density_estimate_s_nearest(t,h));                
                                
    end
    % End loop over points where density is estimated
    
    % We have to normalize the s-nearest estimates
    area(t) = trapz(x_vec, density_estimate_s_nearest(t,:));
    density_estimate_s_nearest_normalized(t,:) = density_estimate_s_nearest(t,:)/area(t);  
    error_s_nearest_normalized(t,:) = abs(cauchypdf(x_vec) - density_estimate_s_nearest_normalized(t,:));
    error_L1_s_nearest(t) = trapz(x_vec, error_s_nearest(t,:));
    error_L2_s_nearest(t) = sqrt(trapz(x_vec, error_s_nearest(t,:).^2));
    
    % The densities from Scott's rule
    [fixed_width_density, bin_centres] = scotts_rule_density_estimation(y, x_vec);
    density_estimate_scott(t,:) = fixed_width_density;
    
    % Error for trial t
    error_scott(t,:) = abs(cauchypdf(x_vec) - density_estimate_scott(t,:));
    error_L1_scott(t) = trapz(x_vec, error_scott(t,:));
    error_L2_scott(t) = sqrt(trapz(x_vec, error_scott(t,:).^2));
                
    % B-G-K Kernel Density Estimation    
    [bandwidth density_estimate_kernel_temp xmesh cdf] = kde(y, N_points, [-3 3]);
    density_estimate_kernel(t,:) = interp1(xmesh,density_estimate_kernel_temp,x_vec,'nearest','extrap');
    clear density_estimate_kernel_temp;
    clear xmesh;
    
    % Error for trial t
    error_kde(t,:) = abs(density_estimate_kernel(t,:) - cauchypdf(x_vec));
    error_L1_kde(t) = trapz(x_vec, error_kde(t,:));
    error_L2_kde(t) = sqrt(trapz(x_vec, error_kde(t,:).^2));
                    
end
% End loop over trials
    
%% Errors and Plots

N_representations
N_total_samples

a

mean(error_L2_s_nearest)
mean(error_L2_scott)
mean(error_L2_kde)

std(error_L2_s_nearest)/sqrt(N_total_samples)
std(error_L2_scott)/sqrt(N_total_samples)
std(error_L2_kde)/sqrt(N_total_samples)

%
h1 = figure(1)
plot(x_vec, mean(error_s_nearest),'r-', 'LineWidth', 2)
hold on
plot(x_vec, mean(error_scott),'g-', 'LineWidth', 2)
hold on

xlabel_legend = xlabel('x');
ylabel_legend = ylabel('Mean Absolute Error');
fig_legend = legend('MLD-DE', 'IQR');
set(xlabel_legend,'FontSize',18);
set(ylabel_legend,'FontSize', 18);
set(fig_legend,'FontSize',18)
orient(figure(1),'portrait')
%

%
h2 = figure(2);
set(h2,'DefaultLineMarkerSize',12)
set(h2,'DefaultLineLineWidth',1)
set(h2,'DefaultAxesFontSize',14)
plot(x_vec, cauchypdf(x_vec),'k-','LineWidth', 2)
hold on
plot(x_vec, mean(density_estimate_s_nearest),'r-','LineWidth', 2)
hold on
plot(x_vec, mean(density_estimate_s_nearest) + std(density_estimate_s_nearest),'b--', 'LineWidth', 2)
hold on
plot(x_vec, mean(density_estimate_s_nearest) - std(density_estimate_s_nearest),'b--', 'LineWidth', 2)

xlabel_legend = xlabel('x');
ylabel_legend = ylabel('Density');
fig_legend = legend('True PDF','mean(MLD-DE)', 'mean \pm sd');
set(xlabel_legend,'FontSize',18);
set(ylabel_legend,'FontSize', 18);
set(fig_legend,'FontSize',18)
orient(figure(2),'portrait')
%

%%
figure(3)
plot(x_vec, cauchypdf(x_vec),'k-','LineWidth', 2)
hold on
plot(x_vec, mean(density_estimate_kernel),'r-','LineWidth', 2)
hold on
plot(x_vec, mean(density_estimate_kernel) + std(density_estimate_kernel),'b--', 'LineWidth', 2)
hold on
plot(x_vec, mean(density_estimate_kernel) - std(density_estimate_kernel),'b--', 'LineWidth', 2)

xlabel_legend = xlabel('x');
ylabel_legend = ylabel('Density');
fig_legend = legend('True PDF','mean(KDE)', 'mean \pm sd');
set(xlabel_legend,'FontSize',18);
set(ylabel_legend,'FontSize', 18);
set(fig_legend,'FontSize',18)
orient(figure(2),'portrait')
%%

%
figure(4)
plot(x_vec, cauchypdf(x_vec),'k-','LineWidth', 2)
hold on
plot(x_vec, mean(density_estimate_scott),'r-','LineWidth', 2)
hold on
plot(x_vec, mean(density_estimate_scott) + std(density_estimate_scott),'b--', 'LineWidth', 2)
hold on
plot(x_vec, mean(density_estimate_scott) - std(density_estimate_scott),'b--', 'LineWidth', 2)

xlabel_legend = xlabel('x');
ylabel_legend = ylabel('Density');
fig_legend = legend('True PDF','mean(IQR)', 'mean \pm sd');
set(xlabel_legend,'FontSize',18);
set(ylabel_legend,'FontSize', 18);
set(fig_legend,'FontSize',18)
orient(figure(2),'portrait')
%

%
figure(5)
plot(x_vec,cauchypdf(x_vec),'k--','LineWidth', 2)
hold on
plot(x_vec, density_estimate_s_nearest(t,:),'r-','LineWidth', 2)
hold on
plot(x_vec, density_estimate_scott(t,:) ,'b-', 'LineWidth', 2)
hold on
%plot(x_vec, density_estimate_kernel(t,:),'b-', 'LineWidth', 2)

xlabel_legend = xlabel('x');
ylabel_legend = ylabel('Density');
fig_legend = legend('True PDF','MLD-DE','Histogram (IQR)');
set(xlabel_legend,'FontSize',18);
set(ylabel_legend,'FontSize', 18);
set(fig_legend,'FontSize',18)
orient(figure(5),'portrait')
%
