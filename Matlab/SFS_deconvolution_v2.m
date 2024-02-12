clear
%-------------------------------------------------Number of cells at T_N
%----------------------(the exact SFS will be normalized to this number)
global                  N_end
N_end                   = 1000;
%------------------------------------Minimum and maximum number of reads
global                  r_min r_max
r_min                   = 0;
r_max                   = 500;
%------------------------------Minimum variant read count to be accepted
global                  min_variant_read
min_variant_read        = 5;
%--------------------------------Minimum total read count to be accepted
global                  min_total_read
min_total_read          = 0;
%---------------------Number of steps to divide SFS frequencies in [0,1]
global                  SFS_totalsteps SFS_totalsteps_base
SFS_totalsteps          = 50;
SFS_totalsteps_base     = 100;
vec_freq                = (1:SFS_totalsteps)/SFS_totalsteps;
%----------------------Choice of ploidy, which changes the binomial rate
%---------------------------------------------Z|R ~ Binomial(R,k/(pi*n))
%                         1 means haploid
%                         2 means diploid
option_ploidy           = 2;
%------------------------------------------------Maximum number of humps
global                  max_hump_count
max_hump_count          = 2;
%------------Choice of sampling coverage distribution and its parameters
% global                  option_dist_coverage
% option_dist_coverage    = 'TCGA';



global                  option_dist_coverage dist_coverage_var_1
option_dist_coverage    = 'binomial';
dist_coverage_var_1     = 100;
%----------------------------------------------------Options for fitting
N_SFS_positions         = 100;
N_fitting_rounds        = 200;
threshold_stop          = 0.95;
threshold_error         = 0.2;
%--------------------------------------------Categories for plotting SFS
SFS_categories          = {'Foreground_1','Foreground_0','Background_1','Truncal'};
SFS_categories_colormap = [0.5 0.5 0.5; 0.2 0.2 0.2; 0.6350 0.0780 0.1840; 0.3010 0.7450 0.9330];




%---------------------------------------------------Input binomial table
global          matrix_binomial_PDF
filename        = ['/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/MK-Cod.Analysis of the SFS/Binomial_tables/Binomial_PDF_' ...
    num2str(N_end) '_' ...
    num2str(r_max) '_' ...
    num2str(min_variant_read) '_' ...
    num2str(min_total_read) '_'...
    num2str(SFS_totalsteps_base) '_' ...
    num2str(option_ploidy) '.mat'];
inputBinomialMatrix = matfile(filename);
matrix_binomial_PDF = inputBinomialMatrix.matrix_binomial_PDF;
%--------------------------------Fit the SFS for each of the simulations
n_simulations           = 100;
table_parameters        = zeros(n_simulations,2*max_hump_count+1);
for n_simulation = 1:n_simulations
    %------------------------------------------Input mutational SFS data
    filename        = strcat('/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/vignettes/TEST_SFS_DECONVOLUTION/SFS_',num2str(n_simulation),'.txt');
    T               = table2cell(readtable(filename,'Delimiter',' ','ReadVariableNames',false));
    vec_refcount    = str2double(string(T(:,1)));
    vec_altcount    = str2double(string(T(:,2)));
    vec_marker      = string(T(:,3));
    vec_totcount    = vec_refcount+vec_altcount;
    vec_marker      = vec_marker(vec_altcount>=5);
    vec_refcount    = vec_refcount(vec_altcount>=5);
    vec_altcount    = vec_altcount(vec_altcount>=5);
    tmp_VAF         = vec_altcount./(vec_altcount+vec_refcount);
    %--------------------Check that all SFS categories are accounted for
    notInCategories = setdiff(unique(vec_marker), SFS_categories);
    if ~isempty(notInCategories)
        disp('The following strings in vec_marker are not in SFS_categories:')
        disp(notInCategories)
    end
    %--------------------------------------Prepare coverage distribution
    prep_distribution_patient(vec_totcount);
    %--------------------------------------Build reference SFS libraries
    vec_SFS_positions               = [1:N_SFS_positions]/N_SFS_positions;
    library_SFS_component           = cell(2,N_SFS_positions);
    %   Build SFS library for neutral component
    vec_para                        = [1];
    vec_SFS                         = SFS_expected(vec_para);
    library_SFS_component{1,1}      = vec_SFS;
    %   Build SFS library for binomial humps
    for i=1:N_SFS_positions
        p                           = vec_SFS_positions(i);
        vec_para                    = [0 p 1];
        vec_SFS                     = SFS_expected(vec_para);
        library_SFS_component{2,i}  = vec_SFS;
    end
    %------------------------------------------Prepare the SFS from data
    %   Initialize the SFS
    vec_SFS_real             = zeros(1,SFS_totalsteps);
    mutation_count      = 0;
    %   Build the SFS
    for j=1:length(vec_refcount)
        no_variant  = vec_altcount(j);
        no_total    = vec_refcount(j)+vec_altcount(j);
        if (no_variant>=min_variant_read)&&(no_total>=min_total_read)
            mutation_count  = mutation_count+1;
            %           Find and record the VAF
            VAF             = no_variant/no_total;
            pos         = find(vec_freq>=VAF,1);
            vec_SFS_real(pos)= vec_SFS_real(pos)+1;
        end
    end
    %--------------------------------------------------------Fit the SFS
    N_humps                             = -1;
    ratio_error                         = 0;
    err_best_final                      = Inf;
    while ratio_error<threshold_stop && err_best_final>threshold_error && N_humps<max_hump_count
        N_humps                         = N_humps+1;
        err_best_current                = Inf;
        vec_para_best_current           = [];
        N_fitting_rounds_current        = factorial(N_humps+1)*N_fitting_rounds;
        for i=1:N_fitting_rounds_current
            [err,vec_para]              = fit_SFS_one_iteration(vec_SFS_real,N_humps,vec_SFS_positions,library_SFS_component);
            if err<err_best_current
                err_best_current        = err;
                vec_para_best_current   = vec_para;
            end
        end
        if N_humps==0
            ratio_error                 = 0;
            err_best_final              = err_best_current;
            vec_para_best_final         = vec_para_best_current;
        else
            ratio_error                 = err_best_current/err_best_final;
            if ratio_error<threshold_stop
                err_best_final          = err_best_current;
                vec_para_best_final     = vec_para_best_current;
            end
        end
        fprintf('%d humps: error = %f; error ratio = %f:      ',N_humps,err_best_current,ratio_error);
        for i=1:length(vec_para_best_current)
            fprintf('%.3f   ',vec_para_best_current(i));
        end
        fprintf('\n');
    end
    table_parameters(n_simulation,1:length(vec_para_best_final))    = vec_para_best_final;
    filename        = strcat('/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/vignettes/TEST_SFS_DECONVOLUTION/TEST_SFS_DECONVOLUTION_deconvolution_parameters_',num2str(n_simulation),'.txt');
    fileID          = fopen(filename,'w');
    fprintf(fileID,'%.3f\t',vec_para_best_final);
    fclose(fileID);
    %---------------------------------------Plot the SFS fitting results
    figure(1);clf;
    %   Prepare the SFS for each category
    vec_SFS_category = zeros(length(SFS_categories), SFS_totalsteps);
    %   Categorize the SFS
    for i=1:length(SFS_categories)
        for j=1:length(vec_refcount)
            no_variant  = vec_altcount(j);
            no_total    = vec_refcount(j)+vec_altcount(j);
            if (no_variant>=min_variant_read)&&(no_total>=min_total_read)&&strcmp(vec_marker{j},SFS_categories{i})
                VAF         = no_variant/no_total;
                pos         = find(vec_freq>=VAF,1);
                vec_SFS_category(i, pos) = vec_SFS_category(i, pos) + 1;
            end
        end
    end
    %   Plot the categorized SFS data
    flag_plot = any(vec_SFS_category' ~= 0);
    SFS_categories_colormap_mini = SFS_categories_colormap(flag_plot, :);
    b = bar(vec_freq,vec_SFS_category(flag_plot,:)', 'stacked');hold on;
    for i = 1:length(b)
        b(i).FaceColor = SFS_categories_colormap(i, :);
    end
    legend(strrep(SFS_categories, '_', ' '), 'Location', 'best');
    %   Plot the total the SFS fit
    vec_A_and_K     = vec_para_best_final(1:2:end);
    vec_p           = vec_para_best_final(2:2:end-1);
    vec_SFS_model   = compute_SFS_one_iteration(vec_A_and_K,vec_p,vec_SFS_positions,library_SFS_component);
    plot(vec_freq, vec_SFS_model, 'r', 'LineWidth', 3, 'DisplayName', 'SFS fit');
    %   Plot the SFS components
    vec_mini_A_and_K    = [vec_A_and_K(1) 0];
    vec_mini_p          = [];
    vec_SFS_model       = compute_SFS_one_iteration(vec_mini_A_and_K,vec_mini_p,vec_SFS_positions,library_SFS_component);
    plot(vec_freq(vec_SFS_model > 0.01), vec_SFS_model(vec_SFS_model > 0.01), 'k:', 'LineWidth', 3, 'DisplayName', 'SFS fit - neutral component');
    for i=1:(length(vec_A_and_K)-1)
        vec_mini_A_and_K    = [0 vec_A_and_K(i+1)];
        vec_mini_p          = [vec_p(i)];
        vec_SFS_model       = compute_SFS_one_iteration(vec_mini_A_and_K,vec_mini_p,vec_SFS_positions,library_SFS_component);
        if (i==1)
            plot(vec_freq(vec_SFS_model > 0.01), vec_SFS_model(vec_SFS_model > 0.01), 'b:', 'LineWidth', 3, 'DisplayName', 'SFS fit - binomial hump(s)');
        else
            plot(vec_freq(vec_SFS_model > 0.01), vec_SFS_model(vec_SFS_model > 0.01), 'b:', 'LineWidth', 3, 'HandleVisibility', 'off');
        end
    end
    %   Save the plot
    filename        = strcat('/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SFS_CNA_deconvolution/vignettes/TEST_SFS_DECONVOLUTION/TEST_SFS_DECONVOLUTION_deconvolution_',num2str(n_simulation),'.png');
    saveas(gcf,filename)
end



%=============================Prepare the sampling coverage distribution
%===============================================for a particular patient
function prep_distribution_patient(vec_totcount)
global TCGA_coverage_values TCGA_coverage_PDF
global TCGA_coverage_values_allTCGA TCGA_coverage_PDF_allTCGA
global r_min r_max
%------------------------------------------Compute coverage distribution
L                       = max(vec_totcount);
TCGA_coverage_values    = (1:L);
TCGA_coverage_PDF       = zeros(1,L);
for i=1:length(vec_totcount)
    pos                     = vec_totcount(i);
    TCGA_coverage_PDF(pos)  = TCGA_coverage_PDF(pos)+1;
end
TCGA_coverage_PDF       = TCGA_coverage_PDF/sum(TCGA_coverage_PDF);
%--------------------------------------If the coverage is over the range
if sum(TCGA_coverage_PDF(1:min(r_max,length(TCGA_coverage_PDF))))<=0
    TCGA_coverage_PDF(r_max)    = 1;
end
end
%===============================================Distribution of coverage
%============================================================phi=Prob(r)
function phi_r = pdf_coverage(r)
global option_dist_coverage dist_coverage_var_1 dist_coverage_var_2
global TCGA_coverage_values TCGA_coverage_PDF
global r_min r_max
global N_end
%   Compute the probability of a read number
%   based on choice of sampling coverage distribution
if strcmp(option_dist_coverage,'uniform')==1
    if (r<r_min)||(r>r_max)||(r_min>r_max)
        phi_r   = 0;
    elseif r_min==r_max
        if (dist_coverage_var_1<=r_min)&&(dist_coverage_var_2>=r_min)
            phi_r   = 1;
        else
            phi_r   = 0;
        end
    else
        phi_r   = 1/(dist_coverage_var_2-dist_coverage_var_1);
    end
elseif strcmp(option_dist_coverage,'binomial')==1
    D       = dist_coverage_var_1;
    if (r<r_min)||(r>r_max)||(r_min>r_max)
        phi_r   = 0;
    elseif D>0
        phi_r   = binopdf(r,N_end,D/N_end);
    else
        phi_r   = -1;
    end
elseif strcmp(option_dist_coverage,'TCGA')==1
    pos = find(TCGA_coverage_values==r);
    if (isempty(pos)==1)||(r<r_min)||(r>r_max)||(r_min>r_max)
        phi_r   = 0;
    else
        phi_r   = TCGA_coverage_PDF(pos);
    end
end
end
%==============================Compute the expected Griffiths-Tavare SFS
function vec_SFS_GT = SFS_Griffiths_Tavare(vec_para)
global N_end
%-----------------------------------------------------Get the parameters
no_hump = (length(vec_para)-1)/2;
para_A  = vec_para(1);
para_K  = zeros(1,no_hump);
para_P  = zeros(1,no_hump);
for i=1:no_hump
    para_P(i)   = vec_para(2*i);
    para_K(i)   = vec_para(2*i+1);
    % para_K(i)   = vec_para(i+1);
    % para_P(i)   = vec_para(no_hump+i+1);
end
%---------------------------------------Compute the Griffiths-Tavare SFS
vec_SFS_GT  = zeros(1,N_end);
for m=2:N_end
    vec_SFS_GT(m)   = para_A*N_end/(m*(m-1));
    for i=1:no_hump
        K   = para_K(i);
        P   = para_P(i);
        vec_SFS_GT(m)   = vec_SFS_GT(m)+K*binopdf(m,N_end,P);
    end
end
end
%=================================Compute the theoretically expected SFS
%================based on the Griffiths-Tavare SFS and the read coverage
function vec_SFS_expected = SFS_expected(vec_para,package_input)
global matrix_binomial_PDF
global SFS_totalsteps SFS_totalsteps_base
global vec_SFS_freq
global vec_SFS_expected
global N_end
global r_min r_max
%---------------------------------------Compute the Griffiths-Tavare SFS
vec_SFS_GT = SFS_Griffiths_Tavare(vec_para);
%---------------------------------Compute the theoretically expected SFS
%   Prepare the vectors
vec_SFS_freq    = [0:1:SFS_totalsteps]/SFS_totalsteps;
vec_SFS_expected= zeros(1,SFS_totalsteps);
%   Compute the theoretically expected SFS
for i=1:SFS_totalsteps
    j_lower = round(SFS_totalsteps_base*vec_SFS_freq(i))+1;
    j_upper = round(SFS_totalsteps_base*vec_SFS_freq(i+1));
    omega   = 0;
    for r=max(r_min,1):r_max
        PDF_coverage= pdf_coverage(r);
        if PDF_coverage>0
            Sum = 0;
            for m=1:N_end
                q_m = vec_SFS_GT(m);
                if q_m>0
                    PDF_success = 0;
                    for j=j_lower:j_upper
                        PDF_success = PDF_success+matrix_binomial_PDF(r,m,j);
                    end
                    Sum         = Sum+q_m*PDF_success;
                end
            end
            omega   = omega+PDF_coverage*Sum;
        end
    end
    vec_SFS_expected(i)   = omega;
end
end
%=======================================================================
function [err,vec_para] = fit_SFS_one_iteration(vec_SFS_real,N_humps,vec_SFS_positions,library_SFS_component)
%------------------------------------Choose the hump locations at random
vec_p                           = sort(randsample(vec_SFS_positions,N_humps));
%------------------------------------------------------Fit for A and K's
func_fit                        = @(vec_A_and_K) error_SFS_one_iteration(vec_A_and_K,vec_p,vec_SFS_positions,library_SFS_component,vec_SFS_real);
vec_A_and_K_initial             = [1 100*ones(1,N_humps)];
options                         = optimset('fminsearch');
options.MaxFunEvals             = 10^6;
options.MaxIter                 = 10^6;
[vec_A_and_K,err]               = fminsearch(func_fit,vec_A_and_K_initial,options);
vec_para(1)                     = vec_A_and_K(1);
for i=1:N_humps
    vec_para(2*i)               = vec_p(i);
    vec_para(2*i+1)             = vec_A_and_K(i+1);
end
end
%=======================================================================
function vec_SFS_model = compute_SFS_one_iteration(vec_A_and_K,vec_p,vec_SFS_positions,library_SFS_component)
%   Add the neutral component
A                               = vec_A_and_K(1);
vec_SFS_model                   = A*library_SFS_component{1,1};
%   Add the binomial humps
for i_hump=1:length(vec_p)
    p                           = vec_p(i_hump);
    loc                         = find(vec_SFS_positions==p);
    K                           = vec_A_and_K(i_hump+1);
    vec_SFS_model               = vec_SFS_model+K*library_SFS_component{2,loc};
end
end
%=======================================================================
function output = error_SFS_one_iteration(vec_A_and_K,vec_p,vec_SFS_positions,library_SFS_component,vec_SFS_real)
%--------------------------------------Check point to make sure A, K > 0
if min(vec_A_and_K)<0
    output                      = Inf;
    return;
end
%-------------------------Compute the expected SFS with input parameters
vec_SFS_model                   = compute_SFS_one_iteration(vec_A_and_K,vec_p,vec_SFS_positions,library_SFS_component);
%-----------------Compute the 1-norm error between expected and real SFS
output                          = norm(vec_SFS_model-vec_SFS_real,1)/norm(vec_SFS_model,1);
end
