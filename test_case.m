
% LRSMA Parameters
Threshold = 0.1;

% Simulation Parameters
nsims=100;

common_frequency= 0.15 ;
indv_frequency=0.80;

window_size=3000;
nsamples=100;


% Common peaks parameters
number_of_common_peaks = 1;
common_peak_avg_center= [2500];
common_peak_sd_center= [1];
common_peak_avg_width= [200];
common_peak_sd_width= [10];
common_peak_avg_height= [10.0];
common_peak_sd_height= [0.05];

% Sample-specific peaks parameters
num_indv_peaks = 1;
indv_min_width = 200;
indv_max_width = 200;
indv_max_height = 10.0;


% Output file name
filename = 'test_sim.txt';


% Summary variables
LRMA_tot_TPR = [];
LRMA_tot_FPR = [];


for iter=1:nsims
    % Data Matrices
    D_wig = zeros(nsamples,window_size);
    D_bed = zeros(nsamples,window_size);
    
    common_peaks_list = [];
    
    i=1;
    for cp=1:number_of_common_peaks
        common_samples = datasample(1:nsamples,floor(common_frequency*nsamples), 'Replace',false);
        
        for sample=common_samples
            common_peak_center= floor(abs(common_peak_sd_center(cp).*randn + common_peak_avg_center(cp)));
            common_peak_width= floor(abs(common_peak_sd_width(cp).*randn + common_peak_avg_width(cp)));
            common_peak_max_height= double(abs(common_peak_sd_height(cp).*randn + common_peak_avg_height(cp)));
            
            lower_bound=max(int16(common_peak_center-common_peak_width/2),1);
            upper_bound=min(int16(common_peak_center+common_peak_width/2), window_size);
            
            scaling=normpdf(double(lower_bound:upper_bound), common_peak_center, common_peak_width./6);
            common_peak_height=scaling./max(scaling).*common_peak_max_height;
            
            common_peaks_list(i,:)= [sample, lower_bound,upper_bound] ;
            
            D_wig(sample, lower_bound:upper_bound)=common_peak_height;
            D_bed(sample, lower_bound:upper_bound)=1.0;
            i=i+1;
        end
    end
    
        
    indv_peaks_list = [];
    i=1;
    indv_samples = datasample(1:nsamples,floor(indv_frequency*nsamples));
    for s=1:num_indv_peaks
        for sample=indv_samples
            indv_peak_center=randi([1000+0.5*indv_max_width,2500-0.5*indv_max_width]);
            indv_peak_width=randi([indv_min_width,indv_max_width]);
            indv_peak_max_height=double(randi([1,indv_max_height]));
            
            lower_bound=int16(indv_peak_center-indv_peak_width/2);
            upper_bound=int16(indv_peak_center+indv_peak_width/2);
            
            scaling=normpdf(double(lower_bound:upper_bound), indv_peak_center, indv_peak_width./6);
            indv_peak_height=scaling./max(scaling).*indv_peak_max_height;
            
            indv_peaks_list(i,:)= [sample, lower_bound,upper_bound] ;
            i=i+1;
            
            %         D_wig(sample, indv_peak_center)=indv_peak_width;
            D_wig(sample, lower_bound:upper_bound)=indv_peak_height;
            D_bed(sample, lower_bound:upper_bound)=1.0;
            
        end
    end
    
    D_wig = double(D_wig);
    
    
    
    % LRSMA to get common peaks
    D_wig=D_wig';
    D_wig(isnan(D_wig)) = 0;
    [m,n] = size(D_wig);
    
    p1=0.75*mean(max(D_wig));
    p2=1.5*mean(var(D_wig));
    
    % estimate sigma
    d = D_wig(abs(D_wig)<p1);
    sig = 10*1.48*mean(abs(d(:)-mean(d(:))));
    
    % Run PLA
    beta1 = sqrt(m)*sig;
    beta2 = 0.01*beta1;
    beta3 = 2*sig;
    tol = 1e-4;
        
    
    tic;
    % RPLA Engine
    [B,E] = LRSMA(D_wig,beta1,beta2,beta3,tol);
    t_pla = toc;
    B=abs(B);
    B2=B./norm(B(:)).*norm(D_wig(:));
    B3=floor(B2);
    F=(sum(B3>0,2)/n);
    
    intervals = makeIntervals(F, Threshold, 5);
    
    gain_merge = sum(D_bed'>0,2);
    F_merge=gain_merge./n;
    
    merge_intervals = makeIntervals(F_merge, Threshold, 5);
    
    
    
    
    
    % ROC Curve
    merge_TPR = [];
    merge_FPR = [];
    i = 1;
    for F_th=linspace(1,0,nsteps)
        merge_TP = 0;
        merge_FN = 0;
        for peak=1:size(common_peaks_list,1)
            if max(F_merge(common_peaks_list(peak,2):common_peaks_list(peak,3))) >= F_th
                merge_TP=merge_TP+1;
            else
                merge_FN=merge_FN+1;
            end
        end
        
        merge_TN = 0;
        merge_FP = 0;
        for peak=1:size(indv_peaks_list,1)
            if max(F_merge(indv_peaks_list(peak,2):indv_peaks_list(peak,3))) >= F_th
                merge_FP=merge_FP+1;
            else
                merge_TN=merge_TN+1;
            end
        end
        merge_TPR(i) = merge_TP / (merge_TP+merge_FN);
        merge_FPR(i) = merge_FP / max((merge_FP+merge_TN),1);
        i=i+1;
    end
    merge_tot_TPR = cat(2,merge_tot_TPR, merge_TPR);
    merge_tot_FPR = cat(2,merge_tot_FPR, merge_FPR);

    
    
    % ROC Curve
    LRMA_TPR = [];
    LRMA_FPR = [];
    i = 1;
    for F_th=linspace(1,0,nsteps)
        LRMA_TP = 0;
        LRMA_FN = 0;
        for peak=1:size(common_peaks_list,1)
            if max(F(common_peaks_list(peak,2):common_peaks_list(peak,3))) >= F_th
                LRMA_TP=LRMA_TP+1;
            else
                LRMA_FN=LRMA_FN+1;
            end
        end
        
        LRMA_TN = 0;
        LRMA_FP = 0;
        for peak=1:size(indv_peaks_list,1)
            if max(F(indv_peaks_list(peak,2):indv_peaks_list(peak,3))) >= F_th
                LRMA_FP=LRMA_FP+1;
            else
                LRMA_TN=LRMA_TN+1;
            end
        end
        LRMA_TPR(i) = LRMA_TP / (LRMA_TP+LRMA_FN);
        LRMA_FPR(i) = LRMA_FP / max((LRMA_FP+LRMA_TN),1);
        i=i+1;
    end
    LRMA_tot_TPR = cat(2,LRMA_tot_TPR, LRMA_TPR);
    LRMA_tot_FPR = cat(2,LRMA_tot_FPR, LRMA_FPR);
    
        
    
end




