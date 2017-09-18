%Mehreen Ali
%mehreen.ali@helsinki.fi

load 'DrugResponse.mat';
view_combination_perf(1:12) = struct(); %

%%%% cytotoxic agents - demo
drug_source = cytotoxic_response;
selected_combinations = [1]; %indices from ViewCombinations.mat

%%%% cytotoxic agents
%drug_source = cytotoxic_response;
%selected_combinations = [1 2 3 4 5 6 7 8 9 10 11 12]; %indices from ViewCombinations.mat

%%%% targeted agents
%drug_source = taregted_response;
%selected_combinations = [1 2 3 4 5 6 8 9 10 11 13]; %indices from ViewCombinations.mat

for k = 1:length(selected_combinations)
    view_index = selected_combinations(k);
    [validation_response_CV,predicted_response_CV] = bemkl_loocv(cell2mat(drug_source), view_index);

    view_combination_perf(view_index).validation_response_CV = validation_response_CV;
    view_combination_perf(view_index).predicted_response_CV = predicted_response_CV; 
end


perf_measures = zeros(length(selected_combinations),3);
for j = 1:length(selected_combinations)
    i = selected_combinations(j);
    % Spearman cor
    perf_measures(i,1) = nanmean(diag(corr(view_combination_perf(i).validation_response_CV, view_combination_perf(i).predicted_response_CV,'type','Spearman', 'rows','pairwise')));
    % RMSE
    perf_measures(i,2) = mean(sqrt(nanmean((view_combination_perf(i).validation_response_CV-view_combination_perf(i).predicted_response_CV).^2)));
    % CI-index
    perf_measures(i,3) = civalue(view_combination_perf(i).validation_response_CV, view_combination_perf(i).predicted_response_CV);
end

