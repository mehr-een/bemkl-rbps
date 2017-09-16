%Mehreen Ali
%mehreen.ali@helsinki.fi

function [validation_response_CV,predicted_response_CV] = bemkl_loocv(drug_source, view_index)

    load 'DataViews.mat';
    load 'ViewCombinations.mat'; % combos with mut_data

    %%%% selected views
    view_names = view_combinations{view_index,1}; 
    view_kernels = view_combinations{view_index,2};

    index_k = 1:1:58; 
    num_drugs = size(drug_source,2);

    predicted_response_CV=zeros(58, num_drugs);
    
    P = length(view_names);
    standardize_response = 0; %0:false 1:true

    %%%% LOO-CV
    for k = 1:58
    
        learning_indices = find(index_k ~= k);
        prediction_indices = k;

        Nlearning = length(learning_indices);
        Nprediction = length(prediction_indices);

        %%%% views dimensions
        learning.X = cell(1, P);
        prediction.X = cell(1, P);
        for m = 1:P
            learning.X{m} = eval(sprintf('%s(learning_indices, :)', view_names{m}));
            prediction.X{m} = eval(sprintf('%s(prediction_indices, :)', view_names{m}));
        end

        learning.Y = eval('drug_source(learning_indices, :)');
        validation.Y = eval('drug_source(prediction_indices, :)');
        Ndrug = size(learning.Y, 2);
    
        seed = 1606;
        rand('state', seed); %#ok<RAND>
        randn('state', seed); %#ok<RAND>
        
        %%%% view normalization
        for m = 1:P
            nan_indices_learning = isnan(learning.X{m});
            nan_indices_prediction = isnan(prediction.X{m});
            if strcmp(view_kernels{m}, 'jaccard') ~= 1
                mas.mea = nanmean(learning.X{m}, 1);
                mas.std = nanstd(learning.X{m}, 0, 1);
                mas.std(mas.std == 0) = 1;
                learning.X{m} = bsxfun(@rdivide, bsxfun(@minus, learning.X{m}, mas.mea), mas.std);
                prediction.X{m} = bsxfun(@rdivide, bsxfun(@minus, prediction.X{m}, mas.mea), mas.std);
            end
            learning.X{m}(nan_indices_learning) = 0;
            prediction.X{m}(nan_indices_prediction) = 0;
        end
    
        %%%% drug response standardization
        if standardize_response == 0
            mas.mea = nanmean(learning.Y, 1);
            learning.Y = bsxfun(@minus, learning.Y, mas.mea);
            validation.Y = bsxfun(@minus,validation.Y , mas.mea);
        else
            mas.mea = nanmean(learning.Y, 1);
            mas.std = nanstd(learning.Y, 0, 1);
            mas.std(mas.std == 0) = 1;
            learning.Y = bsxfun(@rdivide, bsxfun(@minus, learning.Y, mas.mea), mas.std);
            validation.Y = bsxfun(@rdivide, bsxfun(@minus, validation.Y, mas.mea), mas.std);
        end
        validation_response_CV(k,:) = validation.Y;
    
    
        %%%% kernels for views
        Kx_learning = zeros(Nlearning, Nlearning, P);
        Kx_prediction = zeros(Nlearning, Nprediction, P);

        for m = 1:P
            if strcmp(view_kernels{m}, 'gaussian') == 1
                Kx_learning(:, :, m) = exp(-pdist2(learning.X{m}, learning.X{m}).^2 / size(learning.X{m}, 2) / 2);
                Kx_prediction(:, :, m) = exp(-pdist2(learning.X{m}, prediction.X{m}).^2 / size(learning.X{m}, 2) / 2);
                
            elseif strcmp(view_kernels{m}, 'exponential') == 1
                Kx_learning(:, :, m) = exp(-pdist2(learning.X{m}, learning.X{m}) / sqrt(size(learning.X{m}, 2)  ));
                Kx_prediction(:, :, m) = exp(-pdist2(learning.X{m}, prediction.X{m}) / sqrt(size(learning.X{m}, 2)  ));

            elseif strcmp(view_kernels{m}, 'jaccard') == 1
                J = 1 - pdist2(learning.X{m}, learning.X{m}, 'jaccard');
                J(isnan(J)) = 0;
                for j = 1:size(J, 1)
                    J(j, j) = 1;
                end
                Kx_learning(:, :, m) = J;
                J = 1 - pdist2(learning.X{m}, prediction.X{m}, 'jaccard');
                J(isnan(J)) = 0;
                Kx_prediction(:, :, m) = J;
            end
        end
    
        %%%% NCI/DREAM7 parameters 
        parameters.alpha_lambda = 1e-10; %positive
        parameters.beta_lambda = 1e-10; %positive
        parameters.alpha_upsilon = 1; %positive
        parameters.beta_upsilon = 1; %positive
        parameters.alpha_gamma = 1e-10; %positive
        parameters.beta_gamma = 1e-10; %positive
        parameters.alpha_omega = 1e-10; %positive
        parameters.beta_omega = 1e-10; %positive
        parameters.alpha_epsilon = 1; %positive
        parameters.beta_epsilon = 1; %positive
        parameters.iteration = 200;

        parameters.seed = 1606;
   
        %%%% combined kernel
        K_learning = cell(1, Ndrug);
        y_learning = cell(1, Ndrug);
        K_prediction = cell(1, Ndrug);
        for i = 1:Ndrug
            indices = find(isnan(learning.Y(:, i)) == 0);
            K_learning{i} = zeros(length(indices), length(indices), P);
            for m = 1:P
                K_learning{i}(:, :, m) = Kx_learning(indices, indices, m);
            end
            y_learning{i} = learning.Y(indices, i);
            K_prediction{i} = zeros(size(Kx_prediction, 2), length(indices), P);
            for m = 1:P
                K_prediction{i}(:, :, m) = Kx_prediction(indices, :, m)';
            end
        end   
        
        state = bemkl_supervised_multitask_regression_variational_train(K_learning, y_learning, parameters);

        output = bemkl_supervised_multitask_regression_variational_test(K_prediction, state);
        predicted_response=zeros(length(prediction_indices), Ndrug);
       
        for drug = 1:Ndrug
            predicted_response(:, drug) = output.f{drug}.mean;
        end
    
        predicted_response_CV(k,:) = predicted_response;

    end
end
