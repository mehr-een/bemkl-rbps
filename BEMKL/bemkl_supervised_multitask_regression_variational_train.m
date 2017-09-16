function state = bemkl_supervised_multitask_regression_variational_train(Km, y, parameters)
    rand('state', parameters.seed); %#ok<RAND>
    randn('state', parameters.seed); %#ok<RAND>

    T = length(Km);
    N = zeros(T, 1);
    for o = 1:T
        N(o) = size(Km{o}, 2);
    end
    P = size(Km{1}, 3);

    lambda = cell(1, T);
    a = cell(1, T);    
    for o = 1:T
        lambda{o}.shape = (parameters.alpha_lambda + 0.5) * ones(N(o), 1);
        lambda{o}.scale = parameters.beta_lambda * ones(N(o), 1);
        a{o}.mean = randn(N(o), 1);
        a{o}.covariance = eye(N(o), N(o));
    end
    upsilon.shape = (parameters.alpha_upsilon + 0.5 * N * P) .* ones(T, 1);
    upsilon.scale = parameters.beta_upsilon * ones(T, 1);
    G = cell(1, T);
    for o = 1:T
        G{o}.mean = randn(P, N(o));
        G{o}.covariance = eye(P, P);
    end
    gamma.shape = (parameters.alpha_gamma + 0.5) * ones(T, 1);
    gamma.scale = parameters.beta_gamma * ones(T, 1);    
    omega.shape = (parameters.alpha_omega + 0.5) * ones(P, 1);
    omega.scale = parameters.beta_omega * ones(P, 1);
    be.mean = [zeros(T, 1); ones(P, 1)];
    be.covariance = eye(T + P, T + P);
    epsilon.shape = (parameters.alpha_epsilon + 0.5 * N) .* ones(T, 1);
    epsilon.scale = parameters.beta_epsilon * ones(T, 1);

    KmKm = cell(1, T);
    for o = 1:T
        KmKm{o} = zeros(N(o), N(o));
        for m = 1:P
            KmKm{o} = KmKm{o} + Km{o}(:, :, m) * Km{o}(:, :, m)';
        end
        Km{o} = reshape(Km{o}, [N(o), N(o) * P]);
    end

    for iter = 1:parameters.iteration
        if mod(iter, 1) == 0
            fprintf(1, '.');
        end
        if mod(iter, 10) == 0
            fprintf(1, ' %5d\n', iter);
        end

        %%%% update lambda
        for o = 1:T
            lambda{o}.scale = 1 ./ (1 / parameters.beta_lambda + 0.5 * (a{o}.mean.^2 + diag(a{o}.covariance)));
        end
        %%%% update upsilon
        for o = 1:T
            upsilon.scale(o) = 1 / (1 / parameters.beta_upsilon + 0.5 * (sum(diag(G{o}.mean * G{o}.mean' + N(o) * G{o}.covariance)) ...
                                                                         - 2 * sum(sum(reshape(a{o}.mean' * Km{o}, [N(o), P])' .* G{o}.mean)) ...
                                                                         + sum(diag(KmKm{o} * (a{o}.mean * a{o}.mean' + a{o}.covariance)))));
        end
        %%%% update a
        for o = 1:T
            a{o}.covariance = (diag(lambda{o}.shape .* lambda{o}.scale) + upsilon.shape(o) * upsilon.scale(o) * KmKm{o}) \ eye(N(o), N(o));
            a{o}.mean = a{o}.covariance * (upsilon.shape(o) * upsilon.scale(o) * Km{o} * reshape(G{o}.mean', N(o) * P, 1));
        end
        %%%% update G        
        for o = 1:T
            G{o}.covariance = (upsilon.shape(o) * upsilon.scale(o) * eye(P, P) + epsilon.shape(o) * epsilon.scale(o) * (be.mean(T + 1:T + P) * be.mean(T + 1:T + P)' + be.covariance(T + 1:T + P, T + 1:T + P))) \ eye(P, P);
            G{o}.mean = G{o}.covariance * (upsilon.shape(o) * upsilon.scale(o) * reshape(a{o}.mean' * Km{o}, [N(o), P])' + epsilon.shape(o) * epsilon.scale(o) * (be.mean(T + 1:T + P) * y{o}' - repmat(be.mean(T + 1:T + P) * be.mean(o) + be.covariance(T + 1:T + P, o), 1, N(o))));
        end   
        %%%% update gamma
        gamma.scale = 1 ./ (1 / parameters.beta_gamma + 0.5 * (be.mean(1:T).^2 + diag(be.covariance(1:T, 1:T))));
        %%%% update omega
        omega.scale = 1 ./ (1 / parameters.beta_omega + 0.5 * (be.mean(T + 1:T + P).^2 + diag(be.covariance(T + 1:T + P, T + 1:T + P))));
        %%%% update epsilon
        for o = 1:T
            epsilon.scale(o) = 1 / (1 / parameters.beta_epsilon + 0.5 * (y{o}' * y{o} - 2 * y{o}' * [ones(1, N(o)); G{o}.mean]' * be.mean([o, T + 1:T + P]) ...
                                                                         + N(o) * (be.mean(o)^2 + be.covariance(o, o)) ...
                                                                         + sum(diag((G{o}.mean * G{o}.mean' + N(o) * G{o}.covariance) * (be.mean(T + 1:T + P) * be.mean(T + 1:T + P)' + be.covariance(T + 1:T + P, T + 1:T + P)))) ...
                                                                         + 2 * sum(diag(sum(G{o}.mean, 2)' * (be.mean(T + 1:T + P) * be.mean(o) + be.covariance(T + 1:T + P, o))))));
        end
        %%%% update b and e
        be.covariance = [diag(gamma.shape .* gamma.scale) + diag(N .* epsilon.shape .* epsilon.scale), zeros(T, P); ...
                         zeros(P, T), diag(omega.shape .* omega.scale)];
        for o = 1:T
            be.covariance(T + 1:T + P, o) = epsilon.shape(o) * epsilon.scale(o) * sum(G{o}.mean, 2);
            be.covariance(o, T + 1:T + P) = epsilon.shape(o) * epsilon.scale(o) * sum(G{o}.mean, 2)';
            be.covariance(T + 1:T + P, T + 1:T + P) = be.covariance(T + 1:T + P, T + 1:T + P) + epsilon.shape(o) * epsilon.scale(o) * (G{o}.mean * G{o}.mean' + N(o) * G{o}.covariance);
        end
        be.covariance = be.covariance \ eye(T + P, T + P);
        be.mean = zeros(T + P, 1);        
        for o = 1:T
            be.mean(o) = epsilon.shape(o) * epsilon.scale(o) * sum(y{o});
            be.mean(T + 1:T + P) = be.mean(T + 1:T + P) + epsilon.shape(o) * epsilon.scale(o) * G{o}.mean * y{o};
        end
        be.mean = be.covariance * be.mean;
    end

    state.lambda = lambda;
    state.a = a;
    state.upsilon = upsilon;
    state.gamma = gamma;
    state.omega = omega;
    state.be = be;
    state.epsilon = epsilon;
end