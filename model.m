function dx = SIV(t, x, betas, sigmas, k, alphas, Incomes, gammas, lambdas)
    num_compartments = length(x) / 5;
    dx = zeros(5 * num_compartments, 1);
    
    for j = 1:num_compartments
        Sj = x((j-1)*5 + 1);
        Ej = x((j-1)*5 + 2);
        Ij = x((j-1)*5 + 3);
        Rj = x((j-1)*5 + 4);
        Vj = x((j-1)*5 + 5);
        
        Beta = betas * reshape(x(3:5:end), num_compartments, 1);

        dx((j-1)*5 + 1) = -Beta(j) * Sj + sigmas(j) * Rj;
        dx((j-1)*5 + 2) = Beta(j) * Sj - k * Ej;
        dx((j-1)*5 + 3) = k * Ej - alphas(j) * Incomes(j) * Ij;
        dx((j-1)*5 + 4) = alphas(j) * Incomes(j) * Ij - sigmas(j) * Rj;
        dx((j-1)*5 + 5) = gammas * (lambdas(1) * Ej + lambdas(2) * Ij + lambdas(3) * Rj);
    end
end

function betas = calculate_betas(n_values, constant, eta_vax, distances)
    num_compartments = length(n_values);
    beta_star = 1.67391974668301e-08;
    betas = zeros(num_compartments, num_compartments);
    
    % Loop to calculate betas matrix
    for i = 1:num_compartments
        for j = 1:num_compartments
            distance_between_i_and_j = distances(i, j); % Assuming the distance matrix is properly formatted
            
            % Calculate betas based on the given formula for transmission matrix
            betas(i, j) = beta_star / n_values(j) * 1 / (1 + (distance_between_i_and_j) / constant) * eta_vax(j);
        end
    end
end