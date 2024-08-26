
# No population mobility
function solve_model_immob(param, fund, dist, tol = 0.01, verbose = false)

    xtic = time()

    A = fund.A
    L = fund.L
    alpha = param.alpha
    sigma = param.sigma

    # Initialize wages and population shapes
    w_i = ones(length(A))
    
    # Trade costs, as in the numerator of trade share equation
    dd = dist .^ (1-sigma)

    num_iter = 1; converged = 0
    tradesh = 0; income = 0; expend = 0;
    abs_ratio_diff = 0
    while num_iter < 100000
            
        # Trade share equation
        pwmat = (w_i ./ (alpha .* A .* L .^ (alpha-1))) .^ (1-sigma)
        nummat = dd .* pwmat'
        tradesh = nummat ./ sum(nummat, dims = 1)
        # # test;
        # test = sum(tradesh, dims = 1)
        # mean(test) # Should be 1

        # Income equals expenditure equilibrium condition;
        income = w_i .* L
        expend = tradesh * income

        ## Rounding estimates to check convergency subject to tolerance:
        # income_r = round.(income ./ tol)
        # expend_r = round.(expend ./ tol)
        # MAPE = mean(abs.(income - expend) ./ expend)

        # Update loop;
        abs_ratio_diff = mean(abs.(income ./ expend .- 1))
        if abs_ratio_diff < tol # income_r == expend_r
            display(">>>> Convergence Achieved <<<<")
            converged = 1
            break
        else
            if verbose 
                print("  ARD(", num_iter, "): ", round(abs_ratio_diff / (tol/100)) * tol)
            end
            # In the model income = w*L = expend: we take the denominator as given and just adjust the numerator in the trade share eq., but not fully...
            w_e = w_i .* (expend ./ income) .^ (1/(sigma-1))
            # Gradual updating
            w_i = (0.25 .* w_e) + (0.75 .* w_i)
            # Normalization: choose geometric mean wage as numeraire;
            w_i /= geomean(w_i)
            num_iter += 1
        end
    end
    xtic = time() - xtic
    return (w = w_i, L = L, tradesh = tradesh, dtradesh = diag(tradesh), income = income, 
            expend = expend, converged = converged, abs_ratio_diff = abs_ratio_diff, 
            tol = tol, num_iter = num_iter, xtic = xtic)
end

# Full population mobility
function solve_model_free_mob(param, fund, dist, tol = 0.1, verbose = false)

    xtic = time()

    A = fund.A
    L = fund.L
    H = fund.H
    B = fund.B
    alpha = param.alpha
    sigma = param.sigma
    theta = param.theta
    beta = param.beta

    # Initialize population
    L_i = L
    L_tot = sum(L)

    # Initialize wages
    w_i = ones(length(A))
    
    # Trade costs, as in the numerator of trade share equation
    dd = dist .^ (1-sigma)

    num_iter = 1
    converged = 0
    tradesh = 0
    income = 0
    expend = 0
    inc_exp_abs_ratio_diff = 0
    pop_abs_ratio_diff = 0
    # TODO: Nested iteration as in Redding (2016)
    while num_iter < 100000
            
        # Trade share equation
        pwmat = (w_i ./ (alpha .* A .* L_i .^ (alpha-1))) .^ (1-sigma)
        nummat = dd .* pwmat'
        tradesh = nummat ./ sum(nummat, dims = 1)
        # # test;
        # test = sum(tradesh, dims = 1)
        # mean(test) # Should be 1

        # Income equals expenditure equilibrium condition;
        income = w_i .* L_i
        expend = tradesh * income

        # Domestic trade share
        dtradesh = diag(tradesh)

        # Population mobility equilibrium condition;
        num = ((1/B) .* dtradesh .^ (theta/(sigma-1)) .* A .^ (-theta) .* H .^ (theta-1)) .^ (1/(theta*alpha + beta))
        L_e = (num ./ sum(num)) .* L_tot # updated population estimate

        # # Rounding estimates to check convergency subject to tolerance:
        # income_r = round.(income ./ tol)
        # expend_r = round.(expend ./ tol)
        # L_i_r = round.(L_i ./ tol)
        # L_e_r = round.(L_e ./ tol)
        
        # Better do this with absolute ratio differneces: less strict & faster 
        inc_exp_abs_ratio_diff = mean(abs.(income ./ expend .- 1))
        pop_abs_ratio_diff = mean(abs.(L_i ./ L_e .- 1))

        # Update loop;
        if inc_exp_abs_ratio_diff < tol & pop_abs_ratio_diff < tol # (income_r == expend_r) & (L_i_r == L_e_r)
            display(">>>> Convergence Achieved <<<<")
            converged = 1
            break
        else
            if verbose
                print("  ARD(", num_iter, "): ", round(inc_exp_abs_ratio_diff / (tol/100)) * tol, " ", round(pop_abs_ratio_diff / (tol/100)) * tol)
            end
            # In the model income = w*L = expend: we take the denominator as given and just adjust the numerator, but not fully...
            w_e = w_i .* (expend ./ income) .^ (1/(sigma-1))
            # Gradual updating
            w_i = (0.25 .* w_e) + (0.75 .* w_i)
            L_i = (0.25 .* L_e) + (0.75 .* L_i)
            # Normalization: choose geometric mean wage as numeraire;
            w_i /= geomean(w_i)
            num_iter += 1
        end
    end
    xtic = time() - xtic
    return (w = w_i, L = L_i, tradesh = tradesh, dtradesh = dtradesh, income = income, 
            expend = expend, converged = converged, inc_exp_abs_ratio_diff = inc_exp_abs_ratio_diff, 
            pop_abs_ratio_diff = pop_abs_ratio_diff, tol = tol, num_iter = num_iter, xtic = xtic)
end

# Full population mobility within country
function solve_model_free_mob_in_ctry(param, fund, dist, tol = 0.1, verbose = false)

    xtic = time()

    A = fund.A
    L = fund.L
    H = fund.H
    B = fund.B
    alpha = param.alpha
    sigma = param.sigma
    theta = param.theta
    beta = param.beta

    # Initialize population
    L_i = L
    L_ctry = zeros(nctry)
    for c in 1:nctry
        L_ctry[c] = sum(L[cindl[c]])
    end

    # Initialize wages and updated population within each country
    w_i = ones(length(A))
    L_e = zeros(length(A))
    
    # Trade costs, as in the numerator of trade share equation
    dd = dist^(1-sigma)

    num_iter = 1
    converged = 0
    tradesh = 0
    income = 0
    expend = 0
    inc_exp_abs_ratio_diff = 0
    pop_abs_ratio_diff = 0
    while num_iter < 100000
            
        # Trade share equation
        pwmat = (w_i ./ (alpha .* A .* L_i .^ (alpha-1))) .^ (1-sigma)
        nummat = dd .* pwmat'
        tradesh = nummat ./ sum(nummat, dims = 1)
        # # test;
        # test = sum(tradesh, dims = 1)
        # mean(test) # Should be 1

        # Income equals expenditure equilibrium condition;
        income = w_i .* L_i
        expend = tradesh * income

        # Domestic trade share
        dtradesh = diag(tradesh)

        # Population mobility equilibrium condition;
        num = ((1/B) .* dtradesh .^ (theta/(sigma-1)) .* A .^ (-theta) .* H .^ (theta-1)) .^ (1/(theta*alpha + beta))
        # This ensures population immobility across countries (ensuring the countries overall population shares are the same)
        for c in 1:nctry
            ind = cindl[c]
            L_e[ind] = (num[ind] ./ sum(num[ind])) .* L_ctry[c]
        end

        # # Rounding estimates to check convergency subject to tolerance:
        # income_r = round.(income ./ tol)
        # expend_r = round.(expend ./ tol)
        # L_i_r = round.(L_i ./ tol)
        # L_e_r = round.(L_e ./ tol)
        # Better do this with absolute ratio differneces: less strict & faster
        inc_exp_abs_ratio_diff = mean(abs.(income ./ expend .- 1))
        pop_abs_ratio_diff = mean(abs.(L_i ./ L_e .- 1))

        # Update loop;
        if inc_exp_abs_ratio_diff < tol & pop_abs_ratio_diff < tol # (income_r == expend_r) & (L_i_r == L_e_r)
            display(">>>> Convergence Achieved <<<<")
            converged = 1
            break
        else
            if verbose
                print("  ARD(", num_iter, "): ", round(inc_exp_abs_ratio_diff / (tol/100)) * tol, " ", round(pop_abs_ratio_diff / (tol/100)) * tol)
            end
            # In the model income = w*L = expend
            # Try to understand this !!
            w_e = w_i .* (expend ./ income) .^ (1/(sigma-1))
            # Gradual updating
            w_i = (0.25 .* w_e) + (0.75 .* w_i)
            L_i = (0.25 .* L_e) + (0.75 .* L_i)
            # Normalization: geometric mean of country 1 wage as numeraire
            w_i[cindl[1]] /= geomean(w_i[cindl[1]])
            num_iter += 1
        end
    end
    xtic = time() - xtic
    return (w = w_i, L = L_i, tradesh = tradesh, dtradesh = dtradesh, income = income, 
            expend = expend, converged = converged, inc_exp_abs_ratio_diff = inc_exp_abs_ratio_diff, 
            pop_abs_ratio_diff = pop_abs_ratio_diff, tol = tol, num_iter = num_iter, xtic = xtic)
end

function solve_model_ab(param, fund, w, dist, tol = 0.1, verbose = false)
    
    xtic = time()

    L = fund.L
    H = fund.H
    income = w .* L
    L_tot = sum(L)

    alpha = param.alpha
    sigma = param.sigma
    theta = param.theta
    beta = param.beta

    A_i = ones(length(L))
    B_i = ones(length(L))

     # Trade costs, as in the numerator of trade share equation
     dd = dist .^ (1-sigma)

    println(">>>> Start productivity convergence <<<<")

    num_iter_prod = 1
    tradesh = 0
    aconverge = 0
    inc_exp_abs_ratio_diff = 0

    while num_iter_prod < 20000
        # Trade share equation
        pwmat = (w ./ (alpha .* A_i .* L .^ (alpha-1))) .^ (1-sigma)
        nummat = dd .* pwmat'
        tradesh = nummat ./ sum(nummat, dims = 1)
        # test = sum(tradesh, dims = 1)
        # mean(test) # Should be 1
        inc_exp_ratio = income ./ (tradesh * income)
        inc_exp_abs_ratio_diff = mean(abs.(inc_exp_ratio .- 1))
        if inc_exp_abs_ratio_diff < tol
            println(">>>> Productivity Convergence Achieved <<<<")
            aconverge = 1
            # break
        else
            if verbose
                print("  ARD(", num_iter_prod, "): ", round(inc_exp_abs_ratio_diff / (tol/100)) * tol)
            end
            A_e = A_i .* inc_exp_ratio .^ (1/(sigma-1))
            A_i = 0.25 .* A_e + 0.75 .* A_i
            # A_i /= geomean(A_i)
            num_iter_prod += 1
        end
    end

    println(">>>> Start amenities convergence <<<<")
    num_iter_pop = 1
    dtradesh = diag(tradesh)
    bconverge = 0
    pop_abs_ratio_diff = 0
    while num_iter_pop < 20000
        num = ((1 ./ B_i) .* dtradesh .^ (theta/(sigma-1)) .* A_i .^ (-theta) .* H .^ (theta-1)) .^ (1/(theta*alpha + beta))
        L_e = (num ./ sum(num)) .* L_tot # updated population estimate
        pop_ratio = L_i ./ L_e
        pop_abs_ratio_diff = mean(abs.(pop_ratio .- 1))
        if  pop_abs_ratio_diff < tol
            println(">>>> Population Convergence Achieved <<<<")
            bconverge = 1
            break
        else
            if verbose
                print("  ARD(", num_iter_pop, "): ", round(pop_abs_ratio_diff / (tol/100)) * tol)
            end
            B_e = B_i .* pop_ratio
            B_i = 0.25 .* B_e + 0.75 .* B_i
            B_i /= geomean(B_i)
            bconverge = 0
            num_iter_pop += 1
        end
    end

    xtic = time() - xtic
    return (A = A_i, B = B_i, tradesh = tradesh, dtradesh = dtradesh, 
            A_converged = aconverge, inc_exp_abs_ratio_diff = inc_exp_abs_ratio_diff, num_iter_prod = num_iter_prod,
            B_converged = bconverge, pop_abs_ratio_diff = pop_abs_ratio_diff, num_iter_pop = num_iter_pop, 
            tol = tol, xtic = xtic)
end
