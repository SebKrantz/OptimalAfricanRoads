#######################
# Helper Functions
#######################

function landprice(param, fund, res)
    return ((1-param.theta)/param.theta) .* ((res.w .* res.L) ./ fund.H)
end

function pindex(param, fund, res)
    return (1 ./ res.dtradesh) .^ (1/(1-param.sigma)) .* (res.w ./ (param.alpha .* fund.A .* res.L .^ (param.alpha-1)))
end

function consumption(param, fund, res)
    return (res.w .* res.L) ./ pindex(param, fund, res)
end

function output(param, fund, res)
    return fund.A .* res.L .^ param.alpha
end

function realw(param, fund, res)
    return res.w ./ (pindex(param, fund, res) .^ param.theta .* landprice(param, fund, res) .^ (1-param.theta))
end

function realinc(param, fund, res)
    return realw(param, fund, res) ./ param.theta
end

function utility(param, fund, res)
    return fund.B .* res.L .^ (1 + param.beta) .* realinc(param, fund, res)
end

function welfare(param, fund, res) # totut
    return sum(utility(param, fund, res))
end

function utility_per_worker(param, fund, res)
    return utility(param, fund, res) ./ res.L
end

function expectut(param, fund, res)
    return welfare(param, fund, res) / sum(res.L)
end

# This gives the same as welfare
function utility2(param, fund, L, dtradesh) # totut2
    variable = L .^ (param.beta + param.alpha*param.theta) .* dtradesh .^ (param.theta/(1-param.sigma))
    constant = fund.B .* (param.alpha .* fund.A ./ param.theta) .^ param.theta .* (fund.H ./ (1-param.theta)) .^ (1-param.theta)
    return variable .* constant
end

function welfare2(param, fund, L, dtradesh)
    return sum(utility2(param, fund, L, dtradesh))
end

function immobile_welfare_gain(param, fund, L_init, dtradesch_init, dtradesh_new)
    welfare_init = welfare2(param, fund, L_init, dtradesch_init)
    welfare_new = welfare2(param, fund, L_init, dtradesh_new)
    return welfare_new / welfare_init
end

function gains_from_trade(param, fund, L_trade, L_autarky, dtradesch)
    welfare_autarky = welfare2(param, fund, L_autarky, 1.0)
    welfare_trade = welfare2(param, fund, L_trade, dtradesch)
    return welfare_trade / welfare_autarky
end

function immobile_gains_from_trade(param, fund, dtradesch)
    welfare_autarky = welfare2(param, fund, 1.0, 1.0)
    welfare_trade = welfare2(param, fund, 1.0, dtradesch)
    return welfare_trade / welfare_autarky
end

function immobile_gains_from_trade_per_worker(param, fund, dtradesch)
    utility_autarky = utility2(param, fund, 1.0, 1.0)
    utility_trade = utility2(param, fund, 1.0, dtradesch)
    return utility_trade ./ utility_autarky
end

