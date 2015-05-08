clear
clc

// f(3.0, 2.0) = 0.0
// f(-2.805118, 3.131312) = 0.0
// f(-3.779310, -3.283186) = 0.0
// f(3.584428, -1.848126) = 0.0
function y = f(x)
    y = (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;
endfunction

function pop = initPop(N_0, dim, Lb, Ub)
    Lb_array = Lb .* ones(1, dim);
    Ub_array = Ub .* ones(1, dim);
    for i = 1 : N_0
        pop(i, :) = Lb_array + (Ub_array - Lb_array).*rand(1, dim, "uniform");
    end
endfunction

function seeds = reproduction(pop, s_min, s_max, dim)
    N_pop = length(pop) / dim;
    fitness = zeros(N_pop, 1);
    for i = 1 : N_pop
        fitness(i) = f(pop(i, :));
    end
    opposite_kathete = s_max - s_min;
    adjacent_kathete = max(fitness) - min(fitness);
    angular_coefficient = atan(opposite_kathete / adjacent_kathete);
    linear_coefficient = s_min - (angular_coefficient * min(fitness));
    seeds = zeros(N_pop, 1);
    for i = 1 : N_pop
        seeds(i) = round((angular_coefficient * fitness(i)) + linear_coefficient);
    end
    seeds = s_max - seeds;
endfunction

function new_plants = spatialDispersal(dim, sigma_initial, sigma_final, n, it_max, it, seeds, pop)
    N_pop = length(pop) / dim;
    mu = 0;
    sigma_it = ((((it_max - it).^n) / (it_max.^n)) * (sigma_initial - sigma_final)) + sigma_final;
    new_plants = zeros(sum(seeds), dim);
    k = 1;
    for i = 1 : N_pop
        if (seeds(i) > 0)
            r = grand(seeds(i), dim, "nor", mu, sigma_it);
            for j = 1 : seeds(i)
                if (rand(1, 1, "uniform") > 0.5)
                    new_plants(k, :) = pop(i, :) + r(j, :);
                else
                    new_plants(k, :) = pop(i, :) - r(j, :);
                end
                k = k + 1;
            end
        end
    end
endfunction

function new_pop = competitive_exclusion(pop, new_plants, dim, p_max)
    num_pop = (length(pop) / dim) + (length(new_plants) / dim);
    temp_pop = zeros(num_pop, dim);
    j = 1;
    for i = 1 : (length(pop) / dim)
        temp_pop(j, :) = pop(i, :);
        j = j + 1;
    end
    for i = 1 : (length(new_plants) / dim)
        temp_pop(j, :) = new_plants(i, :);
        j = j + 1;
    end
    fitness_temp_pop = zeros(num_pop, 1);
    for i = 1 : num_pop
        fitness_temp_pop(i) = f(temp_pop(i, :));
    end
    new_pop = zeros(p_max, dim);
    k = 0;
    for i = 1 : p_max
        min_temp_fitness = min(fitness_temp_pop);
        for j = 1 : (num_pop - k)
            if (fitness_temp_pop(j) == min_temp_fitness)
                new_pop(i, :) = temp_pop(j, :);
                ind = j;
            end
        end
        fitness_temp_pop(ind) = [];
        temp_pop(ind, :) = [];
        k = k + 1;
    end
endfunction

function [bestWeed, bestFitness] = getBestWeed(pop, dim)
    N_pop = length(pop) / dim;
    fitness = zeros(N_pop, 1);
    for i = 1 : N_pop
        fitness(i) = f(pop(i, :));
    end
    bestFitness = min(fitness);
    for i = 1 : N_pop
        if (fitness(i) == bestFitness)
            bestWeed = pop(i, :);
        end
    end
endfunction

function [bestWeed, bestFitness] = iwo()
    N_0 = 10;
    it_max = 100;
    dim = 2;
    p_max = 15;
    s_max = 5;
    s_min = 0;
    n = 3;
    sigma_initial = 1;
    sigma_final = 0.001;
    Lb = [-6, -6];
    Ub = [6, 6];
    pop = initPop(N_0, dim, Lb, Ub);
    for it = 1 : it_max
        seeds = reproduction(pop, s_min, s_max, dim);
        new_plants = spatialDispersal(dim, sigma_initial, sigma_final, n, it_max, it, seeds, pop);
        pop = competitive_exclusion(pop, new_plants, dim, p_max);
    end
    [bestWeed, bestFitness] = getBestWeed(pop, dim);
endfunction

function contour_plot(Lb, Ub)
    vx = linspace(Lb(1), Ub(1), 1000);
    vy = linspace(Lb(2), Ub(2), 1000);
    vf = [];
    matf = [];
    for i = 1 : length(vx)
        for j = 1 : length(vy)
            vf = [vf, f([vx(i); vy(j)])];
        end
        matf = [matf ; vf];
        vf = [];
    end
    contour2d(vx, vy, matf, 25);
    xlabel("$x$", "fontsize", 4);
    ylabel("$y$", "fontsize", 4);
    xset("fpf",string=" ");
    xset("font", 1, 3);
endfunction

[bestWeed, bestFitness] = iwo();
disp(bestWeed, "Best weed:");
disp(bestFitness, "Weed fitness:");
//Lb = [-6, -6];
//Ub = [6, 6];
//contour_plot(Lb, Ub);
