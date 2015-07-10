clc;
clear;

rand('seed', 1);
// begin initialization

function p = population(n, dim, lb, ub)
    p = rand(n,dim)
    p = lb  + p .* (ub - lb);
endfunction

function zi = z_ideal(pop_fit)
    zi = [];
    s = size(pop_fit)
    for i=1:s(2)
        zi(i) = min(pop_fit(:,i));
    end
    zi = zi';
endfunction

function y = f(x, n)
    s = size(x);
    gm = sum((x(n:s(2))- 0.5).^2);
    y = [];
    c = [1];
    s = [1];
    for i=1:n-1
        c = [c , c(i) * cos(0.5*x(i)*%pi)];
        s = [s , sin(0.5*x(n - i)*%pi)];
    end
    for i=1:n
        y(i) = (1+gm)*s(i)*c(n+1-i);
    end
    y = y';
endfunction

function fit = calc_fit(pop, my_f, f_num)
    fit = [];
    s = size(pop)
    for i=1:s(1)
        fit = [fit; my_f(pop(i,:), f_num)];
    end
endfunction

function m = lambda_vecs(f_num)
    m = eye(f_num, f_num);
    for i=1:f_num
        for j=1:f_num
            m(i,j) = max(m(i,j), 10e-6);
        end
    end
endfunction

function clv = class_vec(pop_fit, zi)
    s = size(pop_fit);
    clv = hypermat([s(1), s(2), s(2)])
    for i=1:s(1)//pop size
        for j=1:s(2)// n func
            if(pop_fit(i,j) == zi(j))//workaround for divsion by 0
                den = 1;
            else
                den = pop_fit(i,j) - zi(j);
            end
            for k=1:s(2)
                if j == k then
                    clv(i,j,k) = 10e-6;
                else
                    clv(i,j,k) = (pop_fit(i,k) - zi(k)) / (den);
                end
            end
        end
    end
endfunction
//end initialization

//begin classification

function c = classify(clv)
    c = [];
    s = size(clv)
    for i=1:s(1)
        j_min = 1;
        k_max_o = -1;
        for j=1:s(2)// starting from 1 to check the first k_max
            k_max_i = 1;
            for k=2:s(3)
                if(clv(i,j,k_max_i) < clv(i,j,k))
                    k_max_i = k;
                end
            end
            if j ~= 1 then
                if clv(i,j_min,k_max_o) > clv(i,j,k_max_i) then
                    j_min = j;
                    k_max_o = k_max_i;
                end
            else
                k_max_o = k_max_i //initializes k_max_o
            end
        end
        c(i,1) = j_min;
        c(i,2) = clv(i,j_min,k_max_o);
    end
endfunction

// balance 

function new_p = balance(c, f_num)
    new_p = c;
    count = zeros(f_num, 1);
    s = size(new_p);
    sp = ceil(s(1) / f_num);
    sm = floor(s(1) / f_num);
    for i=1:s(1)
        count(new_p(i)) = count(new_p(i)) + 1
    end
    lp = list();
    lm = list();
    for i=1:f_num
        if count(i) > sp then
            lp($+1) = i;
        end
        if count(i) < sm then
            lm($+1) = i;
        end
    end
    
    while size(lp) > 0 & size(lm) > 0
        for i=1:s(1)
            if new_p(i,1) == lp(1)
                max_i = i;
            end
        end
        for i=max_i+1:s(1)
            if new_p(i,1) == lp(1) & new_p(i,2) > new_p(max_i,2);
                max_i = i;
            end
        end
        new_p(max_i, 1) = lm(1);
        count(lm(1)) = count(lm(1)) + 1;
        count(lp(1)) = count(lp(1)) - 1;
        if(count(lm(1)) >= sm) 
            lm(1) = null();
        end
        if(count(lp(1)) <= sp)
            lp(1) = null();
        end
    end
    new_p = new_p(:,1);
endfunction

function l_p = split_pop(pop, c, f_num)
    l_p = list();
    for i=1:f_num
        l_p($+1) = list();
    end
    s = size(pop);
    for i=1:s(1)
        l_p(c(i))($+1) = pop(i,:);
    end
endfunction

//end classification


//begin crossover-mutation

function b = cross_fact(cross_rate)
    r = rand();
    b = (2*r)^(1/(cross_rate +1));
    if b > 1 then
        b = 1 / ((2*(1-r))^(1/(cross_rate +1)));
    end
endfunction

function d = mut_fact(mut_rate)
    r = rand();
    if r < 0.5 then
        d = (2*r)^(1 / (mut_rate + 1)) - 1;
    else
        d = 1 - (2*(1-r))^(1/(mut_rate + 1));
    end
endfunction

//Simmulated binary crossover
function [ch_p, ch_m] = crossover(p1, p2, cross_rate, lb, ub)
    ch_p = [];
    ch_m = [];
    for i=1:length(p1)
        b = cross_fact(cross_rate);
        ch_p(i) = 0.5 * ((1-b)*p1(i) + (1+b)*p2(i));
        ch_m(i) = 0.5 * ((1+b)*p1(i) + (1-b)*p2(i));
        while(ch_p(i) > ub) 
            ch_p(i) = ch_p(i) - (ub - lb);
        end
        while(ch_m(i) > ub) 
            ch_m(i) = ch_m(i) - (ub - lb);
        end
        while(ch_p(i) < lb) 
            ch_p(i) = ch_p(i) + (ub - lb);
        end
        while(ch_m(i) < lb) 
            ch_m(i) = ch_m(i) + (ub - lb);
        end
    end
endfunction

//Polynomial mutation
function mut = mutation(ch, mut_rate, p1, p2, lb, ub)
    mut = [];
    for i=1:length(ch)
        d = mut_fact(mut_rate);
        mut(i) = ch(i) + (max(p1(i), p2(i)) - min(p1(i), p2(i))) * d;
        while(mut(i) > ub) 
            mut(i) = mut(i) - (ub - lb);
        end
        while(mut(i) < lb) 
            mut(i) = mut(i) + (ub - lb);
        end
    end
endfunction
 
//end crossover-mutation

//begin evolution

function l_p2 = gen_offsprings(l_p, cross_rate, mut_rate, lb, ub)
    l_p2 = l_p;
    for i=1:length(l_p)
        s = length(l_p(i))
        for j=1:ceil(s/2)
            idxs = grand(1,2,'uin',1,s);
            [c1, c2] = crossover(l_p2(i)(idxs(1)), l_p2(i)(idxs(2)), cross_rate, lb, ub);
            c1 = mutation(c1, mut_rate, l_p2(i)(idxs(1)), l_p2(i)(idxs(2)), lb, ub);
            c2 = mutation(c2, mut_rate, l_p2(i)(idxs(1)), l_p2(i)(idxs(2)), lb, ub);
            l_p2(i)($+1) = c1';
            l_p2(i)($+1) = c2';
        end
        if modulo(length(l_p2(i)), 2) == 1 then
            l_p2(i)($) = null();
        end
    end
endfunction

function l_p2 = gen_offsprings2(l_p, cross_rate, mut_rate, lb, ub)
    l_p2 = l_p;
    for i=1:ceil(length(l_p)/2)
            idxs = grand(1,2,'uin',1,length(l_p));
            [c1, c2] = crossover(l_p2(idxs(1)), l_p2(idxs(2)), cross_rate, lb, ub);
            c1 = mutation(c1, mut_rate, l_p2(idxs(1)), l_p2(idxs(2)), lb, ub);
            c2 = mutation(c2, mut_rate, l_p2(idxs(1)), l_p2(idxs(2)), lb, ub);
            l_p2($+1) = c1';
            l_p2($+1) = c2';
    end
    if modulo(length(l_p2), 2) == 1 then
            l_p2($) = null();
    end
endfunction

function [asf, fit] = calc_asf(l_p2, my_f, f_num)
    asf = list();
    fit = list();
    lbd = lambda_vecs(length(l_p2));
    for i=1:length(l_p2)
        s_fit = [];
        for j=1:length(l_p2(i))
            s_fit = [s_fit ;calc_fit(l_p2(i)(j), my_f, f_num)];
        end
        fit(i) = s_fit;
    end
    
    zi = z_ideal(merge_list_2(fit));
    for i=1:length(l_p2)
        t_asf = [];
        for j=1:length(l_p2(i))
            t_asf(j) =  max((fit(i)(j, :) - zi)./ lbd(i, :));
        end
        asf(i) = t_asf;
    end
endfunction

function [s_lp, s_asf, s_fit] = sort_pop_asf_fit(lp, asf, fit)
    s_lp = lp;
    s_asf = asf;
    s_fit = fit;
    for i=1:length(s_lp)
        for j=1:length(s_lp(i))
            for k=j+1:length(s_lp(i))
                if s_asf(i)(j) > s_asf(i)(k) then
                    temp = s_asf(i)(j);
                    s_asf(i)(j) = s_asf(i)(k);
                    s_asf(i)(k) = temp;
                    
                    temp = s_lp(i)(j);
                    s_lp(i)(j) = s_lp(i)(k);
                    s_lp(i)(k) = temp;
                    
                    temp = s_fit(i)(j);
                    s_fit(i)(j) = s_fit(i)(k);
                    s_fit(i)(k) = temp;
                end
            end
        end
    end
endfunction

function [l_p3, asf2, fit2] = halve_pop_asf_fit(pop, asf, fit)
    l_p3 = pop;
    asf2 = asf;
    fit2 = fit;
    for i=1:length(l_p3)
        s = length(l_p3(i))/2;
        while(length(l_p3(i)) > s) 
            l_p3(i)($) = null();
        end;
        asf2(i) = asf2(i)(1:s, :);
        fit2(i) = fit2(i)(1:s, :)
    end
endfunction


//end evolution

//begin reclassification

function merged = merge_list_1(l)
    merged = [];
    for i=1:length(l)
        for j=1:length(l(i))
            merged = [merged; l(i)(j)];
        end
    end
endfunction

function merged = merge_list_2(l)
    merged = [];
    for i=1:length(l)
        merged = [merged; l(i)];
    end
endfunction

function new_l_p = reclassify(l_pop, l_fit, f_num)
    pop = merge_list_1(l_pop);
    fit = merge_list_2(l_fit);
    zi = z_ideal(fit);
    clv = class_vec(fit, zi);
    c = classify(clv);
    new_c = balance(c, f_num);
    new_l_p = split_pop(pop, new_c, f_num);
endfunction

//end reclassification

//begin target points

function [tp, asf] = target_points(l_pop, my_f, f_num)
    s = length(l_pop);
    //disp(l_pop);
    pop = merge_list_1(l_pop);
    fit = calc_fit(pop, my_f, f_num);
    zi = z_ideal(fit);
    lbd = lambda_vecs(s);
    p_s = size(pop);
    asfs = zeros(p_s(1), s);
    for i=1:p_s(1)
        for j=1:s
            asfs(i,j) = max((fit(i, :) - zi)./ lbd(j, :));
        end
    end
    idxs = zeros(s,1);
    tp = [];
    asf = [];
    for i=1:s
        idxs(i) = 1;
        for j=2:p_s(1)
            if asfs(j,i) < asfs(idxs(i), i) then
                idxs(i) = j;
            end
        end
        tp = [tp; pop(idxs(i), :)];
        asf = [asf; asfs(idxs(i), i)];
    end
endfunction

//end target points

//begin reduction objective space

function [tp, asf_] = reduce_space(pop_size, pop_dim, f_num, lb, ub, max_it, cross_rate, mut_rate, f_obj)
    pop = population(pop_size, pop_dim, lb, ub);
    fit = calc_fit(pop, f_obj, f_num);
    zi = z_ideal(fit);
    clv = class_vec(fit, zi);
    c = classify(clv);
    new_c = balance(c, f_num);
    pop = split_pop(pop, new_c, f_num);
    tp = [];
    asf_ = [];
    for i=1:max_it
        pop = gen_offsprings(pop, cross_rate, mut_rate, lb, ub);
        [asf, fit] = calc_asf(pop, f_obj, f_num);
        [pop, asf, fit] = sort_pop_asf_fit(pop, asf, fit);
        [pop, asf, fit] =  halve_pop_asf_fit(pop, asf, fit);
        pop = reclassify(pop, fit, f_num);
    end
    [tmp1, tmp2] = target_points(pop, f_obj, f_num);
    tp = [tp; tmp1];
    asf_ = [asf_; tmp2];
endfunction

//end reduction objective space

// 2nd phase

//begin diversity op

function f_m = f_max(p_fit)
    f_m = [];
    s = size(p_fit)
    for i=1:s(2)
        f_m(i) = max(p_fit(:,i));
    end
    f_m = f_m';
endfunction

function grid = calc_grid(p_fit, final_size)
    s = size(p_fit);
    grid = []
    zi = z_ideal(p_fit);
    f_m = f_max(p_fit);
    den = f_m - zi;
    for i=1:length(den)
        if(den(i) == 0) then
           den(i) = 10e-6
        end
    end
    for i=1:s(1)
        grid = [grid; final_size * ((p_fit(i,:) - zi)./ den)]
    end
endfunction

function A = calc_dist_grid(grid)
    s = size(grid)
    A = zeros(s(1), s(1))
    for i=1:s(1)-1
        for j=i+1:s(1)
            A(i,j) = norm(grid(i,:)- grid(j,:));
            A(j,i) = A(i,j);
        end
    end
endfunction

function sorted = sort_by_dist(A_dist)
    s = size(A_dist);
    mark = zeros(s(1), 1);
    sorted = zeros(s(1),1);
    for ind=1:s(1)-1
        idx = 0;
        min_ = 10e20;
        for i=1:s(1)
            if mark(i) then
                continue;
            end
            for j=1:s(1)
                if(i == j | mark(j))
                    continue;
                end
                if A_dist(i,j) < min_ then
                    idx = i;
                    min_ = A_dist(i,j);
                end
            end 
        end
        mark(idx) = 1;
        sorted(s(1) - ind + 1) = idx;
    end
    for i=1:s(1)
        if ~mark(i) then
            sorted(1) = i;
            break;
        end
    end
endfunction

function [r_pop, r_fit] = reduce_pop(pop, l_fit, num)
    r_pop = list();
    r_fit = list();
    fit = merge_list_2(l_fit);
    grid = calc_grid(fit, num);
    A = calc_dist_grid(grid);
    sorted = sort_by_dist(A);
    for i=1:num
        r_pop($+1) = pop(sorted(i));
        r_fit($+1) = l_fit(sorted(i));
    end
endfunction

//end diversity op

//begin initialization

function sp = sub_population(tps, num, lb, ub, ra)
    sp = list();
    s = size(tps);
    i = 0;
    while(i < num)
        for j=1:s(1)
             p = (2*ra*rand(1, s(2)) - ra) + tps(j,:);
             for k=1:length(p)
                 while(p(k) > ub) 
                    p(k) = p(k) - (ub - lb);
                end
                while(p(k) < lb) 
                    p(k) = p(k) + (ub - lb);
                end
             end
             sp($+1) = p;
            i = i + 1;
            if i >= num then
                break;
            end
        end
    end
endfunction

function sb = space_bounds(tps, my_f)
    sb = [];
    s = size(tps);
    f_tps = [];
    for i=1:s(1)
        f_tps = [f_tps; calc_fit(tps(i,:), my_f, s(1))];
    end
    for i=1:s(1)
        sb(i) = max(f_tps(:,i));
    end
    sb = sb';
endfunction

//end initialization

//begin selection

function [cb, l_fit] = classify_bounds(l_pop, sb, my_f)
    l_fit = list();
    n_f = length(sb);
    for i=1:length(l_pop)
        l_fit(i) = calc_fit(l_pop(i), my_f, n_f);
    end
    
    cb = list();
    for i=1:length(l_fit)
        cb(i) = 1;
        for k=1:n_f
           if(l_fit(i)(k) > sb(k)) 
              cb(i) = 0;
              break;
           end
        end
    end
endfunction

function s_in = calc_s_in(cb)
    s_in = 0;
    for i=1:length(cb)
        if(cb(i) ~= 0)
            s_in = s_in + 1;
        end
    end
endfunction

function d = dominate(f_p1, f_p2)
    d = 1;
    for i=1:length(f_p1)
        if f_p1(i) > f_p2(i) then
            d = 0;
            break;
        end
    end
endfunction

function fronts = nd_sort(p_fit)
    dom_set = list();
    fronts = list();
    fronts(1) =  list();
    dom_count = zeros(1, length(p_fit));
    for i=1:length(p_fit)
        s_p = [];
        for j=1:length(p_fit)
            if i ~= j then
                if dominate(p_fit(i), p_fit(j)) then
                    s_p = [s_p, j];
                elseif dominate(p_fit(j), p_fit(i)) then
                    dom_count(i) = dom_count(i) + 1;
                end
            end
        end
        dom_set(i) = s_p;
        if dom_count(i) == 0 then
            if length(fronts) > 0 then
                fronts(1)($+1) = i;
            end
        end
    end
    f_idx = 1;
    keep = 1;
    while(keep) 
        nxt_front = list();
        fi_len = length(fronts(f_idx));
        for i=1:fi_len
            dj_len = length(dom_set(fronts(f_idx)(i)))
            for j=1:dj_len
                dom_count(dom_set(fronts(f_idx)(i))(j)) = dom_count(dom_set(fronts(f_idx)(i))(j)) - 1;
                if dom_count(dom_set(fronts(f_idx)(i))(j)) == 0 then
                    nxt_front($+1) = dom_set(fronts(f_idx)(i))(j);
                end
            end
        end
        if length(nxt_front) then
            fronts(f_idx + 1) = nxt_front;
            f_idx = f_idx + 1;
        else
            keep = 0;
        end
       
    end
endfunction

function snd = calc_s_nd(non_dominated, cb)
    snd = 0;
    for i=1:length(non_dominated)
        if cb(non_dominated(i)) then
            snd = snd + 1;
        end
    end
endfunction

function idxs = pop_dist(pop1, pop2)
    s_n = length(pop1);
    s_d = length(pop2);
    dist = zeros(s_n, 1);
    val = zeros(s_d, 1);
    for i=1:s_d
        for j=1:s_n
            dist(j) = norm(pop2(i) - pop1(j));
        end
        val(i) = min(dist);
    end
    [val, idxs] = gsort(val, 'g', 'd');
endfunction




function [s_pop, s_fit] = halve_pop_by_diversity(l_pop, tps, sb, my_f, num)
    [cb, l_fit] = classify_bounds(l_pop, sb, my_f);
    s_in = calc_s_in(cb);
    s_fit = list();

    if s_in > num then

        fronts = nd_sort(l_fit);
        s_nd = calc_s_nd(fronts(1), cb);
        if s_nd > num then
            [s_pop, s_fit] = reduce_pop(l_pop, l_fit, num);
        else
            s_pop = list();
            non_d = fronts(1);
            for i=1:length(non_d)
                if cb(non_d(i)) then
                    s_pop($+1) = l_pop(non_d(i));
                    s_fit($+1) = l_fit(non_d(i));
                end
            end
            
            dom = list();
            for i=2:length(fronts)
                for j=1:length(fronts(i))
                    if cb(fronts(i)(j)) then
                        dom($+1) = l_pop(fronts(i)(j));
                    end
                end
            end
            len = length(s_pop);
            idxs = pop_dist(s_pop, dom);
            i = 1;
            while len < num
                s_pop($+1) = l_pop(idxs(i));
                s_fit($+1) = l_fit(idxs(i));
                len = len + 1;
                i = i + 1;
            end
        end
    else
        s_pop = list();
        out_b = list();
        out_idx = list();
        for i=1:length(cb)
           if cb(i) then
              s_pop($+1) = l_pop(i);
              s_fit($+1) = l_fit(i)
           else
               out_b($+1) = l_pop(i);
               out_idx($+1) = i;
           end
        end
        s = size(tps);
        tp_mid = list();
        for i=1:s(1)
            tp_mid($+1) = tps(i,:);
        end
        for i=1:s(1)-1
            for j=i+1:s(1)
                tp_mid($+1) = (tps(i,:) + tps(j,:)) ./ 2;
            end
        end
        idxs = pop_dist(tp_mid, out_b);
        len = length(s_pop);
        i = length(out_b);
        while len < num
            s_pop($+1) = out_b(i);
            s_fit($+1) = l_fit(out_idx(i));
            len = len + 1;
            i = i - 1;
        end
    end
endfunction

//end selection

//begin solve

function [non_d, pop, fit, tp] = solve_reduced(my_f, p_size, tps, rg, lb, ub, cr, mr, mit)
    pop = sub_population(tps, p_size, lb, ub, rg);
    
    s_bounds = space_bounds(tps, my_f);
    s = size(tps);
    for i=1:mit
        sp1 = size(pop);
        pop = gen_offsprings2(pop, cr, mr, lb, ub);
        [pop, fit] = halve_pop_by_diversity(pop, tps, s_bounds, my_f, p_size);
        //classify
        fit_m = merge_list_2(fit);
        zi = z_ideal(fit_m);
        clv = class_vec(fit_m, zi);
        c = classify(clv);
        new_c = balance(c, s(1));
        pop2 = merge_list_2(pop);
        sp2 = size(pop);
        pop2 = split_pop(pop2, new_c, s(1));
        [tps, asf] = target_points(pop2, my_f, s(1));
        s_bounds = space_bounds(tp, my_f);
    end
    fronts = nd_sort(fit);
    non_d = fronts(1);
    tp = tps;
    disp("fim");
endfunction

//end solve




// run test
pop_size = 30;
num_func = 3;
num_dim = 3;
upper_bounds = 1;
lower_bounds = 0;
cross_rate = 1;
mut_rate = 1 / pop_size;
max_it = 30;
tic();
[tp, asf] = reduce_space(pop_size, num_dim, num_func, lower_bounds, upper_bounds, max_it, cross_rate, mut_rate, f)
s = toc();
disp(s, "tempo reduzindo espaco");
range_f = 0.5;
tic();
[non_d, pop, fit] = solve_reduced(f, pop_size, tp, range_f, lower_bounds, upper_bounds, cross_rate, mut_rate, max_it);
s = toc();
disp(s, "tempo resolvendo");

for i=1:length(fit)
    param3d(fit(i)(1), fit(i)(2), fit(i)(3));
    e=gce();
    e.mark_mode = 'on';
    e.mark_size = 1;
end

