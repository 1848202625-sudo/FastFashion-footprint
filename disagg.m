clc
clear

load('FACTORSS3.mat')
load("production_factor.mat")
load("va_factor.mat")
load("G:\gtap\gtap\gtap2017.mat")

oldSec = 65;
newSec = 70;
nReg = 160;
nBrand = 6;
Nold = oldSec * nReg;
Nnew = newSec * nReg;

if size(production_factor,1) == 6
    production_factor = production_factor';
end

if size(FACTORSS3,1) == 6
    FACTORSS3 = FACTORSS3';
end

production_factor(production_factor < 0) = 0;
FACTORSS3(FACTORSS3 < 0) = 0;

XX1 = X';
a1 = Z ./ repmat(XX1, Nold, 1);
a1(~isfinite(a1)) = 0;

x_apparel = X(28:oldSec:end);

pro_factor = reshape(production_factor', nReg * nBrand, 1);
x_apparelcopy = repelem(x_apparel, nBrand);
x_apparel17 = x_apparelcopy .* pro_factor;

x3 = zeros(Nnew, 1);

for r = 1:nReg

    oldBase = (r-1) * oldSec;
    newBase = (r-1) * newSec;

    x3(newBase+1:newBase+27) = X(oldBase+1:oldBase+27);
    x3(newBase+28:newBase+33) = x_apparel17((r-1)*6+1:r*6);
    x3(newBase+34:newBase+70) = X(oldBase+29:oldBase+65);

end

z3 = zeros(Nnew, Nold);

for r = 1:nReg

    oldBase = (r-1) * oldSec;
    newBase = (r-1) * newSec;

    z3(newBase+1:newBase+27, :) = Z(oldBase+1:oldBase+27, :);
    z3(newBase+28:newBase+33, :) = production_factor(r,:)' * Z(oldBase+28, :);
    z3(newBase+34:newBase+70, :) = Z(oldBase+29:oldBase+65, :);

end

z5 = zeros(Nnew, Nnew);

for r = 1:nReg

    oldBase = (r-1) * oldSec;
    newBase = (r-1) * newSec;

    z5(:, newBase+1:newBase+27) = z3(:, oldBase+1:oldBase+27);
    z5(:, newBase+34:newBase+70) = z3(:, oldBase+29:oldBase+65);

end

a_apparel2 = zeros(Nnew, nReg * nBrand);

for c = 1:nReg

    oldCol = 28 + (c-1) * oldSec;
    cols = (c-1)*6 + 1 : c*6;

    a_col = zeros(Nnew, 1);

    for r = 1:nReg

        oldBase = (r-1) * oldSec;
        newBase = (r-1) * newSec;

        a_col(newBase+1:newBase+27) = a1(oldBase+1:oldBase+27, oldCol);
        a_col(newBase+28:newBase+33) = a1(oldBase+28, oldCol);
        a_col(newBase+34:newBase+70) = a1(oldBase+29:oldBase+65, oldCol);

    end

    a_apparel2(:, cols) = repmat(a_col, 1, 6);

end

p_nonfast = 95.47;
p_fastavg = 16.53;

p_off = 20.00;
p_hm = 9.85;
p_zara = 23.02;
p_gap = 16.07;
p_uniqlo = 15.15;

ratio_nf_fast = p_nonfast / p_fastavg;
share_nonfast = 1 / (ratio_nf_fast + 1);
share_fast = 1 - share_nonfast;

ratio_off_hm = p_off / p_hm;
ratio_zara_hm = p_zara / p_hm;
ratio_gap_hm = p_gap / p_hm;
ratio_uniqlo_hm = p_uniqlo / p_hm;

w_fast_raw = [1/ratio_off_hm, 1, 1/ratio_zara_hm, 1/ratio_gap_hm, 1/ratio_uniqlo_hm];
w_fast = w_fast_raw / sum(w_fast_raw);

w6 = zeros(6,1);
w6(1) = share_nonfast;
w6(2) = share_fast * w_fast(1);
w6(3) = share_fast * w_fast(2);
w6(4) = share_fast * w_fast(3);
w6(5) = share_fast * w_fast(4);
w6(6) = share_fast * w_fast(5);

x_temp = reshape(x_apparel17, 6, nReg);
s_temp = zeros(6, nReg);

for r = 1:nReg

    if sum(x_temp(:,r)) == 0
        s_temp(:,r) = ones(6,1) / 6;
    else
        s_temp(:,r) = x_temp(:,r) / sum(x_temp(:,r));
    end

end

k_temp = zeros(6, nReg);

for r = 1:nReg

    denom = sum(s_temp(:,r) .* w6);

    if denom == 0
        k_temp(:,r) = ones(6,1);
    else
        k_temp(:,r) = w6 / denom;
    end

end

for c = 1:nReg

    rows = 28 + (c-1)*newSec : 33 + (c-1)*newSec;
    cols = (c-1)*6 + 1 : c*6;

    for b = 1:6
        a_apparel2(rows, cols(b)) = a_apparel2(rows, cols(b)) * k_temp(b,c);
    end

end

x5 = repmat(x_apparel17', Nnew, 1);
z_apparel = a_apparel2 .* x5;

z6 = z5;

for c = 1:nReg

    newCols = 28 + (c-1)*newSec : 33 + (c-1)*newSec;
    brandCols = (c-1)*6 + 1 : c*6;

    z6(:, newCols) = z_apparel(:, brandCols);

end

va2 = sum(VA, 1);
va4 = zeros(1, Nnew);

va_factor_use = va_factor(:)';

for r = 1:nReg

    oldBase = (r-1) * oldSec;
    newBase = (r-1) * newSec;

    va4(newBase+1:newBase+27) = va2(oldBase+1:oldBase+27);
    va4(newBase+34:newBase+70) = va2(oldBase+29:oldBase+65);

    va_original_apparel = va2(oldBase+28);
    x_brand = x_apparel17((r-1)*6+1:r*6)';

    va_fast = x_brand(2:6) .* va_factor_use(2:6);
    va_nonfast = va_original_apparel - sum(va_fast);

    va4(newBase+28:newBase+33) = [va_nonfast, va_fast];

end

colcon = x3' - va4;

P = production_factor;
S = FACTORSS3;

for r = 1:nReg

    if sum(P(r,:)) == 0
        P(r,:) = [1 0 0 0 0 0];
    else
        P(r,:) = P(r,:) / sum(P(r,:));
    end

end

for s = 1:nReg

    if sum(S(s,:)) == 0
        S(s,:) = [1 0 0 0 0 0];
    else
        S(s,:) = S(s,:) / sum(S(s,:));
    end

end

f2 = zeros(Nnew, 480);

for r = 1:nReg

    oldBase = (r-1) * oldSec;
    newBase = (r-1) * newSec;

    f2(newBase+1:newBase+27, :) = Y(oldBase+1:oldBase+27, :);
    f2(newBase+34:newBase+70, :) = Y(oldBase+29:oldBase+65, :);

end

f_apparel2 = zeros(nReg * nBrand, 480);

tol = 1e-12;
max_iter = 10000;
conv_tol = 1e-10;

for s = 1:nReg

    cols_s = (s-1)*3 + 1 : s*3;
    S_s = S(s,:);

    for t = 1:3

        col = cols_s(t);

        y_old = zeros(nReg,1);

        for r = 1:nReg
            old_row = 28 + (r-1)*oldSec;
            y_old(r) = Y(old_row, col);
        end

        total_y = sum(y_old);

        if total_y <= tol
            continue
        end

        brand_target = total_y * S_s;

        active = y_old > tol;
        active_id = find(active);

        if isempty(active_id)
            continue
        end

        row_target = y_old(active);
        P_active = P(active,:);

        feasible_brand = sum(P_active > tol, 1) > 0;
        infeasible_brand = (~feasible_brand) & (brand_target > tol);

        if any(infeasible_brand)

            adjust_mass = sum(brand_target(infeasible_brand));
            brand_target(infeasible_brand) = 0;

            if feasible_brand(1)
                brand_target(1) = brand_target(1) + adjust_mass;
            else
                feasible_list = find(feasible_brand);

                if isempty(feasible_list)
                    P_active(:,1) = 1;
                    brand_target(1) = brand_target(1) + adjust_mass;
                else
                    brand_target(feasible_list(1)) = brand_target(feasible_list(1)) + adjust_mass;
                end
            end

        end

        if sum(brand_target) > 0
            brand_target = brand_target / sum(brand_target) * total_y;
        else
            brand_target = [total_y 0 0 0 0 0];
            P_active(:,1) = 1;
        end

        M = P_active .* repmat(row_target, 1, nBrand);

        for b = 1:nBrand
            if brand_target(b) <= tol
                M(:,b) = 0;
            end
        end

        for b = 1:nBrand

            if brand_target(b) > tol && sum(M(:,b)) <= tol

                idx_feasible = P_active(:,b) > tol;

                if any(idx_feasible)
                    M(idx_feasible,b) = row_target(idx_feasible);
                end

            end

        end

        for iter = 1:max_iter

            row_sum = sum(M,2);

            for rr = 1:length(row_sum)

                if row_sum(rr) > tol
                    M(rr,:) = M(rr,:) * (row_target(rr) / row_sum(rr));
                elseif row_target(rr) > tol

                    feasible_cols = find(P_active(rr,:) > tol);

                    if isempty(feasible_cols)
                        feasible_cols = 1;
                    end

                    weights = S_s(feasible_cols);

                    if sum(weights) <= tol
                        weights = ones(size(feasible_cols)) / numel(feasible_cols);
                    else
                        weights = weights / sum(weights);
                    end

                    M(rr, feasible_cols) = row_target(rr) .* weights;

                end

            end

            col_sum = sum(M,1);

            for b = 1:nBrand

                if brand_target(b) > tol

                    if col_sum(b) > tol
                        M(:,b) = M(:,b) * (brand_target(b) / col_sum(b));
                    else
                        idx_feasible = P_active(:,b) > tol;

                        if any(idx_feasible)
                            M(idx_feasible,b) = brand_target(b) / sum(idx_feasible);
                        end
                    end

                else
                    M(:,b) = 0;
                end

            end

            err_row = max(abs(sum(M,2) - row_target));
            err_col = max(abs(sum(M,1) - brand_target));

            if max(err_row, err_col) < conv_tol
                break
            end

        end

        M(abs(M) < tol) = 0;

        row_sum = sum(M,2);

        for rr = 1:length(row_sum)

            if row_sum(rr) > tol
                M(rr,:) = M(rr,:) * (row_target(rr) / row_sum(rr));
            elseif row_target(rr) > tol
                M(rr,1) = row_target(rr);
            end

        end

        for ii = 1:length(active_id)

            r = active_id(ii);

            for b = 1:nBrand
                new_row = (r-1)*nBrand + b;
                f_apparel2(new_row, col) = M(ii,b);
            end

        end

    end

end

f3 = f2;

for r = 1:nReg

    rows_new = 28 + (r-1)*newSec : 33 + (r-1)*newSec;
    rows_brand = (r-1)*nBrand + 1 : r*nBrand;

    f3(rows_new, :) = f_apparel2(rows_brand, :);

end

tol_capacity = 1e-9;

for r = 1:nReg

    rows_apparel = 28 + (r-1)*newSec : 33 + (r-1)*newSec;

    fd_sum = sum(f3(rows_apparel,:), 2);
    x_rows = x3(rows_apparel);
    surplus = fd_sum - x_rows;

    while any(surplus > tol_capacity)

        [excess_amount, bad_local] = max(surplus);
        bad_row = rows_apparel(bad_local);

        available = x_rows - fd_sum;
        available(bad_local) = 0;

        p_row = P(r,:)';
        available(p_row <= tol_capacity) = 0;

        receiver_local = find(available > tol_capacity);
        total_available = sum(available(receiver_local));

        if total_available < excess_amount - tol_capacity
            error('第 %d 个国家服装部门内部可用容量不足。', r)
        end

        donor_total = sum(f3(bad_row,:));

        if donor_total <= 0
            error('第 %d 个国家第 %d 个服装子部门最终需求为0。', r, bad_local)
        end

        transfer_by_col = f3(bad_row,:) .* (excess_amount / donor_total);

        f3(bad_row,:) = f3(bad_row,:) - transfer_by_col;

        receiver_weight = available(receiver_local) ./ total_available;

        for k = 1:length(receiver_local)

            rec_row = rows_apparel(receiver_local(k));
            f3(rec_row,:) = f3(rec_row,:) + transfer_by_col .* receiver_weight(k);

        end

        fd_sum = sum(f3(rows_apparel,:), 2);
        surplus = fd_sum - x_rows;

    end

end

f3(abs(f3) < tol_capacity) = 0;

f4 = sum(f3, 2);
rowcon = x3 - f4;

u = rowcon;
v = colcon;
A = z6;

if isrow(u)
    u = u';
end

if iscolumn(v)
    v = v';
end

eps_gras = 1e-6;
halt_gras = 10000;

Z_balanced = GRAS(u, v, A, eps_gras, halt_gras);

X = Z_balanced;

save('2017-new-d.mat', ...
    'X', 'Z_balanced', 'x3', 'va4', 'f3', ...
    'rowcon', 'colcon', ...
    'P', 'S', 'va_factor_use', ...
    '-v7.3');