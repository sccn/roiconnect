function [mic, mim] = roi_mim2_v2(Cohroi, subinds1, subinds2)

    regu = 0.000001; 
    [~, ~, nfreq] = size(Cohroi);
    
    ipcs = numel(subinds1);
    jpcs = numel(subinds2);
    
    cs_red1 = Cohroi(subinds1,subinds1, :);
    cs_red2 = Cohroi(subinds1,subinds2, :);
    cs_red3 = Cohroi(subinds2,subinds2, :);
    
    [mean_diag_r_csred1, dim1] = get_mean_diag(cs_red1, nfreq);
    caainv = pageinv(real(cs_red1) + regu * ones(dim1, dim1, nfreq) .* eye(ipcs) .* mean_diag_r_csred1);
    cab = imag(cs_red2);
    [mean_diag_r_csred3, dim1] = get_mean_diag(cs_red3, nfreq);
    cbbinv = pageinv(real(cs_red3) + regu * ones(dim1, dim1, nfreq) .* eye(jpcs) .* mean(mean_diag_r_csred3));
    X = pagemtimes(pagemtimes(cab, cbbinv), pagetranspose(cab));

    % MIM Ewald Eq. 14
    r_caainv_X = reshape(pagemtimes(caainv, X), [], nfreq);
    diag_idx = 1:dim1+1:size(r_caainv_X, 1);
    trace_operator = zeros(1, size(r_caainv_X, 1));
    trace_operator(diag_idx) = 1;
    mim = trace_operator * r_caainv_X;
    
    caainvsqrt = zeros(size(caainv, 1), size(caainv, 2), size(caainv, 3));
    for ifq = 1:nfreq
        caainvsqrt(:, :, ifq) = sqrtm(caainv(:, :, ifq));
    end
    Y = pagemtimes(pagemtimes(caainvsqrt, X), caainvsqrt); % Eq. 23
    [~, s, ~] = pagesvd(Y);

    % MIC
    mic = sqrt(s(1, 1, :));
end

function [mean_diag_r_csred, dim1] = get_mean_diag(cs_red, nfreq)
    % equivalent to diag(real(cs_red) in 2D
    dim1 = size(cs_red, 1);
    r_csred = reshape(real(cs_red), [], nfreq);
    diag_idx = 1:dim1+1:size(r_csred, 1);
    diag_r_csred = r_csred(diag_idx, :); 
    mean_diag_r_csred = reshape(mean(diag_r_csred, 1), 1, 1, []);
end