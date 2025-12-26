function [BeamSgnl, timeMethod, z2] = TD_L1_Red_ADMM_ConvFast(par, varargin)
%%name for saving results
par.name = strcat(par.name, "_SubImg=", num2str(par.SubImg));

par.Sgnl   = par.Sgnl./norm(reshape(par.Sgnl,[],1),2);
%%Observations fidelity
if nargin > 1 && ~isempty(varargin{1})
    Aty = varargin{1};
else
    Aty = A_transp_conv2D(par);
end
%%Default return
BeamSgnl = reshape(sqrt(mean(Aty.^2, 2)), par.NPixlsX, par.NPixlsZ)';
%%Precomputations
normAty  = norm(Aty(:),2);

%%Cross-validation
errorCross = zeros(length(par.tau_v), length(par.lambd_v), length(par.mu_v));
timeCross  = zeros(size(errorCross));
minError   = Inf;
timeMethod = 0;
tol_alg    = 1e-3;

for iTau = 1:length(par.tau_v)
    tau = par.tau_v(iTau);
    for cont_lambd = 1:length(par.lambd_v)
        lambd = par.lambd_v(cont_lambd)*tau;
        for iM = 1:length(par.mu_v)
            mu = par.mu_v(iM)*tau;
            tic;
            BB_old = Aty;
            z1     = BB_old;
            z2     = zeros(size(BB_old));
            u1     = zeros(size(z1));
            u2     = zeros(size(z2));

            Afun_handle = @(x) Afun(x, tau, par);
            iterative_error = zeros(par.NIter,1);

            for it = 1:par.NIter
                %%Proximal L1
                aux_z1        = BB_old + u1;
                lambd         = (par.lambd_v(cont_lambd)*tau)*(norm(reshape(aux_z1,[],1),2)/normAty);
                z1            = (sign(aux_z1)).*max(abs(aux_z1) - (lambd/tau),0);
                z1            = (z1./norm(z1(:)))*normAty;
                
                u1            = aux_z1 - z1;

                %%Proximal Red with BM4D
                aux_z2  = BB_old + u2;

                % tic;
                % medMat = mean(aux_z2, 2);              % media de cada variable (1 x 880)
                % X      = aux_z2 - medMat;
                % [U, S, V] = svd(X);
                % PCs = V;                      % 880 x 880
                % eigenvalues = (diag(S).^2) / (size(X,1) - 1);
                % explained_variance = eigenvalues / sum(eigenvalues);
                % cum_explained = cumsum(explained_variance);
                % scores = X * PCs;
                % my_k     = 65;%find(cum_explained >= 0.95, 1);
                % PCs_K    = PCs(:, 1:my_k);       % direcciones principales
                % scores_K = scores(:, 1:my_k);    % datos proyectados
                %
                % aux_cube = reshape(scores_K, par.NPixlsX, par.NPixlsZ, my_k);
                % aux_cube = BM4D(aux_cube, sqrt(mu/tau));
                % z2       = reshape(aux_cube, par.NPixls,my_k);
                % z2       = z2 * PCs_K' + medMat;
                % td = toc;

                tic;

                aux_cube = reshape(aux_z2, par.NPixlsX, par.NPixlsZ, par.NSampls);
                
                % wname = 'db2';
                % level = 1;
                % W = wavedec3(aux_cube, level, wname);
                % hf = [];
                % for i = 2:numel(W.dec)
                %     hf = [hf; W.dec{i}(:)]; %#ok<AGROW>
                % end
                % mu =  ((median(abs(hf)) / 0.6745)^2)*tau;

                mu = par.myMu(it);
                aux_cube = BM4D(aux_cube, sqrt(mu/tau));
                z2 = reshape(aux_cube, par.NPixls, par.NSampls);
                %z2 = (z2./norm(z2(:)))*norm(aux_z2(:));
                z2 = (z2./norm(z2(:)))*normAty;
                td2 = toc;
                %mu   = (par.mu_v(iM))*(norm(reshape(aux_z2,[],1),2)/normAty);

                 u2 = aux_z2 - z2;
                
                
                if mod(it,3)==0
                    sdfdsf = 1;
                end


                %%Update x with PCG
                rhs = Aty + 1*tau*(z1 - u1) + 1*tau*(z2 - u2);
                % tic;
                % [BB_new, ~] = pcg(Afun_handle, rhs(:), 1e-4, 100, [], [], BB_old(:));
                % toc
                %BB_new = reshape(BB_new, par.NPixls, par.NSampls);

                %%Update x with closed solution
                %tic;
                %BB_new2   = solve_tikhonov(par, rhs, tau);
                %toc

                tic;
                BB_new   = solve_tikhonov(par, rhs, tau);
                toc

                %%Check convergence
                iterative_error(it) = norm(BB_new - BB_old, 'fro') / norm(BB_old, 'fro');
                if it >= 2 && (iterative_error(it) < tol_alg || ~isfinite(iterative_error(it)))
                    break;
                end
                BB_old = BB_new;
            end

            taux         = toc;
            currentError = norm(BB_new(:),1);
            errorCross(iTau, cont_lambd, iM) = currentError;
            timeCross(iTau, cont_lambd, iM)  = taux;

            inRangeCur = currentError > par.thr_stp(1) && currentError < par.thr_stp(2);
            inRangeMin = minError     > par.thr_stp(1) && minError     < par.thr_stp(2);

            shouldUpdate = ...
                (inRangeCur && inRangeMin  && currentError < minError) || ...
                (inRangeCur && ~inRangeMin) || ...
                (~inRangeCur && ~inRangeMin && currentError < minError && currentError > 0);

            if shouldUpdate
                minError   = currentError;
                timeMethod = taux;
                BeamSgnl   = sqrt(mean(BB_new.^2, 2));
                save(strcat('Results/res_L1_Red_ADMM_', par.name, '.mat'), ...
                    'BB_new', 'minError', 'tau', 'lambd', 'mu', ...
                    'timeMethod', 'errorCross', 'timeCross', 'iterative_error');
            end
        end
    end
end

save(strcat('Results/res_L1_Red_ADMM_', par.name, '.mat'), 'errorCross', 'timeCross');
end


%% ==================== AUXILIARY FUNCTIONS ====================


function BB_est = solve_tikhonov(par, rhs, rho)
% Solve x_est = (K'K + 2*rho*I)^(-1)(K'y+rho*z1+rho*z2)
rhs     = reshape(rhs, 1, par.NPixls, par.NSampls);
RHS     = fft(rhs, [], 3); % 1 × NPixls × Nt
den_fft = sum(abs(par.phase).^2, 1);   % 1 × NPixls × Nt
BB_fft  = RHS ./ (den_fft + 2*rho);        % 1 × NPixls × Nt
BB_est  = squeeze(real(ifft(BB_fft, [], 3)));   % 1 × NPixls × Nt
end


function BB_est = solve_tikhonov_freq(par, rhs, rho)
% Solve x_est = (K'K + rho*I)^(-1)(K'y)
RHS_fft = fft(rhs, [], 2);               % Nsens × Nt
RHS_fft = permute(RHS_fft, [2 1]);       % Nt × Nsens
BB_fft = zeros(par.NSampls, par.NPixls);
for ff = 1:par.NSampls
    Kf = squeeze(par.phase(:, :,ff));  % Nsens × Npix
    KtK = Kf' * Kf;             % Npix × Npix
    Kty = RHS_fft(ff, :).';     % Npix × 1
    BB_fft(ff, :) = (KtK + 3*rho * eye(par.NPixls)) \ Kty;
end
BB_fft  = permute(BB_fft, [2 1]);        % Npix × Nt
BB_est  = real(ifft(BB_fft, [], 2));     % Npix × Nt
end



function z1_init = initz1(Aty, BeamSgnl, par)
aux_cube   = reshape(Aty, par.NPixlsX, par.NPixlsZ, par.NSampls);
max_vals   = squeeze(max(BeamSgnl, [], 2));
thr_col    = quantile(max_vals, par.thr_sparse(1));
[~, locs]  = findpeaks(max_vals, 'MinPeakHeight', thr_col);

offset     = 3;
intervals  = arrayfun(@(x) (x-offset):(x+offset), locs, 'UniformOutput', false);

z1_init    = zeros(size(aux_cube));

for k = 1:length(locs)
    row_vals     = BeamSgnl(locs(k), :);
    thr_row      = quantile(row_vals, par.thr_sparse(1));
    [~, locs_row]= findpeaks(row_vals, 'MinPeakHeight', thr_row);
    nCols = size(z1_init, 2); % número de columnas máximo
    cols = cell2mat(arrayfun(@(x) max(1, x-offset) : min(nCols, x+offset), locs_row, 'UniformOutput', false));
    cols = unique(cols);
    z1_init(intervals{k}, cols, :) = aux_cube(intervals{k}, cols, :);
end

% Normalize across selected rows
norms = arrayfun(@(i) norm(squeeze(z1_init(locs(i), :, :)), "fro"), 1:length(locs))';
max_norm = max(norms);
z1_init(locs, :, :) = (z1_init(locs, :, :) ./ norms) .* max_norm;

z1_init = z1_init ./ max(abs(z1_init(:)));
end


function A_PGDx = Afun(x_vec, tau, par)
par.BB   = reshape(x_vec, par.NPixls, par.NSampls);
par.Sgnl = A_direct_conv2D(par);
AtAx     = A_transp_conv2D(par);
A_PGDx   = AtAx(:) + 2*tau.*x_vec;
end
