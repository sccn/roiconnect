% ssgc_demo
%
% State space Granger causality demonstration script
%
% Generate random state space model and calculate time- and frequency-domain
% causalities and causal graphs. Please see the README and reference article for
% background.
%
%----------------- Simulation parameters ---------------------------------------

n       = 7;      % observation variable dimension
m       = 20;     % state space dimension

rhoa    = 0.85;   % transition matrix spectral norm (must be < 1)

icfac   = 10;     % innovations correlation factor: a positive integer (as icfac
                  % gets bigger, innovation correlations get smaller), or set to
                  % Inf for zero correlation ( =>  no instantaneous causality)

fres    = 1000;   % frequency resolution for spectral causalities

n1      = 2;      % target variable dimension
n2      = 3;      % source variable dimension

%-------------------------------------------------------------------------------

% Generate random state space parameters

[A,C,K] = iss_rand(n,m,rhoa,icfac);

B = A-K*C;
rhob = specnorm(B);
fprintf('\nrho(A) = %8.6f\n',rhoa);
fprintf('rho(B) = %8.6f\n\n',rhob);

% Generate random (positive-definite) covariance matrix

V = cov_rand(n,icfac);

% Choose random source and target

assert(n1+n2 <= n);
nperm = randperm(n);
i1 = nperm(1:n1);          % target multi-index
i2 = nperm(n1+1:n1+n2);    % source multi-index
i3 = setdiff(1:n,[i1 i2]); % multi-index of remaining (conditioning) variables

% Calculate time-domain GC

F = iss_GC(A,C,K,V,i1,i2);

% (Normalised) frequencies for frequency-domain GCs

omega = ((0:fres)/fres)*pi; % frequencies in range [0,pi]
z = exp(-1i*omega);         % on unit circle in complex plane

% Calculate frequency-domain GC

f = iss_SGC(A,C,K,V,z,i1,i2);

% Integrate frequency-domain GC up to Nyqvist frequency

Fint = trapz(f)/fres;

% Now we check Geweke's GC frequency decomposition: we should find Fint <= F,
% with equality under the condition that a certain spectral matrix is
% nonsingular on the unit disk in the complex plane. See J. Geweke, J. Am. Stat.
% Assoc. 77(378), 1982 and J. Geweke, J. Am. Stat. Assoc. 79(388), 1984 for
% details. (Geweke is rather optimistic that the equality condition is, in his
% experience with empirical data, "usually" met - here you may find that it
% fails occasionally for randomly generated state space models.)

Fdiff = F-Fint;
fprintf('F    = %f\n',F);
if abs(Fdiff) < 1e-8
    fprintf('Fint = %f : Geweke frequency decomposition PASS (|F-Fint| = %g)\n\n',Fint,abs(Fdiff));
else
    fprintf('Fint = %f : Geweke frequency decomposition FAIL (F-Fint = %g)\n\n',Fint,Fdiff);
end

% Plot frequency-domain GCs

figure(1); clf
plot(omega,f);
xlim([0 pi]);
title(sprintf('Causality [%s] -> [%s] conditional on [%s]:\n',num2str(i2),num2str(i1),num2str(i3)));
xlabel('normalised frequency');
set(gca,'XTick',[0,pi/2,pi]);
set(gca,'XTickLabel',{'0','pi/2','pi'});
ylabel('spectral causality');

% Calculate time-domain pairwise-conditional GCs ("causal graph")

G = iss_PWGC(A,C,K,V);

% Plot time-domain causal graph

figure(2); clf
colormap(flipud(bone));
maxG = max(G(:));
imagesc(G,[0 maxG]);
title('Causal graph');
axis('square');
xlabel('source');
ylabel('target');
set(gca,'XTick',1:n);
set(gca,'XTickLabel',1:n);
set(gca,'YTick',1:n);
set(gca,'YTickLabel',1:n);

% Calculate frequency-domain pairwise-conditional GCs

g = iss_SPWGC(A,C,K,V,z);

% Plot frequency-domain causal graph

figure(3); clf
ylims = [min(g(:)) 1.1*max(g(:))];
k = 0;
for i = 1:n
    for j = 1:n
        k = k+1;
        if i ~= j
            subplot(n,n,k);
            plot(omega',squeeze(g(i,j,:)));
            axis('square');
            xlim([0 pi]);
            ylim(ylims);
            xlabel('frequency')
            ylabel(sprintf('%d -> %d',j,i));
            set(gca,'XTick',[0,pi/2,pi]);
            set(gca,'XTickLabel',{'0','pi/2','pi'});
       end
    end
end
