%% Example of using function xmap to perform Convergent Cross Mapping
% Dan Mønster, February 2016

%% Plot convergence of the correlation coefficient as a function of L
N_data = 1000; % length of the time series to be generated
N_trans = 1000; % Number of generations to skip to avoid transient data

% Model parameterss
rx = 3.65;
ry = 3.77;
bxy = 0;
byx = 0.05;

% Parameters for phase space emdbedding
tau = 1; % Time delay
m = 2;   % Embedding dimension


Lmin = 6;
Lmax = N_data-(m-1)*tau; % Max possible Library size
NL = 100;
Lstep = round((Lmax-Lmin)/(NL-1));
LibLen = [Lmin:Lstep:Lmax]';
NL = numel(LibLen);
rho_X = zeros(NL,1);
rho_Y = zeros(NL,1);
p_X = zeros(Lmax-Lmin+1,1);
p_Y = zeros(Lmax-Lmin+1,1);

%
% Iterate for N_trans generations without collecting data (only
% X(1) and Y(1) are updated, so the last data point will be
% used as the new X(1) and Y(1)).
%
if (N_trans > 0)
    initial = rand(2,1);
    X0 = initial(1);
    Y0 = initial(2);
    [Xtrans, Ytrans] = coupled_logistic(X0, Y0, rx, ry, bxy, byx, N_trans);
end
% Iterate the coupled maps for N_data generations
X0 = Xtrans(end);
Y0 = Ytrans(end);
[X, Y] = coupled_logistic(X0, Y0, rx, ry, bxy, byx, N_data);
%
% Perform phase space embedding
%
MX = psembed(X,m,tau);
MY = psembed(Y,m,tau);
%
% Loop over Library length
%
for l=1:NL
    L=LibLen(l);
    %
    % M is the size of the ensemble, i.e. the number of
    % xmaps performed for each L. This only makes sense when the
    % sampling method in xmap is set to 'random'
    %
    M = round(Lmax/L);
    rho_XM = zeros(M,1);
    rho_YM = zeros(M,1);
    for r = 1:M
        [ X_MY, Y_MX, X1, Y1] = ...
            xmap( X, Y, MX, MY, m, tau, L,'random');
        rho_XM(r) = corr(X_MY,X1');
        rho_YM(r) = corr(Y_MX,Y1');
    end
    %
    % Take the mean of the correlation coefficients in the ensemble
    %
    rho_X(l) = mean(rho_XM);
    rho_Y(l) = mean(rho_YM);
end



%% Plot
figure(4)
plot(LibLen,rho_X,'ob')
hold on
plot(LibLen,rho_Y,'or')
xlabel('$L$','interpreter','latex')
ylabel('$\rho$','interpreter','latex')
