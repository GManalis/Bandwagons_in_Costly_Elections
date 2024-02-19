% The present script generates Figures 5,6 and 7 in Section 5.3 of the 
% paper "Bandwagons in costly elections: The role of loss aversion", 2023
% by A. Leontiou, G. Manalis and D. Xefteris. 
clear all 
close all 
clc 
options = optimoptions('fsolve','Display','iter','Algorithm','trust-region', 'FinDiffType', 'central', 'FunctionTolerance', 1.0000e-3)
%% Parametrization
N = 20 % Population of fixed size

% Figure 5.1
eta = 0.5; 
lambda = 1.5; 
x0  = [0.01;0.99];

% Figure 5.2
% Uncomment 3 lines below for Fig 5.2
% eta = 0.5; 
% lambda = 1.9; 
% x0 = [0.01;0.99];

% Figure 6.1
% Uncomment 3 lines below for Fig 6.1
% eta = 1;
% lambda = 1.5;
% x0 = [0.01;0.99];

% Figure 6.2
% Uncomment 3 lines below for Fig 6.2
% eta = 1;
% lambda = 1.9;
% x0  = [0.1;0.9];

% Figure 7.1
% Uncomment 3 lines below for Fig 7.1
% eta = 1;
% lambda = 1.9;
% x0 = [0.01;0.99];


% Figure 7.3
% Uncomment 3 lines below for Fig 7.2
% eta = 1;
% lambda = 1.9;
% x0 = [0.99;0.01];

%% Main Code:
fi_vec = linspace(0.05,0.5,10);
NA = unique(round(N*fi_vec));
NB = N-NA; 
syms a b j k

% Party A
%Initialize the matrices where the sym expressions are stored: 
eq1A = cell(size(NA));
eq2A = cell(size(NA));
eq3A = cell(size(NA));
eq4A = cell(size(NA));
pA = cell(size(NA));
qA = cell(size(NA));
pA_minus_qA = cell(size(NA));
pA_plus_qA = cell(size(NA));


%Below I express symbolically in terms of alpha and beta the equations in
%the draft. Equation 1,2, 3,4

for i = 1:numel(NA)
    eq1A{i} = symsum(nchoosek(NA(i)-1,k)*a^k*(1-a)^(NA(i)-k-1)*symsum(nchoosek(NB(i),k-j)*b^(k-j)*(1-b)^(NB(i)-k+j), j,0,k),k,0,NA(i)-1);
    eq2A{i} = symsum(nchoosek(NA(i)-1,k)*nchoosek(NB(i),k+1)*a^k*(1-a)^(NA(i)-1-k)*b^(k+1)*(1-b)^(NB(i)-k-1),k,0,NA(i)-1);
    eq3A{i} = symsum(nchoosek(NA(i)-1,k)*a^k*(1-a)^(NA(i)-1-k)*symsum(nchoosek(NB(i),k-j)*b^(k-j)*(1-b)^(NB(i)-k+j),j,1,k),k,1,NA(i)-1);
    eq4A{i} = symsum(nchoosek(NA(i)-1,k)*nchoosek(NB(i),k)*a^k*(1-a)^(NA(i)-1-k)*b^k*(1-b)^(NB(i)-k),k,0,NA(i)-1);
%Then I can form pA = 1*eq1 + (1/2)*eq2 and qA = 1*eq3 + (1/2)*eq4
    pA{i} = 1*eq1A{i} + (1/2)*eq2A{i};
    qA{i} = 1*eq3A{i} + (1/2)*eq4A{i};
    pA_minus_qA{i} = pA{i}-qA{i}; 
    pA_plus_qA{i} = pA{i}+qA{i};
end

% Doing the same for PARTY B

% Initialize the matrices where the sym expressions are stored: 
eq1B = cell(size(NA));
eq2B = cell(size(NA));
eq3B = cell(size(NA));
eq4B = cell(size(NA));
pB = cell(size(NA));
qB = cell(size(NA));
pB_minus_qB = cell(size(NA));
pB_plus_qB = cell(size(NA));

for i = 1:numel(NA)
    eq1B{i} = symsum(nchoosek(NB(i)-1,k)*b^k*(1-b)^(NB(i)-k-1)*symsum(nchoosek(NA(i),k-j)*a^(k-j)*(1-a)^(NA(i)-k+j), j,0,k),k,0,NB(i)-1);
    eq2B{i} = symsum(nchoosek(NB(i)-1,k)*nchoosek(NA(i),k+1)*b^k*(1-b)^(NB(i)-1-k)*a^(k+1)*(1-a)^(NA(i)-k-1),k,0,NB(i)-1);
    eq3B{i} = symsum(nchoosek(NB(i)-1,k)*b^k*(1-b)^(NB(i)-1-k)*symsum(nchoosek(NA(i),k-j)*a^(k-j)*(1-a)^(NA(i)-k+j),j,1,k),k,1,NB(i)-1);
    eq4B{i} = symsum(nchoosek(NB(i)-1,k)*nchoosek(NA(i),k)*b^k*(1-b)^(NB(i)-1-k)*a^k*(1-a)^(NA(i)-k),k,0,NB(i)-1);
%Then I can form pA = 1*eq1 + (1/2)*eq2 and qA = 1*eq3 + (1/2)*eq4
    pB{i} = 1*eq1B{i} + (1/2)*eq2B{i};
    qB{i} = 1*eq3B{i} + (1/2)*eq4B{i};
    pB_minus_qB{i} = pB{i}-qB{i}; 
    pB_plus_qB{i} = pB{i}+qB{i};
end

%Now I store equation 1 and equation 2 of the system into a cell array for
%different levels of Φ. 
one = cell(size(NA))
two = cell(size(NA))
sol = zeros(2,numel(NA))
% x0 = [0.01;0.99]
for i = 1:numel(NA)
    one{i} = matlabFunction(pA_minus_qA{i}*(1-eta*(lambda-1)*(1-pA_plus_qA{i}))-a)
    two{i} = matlabFunction(pB_minus_qB{i}*(1-eta*(lambda-1)*(1-pB_plus_qB{i}))-b)
    sys = @(x)[one{i}(x(1),x(2));two{i}(x(1),x(2))]
    sol(:,i) = fsolve(@(x)sys(x), x0, options);
    x0 = sol(:,i); % To change the initial values with the sol of the previous one
%     if NA(i) > 4 
%         x0 = [0.99;0.01]
%     end
end


fi_vec = NA./N

a_star = sol(1,:)
b_star = sol(2,:)
share_astar = fi_vec.*a_star
share_bstar = (1-fi_vec).*b_star

%% Plot
txt = ['$\eta =$ ' num2str(eta) ' , $\lambda=$ ' num2str(lambda) ' and $N =$ ' num2str(N)];
f = figure
subplot(1,3,[1 2])
 plot(fi_vec, a_star,'-o','Linewidth', 1.5, 'color', 'b', 'DisplayName', ['\alpha^*']);
    title('Equilibrium $\alpha^*$ and $\beta^*$', 'interpreter','latex')
    subtitle(txt,'interpreter','latex')
    xlabel('$\phi$', 'interpreter', 'latex')
    hold on 
    grid on
    plot(fi_vec, b_star,'--*','Linewidth', 1.5, 'color', 'r', 'DisplayName', ['\beta^*'])
    ylabel('$\alpha^*, \beta^*$', 'interpreter', 'latex')
    legend('show', 'Location', 'southeast')
hold on
 subplot(1,3,3)
 plot(fi_vec, share_astar,'-o','Linewidth', 1.5, 'color', 'b', 'DisplayName', ['S_A']);
    title('$S_A = \phi\cdot\alpha^*$ and $S_B = (1-\phi)\cdot\beta^*$', 'interpreter','latex')
    xlabel('$\phi$', 'interpreter', 'latex')
    hold on 
    grid on 
    plot(fi_vec, share_bstar,'--*','Linewidth', 1.5, 'color', 'r', 'DisplayName', ['S_B'])
    ylabel('$S_\alpha^*, S_\beta^*$', 'interpreter', 'latex')
    legend('show', 'Location', 'southeast')

f.Position = [100 100 700 550];
fig_name = strcat('LASTeta', num2str(eta), 'lambda', num2str(lambda), 'N', num2str(N),'.png');
exportgraphics(f, fig_name)