function [shortestLength, shortestRoute] = acotsp(cityPosition, isPlot, ...
    numAnt, maxIter, alpha, beta, rho, q)
% acotsp use ACO(Ant Colony Algorithm) to solve TSP(Travelling Salesman Problem)
%	[SHORTESTLENGTH, SHORTESTROUTE] = ACOTSP(CITYPOSITION, ISPLOT, ...
%   NUMANT, MAXITER, ALPHA, BETA, RHO, Q)
%
%	Solve TSP(Travelling Salesman Problem) using ACO(Ant Colony Algorithm)
%
%   Inputs
%		CITYPOSITION : cities to visit
%		ISPLOT : plot or not
%   	NUMANT : number of the ant
%   	MAXITER : max iteration
%		ALPHA : a parameter to control the influence of pheromone deposited(tau(i, j))
%		BETA : a parameter to control the influence of eta(i, j)
%		RHO : pheromone evaporation cofficient
%		Q :	constant used to update pheromone
%   Outputs
%   	SHORTESTLENGTH : shortest length found
%   	SHORTESTROUTE : shortest Route found

%   Author:		Yan
%   Email:		myncepu@gmail.com
%   References:	WenZheng, Proficient in MATLAB Intelligent Algorithm
%				https://en.wikipedia.org/wiki/Ant_colony_optimization_algorithms

%	Important parameters
%		tau(i, j) : the amount of pheromone deposited for a state transition i-j
%		eta(i, j) : the desirability of state transition i-j, typically 1/distance(i,j)
%% default arguments
if nargin < 8
    q = 1e6;
    if nargin < 7
        rho = 0.15;
        if nargin < 6
            beta = 2.2;
            if nargin < 5
                alpha = 1.4;
                if nargin < 4
                    maxIter = 100;
                    if nargin < 3
                        numAnt = 50;
                        if nargin < 2
                            isPlot = 1;
                            if nargin < 1
                                x = [41 37 54 25 7 2 68 71 54 83 64 18 22 83 ...
                                    91 25 24 58 71 74 87 18 13 82 62 58 45 41 4 ...
                                    4 4]';
                                y = [94 84 67 62 64 99 58 44 62 69 60 54 60 4 ...
                                    6 38 38 42 69 71 78 76 40 40 7 32 35 21 26 ...
                                    35 50]';
                                cityPosition = [x, y];
                            end
                        end
                    end
                end
            end
        end
    end
end

%% calculate distance between cities
numCity = size(cityPosition, 1);
cityDistance = sqrt((bsxfun(@minus, cityPosition(:, 1)', cityPosition(:, 1))).^2 + ...
    (bsxfun(@minus, cityPosition(:, 2)', cityPosition(:, 2))).^2);
cityDistance(find(cityDistance) == 0) = eps;

%%
eta = 1 ./ cityDistance;
tau = ones(numCity, numCity);
antPosition = zeros(numAnt, numCity);
R_best = zeros(maxIter, numCity);
L_best = inf .* ones(maxIter, 1);
L_ave = zeros(maxIter, 1);

%%
for iter = 1 : maxIter
    initPosition = [];
    for i = 1:ceil(numAnt/numCity)
        initPosition = [initPosition, randperm(numCity)];
    end
    % init position of ants
    antPosition(:, 1) = initPosition(1:numAnt)';
    for iCity = 2 : numCity
        for iAnt = 1 : numAnt
            visitedCity = antPosition(iAnt, 1 : (iCity - 1));
            toVisitCity = setdiff(1 : numCity, visitedCity);
            % the iter_th ant moves for city i to city j with probability
            probability = zeros(1, (numCity - iCity + 1));
            for k = 1:length(toVisitCity)
                probability(k) = (tau(visitedCity(end), toVisitCity(k)) ^ alpha) * ...
                    (eta(visitedCity(end), toVisitCity(k)) ^ beta);
            end
            probability = probability / sum(probability);
            Pcum = cumsum(probability);
            select = find(Pcum >= rand);
            toVisit = toVisitCity(select(1));
            antPosition(iAnt, iCity) = toVisit;
        end
    end
    if iter >= 2
        antPosition(1, :) = R_best(iter - 1, :);
    end
    
    L = zeros(numAnt, 1);
    for iAnt = 1:numAnt
        R = antPosition(iAnt, :);
        for iCity = 1:(numCity - 1)
            L(iAnt) = L(iAnt) + cityDistance(R(iCity), R(iCity + 1));
        end
        L(iAnt) = L(iAnt) + cityDistance(R(1), R(numCity));
    end
    [L_best(iter), pos] = min(L);
    R_best(iter, :) = antPosition(pos(1), :);
    L_ave(iter) = mean(L);
    
	% update pheromone 
    deltaTau = zeros(numCity, numCity);
    for iAnt = 1 : numAnt
        for iCity = 1 : (numCity - 1)
            deltaTau(antPosition(iAnt, iCity), antPosition(iAnt, iCity + 1)) = ...
                deltaTau(antPosition(iAnt, iCity), antPosition(iAnt, iCity + 1)) + q / L(iAnt);
        end
        deltaTau(antPosition(iAnt, numCity), antPosition(iAnt, 1)) = ...
            deltaTau(antPosition(iAnt, numCity), antPosition(iAnt, 1)) + q / L(iAnt);
    end
    tau = (1 - rho) .* tau + deltaTau;
	
    antPosition = zeros(numAnt, numCity);
end
[shortestLength, idx] = min(L_best);
shortestRoute = R_best(idx, :);

%%
if isPlot
    figure;
    subplot(1, 2, 1);
    drawRoute(cityPosition, shortestRoute);
    subplot(1, 2, 2);
    plot(L_best);
    hold on;
    plot(L_ave, 'r');
    legend('Shortest length', 'Average length');
    title('Average length and shortest length');
end
end

%%
function drawRoute(cityPosition, shortestRoute)
N = length(shortestRoute);
scatter(cityPosition(:, 1), cityPosition(:, 2));
hold on;
plot([cityPosition(shortestRoute(1), 1), cityPosition(shortestRoute(N), 1)], ...
    [cityPosition(shortestRoute(1), 2), cityPosition(shortestRoute(N), 2)], 'g');
hold on;
for ii = 2:N
    plot([cityPosition(shortestRoute(ii - 1), 1), cityPosition(shortestRoute(ii), 1)], ...
        [cityPosition(shortestRoute(ii - 1), 2), cityPosition(shortestRoute(ii), 2)], 'g');
    hold on
end
title('Final result of TSP problem');
hold off
end
