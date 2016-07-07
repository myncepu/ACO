function maximum = simpleACO(objectFunction, isPlot, numAnt, maxIter, rou, ...
    probability, step, upperBound, lowerBound)
% simpleACO use ACO(Ant Colony  Algorithm) to find the maximum of a function
%   MAXIMUM = SIMPLEACO(OBJECTFUNCTION, ISPLOT, NUMANT, MAXITER, ROU, ...
%   PROBABILITY, STEP, UPPERBOUND, LOWERBOUND)
%
%   Find maximum of a function with two variables.
%
%   Inputs
%   	OBJECTFUN : handle of object function
%   	ISPLOT : plot or not
%   	NUMANT : number of the ant
%   	MAXITER : max iteration
%   	ROU : pheromone evaportion coefficient
%   	PROBABILITY : probability
%   	STEP : step
%   	UPPERBOUND : upper bounds
%   		upperbound = [xCordinateLowerBound, yCordinateUpperBound]
%   	LOWERBOUND : lower bounds
%		lowerbound = [xCordinateLowerBound, yCordinateUpperBound]
%   Outputs
%   	MAXIMUM : maximum value of the objectFun found

%   Author:		Yan
%   Email:		myncepu@gmail.com
%   References:	WenZheng, Proficient in MATLAB Intelligent Algorithm

%%----- initialize -----
%% initialize ant
position = zeros(numAnt, 2);
fitness = zeros(numAnt, 1);
for iAnt = 1:numAnt
    position(iAnt, :) = (upperBound - lowerBound) .* rand(1, 2) + lowerBound;
    fitness(iAnt) = objectFunction(position(iAnt, 1), position(iAnt, 2));
end

%% plot initial position
if isPlot
    [X, Y] = meshgrid(lowerBound(1):step:upperBound(1), ...
        lowerBound(2):step:upperBound(2));
    Z = objectFunction(X, Y);
    figure;
    subplot(1, 2, 1);
    mesh(X, Y, Z);
    hold on;
    plot3(position(:, 1), position(:, 2), fitness, 'k*');
    text(0.1, 0.8, -0.1, 'Initial position');
    xlabel('x'); ylabel('y'); zlabel('z');
	hold off;
end

%% iteration
for iter = 1:maxIter
	lambda = 1/iter;
	[maxValue(iter), maxIdx(iter)] = max(fitness);
	for iAnt = 1:numAnt
		prob(iter, iAnt) = (fitness(maxIdx(iter)) - fitness(iAnt)) / maxValue(iter);
		if prob(iter, iAnt) < probability
			temp = position(iter, :) + (2 * rand(1, 2) - 1) * lambda;
		else
			temp = position(iter, :) + (upperBound - lowerBound) * rand(1, 2) - 0.5;
		end
		temp = max(temp, lowerBound);
		temp = min(temp, upperBound);
		if objectFunction(temp(1), temp(2)) > objectFunction(position(1), position(2))
			positon = temp;
		end
% 		rou(iAnt) = (1 - rou) .* rou + objectFunction(position(1), position(2));
	end
end

%% plot final position
if isPlot
    [X, Y] = meshgrid(lowerBound(1):step:upperBound(1), ...
        lowerBound(2):step:upperBound(2));
    Z = objectFunction(X, Y);
    subplot(1, 2, 2);
    mesh(X, Y, Z);
    hold on;
    plot3(position(:, 1), position(:, 2), fitness, 'k*');
    text(0.1, 0.8, -0.1, 'Final position');
    xlabel('x'); ylabel('y'); zlabel('z');
	hold off;
end

%% output
maximum = max(fitness);
