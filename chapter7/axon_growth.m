% Author: Marcus Kaiser    Date: 19 October 2016

%% initialisation

cone_range = 0:20:60; % axon direction can vary by +- cone_range for
                 % each growth step
                 
trials = 10000;     % number of trials for each configuration
neuron_density = 0.05; % 5 percent of all unit spaces are filled
side_length = 31; % number of unit spaces on the side of embedding square
centre_neuron = ceil(side_length / 2)
steps = 100     % number of growth steps


% xy position of growth cone after each step
C = zeros(length(cone_range), trials, steps, 2);

% matrix of neuron position unit spaces
A = zeros(length(cone_range), trials, side_length, side_length);

% number of steps before hitting a neuron 
S = zeros(length(cone_range), trials);

% variation in direction before hitting a neuron
D = zeros(length(cone_range), trials);


% straight growth
cc  = 0;

matrix = double(rand(side_length) <= neuron_density);
matrix(centre_neuron,centre_neuron) = 2; % growth always starts from the neuron at the centre 


for c = cone_range

    cc = cc + 1;
    tt = 0;

    for t = 1:trials
        tt = tt + 1;

        direction = 360 * rand(1);
        A(cc,tt,:,:) = matrix;
        i = 1;
        C(cc,tt,i,:) = [centre_neuron + 0.5, centre_neuron + 0.5]; % start in centre neuron
        directions = direction;
        
        while (i < steps) ... 
                && (floor(C(cc,tt,i,1)) > 3) && ceil(C(cc,tt,i,1)) < side_length-3 ...
                && (floor(C(cc,tt,i,2)) > 3) && ceil(C(cc,tt,i,2)) < side_length-3 
                
            i = i + 1;
            C(cc,tt,i,1) = C(cc,tt,i-1,1) + 0.5 * cosd(direction); 
            C(cc,tt,i,2) = C(cc,tt,i-1,2) + 0.5 * sind(direction);
            if matrix( floor(C(cc,tt,i,1)), floor(C(cc,tt,i,2)) ) ~= 0 ...
                    && floor(C(cc,tt,i,1)) ~= centre_neuron ...
                    && floor(C(cc,tt,i,2)) ~= centre_neuron,
                S(cc,tt) = i;
                % A(cc,tt, round(C(cc,tt,i,1)), round(C(cc,tt,i,2)) ) = 3;
                % C(cc,t,i,1) = C(cc,t,i,1) + 0.5 * cosd(direction);
                % C(cc,t,i,2) = C(cc,t,i,2) + 0.5 * sind(direction);
                break
            end
            direction = direction + 2 * c * rand(1) - c; % new direction
            directions = [directions direction];
        end;
        if S(cc,tt) > 2
            % D(cc,tt) = std(directions);
            D(cc,tt) = mean(abs(directions(2:end) - directions(1:end-1)));
        else
            S(cc,tt) = 0;
            tt = tt - 1;
        end;
    end
end
    
nnz(S)

%%
% visualisation
figure(1);
clf;
set(gcf,'Color','w')

c = colormap(lines(length(cone_range)));
colormap(parula);

subplot(2,2,1); % examples for straight growth
pcolor(squeeze(A(1,1,:,:))');
colormap('gray');
hold on;

for i = 1:length(cone_range)
    % M = [zeros(1,S(i,1)); eye(S(i,1)-1,S(i,1))];
    % gplot(M, squeeze(C(i,1,1:S(i,1),:)));    
    xy = squeeze(C(i,1,1:S(i,1),:));
    plot(xy(:,1),xy(:,2),'Color',c(i,:));
    hold on;
end;


subplot(2,2,2); % summary plot

% plot(cone_range,mean(S'));
Sm = mean(S');
hold on;
%colormap(parula);
for i = 1:length(cone_range)
    % plot(cone_range(i), Sm(i), 'o');
    bar(cone_range(i), Sm(i), 5, 'FaceColor',c(i,:));
    hold on;
end;
    
xlabel('Maximum cone direction change [deg]');
ylabel('Axon length [mm]');
axis([-5 max(cone_range)+5 0 max(Sm)])
set(gcf,'Color','w')


subplot(2,2,4); % summary plot

for i = 1:length(cone_range)
    [y,x] = hist(S(i,:),4);
    plot(x,y,'Color',c(i,:));
    hold on;
end;

ylabel('Number of occurrences');
xlabel('Axon length [mm]');


% plot(mean(D'),mean(S'));
% xlabel('Average cone direction change [deg]');
% ylabel('Axon length [mm]');


% subplot(3,2,5); % summary plot
% hist(S(1,:),5)
% ylabel('Number of occurrences');
% xlabel('Axon length [mm]');
% 
% subplot(3,2,6); % summary plot
% hist(S(length(cone_range),:),5)
% ylabel('Number of occurrences');
% xlabel('Axon length [mm]');

[cone_range; mean(S')]
