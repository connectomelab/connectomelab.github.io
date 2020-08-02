% Program:     sng - spatial network growth
% Author:      Marcus Kaiser, International University Bremen
%              http://www.biological-networks.org/
% Date:        10.11.2002


% constants
NODES = 500;  % final number of nodes in the network
INISIZE = 1;  % number of initial nodes
TRIALS = 1;   % number of networks to be generated


% parameters
alpha = 50;   % distance dependence

beta = 0.001;    % 'density' parameter




for tr = 1 : TRIALS
  matrix = zeros(NODES,NODES);   % connectivity matrix (no distances!)
  position = zeros(NODES,2);    % (x,y) positions of the nodes
  distance = zeros(NODES,1);    % distances of new node to existing nodes



  % inimatrix
  position(1,:) = [0.5 0.5]; 


  n = INISIZE + 1;

  while n <= NODES
      position(n,:) = rand(1,2);    % random position for candidate node
      for i=1:n-1                   % distances to node n
          distance(i) = sqrt( (position(n,1)-position(i,1))^2 + (position(n,2)-position(i,2))^2 );        
          prob = beta * exp( -alpha * distance(i) );   % spatial contraint
          if rand(1) <= prob   % bidirectional edge will be established
              matrix(i,n) = 1;
              matrix(n,i) = 1;
          end; % if
      end; % for
      if sum(matrix(n,:))+sum(matrix(:,n)) > 0 
          % new node has at least established one connection
          n = n + 1;
      end; % if
      
  end; % while n


  % calculations of network properties
  % [...]


end; % trials





