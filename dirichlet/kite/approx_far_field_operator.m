%main program
global partition
partition = 24;
k = 8.75;
gridsize = 0.1;
F = zeros(partition,partition);
% noise = 0.01; % noise percentage

angle = 2*pi/partition;
[far,XY,utotal] = incident(k,angle);
F(1,:) = far;
clearvars far
utotal_tensor = zeros(length(transpose(utotal)),partition);
utotal_tensor(:,1) = utotal;

for count=2:partition
    angle = 2*pi*count/partition;
    [far,~,utotal] = incident(k,angle);
    F(count,:) = far;
    clearvars far
    utotal_tensor(:,count) = utotal;
end

% % generate noise
% for count = 1:partition
%     for count1 = 1:partition
%         F(count,count1) = F(count,count1) + normrnd(0,abs(F(count,count1))*noise);
%     end
% end

[U,S,V] = svd(F);
aux_model = createpde(1);
g = @grid;
geometryFromEdges(aux_model,g);
generateMesh(aux_model, 'Hmax',gridsize);
points = aux_model.Mesh.Nodes;
r = zeros(partition,length(points));

for count=1:partition
    angle = 2*pi*count/partition;
    G = scatteredInterpolant(transpose(XY(1,:)),transpose(XY(2,:)),utotal_tensor(:,count),'natural');
    for num=1:length(points)
        r(count,num) = G(points(1,num), points(2,num));
    end
    clearvars G
end

P = (transpose(V))*r;
W0 = zeros(length(points),1);
for num=1:length(points)
    for count=1:partition
        W0(num,1) = W0(num,1) + ((abs(P(count,num)))^2)/S(count,count);
    end
end
W = 1./W0;

figure
pdeplot(aux_model,'XYData',W,'Mesh','off')
colormap(jet)
xlabel 'x'
ylabel 'y'
hold on

%plot the actual obstacle
t = linspace(0,2*pi,100);
x = 0.5*(cos(t) + 0.65*cos(2*t) - 0.65);
y = 0.5*1.5*sin(t);
plot(x,y,'color','k')

% -------------------- Nested functions --------------------
% no obstacle
    function [x,y] = grid(bs,s) 
    if nargin == 0  
        x = 4; % 8 segments
        return 
    end
    if nargin == 1
        % circle with radius 1
        dl = [0      pi/2   pi       3*pi/2
              pi/2    pi     3*pi/2   2*pi
              1       1      1        1 % region label to left (anticlockwise)
              0       0      0        0]; % region label to right (anticlockwise)
        x = dl(:,bs);   
        return 
    end 
    x = zeros(size(s)); 
    y = zeros(size(s)); 
    if numel(bs) == 1 % Does bs need scalar expansion?
        bs = bs*ones(size(s)); % Expand bs
    end
    cbs = find(bs <= 4); 
    x(cbs) = 2*cos(s(cbs));
    y(cbs) = 2*sin(s(cbs));
    end