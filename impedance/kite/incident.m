function [far,XYData_ref,utotal] = incident(k,theta)

FEM_mesh_size = 0.1;
R = 5;
global partition

model_ref = createpde(1); % Create a PDE model with a single equation
g = @circleRfunction;
geometryFromEdges(model_ref,g); % Specify geometry from boundary

a = @(location,state) (location.x.^2+location.y.^2<1).*(-(k^2).*(1+exp(-1./(1+location.x.^2+location.y.^2)))) + (location.x.^2+location.y.^2>=1).*(-k^2);
f = @(location,state) complex((location.x.^2+location.y.^2<1).*((k^2).*(exp(-1./(1+location.x.^2+location.y.^2))).*exp(1i*k*(location.x.*(-cos(theta))+location.y.*(-sin(theta))))));
specifyCoefficients(model_ref, 'm', 0, 'd', 0, 'c', 1, 'a', a, 'f', f); % specify coefficiens
applyBoundaryCondition(model_ref,'neumann','Edge',(1:4),'g',0,'q',-k*1i); % boundary condition approximated to Sommerfeld radiation condition
generateMesh(model_ref,'Hmax',FEM_mesh_size); % generate a mesh ;

% Solve the equation
result = solvepde(model_ref);
uscattered = result.NodalSolution;

% Express the incident field exp(1i*k*(x*cos(theta) + y*sin(theta))) by using the mesh
uincident = zeros(length(model_ref.Mesh.Nodes),1);
for n=1:length(model_ref.Mesh.Nodes)
    uincident(n,1) = exp(1i*k*((model_ref.Mesh.Nodes(1,n)*(-cos(theta)) + (model_ref.Mesh.Nodes(2,n)*(-sin(theta))))));
end
utotal = uincident + uscattered;
conj_utotal = conj(utotal);
model = createpde(1);
g1 = @circleobstacle;
geometryFromEdges(model,g1);
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', a, 'f', 0);
applyBoundaryCondition(model,'neumann','Edge',(1:4),'g',0,'q',-k*1i);

lambda = @(location,state) ((location.x)*((location.y)^2)<0).*(5 + 10*cos(((location.x)^2)*location.y)*1i) + ((location.x)*((location.y)^2)>=0).*(2*sin((location.x)^3)-5*1i);
XYData_ref = model_ref.Mesh.Nodes;
F = scatteredInterpolant(transpose(XYData_ref(1,:)),transpose(XYData_ref(2,:)),conj_utotal,'natural');
trace_utotal = @(location,state) F(location.x,location.y);
applyBoundaryCondition(model,'neumann','Edge',(5:8),'g',trace_utotal,'q',lambda);

generateMesh(model,'Hmax',FEM_mesh_size); %generate a mesh 

% Solve the equation
result2 = solvepde(model);
uscattered2 = result2.NodalSolution;

% % Plot the incident plane field
% figure
% subplot(2,2,1)
% pdeplot(model_ref,'XYData',real(uincident),'Mesh','off');
% colormap(jet)
% xlabel 'x'
% ylabel 'y'
% title 'Re(incident plane wave)'
% axis equal
% 
% % Plot the solution (scattered field)
% subplot(2,2,2)
% pdeplot(model_ref,'XYData',real(uscattered),'Mesh','off');
% colormap(jet)
% xlabel 'x'
% ylabel 'y'
% title 'Re(scattered field)'
% axis equal
% 
% % Plot the total field
% subplot(2,2,3)
% pdeplot(model_ref,'XYData',real(utotal),'Mesh','off');
% colormap(jet)
% xlabel 'x'
% ylabel 'y'
% title 'Re(total field) (total field will serve as incident field to detect obstacle)'
% axis equal
% 
% % Plot the solution (scattered field with obstacle)
% subplot(2,2,4)
% pdeplot(model,'XYData',real(uscattered2),'Mesh','off');
% colormap(jet)
% xlabel 'x'
% ylabel 'y'
% title 'Re(scattered field), with the computed incident field'
% axis equal

% compute far-field (estimate by the values at |x| = 3 = 0.75*R)
Rfar = 0.75*R;

XYData = model.Mesh.Nodes;
G = scatteredInterpolant(transpose(XYData(1,:)),transpose(XYData(2,:)),uscattered2,'natural');
farcoef = exp(1i*pi/4 + 1i*k*Rfar)/(sqrt(8*pi*k*Rfar));
far = zeros(1,partition);
for count=1:partition
    angle = 2*pi*count/partition;
    xfar = Rfar*cos(angle);
    yfar = Rfar*sin(angle);
    far(1,count) = G(xfar,yfar)./farcoef;
end

% -------------------- Nested functions --------------------
% no obstacle
    function [x,y] = circleRfunction(bs,s) 
    if nargin == 0  
        x = 8; % 8 segments
        return 
    end
    if nargin == 1
        % Outer circle with radius R
        dl1 = [0      pi/2   pi       3*pi/2
              pi/2    pi     3*pi/2   2*pi
              1       1      1        1 % region label to left (anticlockwise)
              0       0      0        0]; % region label to right (anticlockwise)
        % Inner circular obstacle with radius 3/4
        dl2 = [0      pi/2   pi       3*pi/2
               pi/2   pi     3*pi/2   2*pi
               2      2      2        2 % region label to left (anticlockwise)
               1      1      1        1]; % region label to right (anticlockwise)
        % Combine the two edge matrices
        dl = [dl1,dl2];
        x = dl(:,bs);   
        return 
    end 
    x = zeros(size(s)); 
    y = zeros(size(s)); 
    if numel(bs) == 1 % Does bs need scalar expansion?
        bs = bs*ones(size(s)); % Expand bs
    end
    cbs = find(bs <= 4); % Outer circle with radius R
    x(cbs) = R*cos(s(cbs));
    y(cbs) = R*sin(s(cbs));
    cbs = find(bs > 4); % Kite shaped
    x(cbs) = 0.5*(cos(s(cbs)) + 0.65*cos(2*(s(cbs))) - 0.65);
    y(cbs) = 0.5*1.5*sin(s(cbs));
    end

% with obstacle
    function [x,y] = circleobstacle(bs,s) 
    if nargin == 0  
        x = 8; % 8 segments
        return 
    end
    if nargin == 1
        % Outer circle with radius R
        dl1 = [0      pi/2   pi       3*pi/2
              pi/2    pi     3*pi/2   2*pi
              1       1      1        1   % region label to left (anticlockwise)
              0       0      0        0]; % region label to right (anticlockwise)
        % Inner circular obstacle with radius 1/2
        dl2 = [0      pi/2   pi       3*pi/2
               pi/2   pi     3*pi/2   2*pi
               0      0      0        0   % region label to left (anticlockwise)
               1      1      1        1]; % region label to right (anticlockwise)
        % Combine the two edge matrices
        dl = [dl1,dl2];
        x = dl(:,bs);   
        return 
    end 
    x = zeros(size(s)); 
    y = zeros(size(s)); 
    if numel(bs) == 1 % Does bs need scalar expansion?
        bs = bs*ones(size(s)); % Expand bs
    end
    cbs = find(bs <= 4); % Outer circle with radius R
    x(cbs) = R*cos(s(cbs));
    y(cbs) = R*sin(s(cbs));
    cbs = find(bs > 4); % Kite shaped
    x(cbs) = 0.5*(cos(s(cbs)) + 0.65*cos(2*(s(cbs))) - 0.65);
    y(cbs) = 0.5*1.5*sin(s(cbs));
    end

end