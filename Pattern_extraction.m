function PP = Pattern_extraction(xy,NN,XY,R_data,coreind)
% PP = Pattern_extraction(xy,NN)
%
% It applies for extracting spatial patterns of regions of interest (ROIs)
%
% Inputs of the function are: 
%     - xy: boudary of ROI (x,y; matrix in two columns)
%     - NN: the number of grids inside ROI.
%     - XY: coordinates of grids inside ROI (X,Y; matrix in two columns).
%     - R_data: hydroclimatic variable from datasets.
%
% The ouput is a structure with fields containing:
%     - PP.Area: the size of ROI
%     - PP.XCentroid: the location of ROI demonstrated by ROI's geometric
%     centre (x-coordinate)
%     - PP.YCentroid: the location of ROI demonstrated by ROI's geometric
%     centre (y-coordinate)
%     - PP.s: the tan of rotated angle of the princpal axes
%     - PP.MajorAxesAngle = the orientation of ROI with respect to the North/Upward direction;
%     - PP.MaxMajorAxes = the longest side of outside rectangle of ROI;
%     - PP.MinMajorAxes = the shortest side of outside rectangle of ROI;
%
% Code written by: H. Wang, <zzwang0924@gmail.com> 
% Partly redistribute the code for calculating moments of polygon (Software written by: Ciro A. Soto, <ciro@kavyata.com>
% Created: 2015-11-22)
% 
% Created: 2019-11-27
% Copyright (c) 2019, H. Wang
% All rights reserved.
% 
% Disclaimer:
% This program (hereafter, software) is designed for instructional, educational and research use only.
% Commercial use is prohibited. The software is provided 'as is' without warranty
% of any kind, either express or implied. The software could include technical or other mistakes,
% inaccuracies or typographical errors. The use of the software is done at your own discretion and
% risk and with agreement that you will be solely responsible for any damage and that the authors
% and their affiliate institutions accept no responsibility for errors or omissions in the software
% or documentation. In no event shall the authors or their affiliate institutions be liable to you or
% any third parties for any special, indirect or consequential damages of any kind, or any damages whatsoever.
%
%==========================================================================
N=size(xy,1);
x=xy(:,1);
y=xy(:,2);
x=[x;x(1)];
y=[y;y(1)];
PP.boundary = [x,y];

% Compute common values and Area Centroid
wi=x(1:N).*y(2:N+1)-x(2:N+1).*y(1:N);

% area of the polygon
Area=0.5*sum(wi);
PP.Area=NN;

% area moments
MAx=1/6*sum(wi.*(y(1:N) + y(2:N+1)));
MAy=1/6*sum(wi.*(x(1:N) + x(2:N+1)));

% Area moments of second degree
Ixx=1/12*sum(wi.*((y(1:N)+y(2:N+1)).^2 - y(1:N).*y(2:N+1)));
Iyy=1/12*sum(wi.*((x(1:N)+x(2:N+1)).^2 - x(1:N).*x(2:N+1)));
Ixy=1/24*sum(wi.*((x(1:N)+x(2:N+1)).*(y(1:N)+y(2:N+1)) + x(1:N).*y(1:N) + x(2:N+1).*y(2:N+1)));


% coordinates of the area centroid
ACx=MAy/Area;
ACy=MAx/Area;
PP.XCentriod=ACx;
PP.YCentriod=ACy;
% area moments in the area centroid
IxxAC=Ixx-ACy^2*Area;
IyyAC=Iyy-ACx^2*Area;
IxyAC=Ixy-ACx*ACy*Area;
PrincAxesRotationDeg=0.5*atan(2*IxyAC/(IyyAC-IxxAC))*180/pi;

% plot area in gray
% patch(x,y,0.7*[1 1 1],'LineStyle','none')

a=PrincAxesRotationDeg*pi/180;
xcg=ACx;
ycg=ACy;
s=tan(a);
PP.s = s;
xx=[min(x) max(x)];
yy=s*(xx-xcg)+ycg;
% plot principal axes
% plot(xx,yy,'-b','linewidth',2)
yy1=[min(y) max(y)];
xx1 = s*(ycg-yy1)+xcg;

% plot(xx1,yy1,'-b','linewidth',2)
leng1 = sqrt((xx(1)-xx(2))^2+(yy(1)-yy(2)).^2);
leng2 = sqrt((xx1(1)-xx1(2))^2+(yy1(1)-yy1(2)).^2);
id=[];flag=[];
if leng1>leng2
    id = find(yy == min(yy));
    id1 = find(yy == max(yy));
    v_1 = [xx(id1),yy(id1),0] - [xx(id),yy(id),0];
    v_2 = [xx(id),yy(id)+100,0] - [xx(id),yy(id),0];
    if xx(id1)<=xx(id)
        flag = 1;
    else
        flag = 0;
    end
else
    id = find(yy1 == min(yy1));
    id1 = find(yy1 == max(yy1));
    v_1 = [xx1(id1),yy1(id1),0] - [xx1(id),yy1(id),0];
    v_2 = [xx1(id),yy1(id)+100,0] - [xx1(id),yy1(id),0];
    if xx1(id1)<=xx1(id)
        flag = 1;
    else
        flag = 0;
    end
end
if length(id)~=1 || isempty(id)
    if leng1>leng2
        Theta = 90;
    elseif leng1==leng2
        Theta = NaN;
    else
        Theta = 0;
    end
else
    Theta = atan2d(norm(cross(v_1, v_2)), dot(v_1, v_2))+360*(norm(cross(v_1, v_2))<0);
end
if flag == 1
    Theta = -Theta;
end
if Theta<0
    Theta = -90-Theta;
elseif Theta>0
    Theta = 90-Theta;
end
PP.MajorAxesAngle = -Theta;
[recty,rectx,~,~] =  minboundrect(x, y);
[wei,hei] = minboxing(rectx(1:end-1),recty(1:end-1));
PP.MaxMajorAxes = hei;
PP.MinMajorAxes = wei;
PP.sp = PP.MaxMajorAxes/PP.MinMajorAxes;
[output,cidx,indn,ind1,cmeans] = core(XY,R_data,coreind);
PP.core = output;
PP.cidx = cidx;
PP.indn = indn;
PP.ind1 = ind1;
PP.cmeans =cmeans;
end

function [ wei,hei ] = minboxing( x, y )

xy = [x, y];
xy1 = xy([4 1 2 3],:);
side = sqrt(sum((xy-xy1).^2,2));
wei = min(side(1:2));
hei = max(side(1:2));

end


function [rectx,recty,area,perimeter] = minboundrect(x,y,metric)

if (nargin<3) || isempty(metric)
    metric = 'a';
elseif ~ischar(metric)
    error 'metric must be a character flag if it is supplied.'
else
    % check for 'a' or 'p'
    metric = lower(metric(:)');
    ind = strmatch(metric,{'area','perimeter'});
    if isempty(ind)
        error 'metric does not match either ''area'' or ''perimeter'''
    end
    
    % just keep the first letter.
    metric = metric(1);
end

% preprocess data
x=x(:);
y=y(:);

% not many error checks to worry about
n = length(x);
if n~=length(y)
    error 'x and y must be the same sizes'
end

if n>3
    if (var(x)== 0|| var(y)==0)
        if var(x)== 0
            x = [x-1;x(1); x+1 ];
            y = [y ;y(1);y];
            flag = 1;
        else
            y = [y-1;y(1); y+1 ];
            x = [x ;x(1);x];
            flag = 1;
        end
    else
        flag = 0;
        edges = convhull(x,y);
    end
    
    if flag == 0
        x = x(edges);
        y = y(edges);
        
    end
    nedges = length(x) - 1;
elseif n>1
    % n must be 2 or 3
    nedges = n;
    x(end+1) = x(1);
    y(end+1) = y(1);
else
    % n must be 0 or 1
    nedges = n;
end

switch nedges
    case 0
        rectx = [];
        recty = [];
        area = [];
        perimeter = [];
        return
    case 1
        rectx = repmat(x,1,5);
        recty = repmat(y,1,5);
        area = 0;
        perimeter = 0;
        return
    case 2
        rectx = x([1 2 2 1 1]);
        recty = y([1 2 2 1 1]);
        area = 0;
        perimeter = 2*sqrt(diff(x).^2 + diff(y).^2);
        return
end
% 3 or more points.
Rmat = @(theta) [cos(theta) sin(theta);-sin(theta) cos(theta)];

% get the angle of each edge of the hull polygon.
ind = 1:(length(x)-1);
edgeangles = atan2(y(ind+1) - y(ind),x(ind+1) - x(ind));
% move the angle into the first quadrant.
edgeangles = unique(mod(edgeangles,pi/2));

% now just check each edge of the hull
nang = length(edgeangles);
area = inf;
perimeter = inf;
met = inf;
xy = [x,y];
for i = 1:nang
    % rotate the data through -theta
    rot = Rmat(-edgeangles(i));
    xyr = xy*rot;
    xymin = min(xyr,[],1);
    xymax = max(xyr,[],1);
    % The area is simple, as is the perimeter
    A_i = prod(xymax - xymin);
    P_i = 2*sum(xymax-xymin);
    if metric=='a'
        M_i = A_i;
    else
        M_i = P_i;
    end
    if M_i<met
        % keep this one
        met = M_i;
        area = A_i;
        perimeter = P_i;
        rect = [xymin;[xymax(1),xymin(2)];xymax;[xymin(1),xymax(2)];xymin];
        rect = rect*rot';
        rectx = rect(:,1);
        recty = rect(:,2);
    end
end

end

function [output,cidx,indn,ind1,cmeans] = core(XY,R_data,coreind)
XY1 = XY;
RP = R_data;

% find the core
meas = [XY1,RP];
if length(XY1)<500
    klist=2:10;%the number of clusters
    eva.OptimalK = klist(end-1);
elseif length(XY1)<1000 && length(XY1)>500
    klist=2:20;%the number of clusters
    eva.OptimalK = klist(end-1);
else
    klist=2:30;%the number of clusters
    eva.OptimalK = klist(end-1);
end
myfunc = @(X,K)(kmeans(X, K));
while eva.OptimalK == klist(end-1)
    eva = evalclusters(meas,myfunc,'CalinskiHarabasz','klist',klist);
    klist(end+1) = klist(end)+1;
end
[cidx,cmeans] = kmeans(meas,eva.OptimalK,'dist','sqeuclidean','replicates',10);

if coreind == 1
[~,ind1] = sort(cmeans(:,3),'descend');
deta = max(RP)-cmeans(ind1(1),3);
deta = deta*1.2;
indn = cmeans(:,3)>=max(RP)-deta;
else
    [~,ind1] = sort(cmeans(:,3),'ascend');
    deta = min(RP)-cmeans(ind1(1),3);
deta = deta*1.2;
indn = cmeans(:,3)<=min(RP)-deta;
end

output = [cmeans(ind1(1),2),cmeans(ind1(1),1),cmeans(ind1(1),3)];


end


