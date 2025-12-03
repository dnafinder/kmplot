function varargout = kmplot(x,varargin)
% KMPLOT Plot the Kaplan-Meier estimation of the survival function
% Survival times are data that measure follow-up time from a defined
% starting point to the occurrence of a given event, for example the time
% from the beginning to the end of a remission period or the time from the
% diagnosis of a disease to death. Standard statistical techniques cannot
% usually be applied because the underlying distribution is rarely Normal
% and the data are often "censored". A survival time is described as
% censored when there is a follow-up time but the event has not yet
% occurred or is not known to have occurred. For example, if remission time
% is being studied and the patient is still in remission at the end of the
% study, then that patient's remission time would be censored. If a patient
% for some reason drops out of a study before the end of the study period,
% then that patient's follow-up time would also be considered to be
% censored. The survival function S(t) is defined as the probability of
% surviving at least to time t. The graph of S(t) against t is called the
% survival curve. The Kaplan-Meier method can be used to estimate this
% curve from the observed survival times without the assumption of an
% underlying probability distribution.
%
% Syntax: 	kmplot(x,alpha,censflag)
%      
%     Inputs:
%           X (mandatory)- N-by-2 data matrix:
%                          (X(:,1)) = survival time of the i-th subject
%                          (X(:,2)) = censored flag 
%                                     (0 if not censored; 1 if censored)
%           Note: if X is a vector, it will be internally converted to
%                 [X(:) zeros(N,1)], i.e. all observations are treated as
%                 uncensored (flag = 0).
%
%           ALPHA (optional)    - significance level (default 0.05) 
%           CENSFLAG (optional) - censored plot flag (default 0). If 0
%                                 censored data will be plotted spread over
%                                 the horizontal segment; if 1 they will be
%                                 plotted at the given time of censoring.
%           FLAG (optional)     - internal flag for LOGRANK (default 1).
%                                 If 1, kmplot performs full plotting and
%                                 confidence intervals. If 0, it skips
%                                 plotting and some computations, and is
%                                 meant to be called internally by LOGRANK.
%
%     Outputs:
%           Kaplan-Meier plot (if FLAG == 1) and, if requested:
%           table1  - table of event times, counts and numbers at risk
%           table12 - table of censored times (from TABULATE)
%           t1      - time points used for the step function
%           T1      - Kaplan-Meier survival estimates at t1
%           xcg     - x-coordinates of censored marks (if any)
%           ycg     - y-coordinates of censored marks (if any)
%           lambda  - average hazard rate estimate
%
%      Example: (+ indicate that patient is censored)
%  
%                   ---------------------
%                   Patient     Survival
%                               time       
%                   ---------------------
%                      1        2+    
%                      2        2+   
%                      3        2+    
%                      4        3
%                      5        4
%                      6        5    
%                      7        5+    
%                      8        6
%                      9        7
%                     10        8    
%                     11        9+    
%                     12        10
%                     13        14    
%                     14        14+
%                     15        18+
%                     16        19+  
%                     17        20   
%                     18        22
%                     19        22+
%                     20        31+  
%                     21        33
%                     22        39
%                     23        195+    
%                   ----------------------
%    x=[2 1; 2 1; 2 1; 3 0; 4 0; 5 0; 5 1; 6 0; 7 0; 8 0; 9 1; 10 0;...
%    14 0; 14 1; 18 1; 19 1; 20 0; 22 0; 22 1; 31 1; 33 0; 39 0; 195 1];
%
%    Calling in Matlab the function: kmplot(X) the function will plot the
%    Kaplan-Meier estimation of the survival function
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo.75@gmail.com
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008). KMPLOT: Kaplan-Meier estimation of the survival
% function.
% http://www.mathworks.com/matlabcentral/fileexchange/22293

% Input Error handling
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'}, ...
    {'2d','real','finite','nonnan','nonempty'}));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'}, ...
    {'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p,'cflag',0, @(x) isnumeric(x) && isreal(x) && ...
    isfinite(x) && isscalar(x) && (x==0 || x==1));
addOptional(p,'flag',1, @(x) isnumeric(x) && isreal(x) && ...
    isfinite(x) && isscalar(x) && (x==0 || x==1));
parse(p,x,varargin{:});
alpha = p.Results.alpha;
cflag = p.Results.cflag;
flag  = p.Results.flag;
x     = p.Results.x;
clear p

% Allow x as either a vector (all uncensored) or an N-by-2 matrix
if isvector(x)
    % Ensure column vector and create censor flag column (all zeros)
    x = x(:);
    x = [x zeros(size(x,1),1)];
elseif size(x,2) == 2
    % Already in [time, censorFlag] form
else
    error('kmplot:InvalidXSize', ...
          'Input x must be either a vector or an N-by-2 matrix.');
end

% Check censoring flags
assert(all(x(:,2)==0 | x(:,2)==1), ...
    'kmplot:InvalidCensorFlag', ...
    'All x(:,2) values must be 0 (uncensored) or 1 (censored).');

% String for LEGEND function
str1 = [num2str((1-alpha)*100) '% confidence interval']; 

% Sort data by survival time
x = sortrows(x,1);

% Table of patients observed for each survival time
% table1 = [time count percent(on total)]
table1 = [0 size(x,1) 1; tabulate(x(:,1))];
% Remove not observed times added by TABULATE function
table1(table1(:,3)==0,:) = [];

% Table of censored data
table12 = tabulate(x(x(:,2)==1));
if ~isempty(table12)
    % Remove not observed times added by TABULATE function
    table12(table12(:,3)==0,:) = [];
    % Setup the vector of the censored data
    [cens,loc] = ismember(table1(:,1),table12(:,1)); % find censored data
end    

% The percents stored in the third column are unuseful;
% place in the third column how many subjects are still alive at the
% beginning of the i-th interval.
a1 = [table1(1,2); -1.*table1(2:end,2)];
table1(:,3) = cumsum(a1);
table1(2:end,3) = table1(1:end-1,3);

% Number of deaths in the intervals (don't take into account censored data)
if ~isempty(table12)
    table1(cens,2) = table1(cens,2) - table12(loc(cens),2);
end

% Finally, delete the first row that is now useless
table1(1,:) = [];

t1 = [0; table1(:,1)]; % this is the x variable (time)
% This is the y variable (survival function)
T1 = [1; cumprod(1-(table1(:,2)./table1(:,3)))];

if flag % if this function was not called by LOGRANK function
    % Compute the standard error of the survival function (Greenwood)
    SE = [0; T1(2:end).*sqrt(cumsum(table1(:,2)./(table1(:,3).* ...
        (table1(:,3)-table1(:,2)))))];
end

% Censored data plotting
if ~isempty(table12) 
    % If there are censored data after max(t1), add a new cell into the t1,
    % T1 and SE arrays
    if table12(end,1) >= t1(end,1)
        t1(end+1,1) = table12(end,1)+1;
        T1(end+1,1) = T1(end,1);
        if flag % if this function was not called by LOGRANK function
            SE(end+1,1) = SE(end,1);
        end
    end
    if ~cflag
        % Vectors preallocation
        xcg = zeros(1,sum(table12(:,2)));
        ycg = xcg;
        J   = 1;
        % For each censored data into the i-th time interval...
        for I = 1:size(table12,1)
            % Compute how many positions into the array they must occupy
            JJ = J + table12(I,2) - 1;
            % Find the correct time interval in which censored data must be
            % placed
            A = find(t1<=table12(I,1),1,'last');
            B = find(t1>table12(I,1),1,'first');
            % Equally divide this interval
            int = linspace(table12(I,1),t1(B,1),table12(I,2)+2);
            % Put all in the vectors of the plotting variables
            xcg(J:JJ) = int(2:end-1);
            ycg(J:JJ) = T1(A);
            % Update the counter
            J = JJ+1;
        end
    else
        xcg = table1(table1(:,2)==0,1);
        T1tmp = T1;
        ycg = T1tmp(table1(:,2)==0);
    end
else
    if ~flag % if this function was called by LOGRANK function
        xcg = [];
        ycg = [];
    end
end

% Compute the hazard rate (average hazard from log-survival)
c1 = T1.*numel(x);
c2 = -(diff(log(c1(1:end-1)))./diff(t1(1:end-1)));
lambda = mean(c2(c2~=0));

if flag % if this function was not called by LOGRANK function
    % Compute the (1-alpha)*100% confidence interval curves
    cv = realsqrt(2)*erfcinv(alpha); % critical value
    % Lower curve (remember that: the lower curve values can't be negative)
    lowc = max(0,T1-SE.*cv);
    % If the lower curve reaches 0 earlier than survival function, trim the data
    if isequal(lowc(end-1:end),[0; 0])
        lowcend = find(lowc==0,1,'first');
    else
        lowcend = length(lowc);
    end
    % Upper curve (remember that the upper curve values can't be >1)
    upc = min(1,T1+SE.*cv);
    % Eventually, correct the data
    if isequal(upc(end),1) 
        cupend = find(upc<1,1,'last');
        upc(cupend:end) = upc(cupend);
    end

    % Compute the median survival time (if it exists...)
    if isempty(T1(T1==0.5)) % if there is not a point where T=0.5...
        I = find(T1>0.5,1,'last'); % first point where T>0.5
        J = find(T1<0.5,1,'first'); % first point where T<0.5
        if isempty(J) % if all points are >0.5...
            mt = 0; % there is no median time
        else 
            % Compute the median time by linear interpolation
            p  = polyfit([t1(I) t1(J)],[T1(I) T1(J)],1);
            mt = (0.5-p(2))/p(1);
            str2 = ['Median time ' num2str(mt)]; % string for LEGEND
        end
    else
        mt   = t1(T1==0.5);
        str2 = ['Median time ' num2str(mt)]; % string for LEGEND
    end
 
    % Plot all the data
    hFig = figure;
    units = get(hFig,'units');
    set(hFig,'units','normalized','outerposition',[0 0 1 1], 'Color', 'white');
    set(hFig,'units',units);
    hold on
    S2 = stairs(t1(1:lowcend),lowc(1:lowcend),'g--'); % lower CI curve
    stairs(t1,upc,'g--');                             % upper CI curve
    S1 = stairs(t1,T1,'b');                           % KM survival function
    grid on
    if mt>0 % if there exists a median time...
        S3 = plot([0 mt mt],[0.5 0.5 0],'k:'); 
    else
        S3 = [];
    end
    if ~isempty(table12) % if there are censored data...
        S4 = plot(xcg,ycg,'r+');
    else
        S4 = [];
    end
    hold off

    % Set the axis properly
    xmax = max(t1)+1;
    axis([0 xmax 0 1.2]);
    axis square

    % Add labels and legend
    txt = sprintf('Kaplan-Meier estimate of survival function (hazard rate: %0.4f)\n',lambda);
    title(txt,'FontName','Arial','FontSize',14,'FontWeight','Bold'); 
    ylabel('Estimated survival function','FontName','Arial','FontSize',14,'FontWeight','Bold'); 
    xlabel('Time','FontName','Arial','FontSize',14,'FontWeight','Bold'); 

    if mt
        if isempty(S4)
            legend([S1 S2 S3],'Data',str1,str2)
        else
            legend([S1 S2 S3 S4],'Data',str1,str2,'Censored')
        end
    else
        if isempty(S4)
            legend([S1 S2],'Data',str1)
        else
            legend([S1 S2 S4],'Data',str1,'Censored')
        end
    end
end

% Outputs
if nargout
    out = {table1, table12, t1, T1, xcg, ycg, lambda};
    varargout = out(1:nargout);
end
end
