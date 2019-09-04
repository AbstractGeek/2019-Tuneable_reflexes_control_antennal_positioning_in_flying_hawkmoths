function [projAngles,debugvars] = ...
    calculatePerturbedAntennalAngleUsingHeadVect(digitizedData,filterflag)

% 
N = size(digitizedData,1);
projAngles = NaN(N,2);
defaultvars = who;
defaultvars = [defaultvars;{'defaultvars';'currentvars';'debugvars'}];

% Extract coordinates for easy reference
righttip = digitizedData(:,1:3);
lefttip = digitizedData(:,4:6);
rightbase = digitizedData(:,7:9);
leftbase = digitizedData(:,10:12);
head = digitizedData(:,13:15);

% Save data as the debug variable
currentvars = who;
currentvars(ismember(currentvars,defaultvars)) = [];
debugvars.rawData = cell2struct(...
    cellfun(@eval,currentvars,'UniformOutput',false),currentvars,1);

% Get unfiltered head and base vector.
unfiltbasemidpt = (leftbase + rightbase)./2;
unfiltheadvect = unitizeVect(head - unfiltbasemidpt);    % Using filtered points to make it noise invariant
unfiltbasevect = unitizeVect(rightbase - leftbase);

% Filter head and base points
[filtrightbase,filtleftbase,filthead] = ...
   smartFilterBaseHeadPoints(rightbase,leftbase,head);

% Calculate mid points and obtain filtered head and based vectors
filtbasemidpt = (filtleftbase+filtrightbase)./2;
filtheadvect = unitizeVect(filthead - filtbasemidpt);    % Using filtered points to make it noise invariant
filtbasevect = unitizeVect(filtrightbase - filtleftbase);

% Create filtered antenna vectors
filtrighttipvect = unitizeVect(righttip - filtrightbase); 
filtlefttipvect = unitizeVect(lefttip - filtleftbase);

% Create vectors
righttipvect = unitizeVect(righttip - rightbase); 
lefttipvect = unitizeVect(lefttip - leftbase);

if filterflag
    headvect = filtheadvect;
    basevect = filtbasevect;
else
    headvect = unfiltheadvect;
    basevect = unfiltbasevect;
end

% Obtain antennal angles
projAngles(:,1) = atan2(matrixNorm(cross(righttipvect,headvect,2)),...
    dot(righttipvect,headvect,2)).*180/pi;
projAngles(:,2) = atan2(matrixNorm(cross(lefttipvect,headvect,2)),...
    dot(lefttipvect,headvect,2)).*180/pi;

% Save data as the debug variable
currentvars = who;
currentvars(ismember(currentvars,defaultvars)) = [];
currentvars(ismember(currentvars,fieldnames(debugvars.rawData))) = [];
debugvars.processedData = cell2struct(...
    cellfun(@eval,currentvars,'UniformOutput',false),currentvars,1);

% get xyz points for debugging
xyzPoints = NaN(N,3*6);
xyzPoints(:,1:15) = digitizedData;
xyzPoints(:,16:18) = filtbasemidpt;

debugvars.xyzPoints = xyzPoints;

end