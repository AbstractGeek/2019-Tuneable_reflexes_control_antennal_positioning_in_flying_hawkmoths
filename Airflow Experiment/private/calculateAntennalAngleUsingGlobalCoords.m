function [antAngles] = calculateAntennalAngleUsingGlobalCoords(digitizedData)
% function [] = calculateAntennalAngle()
% Details about the function goes here
% 
% 

% Extract coordinates for easy reference
righttip = digitizedData(:,1:3);
lefttip = digitizedData(:,4:6);
rightbase = digitizedData(:,7:9);
leftbase = digitizedData(:,10:12);

rowNum = size(righttip,1);
% From Calibration object - Coordinate specifications
% X axis along wind direction.
% Y axis, perpendicular to the wind direction, the width of the wind
% tunnel.
% Z axis, the height of the wind tunnel.

% Base vector - gives the head movements
% Should I nullify the base movements by fixing their coordinates?
% The pivot point of the base-vector is the neck. So rotation might not be
% good.
% 
% Interantennal angle will be independent of head rotation
% Individual antennal angles will definitely depend on head rotation.
% How do i remove this out?
% 
% Plot out yaw, pitch and roll of the head, using the base vector?
% 

% Global coordinates vector definitions
i = repmat([1 0 0],rowNum,1);
j = repmat([0 1 0],rowNum,1);
k = repmat([0 0 1],rowNum,1);

% % Define relevant vectors
basevect = leftbase - rightbase;
rightantvect = righttip - rightbase; 
leftantvect = lefttip - leftbase;
basemidvect = (leftbase+rightbase)/2;

%%% Get vector interantennal angles.
iaaVect = atan2(matrixNorm(cross(rightantvect,leftantvect,2)),dot(rightantvect,leftantvect,2)).*180/pi;

% Define zero vectors (with NaN) to make combination easier.
zeroVect = NaN(size(rightantvect));
zeroVect(~isnan(rightantvect)) = 0;

% Find xy component of interantennal angle.
rightantXYProj = [rightantvect(:,1:2),zeroVect(:,3)];
leftantXYProj = [leftantvect(:,1:2),zeroVect(:,3)];
iaaXY = atan2(matrixNorm(cross(rightantXYProj,leftantXYProj,2)),dot(rightantXYProj,leftantXYProj,2)).*180/pi;
laaXY = atan2(matrixNorm(cross(leftantXYProj,i,2)),dot(leftantXYProj,i,2)).*180/pi;
raaXY = atan2(matrixNorm(cross(rightantXYProj,i,2)),dot(rightantXYProj,i,2)).*180/pi;

% Find xz component of interantennal angle
rightantXZProj = [rightantvect(:,1),zeroVect(:,2),rightantvect(:,3)];
leftantXZProj = [leftantvect(:,1),zeroVect(:,2),leftantvect(:,3)];
iaaXZ = atan2(matrixNorm(cross(rightantXZProj,leftantXZProj,2)),dot(rightantXZProj,leftantXZProj,2)).*180/pi;
laaXZ = atan2(matrixNorm(cross(leftantXZProj,i,2)),dot(leftantXZProj,i,2)).*180/pi;
raaXZ = atan2(matrixNorm(cross(rightantXZProj,i,2)),dot(rightantXZProj,i,2)).*180/pi;

% Find yz component of interantennal angle
% Probably (almost definitely) redundant. Let's calculate it anyway and see.
rightantYZProj = [zeroVect(:,1),rightantvect(:,2:3)];
leftantYZProj = [zeroVect(:,1),leftantvect(:,2:3)];
iaaYZ = atan2(matrixNorm(cross(rightantYZProj,leftantYZProj,2)),dot(rightantYZProj,leftantYZProj,2)).*180/pi;
laaYZ = atan2(matrixNorm(cross(leftantYZProj,j,2)),dot(leftantYZProj,j,2)).*180/pi;
raaYZ = atan2(matrixNorm(cross(rightantYZProj,j,2)),dot(rightantYZProj,j,2)).*180/pi;

% The above method does not give proper sign to the angles, with respect to
% the comparision axis. Using atan2 should solve the problem. Lets try it
% out.
% XY plane
raaXYws = atan2(rightantvect(:,2),rightantvect(:,1));       % ws = with sign
laaXYws = atan2(leftantvect(:,2),leftantvect(:,1));
% XZ plane
raaXZws = atan2(rightantvect(:,3),rightantvect(:,1));
laaXZws = atan2(leftantvect(:,3),leftantvect(:,1));
% YZ plane
raaYZws = atan2(rightantvect(:,3),rightantvect(:,2));
laaYZws = atan2(leftantvect(:,3),leftantvect(:,2));

% Head rotation calculated using base-vectors
% Head Yaw is the angle of the base-vector projection in the xy plane. (with respect to y)
basevectXYProj= [basevect(:,1:2),zeroVect(:,3)];
headYaw = atan2(matrixNorm(cross(basevectXYProj,j,2)),dot(basevectXYProj,j,2)).*180/pi;
% Head roll is the angle of base-vector projection in the yz plane (with respect to y)
basevectYZProj= [zeroVect(:,1),basevect(:,2:3)];
headRoll = atan2(matrixNorm(cross(basevectYZProj,j,2)),dot(basevectYZProj,j,2)).*180/pi;
% Head pitch is a bit harder to calculate. Using the midpoint of the
% basevect as a simple proxy to calculate the gross pitch angle of the
% head.
basemidXZProj = [basemidvect(:,1),zeroVect(:,2),basemidvect(:,3)];
headPitch = atan2(matrixNorm(cross(basemidXZProj,i,2)),dot(basemidXZProj,i,2)).*180/pi;

% Calculations Done!
% Add the Custom Plane and XY plane data to the array.
antAngles = calculateAntennalAngleUsingCustomPlane(digitizedData);

% Sort into output structure
antAngles.xyangle.iaa = iaaXY;
antAngles.xyangle.laa = correctAngle(laaXY);
antAngles.xyangle.raa = correctAngle(raaXY);

antAngles.xzangle.iaa = iaaXZ;
antAngles.xzangle.laa = correctAngle(laaXZ);
antAngles.xzangle.raa = correctAngle(raaXZ);

antAngles.yzangle.iaa = iaaYZ;
antAngles.yzangle.laa = correctAngle(laaYZ);
antAngles.yzangle.raa = correctAngle(raaYZ);

antAngles.xysignedangle.iaa = iaaXY;
antAngles.xysignedangle.laa = laaXYws;
antAngles.xysignedangle.raa = raaXYws;

antAngles.xzsignedangle.iaa = iaaXZ;
antAngles.xzsignedangle.laa = laaXZws;
antAngles.xzsignedangle.raa = raaXZws;

antAngles.yzsignedangle.iaa = iaaYZ;
antAngles.yzsignedangle.laa = laaYZws;
antAngles.yzsignedangle.raa = raaYZws;

antAngles.headangle.Roll = headRoll;
antAngles.headangle.Yaw = headYaw;
antAngles.headangle.Pitch = headPitch;

% Add the Head point data to the array
if size(digitizedData,2)==15 && sum(isnan(digitizedData(:,15)))<=90 
    antAngles.headplane = calculateAntennalAngleUsingHeadPoint(digitizedData);
end

antAngles.VectAngle.iaa = iaaVect;
antAngles.VectAngle.laa = NaN(size(iaaVect));
antAngles.VectAngle.raa = NaN(size(iaaVect));

end

function [unitVect] = unitizeVect(inVect)
% function [unitVect] = unitizeVect(inVect)
% Converts the input vector into a unit vector
% 
% Dinesh Natesan, 10.10.14
unitVect = inVect./repmat(sqrt(sum(inVect.^2,2)),1,3);
end

function [corrAngle] = correctAngle(angle)
incorr = angle>90;
corrAngle = incorr.*180 + (incorr.*-1 + not(incorr)).*angle;
end