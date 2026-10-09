% NEPAL GLOF ORBIT OPTIMIZATION
% -------------------------------------------------------------------------
% Objective
%   Optimize a 550 km circular LEO orbit to balance:
%     1) Nepal swath coverage,
%     2) number of Nepal swath passes per day,
%     3) Sagarmatha ground-station contact time,
%     4) revisit time and observation-to-contact latency,
%     5) coverage of six GLOF-vulnerable glacial lakes,
%     6) lake revisit performance and minimum lake passes/day.
%
% Workflow
%   A. Define mission/orbit/constellation configuration.
%   B. Load and prepare the Nepal polygon once.
%   C. Build fixed geometry, scoring, and optimization settings.
%   D. Generate the initial Latin-hypercube design.
%   E. Fit a GP surrogate and maximize Expected Improvement.
%   F. Select and report the best evaluated orbit.
%   G. Save all optimization results.
%   H. Re-simulate the best orbit at high temporal resolution, generate the
%      static 2-D map, 3-D scenario view, animated pass build-up, overlap
%      geometry, and mission summaries.
%
% Core geometry
%   Sensor half-angle: theta = atan(r / h), so r = h*tan(theta).
%   Local footprint ellipse in degree space:
%       rLat = r / 111.32
%       rLon = r / (111.32*cosd(latitude))
%   Point inclusion test:
%       ((lonTarget-lonSat)/rLon)^2 + ((latTarget-latSat)/rLat)^2 <= 1
%
% Coverage/overlap
%   A_i = area of pass i inside Nepal.
%   S   = sum_i A_i              (overlaps counted repeatedly)
%   U   = area(union_i pass_i)   (unique coverage)
%   S-U = algebraic overlap redundancy.
%   A>=2 = union of all pairwise intersections.
%   Consistency check used here:
%       S-U = A>=2 + A_triple
%
% Optimization
%   x = [inclination, RAAN, true anomaly].
%   f(x) is the weighted mission score. A GP models f(x), and Expected
%   Improvement is maximized to select the next expensive simulation.
%
% Constellation readiness
%   The mission is represented as arrays of orbital elements, satellites, and
%   ground stations. Current production mode remains 1 plane x 1 satellite x
%   1 ground station. Coverage, revisit, contact, lake detection, union, and
%   per-satellite/per-plane metrics are computed from the same array-based
%   path, so later constellation growth does not require a second model.
%
% Accuracy range
%   The nominal area uses a 16-sided footprint approximation. The final
%   best-orbit union is independently recomputed with 64 footprint vertices.
%   The reported numerical geometry sensitivity is
%       margin = |A_64 - A_16|
%   and the displayed range is A_16 +/- margin. This is a numerical-model
%   sensitivity band, not a claim of absolute geodetic or sensor accuracy.
% -------------------------------------------------------------------------

clc;
clear;
close all;

%% =========================================================
% A. FIXED PARAMETERS
%% =========================================================

% Earth/orbit/sensor constants.
earthRadius     = 6378e3;
perigeeAltitude = 550e3;
eccentricity    = 0.0;
semiMajorAxis   = (earthRadius + perigeeAltitude) / (1 - eccentricity);
sensorHalfAngle = 10;
theta_half      = deg2rad(sensorHalfAngle);
startTime       = datetime(2025,5,1,0,0,0);
stopTime        = startTime + hours(24);
baseSampleTime  = 2;       % high-resolution output/plots
optimizationSampleTime = 60; % faster objective evaluations

% Fixed circular-footprint geometry constants.
kmPerDegLat = 111.32;
tanHalfAngle = tan(theta_half);
nCirclePts = 16;
thetaCircle = linspace(0,2*pi,nCirclePts+1);
thetaCircle(end) = [];
circleCos = cos(thetaCircle);
circleSin = sin(thetaCircle);

% Ground-station configuration. Keep one station for the current run; the
% array-based representation is ready for additional stations.
groundStations = struct( ...
    'name','Sagarmatha Satellite E/S', ...
    'lat',27.6889797780441, ...
    'lon',85.23841620362452);

% Constellation configuration. Keep 1 plane x 1 satellite for the current
% numerical case. Every computation downstream operates on these arrays, so
% adding planes/satellites later does not require changing the mission API.
constellation = struct( ...
    'numOrbits',1, ...
    'satsPerOrbit',1, ...
    'raanSpacingDeg',[], ...
    'trueAnomalySpacingDeg',[], ...
    'raanOffsetsDeg',[], ...
    'trueAnomalyOffsetsDeg',[], ...
    'inclinationOffsetsDeg',[]);
constellation.raanSpacingDeg = 360/max(constellation.numOrbits,1);
constellation.trueAnomalySpacingDeg = 360/max(constellation.satsPerOrbit,1);

% Accuracy/sensitivity calculation. The nominal model retains 16 footprint
% vertices for speed. A final 64-vertex recomputation estimates the numerical
% geometry-resolution sensitivity of the reported Nepal union area.
accuracy = struct( ...
    'enabled',true, ...
    'referenceCirclePts',64);

% Six GLOF-vulnerable glacial lakes.
% DMS conversion: decimal degrees = deg + min/60 + sec/3600.
glacierNames = {'Imja Tsho','Tsho Rolpa','Thulagi','Chamlang South', ...
    'Barun Tsho','Lumding Tsho'};
glacierLat = [ ...
    27 + 53/60 + 56.53/3600, ...
    27 + 51/60 + 40.48/3600, ...
    28 + 29/60, ...
    27 + 45/60, ...
    27 + 47/60 + 51/3600, ...
    27 + 46/60 + 48/3600];
glacierLon = [ ...
    86 + 55/60 + 20.13/3600, ...
    86 + 28/60 + 32.25/3600, ...
    84 + 29/60, ...
    86 + 57/60, ...
    87 + 5/60 + 26/3600, ...
    86 + 36/60 + 48/3600];
glacierElevM   = [5003.13,4549.91,4003,4924,4530,4831];
glacierAreaKm2 = [1.741,1.712,0.94,0.86,2.193,1.475];

% Backward-compatible aliases for the current single ground station.
gsLat = groundStations(1).lat;
gsLon = groundStations(1).lon;

%% =========================================================
% B. LOAD NEPAL BORDER
%% =========================================================
% Process: load the Nepal administrative polygon, select the Nepal feature,
% remove degenerate points, create a polyshape, and cache its area and local
% degree-to-kilometre conversion factors. This work is done once because it
% is identical for every orbit evaluation.

shpFile = 'ne_10m_admin_0_countries_nep.shp';
nepalBorderLon = [];
nepalBorderLat = [];

try
    if contains(shpFile,'nep.shp')
        nepalShape = shaperead(shpFile,'UseGeoCoords',true);
    else
        nepalShape = shaperead(shpFile,'UseGeoCoords',true, ...
            'Selector',{@(name) strcmpi(name,'Nepal'),'NAME'});
    end

    if isempty(nepalShape)
        error('Nepal feature not found in shapefile.');
    end

    if numel(nepalShape) > 1
        selectedIdx = [];
        nameFields = {'NAME','ADMIN','name','Name','ADMIN_NAME','NAME_EN','NAME_0','Country'};
        for f = 1:numel(nameFields)
            field = nameFields{f};
            if isfield(nepalShape,field)
                vals = {nepalShape.(field)};
                idx = find(contains(vals,'Nepal','IgnoreCase',true),1);
                if ~isempty(idx)
                    selectedIdx = idx;
                    break;
                end
            end
        end

        if isempty(selectedIdx)
            for i = 1:numel(nepalShape)
                bbox = nepalShape(i).BoundingBox;
                lonCenter = mean(bbox(:,1));
                latCenter = mean(bbox(:,2));
                if lonCenter >= 80 && lonCenter <= 88 && ...
                        latCenter >= 26 && latCenter <= 30
                    selectedIdx = i;
                    break;
                end
            end
        end

        if isempty(selectedIdx)
            selectedIdx = 1;
        end
        nepalShape = nepalShape(selectedIdx);
    end

    if iscell(nepalShape.Lon)
        rings = nepalShape.Lon;
        latRings = nepalShape.Lat;
        for r = 1:numel(rings)
            if isempty(nepalBorderLon)
                nepalBorderLon = rings{r};
                nepalBorderLat = latRings{r};
            else
                nepalBorderLon = [nepalBorderLon,NaN,rings{r}]; %#ok<AGROW>
                nepalBorderLat = [nepalBorderLat,NaN,latRings{r}]; %#ok<AGROW>
            end
        end
    else
        nepalBorderLon = nepalShape.Lon;
        nepalBorderLat = nepalShape.Lat;
    end

    if isempty(nepalBorderLon) || isempty(nepalBorderLat)
        error('Extracted coordinates are empty.');
    end

    fprintf('Nepal border loaded successfully from %s\n',shpFile);
catch ME
    error('Failed to load shapefile: %s',ME.message);
end

% Keep the longest non-NaN ring as the mainland outline.
segBStarts = [1,find(isnan(nepalBorderLon))+1];
segBEnds   = [find(isnan(nepalBorderLon))-1,length(nepalBorderLon)];
segmentLengths = segBEnds-segBStarts+1;
[~,mainIdx] = max(segmentLengths);
borderLonMain = nepalBorderLon(segBStarts(mainIdx):segBEnds(mainIdx));
borderLatMain = nepalBorderLat(segBStarts(mainIdx):segBEnds(mainIdx));

% Remove adjacent duplicate points before polygon construction.
tol = 1e-6;
keepBorder = [true,sqrt(diff(borderLonMain).^2 + diff(borderLatMain).^2) > tol];
borderLonMain = borderLonMain(keepBorder);
borderLatMain = borderLatMain(keepBorder);

warnBorder = warning('off','MATLAB:polyshape:repairedBySimplify');
nepalPoly = polyshape(borderLonMain,borderLatMain,'Simplify',false);
warning(warnBorder);

meanLatNepal = mean([min(borderLatMain);max(borderLatMain)]);
kmPerDegLonMean = kmPerDegLat*cosd(meanLatNepal);
nepalAreaKm2 = area(nepalPoly)*kmPerDegLat*kmPerDegLonMean;

%% =========================================================
% C. OPTIMIZATION SETUP
%% =========================================================
% The optimization variables are x = [i, RAAN, nu]. The expensive orbit
% simulation uses a 60 s sampling interval; the final best-orbit analysis
% is rerun at 2 s to recover the high-resolution output and plots.

fixed = struct();
fixed.semiMajorAxis = semiMajorAxis;
fixed.eccentricity = eccentricity;
fixed.theta_half = theta_half;
fixed.tanHalfAngle = tanHalfAngle;
fixed.gsLat = gsLat;
fixed.gsLon = gsLon;
fixed.groundStations = groundStations;
fixed.accuracy = accuracy;
fixed.startTime = startTime;
fixed.stopTime = stopTime;
fixed.sampleTime = optimizationSampleTime;
fixed.nepalPoly = nepalPoly;
fixed.nepalAreaKm2 = nepalAreaKm2;
fixed.degPerKmLat = 1/kmPerDegLat;
fixed.kmPerDegLat = kmPerDegLat;
fixed.kmPerDegLonMean = kmPerDegLonMean;
fixed.circleCos = circleCos;
fixed.circleSin = circleSin;
fixed.nepalBuffered = polybuffer(nepalPoly, ...
    perigeeAltitude/1000*tanHalfAngle*(1/kmPerDegLat));
fixed.lonMin = 60;
fixed.lonMax = 105;
fixed.latMin = 10;
fixed.latMax = 45;
fixed.glacierNames = glacierNames;
fixed.glacierLat = glacierLat;
fixed.glacierLon = glacierLon;
fixed.glacierElevM = glacierElevM;
fixed.glacierAreaKm2 = glacierAreaKm2;
fixed.constellation = constellation;

lb = [0,0,0];
ub = [90,360,360];
varNames = {'inclination','raan','trueAnomaly'};

% Weighted objective. Positive weights reward a metric; negative weights
% penalize it after normalization by the corresponding scale.
scoreWeights = struct( ...
    'coverageFraction',0.35, ...
    'passesPerDay',0.20, ...
    'sgsTotalMin',0.15, ...
    'avgRevisitMin',-0.20, ...
    'avgLatencyMin',-0.10, ...
    'numLakesCovered',0.20, ...
    'worstRevisitMin',-0.15, ...
    'minPassesPerDayLakes',0.10);

scoreScales = struct( ...
    'coverageFraction',1.0, ...
    'passesPerDay',6, ...
    'sgsTotalMin',60, ...
    'avgRevisitMin',180, ...
    'avgLatencyMin',120, ...
    'numLakesCovered',6, ...
    'worstRevisitMin',720, ...
    'minPassesPerDayLakes',3);

fprintf('\n================ OPTIMIZATION SETUP ================\n');
fprintf('Search bounds:\n');
for v = 1:numel(varNames)
    fprintf('  %-14s [%.1f, %.1f] deg\n',varNames{v},lb(v),ub(v));
end
fprintf('Constellation configuration: %d orbit(s) x %d satellite(s)/orbit\n', ...
    constellation.numOrbits,constellation.satsPerOrbit);
fprintf('====================================================\n\n');

%% =========================================================
% D. INITIAL LATIN HYPERCUBE SAMPLING
%% =========================================================
% LHS spreads the first expensive evaluations across the bounded search
% volume rather than clustering them randomly.

rng(42,'twister');
nInitial = 100;
lhsUnit = lhsdesign(nInitial,3);
X0 = lb + lhsUnit.*(ub-lb);

y0 = nan(nInitial,1);
metrics0 = cell(nInitial,1);

fprintf('Evaluating %d initial Latin Hypercube samples...\n',nInitial);

useParallel = ~isempty(ver('parallel'));
if useParallel
    try
        if isempty(gcp('nocreate'))
            parpool;
        end
    catch
        useParallel = false;
    end
end

if useParallel
    parfor i = 1:nInitial
        m = runSim(X0(i,:),fixed);
        y0(i) = computeScore(m,scoreWeights,scoreScales);
        metrics0{i} = m;
    end
else
    for i = 1:nInitial
        m = runSim(X0(i,:),fixed);
        y0(i) = computeScore(m,scoreWeights,scoreScales);
        metrics0{i} = m;
        if mod(i,10) == 0
            fprintf('  ... %d / %d initial samples done\n',i,nInitial);
        end
    end
end

fprintf('Initial sampling complete. Best initial score: %.4f\n\n',max(y0));

%% =========================================================
% E. SURROGATE MODEL + BAYESIAN OPTIMIZATION LOOP
%% =========================================================
% GP model:
%   f(x) ~ GP(m(x),k(x,x'))
% with ARD squared-exponential kernel. Expected Improvement is
%   EI(x) = (mu-fbest)Phi(z) + sigma*phi(z),
%   z = (mu-fbest)/sigma.
%
% The candidate set uses LHS; fmincon locally refines its best EI point.

X = zeros(nInitial+100,3);
y = nan(nInitial+100,1);
X(1:nInitial,:) = X0;
y(1:nInitial) = y0;
allMetrics = cell(nInitial+100,1);
allMetrics(1:nInitial) = metrics0;

maxEvals = 100;
eiStallLimit = 10;
eiStallFrac = 0.005;
stallCount = 0;
nCand = 2000;
optsFmin = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',50);

bestScoreHistory = nan(maxEvals,1);
fprintf('%-6s %-12s %-8s %-14s %-10s %-18s %-10s\n', ...
    'Iter','inclination','RAAN','true anomaly','Score','Best score so far','EI');

nEval = nInitial;
for iter = 1:maxEvals

    Xfit = X(1:nEval,:);
    yfit = y(1:nEval);

    gprMdl = fitrgp(Xfit,yfit, ...
        'KernelFunction','ardsquaredexponential', ...
        'Standardize',true, ...
        'BasisFunction','constant');

    bestSoFar = max(yfit);

    candUnit = lhsdesign(nCand,3);
    Xcand = lb + candUnit.*(ub-lb);
    [muCand,sdCand] = predict(gprMdl,Xcand);
    eiCand = expectedImprovement(muCand,sdCand,bestSoFar);

    [~,bestCandIdx] = max(eiCand);
    xNext0 = Xcand(bestCandIdx,:);

    negEIfun = @(x) negativeExpectedImprovement(gprMdl,x,bestSoFar);
    try
        xNext = fmincon(negEIfun,xNext0,[],[],[],[],lb,ub,[],optsFmin);
    catch
        xNext = xNext0;
    end
    xNext = min(max(xNext,lb),ub);

    [muNext,sdNext] = predict(gprMdl,xNext);
    eiNext = expectedImprovement(muNext,sdNext,bestSoFar);

    mNext = runSim(xNext,fixed);
    yNext = computeScore(mNext,scoreWeights,scoreScales);

    nEval = nEval + 1;
    X(nEval,:) = xNext;
    y(nEval) = yNext;
    allMetrics{nEval} = mNext;

    bestScoreHistory(iter) = max(y(1:nEval));

    fprintf('%-6d %-12.3f %-8.3f %-14.3f %-10.4f %-18.4f %-10.6f\n', ...
        iter,xNext(1),xNext(2),xNext(3),yNext,bestScoreHistory(iter),eiNext);

    if eiNext < eiStallFrac*max(abs(bestSoFar),eps)
        stallCount = stallCount + 1;
    else
        stallCount = 0;
    end

    if stallCount >= eiStallLimit
        fprintf(['\nEarly stopping: EI stayed below %.2f%% of best score ' ...
            'for %d consecutive iterations.\n'], ...
            eiStallFrac*100,eiStallLimit);
        break;
    end
end

% Trim preallocation to evaluated rows only.
X = X(1:nEval,:);
y = y(1:nEval);
allMetrics = allMetrics(1:nEval);

fprintf('\nOptimization loop finished after %d additional evaluations.\n',iter);
fprintf('Total evaluations (initial + optimization): %d\n\n',numel(y));

%% =========================================================
% F. POST-PROCESSING
%% =========================================================
% Sort all evaluated candidates by score, select the maximum, then report
% the best parameter set and the best observed value of every metric.

[sortedScores, sortIdx] = sort(y, 'descend');
Xsorted = X(sortIdx, :);
metricsSorted = allMetrics(sortIdx);

bestIdx = sortIdx(1);
bestParams  = X(bestIdx, :);
bestScore   = y(bestIdx);
bestMetrics = allMetrics{bestIdx};

fprintf('================ BEST OVERALL ORBIT ================\n');
fprintf('  Inclination:   %.4f deg\n', bestParams(1));
fprintf('  RAAN:          %.4f deg\n', bestParams(2));
fprintf('  True anomaly:  %.4f deg\n', bestParams(3));
fprintf('  Overall score: %.4f\n', bestScore);
fprintf('======================================================\n\n');

metricNames = {'coverageFraction','passesPerDay','sgsTotalMin', ...
    'avgRevisitMin','avgLatencyMin','numLakesCovered','worstRevisitMin', ...
    'minPassesPerDayLakes'};
fprintf('================ BEST PER METRIC ================\n');
for m = 1:numel(metricNames)
    field = metricNames{m};
    vals = cellfun(@(s) s.(field), allMetrics);
    if strcmp(field, 'avgRevisitMin') || strcmp(field, 'avgLatencyMin') || ...
            strcmp(field, 'worstRevisitMin')
        vals(isnan(vals)) = Inf;
        if all(isinf(vals))
            [bestVal,idxM] = min(vals);
        else
            finiteIdx = find(isfinite(vals));
            [bestFinite,localIdx] = min(vals(finiteIdx));
            bestVal = bestFinite;
            idxM = finiteIdx(localIdx);
        end
    else
        vals(isnan(vals)) = -Inf;
        [bestVal, idxM] = max(vals);
    end
    fprintf('  %-18s = %10.4f  at  incl=%.2f, raan=%.2f, nu=%.2f\n', ...
        field, bestVal, X(idxM,1), X(idxM,2), X(idxM,3));
end
fprintf('===================================================\n\n');

%% =========================================================
% G. SAVE RESULTS
%% =========================================================
% Save the complete evaluated design matrix, objective values, metric
% structs, best solution, weights/scales, and constellation configuration.

optResults = struct();
optResults.varNames      = varNames;
optResults.lb            = lb;
optResults.ub            = ub;
optResults.X             = X;
optResults.y             = y;
optResults.allMetrics    = {allMetrics{:}}; %#ok<CCAT1>
optResults.sortedScores  = sortedScores;
optResults.Xsorted       = Xsorted;
optResults.metricsSorted = {metricsSorted{:}}; %#ok<CCAT1>
optResults.bestParams    = bestParams;
optResults.bestScore     = bestScore;
optResults.bestMetrics   = bestMetrics;
optResults.scoreWeights  = scoreWeights;
optResults.scoreScales   = scoreScales;
optResults.constellation = constellation;
optResults.groundStations = groundStations;
optResults.accuracy = accuracy;

saveFile = sprintf('cubesat_orbit_optimization_%s.mat', ...
    datestr(now, 'yyyymmdd_HHMMSS'));
save(saveFile, 'optResults');
fprintf('Saved optimization results to: %s\n\n', saveFile);

%% =========================================================
% H. VISUALIZE THE BEST OVERALL ORBIT (HIGH-RESOLUTION)
%% =========================================================
% Process: rebuild every configured satellite and ground station at the best
% solution, collect all tracks/access intervals, form the constellation union,
% and reuse those exact geometries for the static map, metrics, and H2.
% With the present 1 x 1 x 1 configuration, this is numerically identical to
% the original single-satellite result.

bestInclination = bestParams(1);
bestRaan = bestParams(2);
bestTrueAnomaly = bestParams(3);

scBest = satelliteScenario(startTime,stopTime,baseSampleTime);
constellationElements = expandConstellationElements(bestParams,constellation);
nSatBest = numel(constellationElements);

satBest = cell(1,nSatBest);
for si = 1:nSatBest
    e = constellationElements(si);
    satName = sprintf('NepalSat_P%02d_S%02d',e.orbitIndex,e.satIndex);
    satBest{si} = satellite(scBest,semiMajorAxis,eccentricity, ...
        e.inclination,e.raan,0,e.trueAnomaly,Name=satName);
    satBest{si}.MarkerColor = [0 1 0];
    satBest{si}.MarkerSize = 8;
    satBest{si}.Orbit.LineWidth = 2;
end

gsBest = cell(1,numel(groundStations));
for gi = 1:numel(groundStations)
    gsBest{gi} = groundStation(scBest,groundStations(gi).lat, ...
        groundStations(gi).lon,Name=groundStations(gi).name);
end

% Add the six GLOF-vulnerable lakes as visual markers in the 3-D viewer.
% These are groundStation objects used purely for display; no access()
% links are created for them, so they do not affect SGS contact metrics.
glacierMarkersBest = cell(1,numel(glacierNames));
for gi = 1:numel(glacierNames)
    glacierMarkersBest{gi} = groundStation(scBest, ...
        glacierLat(gi), glacierLon(gi), ...
        Name=glacierNames{gi});
    glacierMarkersBest{gi}.MarkerColor = [0 1 1];   % cyan
    glacierMarkersBest{gi}.MarkerSize  = 6;
end

% Build all satellite-ground-station access links. In 1x1 there is exactly
% one access object, preserving the current 3-D view and contact calculation.
accBest = cell(nSatBest,numel(gsBest));
for si = 1:nSatBest
    for gi = 1:numel(gsBest)
        accBest{si,gi} = access(satBest{si},gsBest{gi});
        accBest{si,gi}.LineColor = [0 1 0];
        accBest{si,gi}.LineWidth = 2;
    end
end
accBestGS = accBest{1,1};

% Create a sensor for each satellite; only the first is used for the current
% 3-D field-of-view display, while all satellites remain available for scale-up.
sensorBest = cell(1,nSatBest);
for si = 1:nSatBest
    sensorBest{si} = conicalSensor(satBest{si},'MaxViewAngle',sensorHalfAngle, ...
        'Name',sprintf('Best Orbit Swath P%02d S%02d', ...
        constellationElements(si).orbitIndex,constellationElements(si).satIndex));
end

viewerBest = satelliteScenarioViewer(scBest,'Dimension','3D', ...
    'Name','Best Orbit - Constant Footprint Viewer');
for si = 1:nSatBest
    fieldOfView(sensorBest{si});
end
playScenarioSafely(scBest,sensorBest,'Best Orbit viewer');

% Collect high-resolution geographic states for every satellite.
bestTracks = cell(1,nSatBest);
for si = 1:nSatBest
    [posBest,~,~] = states(satBest{si},'CoordinateFrame','geographic');
    latRaw = posBest(1,:);
    lonRaw = posBest(2,:);
    altRaw = posBest(3,:)/1000;
    valid = isfinite(latRaw) & isfinite(lonRaw) & isfinite(altRaw);
    latRaw = latRaw(valid);
    lonRaw = lonRaw(valid);
    altRaw = altRaw(valid);
    [latSeg,lonSeg,altSeg] = segmentGroundTrack(latRaw,lonRaw,altRaw,10);
    bestTracks{si} = struct('lat',latRaw,'lon',lonRaw,'alt',altRaw, ...
        'latSeg',latSeg,'lonSeg',lonSeg,'altSeg',altSeg);
end

% Current 1x1 aliases preserve the established H2 variable names.
latBestRaw = bestTracks{1}.lat;
lonBestRaw = bestTracks{1}.lon;
altBestRaw = bestTracks{1}.alt;
latBestSeg = bestTracks{1}.latSeg;
lonBestSeg = bestTracks{1}.lonSeg;
altBestSeg = bestTracks{1}.altSeg;

% Authoritative constellation swath geometry: union all satellite swaths.
satSwathPolysBest = cell(1,nSatBest);
planeSwathPolysBest = cell(1,constellation.numOrbits);
for oi = 1:constellation.numOrbits
    planeSwathPolysBest{oi} = polyshape();
end
unionBest = polyshape();

for si = 1:nSatBest
    satPoly = computeSwathPolygon(bestTracks{si}.latSeg,bestTracks{si}.lonSeg, ...
        bestTracks{si}.altSeg,fixed);
    satSwathPolysBest{si} = satPoly;
    if ~isempty(satPoly.Vertices)
        unionBest = union(unionBest,satPoly);
        oi = constellationElements(si).orbitIndex;
        planeSwathPolysBest{oi} = union(planeSwathPolysBest{oi},satPoly);
    end
end

if isempty(unionBest.Vertices)
    areaKm2Best = 0;
else
    areaKm2Best = area(unionBest)*kmPerDegLat*kmPerDegLonMean;
end

% Primary-satellite passes remain the animation input for the current 1x1
% case. All constellation members are nevertheless already included in the
% metrics and union calculations above.
[passShapesIdx,passAreasIdx,passIndexRanges] = computePerPassSwathShapesIdx( ...
    latBestSeg,lonBestSeg,altBestSeg,nepalPoly,fixed);
perPassAreasBest = passAreasIdx;

fprintf('\nContinuous swath coverage inside Nepal (best orbit/constellation): %.2f km²\n',areaKm2Best);

%% ---------------------------------------------------------
% NUMERICAL AREA SENSITIVITY / ACCURACY BAND
%% ---------------------------------------------------------
if accuracy.enabled
    referenceAreaKm2 = computeConstellationSwathArea(bestTracks,fixed, ...
        accuracy.referenceCirclePts);
    accuracyMarginKm2 = abs(referenceAreaKm2-areaKm2Best);
    accuracyLowerKm2 = areaKm2Best-accuracyMarginKm2;
    accuracyUpperKm2 = areaKm2Best+accuracyMarginKm2;
    accuracyPct = 100*accuracyMarginKm2/max(abs(areaKm2Best),eps);
else
    referenceAreaKm2 = areaKm2Best;
    accuracyMarginKm2 = 0;
    accuracyLowerKm2 = areaKm2Best;
    accuracyUpperKm2 = areaKm2Best;
    accuracyPct = 0;
end

fprintf('Numerical geometry reference area (%d-gon): %.2f km²\n', ...
    accuracy.referenceCirclePts,referenceAreaKm2);
fprintf('Estimated numerical area margin:              ±%.2f km² (%.3f%%)\n', ...
    accuracyMarginKm2,accuracyPct);
fprintf('Estimated numerical area range:               %.2f to %.2f km²\n', ...
    accuracyLowerKm2,accuracyUpperKm2);

%% ---------------------------------------------------------
% STATIC WORLD 2-D MAP
%% ---------------------------------------------------------
figure('Name','Best Orbit - Full World Ground Track with Continuous Swath', ...
    'Position',[100,100,1200,700]);
set(gcf,'Color',[0.1 0.1 0.1]);
axBest = axes;
hold(axBest,'on');
set(axBest,'Color',[0 0 0]);
grid(axBest,'on');
axBest.GridColor = [0.3 0.3 0.3];
axBest.XColor = 'w';
axBest.YColor = 'w';
box(axBest,'on');
load coastlines;
plot(axBest,coastlon,coastlat,'Color',[0.6 0.6 0.6], ...
    'LineWidth',0.5,'DisplayName','Coastlines');
plot(axBest,nepalBorderLon,nepalBorderLat,'Color',[0.8 0.6 0.4], ...
    'LineWidth',1.5,'DisplayName','Nepal Border');

for si = 1:nSatBest
    if si == 1
        dispName = 'Best Orbit Ground Track';
    else
        dispName = sprintf('Ground Track P%02d S%02d', ...
            constellationElements(si).orbitIndex,constellationElements(si).satIndex);
    end
    plot(axBest,bestTracks{si}.lonSeg,bestTracks{si}.latSeg, ...
        'Color',[0.3 0.85 0.3],'LineWidth',1.0,'DisplayName',dispName);
end

if ~isempty(unionBest.Vertices)
    plot(axBest,unionBest,'FaceColor',[1 1 1],'FaceAlpha',0.35, ...
        'EdgeColor','none','DisplayName','Swath Inside Nepal');
end

for gi = 1:numel(groundStations)
    plot(axBest,groundStations(gi).lon,groundStations(gi).lat,'p', ...
        'MarkerSize',12,'MarkerFaceColor','yellow','MarkerEdgeColor','red', ...
        'DisplayName',groundStations(gi).name);
end

plot(axBest,glacierLon,glacierLat,'d','MarkerSize',8, ...
    'MarkerFaceColor','cyan','MarkerEdgeColor','k','LineStyle','none', ...
    'DisplayName','GLOF Lakes');
for gi = 1:numel(glacierLat)
    text(axBest,glacierLon(gi)+0.15,glacierLat(gi),glacierNames{gi}, ...
        'Color','w','FontSize',7,'Interpreter','none','HandleVisibility','off');
end
xlabel(axBest,'Longitude (°)','Color','w');
ylabel(axBest,'Latitude (°)','Color','w');
title(axBest,sprintf('Best Orbit: incl=%.2f°, RAAN=%.2f°, \\nu_0=%.2f°', ...
    bestInclination,bestRaan,bestTrueAnomaly),'Color','w');
xlim(axBest,[-180 180]);
ylim(axBest,[-90 90]);
legend(axBest,'Location','southoutside','Orientation','horizontal', ...
    'TextColor','w','Color',[0.1 0.1 0.1]);

% Constellation-wide final metrics combine all satellite tracks and all GS
% access links. With one satellite/one GS this reduces exactly to the current
% single-object metrics.
[bestHighResMetrics,constellationDetails] = computeConstellationMetrics( ...
    bestTracks,accBest,groundStations,startTime,baseSampleTime, ...
    nepalAreaKm2,areaKm2Best,fixed.nepalBuffered,constellationElements);

fprintf('\n========= COVERAGE METRICS OVER NEPAL (BEST ORBIT, HIGH-RES) =========\n');
printMetricsBlock(bestHighResMetrics,perPassAreasBest);
fprintf('Unique union accuracy range: %.2f ± %.2f km²  (%.2f to %.2f km²)\n', ...
    areaKm2Best,accuracyMarginKm2,accuracyLowerKm2,accuracyUpperKm2);
fprintf('========================================================================\n\n');

[bestGlacierMetrics,glacierHitBySat] = computeConstellationGlacierMetrics( ...
    bestTracks,glacierNames,glacierLat,glacierLon,theta_half, ...
    baseSampleTime,fixed);
%% =========================================================
% H2. PER-PASS SWATH ADDITION & OVERLAP — ANIMATED VERSION
%% =========================================================
% Reuses precomputed pass polygons. Active passes are rendered as polyshape
% objects so multi-boundary geometries remain valid during animation.
% Growth state is the cumulative union of instantaneous footprints; completed
% passes remain frozen. Final overlap is computed from pairwise intersections.

nPassesAnim = numel(passShapesIdx);

if nPassesAnim == 0

    warning('No qualifying passes found — skipping the H2 animation.');

else

    %% ---------------------------------------------------------
    % PASS COLORS
    %% ---------------------------------------------------------
    passColorTable = [ ...
        0.95 0.45 0.10;   % Pass 1 - orange
        0.15 0.65 0.95;   % Pass 2 - sky blue
        0.85 0.20 0.65;   % Pass 3 - magenta
        0.40 0.85 0.30;   % Pass 4 - green
        0.95 0.85 0.20];  % Pass 5 - yellow

    if nPassesAnim > size(passColorTable,1)
        passColorTable = [ ...
            passColorTable;
            lines(nPassesAnim - size(passColorTable,1))];
    end

    passColors = passColorTable(1:nPassesAnim, :);

    kmPerDegLatA = kmPerDegLat;
    kmPerDegLonA = kmPerDegLonMean;

    %% ---------------------------------------------------------
    % PASS NEAR-RANGES
    %
    % These are used only to decide when a pass is being animated.
    %% ---------------------------------------------------------
    nTotal = numel(latBestSeg);

    passNearRange = zeros(nPassesAnim,2);

    for p = 1:nPassesAnim

        r = passIndexRanges(p,:);
        seg = r(1):r(2);

        inBox = ...
            lonBestSeg(seg) >= fixed.lonMin & ...
            lonBestSeg(seg) <= fixed.lonMax & ...
            latBestSeg(seg) >= fixed.latMin & ...
            latBestSeg(seg) <= fixed.latMax & ...
            ~isnan(lonBestSeg(seg));

        firstNR = seg(find(inBox,1,'first'));
        lastNR  = seg(find(inBox,1,'last'));

        if isempty(firstNR)
            firstNR = r(1);
        end

        if isempty(lastNR)
            lastNR = r(2);
        end

        passNearRange(p,:) = [firstNR,lastNR];
    end

    %% ---------------------------------------------------------
    % PASS FLAG FOR EACH TRACK SAMPLE
    %% ---------------------------------------------------------
    inPassFlag = false(1,nTotal);
    passOfIdx  = zeros(1,nTotal);

    for p = 1:nPassesAnim

        r = passNearRange(p,:);

        if r(2) >= r(1)
            inPassFlag(r(1):r(2)) = true;
            passOfIdx(r(1):r(2)) = p;
        end
    end

    %% ---------------------------------------------------------
    % FIGURE / BASE MAP
    %% ---------------------------------------------------------
    figAnim = figure( ...
        'Name','Best Orbit - Animated Per-Pass Swath Build-up', ...
        'Position',[200,100,1200,750], ...
        'Interruptible','on', ...
        'BusyAction','queue');

    set(figAnim,'Color',[0.1 0.1 0.1]);

    axAnim = axes(figAnim);
    hold(axAnim,'on');

    set(axAnim,'Color',[0 0 0]);
    grid(axAnim,'on');

    axAnim.GridColor = [0.3 0.3 0.3];
    axAnim.XColor = 'w';
    axAnim.YColor = 'w';

    box(axAnim,'on');

    xlim(axAnim,[-180 180]);
    ylim(axAnim,[-90 90]);

    axAnim.XLimMode = 'manual';
    axAnim.YLimMode = 'manual';

    xlabel(axAnim,'Longitude (°)','Color','w');
    ylabel(axAnim,'Latitude (°)','Color','w');

    load coastlines;

    plot(axAnim,coastlon,coastlat, ...
        'Color',[0.6 0.6 0.6], ...
        'LineWidth',0.5, ...
        'DisplayName','Coastlines');

    plot(axAnim,nepalBorderLon,nepalBorderLat, ...
        'Color',[0.8 0.6 0.4], ...
        'LineWidth',1.5, ...
        'DisplayName','Nepal Border');

    for gi = 1:numel(groundStations)
        plot(axAnim,groundStations(gi).lon,groundStations(gi).lat,'p', ...
            'MarkerSize',12, ...
            'MarkerFaceColor','yellow', ...
            'MarkerEdgeColor','red', ...
            'DisplayName',groundStations(gi).name);
    end

    %% ---------------------------------------------------------
    % GLOF LAKES
    %% ---------------------------------------------------------
    plot(axAnim,glacierLon,glacierLat,'d', ...
        'MarkerSize',8, ...
        'MarkerFaceColor','cyan', ...
        'MarkerEdgeColor','k', ...
        'LineStyle','none', ...
        'DisplayName','GLOF Lakes');

    for gi = 1:numel(glacierLat)

        text(axAnim, ...
            glacierLon(gi)+0.15, ...
            glacierLat(gi), ...
            glacierNames{gi}, ...
            'Color','w', ...
            'FontSize',7, ...
            'Interpreter','none', ...
            'HandleVisibility','off');

    end

    %% ---------------------------------------------------------
    % TRACK + SATELLITE
    %% ---------------------------------------------------------
    titleAnim = title(axAnim, ...
        'Best Orbit Ground Track  |  t = 00:00:00', ...
        'Color','w');

    hTrack = plot(axAnim,NaN,NaN, ...
        'Color',[0.3 0.85 0.3], ...
        'LineWidth',1.2, ...
        'DisplayName','Ground Track');

    hSatMarker = plot(axAnim,NaN,NaN,'o', ...
        'MarkerSize',7, ...
        'MarkerFaceColor',[0.3 1 0.3], ...
        'MarkerEdgeColor','w', ...
        'HandleVisibility','off');

    %% ---------------------------------------------------------
    % IMPORTANT:
    % DO NOT create fill(NaN,NaN,...) patches here.
    %
    % hPassFill(p) is created only when a real polyshape exists.
    %% ---------------------------------------------------------
    hPassFill = gobjects(1,nPassesAnim);

    for p = 1:nPassesAnim
        hPassFill(p) = gobjects(1);
    end

    passLabelText = text(axAnim,0.02,0.06,'', ...
        'Units','normalized', ...
        'Color','w', ...
        'FontSize',11, ...
        'FontWeight','bold', ...
        'VerticalAlignment','bottom');

    %% ---------------------------------------------------------
    % INITIAL LEGEND
    %% ---------------------------------------------------------
    legend(axAnim, ...
        'Location','southoutside', ...
        'Orientation','horizontal', ...
        'NumColumns',4, ...
        'TextColor','w', ...
        'Color',[0.1 0.1 0.1]);

    %% ---------------------------------------------------------
    % ANIMATION TIMING
    %% ---------------------------------------------------------
    strideOutside = 80;
    strideInside  = 5;
    framePause    = 0.01;

    %% ---------------------------------------------------------
    % TRACK SEGMENTS
    %% ---------------------------------------------------------
    segStartsAll = [1, find(isnan(latBestSeg))+1];
    segEndsAll   = [find(isnan(latBestSeg))-1, length(latBestSeg)];

    %% ---------------------------------------------------------
    % BUILD FRAME INDEX LIST
    %% ---------------------------------------------------------
    frameIdxList = [];

    i = 1;

    while i <= nTotal

        frameIdxList(end+1) = i; %#ok<AGROW>

        if inPassFlag(i)
            step = strideInside;
        else
            step = strideOutside;
        end

        i = i + step;

    end

    % Force dense sampling through every pass.
    forcedIdx = [segStartsAll, segEndsAll];

    for p = 1:nPassesAnim

        r = passNearRange(p,:);

        forcedIdx = [ ...
            forcedIdx, ...
            r(1):strideInside:r(2), ...
            r(2)]; %#ok<AGROW>

    end

    forcedIdx = unique(forcedIdx);

    frameIdxList = unique([frameIdxList,forcedIdx]);

    frameIdxList = frameIdxList( ...
        frameIdxList >= 1 & frameIdxList <= nTotal);

    frameIdxList = sort(frameIdxList);

    %% ---------------------------------------------------------
    % TIME INDEX
    %% ---------------------------------------------------------
    validMaskSeg = ~isnan(latBestSeg);
    cumValidCount = cumsum(validMaskSeg);

    %% ---------------------------------------------------------
    % PRECOMPUTE GROWING POLYGONS
    %% ---------------------------------------------------------
    nFrames = numel(frameIdxList);
    warnStateH2 = warning('off','MATLAB:polyshape:repairedBySimplify');

    growingPolyByFrame = cell(nFrames,1);

    currentAccum = polyshape();
    prevPassIdx = 0;

    for f = 1:nFrames

        idx = frameIdxList(f);

        % Outside all passes.
        if ~inPassFlag(idx)

            prevPassIdx = 0;
            continue;

        end

        p = passOfIdx(idx);

        % New pass: start accumulating from zero.
        if p ~= prevPassIdx

            currentAccum = polyshape();
            prevPassIdx = p;

        end

        % Add the instantaneous footprint.
        if ...
                ~isnan(latBestSeg(idx)) && ...
                ~isnan(lonBestSeg(idx)) && ...
                lonBestSeg(idx) >= fixed.lonMin && ...
                lonBestSeg(idx) <= fixed.lonMax && ...
                latBestSeg(idx) >= fixed.latMin && ...
                latBestSeg(idx) <= fixed.latMax

            R_local_km = altBestSeg(idx)*tanHalfAngle;

            rLat = R_local_km*fixed.degPerKmLat;
            rLon = R_local_km/(fixed.kmPerDegLat*cosd(latBestSeg(idx)));

            circLon = lonBestSeg(idx) + rLon*fixed.circleCos;
            circLat = latBestSeg(idx) + rLat*fixed.circleSin;

            newCirc = polyshape( ...
                circLon, ...
                circLat, ...
                'Simplify',false);

            currentAccum = union(currentAccum,newCirc);

        end

        % Clip to Nepal.
        if ~isempty(currentAccum.Vertices)

            growingPolyByFrame{f} = ...
                intersect(currentAccum,nepalPoly);

        else

            growingPolyByFrame{f} = polyshape();

        end

    end

    %% ---------------------------------------------------------
    % PASS TIME RANGES
    %% ---------------------------------------------------------
    passTimeRanges = nan(nPassesAnim,2);

    for p = 1:nPassesAnim

        segIdx = passIndexRanges(p,1):passIndexRanges(p,2);

        hits = false(size(segIdx));

        for kk = 1:numel(segIdx)

            k = segIdx(kk);

            if isnan(latBestSeg(k)) || isnan(lonBestSeg(k))
                continue;
            end

            R_local_km = altBestSeg(k)*tanHalfAngle;

            rLat = R_local_km*fixed.degPerKmLat;
            rLon = R_local_km/(fixed.kmPerDegLat*cosd(latBestSeg(k)));

            circLon = lonBestSeg(k) + rLon*fixed.circleCos;
            circLat = latBestSeg(k) + rLat*fixed.circleSin;

            fp = polyshape( ...
                circLon, ...
                circLat, ...
                'Simplify',false);

            tmpIntersection = intersect(fp,nepalPoly);

            if ~isempty(tmpIntersection.Vertices)
                hits(kk) = true;
            end

        end

        hitIdx = segIdx(hits);

        if ~isempty(hitIdx)

            passTimeRanges(p,1) = ...
                (cumValidCount(hitIdx(1))-1) * baseSampleTime;

            passTimeRanges(p,2) = ...
                (cumValidCount(hitIdx(end))-1) * baseSampleTime;

        end

    end

    warning(warnStateH2);

    %% ---------------------------------------------------------
    % ANIMATION STATE
    %% ---------------------------------------------------------
    completedPasses = false(1,nPassesAnim);

    %% ---------------------------------------------------------
    % MAIN ANIMATION LOOP
    %% ---------------------------------------------------------
    for f = 1:nFrames

        idx = frameIdxList(f);

        %% Ground track
        set(hTrack, ...
            'XData',lonBestSeg(1:idx), ...
            'YData',latBestSeg(1:idx));

        %% Satellite marker
        if ~isnan(lonBestSeg(idx)) && ...
                ~isnan(latBestSeg(idx))

            set(hSatMarker, ...
                'XData',lonBestSeg(idx), ...
                'YData',latBestSeg(idx));

        end

        %% Elapsed simulation time
        elapsedSec = ...
            (cumValidCount(idx)-1)*baseSampleTime;

        hh = floor(elapsedSec/3600);
        mm = floor(mod(elapsedSec,3600)/60);
        ss = floor(mod(elapsedSec,60));

        %% -----------------------------------------------------
        % ACTIVE PASS
        %% -----------------------------------------------------
        if inPassFlag(idx)

            p = passOfIdx(idx);

            growingInside = growingPolyByFrame{f};

            %% -----------------------------------------------
            % ROBUST POLYSHAPE RENDERING
            %
            % Delete old patch and plot the polyshape directly.
            % This handles multiple boundaries / NaN-separated
            % regions correctly.
            %% -----------------------------------------------
            if ~isempty(growingInside) && ...
                    ~isempty(growingInside.Vertices)

                if isgraphics(hPassFill(p))
                    delete(hPassFill(p));
                end

                hPassFill(p) = plot( ...
                    axAnim, ...
                    growingInside, ...
                    'FaceColor',passColors(p,:), ...
                    'FaceAlpha',0.55, ...
                    'EdgeColor','none', ...
                    'DisplayName',sprintf( ...
                        'Pass %d (growing...)',p));

            end

            set(titleAnim,'String',sprintf( ...
                'Flying Pass %d of %d  |  t = %02d:%02d:%02d', ...
                p,nPassesAnim,hh,mm,ss));

            set(passLabelText,'String',sprintf( ...
                'Pass %d in progress — %d / %d complete', ...
                p,sum(completedPasses),nPassesAnim));

            %% -----------------------------------------------
            % PASS COMPLETION
            %% -----------------------------------------------
            if ~completedPasses(p) && ...
                    idx >= passNearRange(p,2)

                finalShape = passShapesIdx{p};

                if isgraphics(hPassFill(p))
                    delete(hPassFill(p));
                end

                if ~isempty(finalShape.Vertices)

                    hPassFill(p) = plot( ...
                        axAnim, ...
                        finalShape, ...
                        'FaceColor',passColors(p,:), ...
                        'FaceAlpha',0.55, ...
                        'EdgeColor','none', ...
                        'DisplayName',sprintf( ...
                            'Pass %d (%.0f km²)', ...
                            p,passAreasIdx(p)));

                end

                completedPasses(p) = true;

            end

        else

            set(titleAnim,'String',sprintf( ...
                'Best Orbit Ground Track  |  t = %02d:%02d:%02d', ...
                hh,mm,ss));

            set(passLabelText,'String',sprintf( ...
                '%d / %d passes complete', ...
                sum(completedPasses),nPassesAnim));

        end

        drawnow;

        if framePause > 0
            pause(framePause);
        end

    end

    %% ---------------------------------------------------------
    % GUARANTEE ALL FINAL PASS POLYGONS ARE PRESENT
    %% ---------------------------------------------------------
    for p = 1:nPassesAnim

        if isgraphics(hPassFill(p))
            delete(hPassFill(p));
        end

        finalShape = passShapesIdx{p};

        displayName = sprintf( ...
            'Pass %d (%.0f km²)', ...
            p,passAreasIdx(p));

        if ~isempty(finalShape.Vertices)

            hPassFill(p) = plot( ...
                axAnim, ...
                finalShape, ...
                'FaceColor',passColors(p,:), ...
                'FaceAlpha',0.55, ...
                'EdgeColor','none', ...
                'DisplayName',displayName);

        end

    end

    completedPasses(:) = true;

    drawnow;

    %% =========================================================
    % FINAL FRAME — INTERSECTION OUTLINES ONLY
    %% =========================================================
    pairInfoList = {};

    pairStyles = {':','--','-.'};

    pairStyleIdx = 0;
    sumPairwise = 0;

    for a = 1:nPassesAnim

        for b = a+1:nPassesAnim

            pairStyleIdx = pairStyleIdx + 1;

            interAB = intersect( ...
                passShapesIdx{a}, ...
                passShapesIdx{b});

            if isempty(interAB.Vertices)
                continue;
            end

            interAB = intersect(interAB,nepalPoly);

            if isempty(interAB.Vertices)
                continue;
            end

            areaAB = area(interAB) * ...
                kmPerDegLatA * kmPerDegLonA;

            if areaAB < 1
                continue;
            end

            [xb,yb] = boundary(interAB);

            plot(axAnim,xb,yb, ...
                'Color',[1 1 1], ...
                'LineWidth',1.8, ...
                'LineStyle',pairStyles{ ...
                    mod(pairStyleIdx-1, ...
                    numel(pairStyles))+1}, ...
                'DisplayName',sprintf( ...
                    'P%d ∩ P%d (%.0f km²)', ...
                    a,b,areaAB));

            pairInfoList{end+1} = ...
                {a,b,areaAB}; %#ok<AGROW>

            sumPairwise = sumPairwise + areaAB;

        end

    end

    %% ---------------------------------------------------------
    % TRIPLE OVERLAP
    %% ---------------------------------------------------------
    tripleAreaKm2Anim = 0;

    if nPassesAnim >= 3

        triplePoly = passShapesIdx{1};

        for p = 2:nPassesAnim

            triplePoly = intersect( ...
                triplePoly, ...
                passShapesIdx{p});

        end

        if ~isempty(triplePoly.Vertices)

            triplePoly = intersect( ...
                triplePoly, ...
                nepalPoly);

            if ~isempty(triplePoly.Vertices)

                tripleAreaKm2Anim = ...
                    area(triplePoly) * ...
                    kmPerDegLatA * ...
                    kmPerDegLonA;

                if tripleAreaKm2Anim >= 1

                    [xb,yb] = boundary(triplePoly);

                    plot(axAnim,xb,yb, ...
                        'Color',[1 1 1], ...
                        'LineWidth',3.0, ...
                        'LineStyle','-', ...
                        'DisplayName',sprintf( ...
                            'P1∩P2∩P3 (%.0f km²)', ...
                            tripleAreaKm2Anim));

                end

            else

                tripleAreaKm2Anim = 0;

            end

        end

    end

    %% ---------------------------------------------------------
    % UNION AREA
    %% ---------------------------------------------------------
    unionAllPasses = polyshape();

    for p = 1:nPassesAnim

        unionAllPasses = union( ...
            unionAllPasses, ...
            passShapesIdx{p});

    end

    unionAreaKm2Anim = 0;

    if ~isempty(unionAllPasses.Vertices)

        unionAreaKm2Anim = ...
            area(unionAllPasses) * ...
            kmPerDegLatA * ...
            kmPerDegLonA;

    end

    sumPassAreas = sum(passAreasIdx);

    redundancyAlgebraic = ...
        sumPassAreas - unionAreaKm2Anim;

    %% ---------------------------------------------------------
    % AT-LEAST-TWO PASSES
    % A>=2 = union of every pairwise intersection. With triple overlap T,
    % algebraic redundancy satisfies S-U = A>=2 + T.
    %% ---------------------------------------------------------
    atLeast2Poly = polyshape();

    for a = 1:nPassesAnim

        for b = a+1:nPassesAnim

            atLeast2Poly = union( ...
                atLeast2Poly, ...
                intersect( ...
                    passShapesIdx{a}, ...
                    passShapesIdx{b}));

        end

    end

    atLeast2AreaKm2 = 0;

    if ~isempty(atLeast2Poly.Vertices)

        atLeast2AreaKm2 = ...
            area(atLeast2Poly) * ...
            kmPerDegLatA * ...
            kmPerDegLonA;

    end

    %% ---------------------------------------------------------
    % FINAL TITLE / LABEL / LEGEND
    %% ---------------------------------------------------------
    set(titleAnim,'String',sprintf([ ...
        'All %d passes complete  |  ' ...
        '\\Sigma passes = %.0f km²,  ' ...
        'Union = %.0f km²,  ' ...
        'S−U = %.0f km²'], ...
        nPassesAnim, ...
        sumPassAreas, ...
        unionAreaKm2Anim, ...
        redundancyAlgebraic));

    set(passLabelText,'String',sprintf( ...
        '%d / %d passes complete', ...
        nPassesAnim,nPassesAnim));

    legend(axAnim, ...
        'Location','southoutside', ...
        'Orientation','horizontal', ...
        'NumColumns',3, ...
        'TextColor','w', ...
        'Color',[0.1 0.1 0.1]);

    drawnow;

    %% ---------------------------------------------------------
    % PASS TIME RANGES
    %% ---------------------------------------------------------
    fprintf('\n========= PASS TIME RANGES (swath within Nepal) =========\n');

    for p = 1:nPassesAnim

        if isnan(passTimeRanges(p,1))

            fprintf( ...
                'Pass %d: swath never intersects Nepal border.\n', ...
                p);

        else

            fprintf( ...
                'Pass %d: %s  to  %s\n', ...
                p, ...
                formatHMS(passTimeRanges(p,1)), ...
                formatHMS(passTimeRanges(p,2)));

        end

    end

    fprintf('==========================================================\n');

    %% ---------------------------------------------------------
    % PER-GLACIER SUMMARY
    %% ---------------------------------------------------------
    numLakesCoveredBest = ...
        sum([bestGlacierMetrics.covered]);

    if all([bestGlacierMetrics.covered])

        worstRevisitBest = ...
            max([bestGlacierMetrics.meanRevisitMin]);

    else

        worstRevisitBest = NaN;

    end

    minPassesPerDayBest = ...
        min([bestGlacierMetrics.passesPerDay]);

    meanPassesPerDayBest = ...
        mean([bestGlacierMetrics.passesPerDay]);

    fprintf('\n========= PER-GLACIER COVERAGE SUMMARY (BEST ORBIT) =========\n');

    fprintf('%-16s %8s %14s %14s %12s %10s\n', ...
        'Lake','Passes','MeanRevisit','MaxRevisit', ...
        'Dwell(min)','Covered');

    for g = 1:numel(bestGlacierMetrics)

        gm = bestGlacierMetrics(g);

        if isnan(gm.meanRevisitMin)
            meanRevisitStr = 'N/A';
        else
            meanRevisitStr = sprintf( ...
                '%.2f',gm.meanRevisitMin);
        end

        if isnan(gm.maxRevisitMin)
            maxRevisitStr = 'N/A';
        else
            maxRevisitStr = sprintf( ...
                '%.2f',gm.maxRevisitMin);
        end

        fprintf( ...
            '%-16s %8d %14s %14s %12.2f %10s\n', ...
            gm.name, ...
            gm.passesPerDay, ...
            meanRevisitStr, ...
            maxRevisitStr, ...
            gm.totalDwellMin, ...
            mat2str(gm.covered));

    end

    fprintf('---------------------------------------------------------------\n');

    fprintf('Lakes covered:            %d / %d\n', ...
        numLakesCoveredBest, ...
        numel(bestGlacierMetrics));

    if isnan(worstRevisitBest)

        fprintf( ...
            'Worst revisit time:       N/A (at least one lake never visited)\n');

    else

        fprintf( ...
            'Worst revisit time:       %.2f min\n', ...
            worstRevisitBest);

    end

    fprintf( ...
        'Minimum passes/day:       %d\n', ...
        minPassesPerDayBest);

    fprintf( ...
        'Mean passes/day:          %.2f\n', ...
        meanPassesPerDayBest);

    fprintf('================================================================\n\n');

    %% ---------------------------------------------------------
    % OVERLAP SUMMARY
    %% ---------------------------------------------------------
    fprintf('\n========= ANIMATED FIGURE - OVERLAP SUMMARY =========\n');

    fprintf( ...
        'Sum of per-pass areas (S):               %.2f km^2\n', ...
        sumPassAreas);

    fprintf( ...
        'Union of all passes (U, measured):       %.2f km^2\n', ...
        unionAreaKm2Anim);

    fprintf( ...
        'S − U (algebraic redundancy):            %.2f km^2\n', ...
        redundancyAlgebraic);

    fprintf( ...
        'At-least-2-passes polygon (measured):    %.2f km^2\n', ...
        atLeast2AreaKm2);

    fprintf( ...
        'Triple-overlap area (all %d passes):      %.2f km^2\n', ...
        nPassesAnim, ...
        tripleAreaKm2Anim);

    if ~isempty(pairInfoList)

        fprintf('Pairwise intersection breakdown:\n');

        for k = 1:numel(pairInfoList)

            pInfo = pairInfoList{k};

            fprintf( ...
                '   P%d ∩ P%d = %.2f km^2\n', ...
                pInfo{1}, ...
                pInfo{2}, ...
                pInfo{3});

        end

    end

    fprintf('Consistency check:  S − U  ==  (at-least-2 area) + (triple area)\n');

    fprintf( ...
        '   S − U                        = %.2f km^2\n', ...
        redundancyAlgebraic);

    fprintf( ...
        '   at-least-2 + triple          = %.2f km^2\n', ...
        atLeast2AreaKm2 + tripleAreaKm2Anim);

    fprintf('=====================================================\n\n');

end % if nPassesAnim == 0

%% =========================================================
% CONSTELLATION / ACCURACY FINAL REPORT
%% =========================================================
% This block exposes the array-based results even though the current mission
% remains 1 plane x 1 satellite x 1 ground station. For future scaling, the
% same fields are populated for every satellite, plane, and ground station.

fprintf('\n========= CONSTELLATION BREAKDOWN =========\n');
for si = 1:numel(bestMetrics.satelliteMetrics)
    sm = bestMetrics.satelliteMetrics(si);
    fprintf('Satellite P%02d-S%02d: coverage=%.2f km², passes/day=%d, SGS contact=%.2f min\n', ...
        sm.orbitIndex,sm.satelliteIndex,sm.accumAreaKm2,sm.passesPerDay,sm.sgsTotalMin);
end
for oi = 1:numel(bestMetrics.planeMetrics)
    pm = bestMetrics.planeMetrics(oi);
    fprintf('Orbit/Plane %02d: coverage=%.2f km², passes/day=%d, SGS contact=%.2f min\n', ...
        pm.orbitIndex,pm.accumAreaKm2,pm.passesPerDay,pm.sgsTotalMin);
end
for gi = 1:numel(bestMetrics.groundStationMetrics)
    gm = bestMetrics.groundStationMetrics(gi);
    fprintf('Ground station %-24s: contacts=%d, contact=%.2f min\n', ...
        gm.name,gm.numContacts,gm.totalContactMin);
end
fprintf('Constellation unique coverage: %.2f km²\n',areaKm2Best);
fprintf('Geometry-resolution margin:     ±%.2f km²\n',accuracyMarginKm2);
fprintf('Numerical range:                %.2f to %.2f km²\n', ...
    accuracyLowerKm2,accuracyUpperKm2);
fprintf('===========================================\n\n');

% Persist the final high-resolution accuracy and constellation metrics in the
% same MAT-file as the optimization history.
optResults.finalCoverageAreaKm2 = areaKm2Best;
optResults.referenceCoverageAreaKm2 = referenceAreaKm2;
optResults.coverageAccuracyMarginKm2 = accuracyMarginKm2;
optResults.coverageAccuracyRangeKm2 = [accuracyLowerKm2,accuracyUpperKm2];
optResults.coverageAccuracyPercent = accuracyPct;
optResults.finalConstellationMetrics = bestMetrics;
optResults.finalConstellationDetails = constellationDetails;
optResults.groundStations = groundStations;
optResults.constellation = constellation;
try
    save(saveFile,'optResults');
catch ME
    warning('Final result update could not be appended to %s: %s',saveFile,ME.message);
end

%% =========================================================
% RESULTS / MODEL SUMMARY
%% =========================================================
% This script deliberately keeps the present single-orbit numerical case:
%   - 1 orbital plane
%   - 1 satellite per plane
%   - 550 km circular orbit
%   - 10 deg sensor half-angle
%   - 100 LHS initial evaluations + up to 100 BO evaluations
%   - 60 s optimization sampling
%   - 2 s high-resolution final simulation
%
% Latest validated 1-orbit/1-satellite/1-GS run used as a regression reference:
%   Best orbit: i = 149.9500 deg, RAAN = 275.8911 deg, nu = 334.3402 deg
%   Best score: 0.0910
%   Nepal continuous swath union: 133186.10 km^2
%   Per-pass areas: P1 = 60939.05, P2 = 29059.51, P3 = 112869.11 km^2
%   Accumulated S = 202867.67 km^2; union U = 133186.10 km^2
%   Redundancy S-U = 69681.57 km^2; triple overlap = 0.00 km^2
%   Pairwise: P1∩P3 = 60270.33 km^2; P2∩P3 = 9411.24 km^2
%   Six lakes covered: 6/6; minimum lake passes/day = 2; worst = 180.27 min
%   These are reference values for regression checks. The 16-gon nominal area
%   is now accompanied by a 64-gon geometry-resolution sensitivity range at runtime.
%
% Reported outputs preserved by this optimized version:
%   - Best overall inclination / RAAN / true anomaly and score.
%   - Best-per-metric parameter locations.
%   - Nepal continuous swath coverage and coverage fraction.
%   - Number of Nepal swath passes/day and revisit statistics.
%   - SGS total/average contact time and number of contacts.
%   - Observation-to-next-SGS latency.
%   - Per-pass Nepal coverage areas.
%   - Accumulated area S, unique union U, and S-U redundancy.
%   - Pairwise overlap areas, triple overlap, and S-U consistency check.
%   - Six-lake pass count, dwell time, revisit, and coverage flags.
%   - Static world 2-D map, 3-D satellite-scenario view, and interactive
%     animated per-pass swath build-up.
%
% Optimization changes made for speed without changing the present 1 x 1
% mission result definition:
%   1) Cache all static geometry/trigonometry once.
%   2) Cache the buffered Nepal polygon instead of rebuilding it per trial.
%   3) Remove duplicate pass-area/swath construction in Section H.
%   4) Vectorize the six-lake footprint tests across all samples.
%   5) Disable/re-enable polyshape repair warnings once per polygon batch,
%      rather than once per footprint.
%   6) Preallocate Bayesian-optimization history arrays.
%   7) Avoid the redundant second GP predict() inside fmincon's EI objective.
%   8) Remove redundant local functions and the unused per-pass-area helper.
%   9) Keep final plotting and animation resolution unchanged.
%
% Scalability status (current run remains 1 x 1 x 1):
%   ✓ 1 orbit / 1 satellite                  -> implemented
%   ✓ Change orbital parameters              -> implemented
%   ✓ Add GLOF targets                       -> implemented as arrays
%   ✓ Add multiple ground stations           -> implemented as an array
%   ✓ Multiple satellites per orbit          -> computation path implemented
%   ✓ Multiple orbital planes                -> computation path implemented
%   ✓ Constellation-wide union coverage      -> implemented
%   ✓ Constellation-wide revisit             -> implemented by interval union
%   ✓ Satellite-to-GS access across fleet   -> implemented for every sat/GS pair
%   ✓ Per-satellite metrics                  -> implemented
%   ✓ Per-plane metrics                     -> implemented
%   ✓ Constellation-wide GLOF metrics        -> implemented
%   The optimizer still exposes one base [i,RAAN,nu] orbit and expands it
%   deterministically to the configured constellation. Current deployment is
%   deliberately one plane, one satellite, and one ground station.
%
% Mathematical consistency targets:
%   S = Σ_i A_i
%   U = area(⋃_i P_i)
%   redundancy = S-U
%   A>=2 = area(⋃_{i<j}(P_i∩P_j))
%   consistency: S-U = A>=2 + A_triple
%
% Numerical area sensitivity added in this version:
%   A_nominal = union area using 16-gon instantaneous footprints.
%   A_ref     = union area recomputed using the configured reference polygon.
%   margin    = |A_ref - A_nominal|.
%   reported range = [A_nominal-margin, A_nominal+margin].
% This is a numerical geometry-resolution sensitivity band, not an absolute
% sensor/geodetic uncertainty model.
%
% Constellation-ready architecture now supports, at the computation layer:
%   - arbitrary orbital planes and satellites/plane;
%   - arbitrary ground-station arrays;
%   - constellation-wide union coverage and revisit intervals;
%   - constellation-wide satellite-to-ground-station access;
%   - per-satellite metrics;
%   - per-plane metrics;
%   - constellation-wide lake coverage metrics.
% The active mission remains exactly 1 plane x 1 satellite x 1 ground station.
%% =========================================================

%% =========================================================
% I. INCREMENTAL CONSTELLATION GROWTH (INNER/OUTER BAYESIAN OPTIMIZATION)
%% =========================================================
% Grows the validated 1x1x1 baseline above (BESTPARAMS) into a larger
% constellation using the two-layer incremental design strategy:
%
%   Inner BO - try adding one satellite to the currently active orbital
%              plane, optimizing its true-anomaly phasing and a
%              collision-safe altitude offset (kept at least
%              COLLISIONOFFSETKM from every satellite already sharing
%              that plane).
%   Outer BO - if same-orbit growth no longer clears the improvement
%              thresholds, optimize a brand-new orbital plane
%              (inclination, RAAN) and place its first satellite with a
%              lighter-weight nested search (no prior satellite in that
%              plane to avoid, so a small LHS sweep over phasing/altitude
%              is sufficient).
%   Decision - accept a same-orbit candidate whenever it clears a
%              user-defined gain threshold on coverage, revisit, latency,
%              or the largest daily pass-time gap (see below); otherwise
%              fall back to the outer loop. A new orbital plane, once
%              found, is adopted (it is the only remaining growth action
%              at that point) and its resulting gains are still reported.
%              Both branches are subject to a hard veto first: any
%              candidate that would let the worst-case (per-lake) revisit
%              regress by more than WORSTREVISITREGRESSIONTOLERANCEMIN is
%              excluded from selection outright, regardless of how well it
%              scores on coverage/revisit/latency/Wasserstein/gap. If
%              every candidate in a step's pool trips that veto, the
%              growth loop stops rather than accept a regression.
%   Objective - both loops maximize the existing mission score
%              (COMPUTESCORE) plus (a) a Wasserstein (earth-mover) term
%              that rewards new pass-times that are temporally
%              complementary to the existing constellation's Nepal
%              pass-time distribution, and (b) a direct largest-gap term
%              that rewards shrinking the single largest empty stretch (in
%              hours, wraparound-aware) in the combined constellation's
%              daily Nepal pass schedule. Wasserstein rewards overall
%              distributional spread but does not specifically target that
%              one largest empty gap, so the gap term is the more surgical
%              fix for gaps such as the pre-dawn/midday/evening blank
%              stretches otherwise visible in the pass-time histogram.
%
% A step-wise table (candidate parameters, new-satellite-only metrics,
% combined constellation metrics, Wasserstein score, largest-gap-before/
% after, decision) is printed and saved after every accepted step,
% mirroring the report format used throughout this script.
%
% Sections A-H above (the validated 1x1x1 baseline) are completely
% untouched by this section. EVALUATEELEMENTSLIST is a self-contained
% generalization of RUNSIM's computation path that also accepts a
% per-satellite altitude, so growth steps can give every added
% satellite/orbit its own independently optimized altitude; RUNSIM itself
% is not modified, so it remains the regression reference.

growthCfg = struct( ...
    'enabled',                 true, ...   % set false to skip this section
    'maxSteps',                6, ...      % hard cap on growth iterations
    'maxOrbits',               3, ...      % hard cap on number of planes
    'maxSatsPerOrbit',         3, ...      % hard cap on satellites/plane
    'collisionOffsetKm',       15, ...     % min |altitude offset| within a shared plane
    'altitudeRangeKm',         40, ...     % +/- search window for altitude offset
    'coverageGainThreshold',   0.03, ...   % min absolute coverage-fraction gain to accept
    'revisitGainThresholdMin', 5, ...      % min avg-revisit reduction (minutes) to accept
    'latencyGainThresholdMin', 3, ...      % min avg-latency reduction (minutes) to accept
    'worstRevisitRegressionToleranceMin', 15, ... % max allowed regression in the worst-case (per-lake) revisit before a candidate is vetoed outright
    'gapGainThresholdHr',      1.0, ...    % min largest-schedule-gap reduction (hours) that alone justifies accepting a candidate
    'nInitialInner',           8, ...      % inner-BO initial LHS samples
    'nBOInner',                10, ...     % inner-BO Expected-Improvement iterations
    'nInitialOuter',           6, ...      % outer-BO initial LHS samples
    'nBOOuter',                8, ...      % outer-BO Expected-Improvement iterations
    'nNestedInnerLHS',         6, ...      % nested LHS trials for a new orbit's 1st sat
    'wassersteinWeight',       0.10, ...   % wW in the combined growth objective
    'gapWeight',               0.10);      % wG: direct reward for shrinking the largest pass-time gap

if growthCfg.enabled

fprintf('\n================ INCREMENTAL CONSTELLATION GROWTH ================\n');
fprintf('Seeding growth from the validated baseline orbit (incl=%.2f, RAAN=%.2f, nu=%.2f).\n', ...
    bestParams(1),bestParams(2),bestParams(3));

periodSec = seconds(stopTime-startTime);

% Seed the growth loop with the validated single-satellite baseline as
% orbit 1 / satellite 1. All later satellites/orbits carry their own
% independently optimized SEMIMAJORAXIS and ALTITUDEOFFSETKM.
currentElements = struct( ...
    'inclination',bestParams(1), ...
    'raan',bestParams(2), ...
    'trueAnomaly',bestParams(3), ...
    'semiMajorAxis',semiMajorAxis, ...
    'orbitIndex',1, ...
    'satIndex',1, ...
    'altitudeOffsetKm',0);

currentMetrics = evaluateElementsList(currentElements,fixed);

% DECIDEGROWTHSTEP needs the schedule period to evaluate the wraparound-
% aware largest-gap criterion; carrying it on GROWTHCFG avoids changing
% that function's call signature at every call site.
growthCfg.periodSec = periodSec;

baselineGapHr = largestScheduleGapHours(currentMetrics.observationIntervals,periodSec);
fprintf('Baseline largest daily pass-time gap: %.2f hr.\n',baselineGapHr);

% Extend the mission score with a Wasserstein complementarity term and a
% direct largest-gap-reduction term for the growth objective only; the
% baseline SCOREWEIGHTS/SCORESCALES (Section C) are left untouched.
growthWeights = scoreWeights;
growthWeights.wasserstein = growthCfg.wassersteinWeight;
growthWeights.gap = growthCfg.gapWeight;
growthScales = scoreScales;
growthScales.wasserstein = 1.0;
growthScales.gap = periodSec/3600;   % normalize by the schedule length (24 h)

optsFminGrowth = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',30);

growthLog = struct('step',{},'action',{},'inclination',{},'raan',{}, ...
    'trueAnomalyDeg',{},'altitudeOffsetKm',{},'newSatMetrics',{}, ...
    'combinedMetrics',{},'wasserstein',{},'oldGapHr',{},'newGapHr',{}, ...
    'gapReductionHr',{},'decision',{});

activeOrbit = 1;
orbitCount = 1;
growthStep = 0;
stopReason = 'maximum number of growth steps reached';

while growthStep < growthCfg.maxSteps

    growthStep = growthStep+1;
    membersActive = find([currentElements.orbitIndex] == activeOrbit);
    acceptInner = false;

    %% ---- INNER BO: try growing the currently active orbital plane. ----
    if numel(membersActive) < growthCfg.maxSatsPerOrbit

        fprintf('\n[Step %d] Inner BO: trying to add a satellite to orbit %d ...\n', ...
            growthStep,activeOrbit);

        lbInner = [0,-growthCfg.altitudeRangeKm];
        ubInner = [360, growthCfg.altitudeRangeKm];

        rng(2000+growthStep,'twister');
        lhsU = lhsdesign(growthCfg.nInitialInner,2);
        Xin = lbInner+lhsU.*(ubInner-lbInner);
        Jin = nan(growthCfg.nInitialInner,1);
        resIn = cell(growthCfg.nInitialInner,1);
        for k = 1:growthCfg.nInitialInner
            [Jin(k),resIn{k}] = evaluateSameOrbitCandidate(Xin(k,:),currentElements, ...
                activeOrbit,fixed,growthWeights,growthScales,periodSec, ...
                growthCfg.collisionOffsetKm,growthCfg.altitudeRangeKm,semiMajorAxis);
        end

        for b = 1:growthCfg.nBOInner
            gprIn = fitrgp(Xin,Jin,'KernelFunction','ardsquaredexponential', ...
                'Standardize',true,'BasisFunction','constant');
            bestSoFarIn = max(Jin);
            candU = lhsdesign(300,2);
            Xc = lbInner+candU.*(ubInner-lbInner);
            [muC,sdC] = predict(gprIn,Xc);
            eiC = expectedImprovement(muC,sdC,bestSoFarIn);
            [~,ixBest] = max(eiC);
            x0 = Xc(ixBest,:);
            try
                xNextIn = fmincon(@(xx) negativeExpectedImprovement(gprIn,xx,bestSoFarIn), ...
                    x0,[],[],[],[],lbInner,ubInner,[],optsFminGrowth);
            catch
                xNextIn = x0;
            end
            xNextIn = min(max(xNextIn,lbInner),ubInner);
            [Jnext,resNext] = evaluateSameOrbitCandidate(xNextIn,currentElements, ...
                activeOrbit,fixed,growthWeights,growthScales,periodSec, ...
                growthCfg.collisionOffsetKm,growthCfg.altitudeRangeKm,semiMajorAxis);
            Xin(end+1,:) = xNextIn; %#ok<AGROW>
            Jin(end+1) = Jnext; %#ok<AGROW>
            resIn{end+1} = resNext; %#ok<AGROW>
        end

        % Hard veto first: exclude any candidate that would let the
        % worst-case (per-lake) revisit regress beyond tolerance, then
        % pick the best-scoring survivor. If nothing in the pool survives
        % the veto, same-orbit growth is treated as unavailable this step
        % (falls through to outer BO below) rather than accepting a
        % regression.
        guardMaskIn = false(numel(resIn),1);
        for gk = 1:numel(resIn)
            guardMaskIn(gk) = passesWorstRevisitGuard(currentMetrics, ...
                resIn{gk}.combinedMetrics,growthCfg);
        end

        if any(guardMaskIn)
            JinGuarded = Jin;
            JinGuarded(~guardMaskIn) = -Inf;
            [~,bestInIdx] = max(JinGuarded);
            innerResult = resIn{bestInIdx};

            [acceptInner,reasonInner,gainsInner] = decideGrowthStep(currentMetrics, ...
                innerResult.combinedMetrics,growthCfg);

            growthLog(end+1) = struct( ...
                'step',growthStep,'action',innerResult.action, ...
                'inclination',innerResult.inclination,'raan',innerResult.raan, ...
                'trueAnomalyDeg',innerResult.trueAnomalyDeg, ...
                'altitudeOffsetKm',innerResult.altitudeOffsetKm, ...
                'newSatMetrics',innerResult.newSatMetrics, ...
                'combinedMetrics',innerResult.combinedMetrics, ...
                'wasserstein',innerResult.wasserstein, ...
                'oldGapHr',gainsInner.oldGapHr,'newGapHr',gainsInner.newGapHr, ...
                'gapReductionHr',gainsInner.gapReductionHr, ...
                'decision',reasonInner); %#ok<AGROW>
            printGrowthStepRow(growthStep,innerResult,reasonInner);

            if acceptInner
                currentElements = innerResult.combinedElements;
                currentMetrics = innerResult.combinedMetrics;
                continue;
            end
        else
            fprintf(['\n[Step %d] Every same-orbit candidate for orbit %d would regress the ' ...
                'worst-case lake revisit beyond the %.0f-min tolerance; skipping inner BO.\n'], ...
                growthStep,activeOrbit,growthCfg.worstRevisitRegressionToleranceMin);
        end
    else
        fprintf('\n[Step %d] Orbit %d already holds the maximum %d satellites; skipping inner BO.\n', ...
            growthStep,activeOrbit,growthCfg.maxSatsPerOrbit);
    end

    %% ---- Same-orbit growth was skipped or rejected: try a new orbit. ----
    if orbitCount >= growthCfg.maxOrbits
        stopReason = 'same-orbit growth is exhausted and the maximum orbit count was reached';
        fprintf('\n[Step %d] Maximum orbit count reached and same-orbit growth is exhausted. Stopping.\n', ...
            growthStep);
        break;
    end

    fprintf('\n[Step %d] Outer BO: searching for a new orbital plane ...\n',growthStep);

    newOrbitIndex = orbitCount+1;
    lbOuter = [0,0];
    ubOuter = [90,360];   % inclination capped at 90 deg for prograde-only orbits

    rng(3000+growthStep,'twister');
    lhsUOut = lhsdesign(growthCfg.nInitialOuter,2);
    Xout = lbOuter+lhsUOut.*(ubOuter-lbOuter);
    Jout = nan(growthCfg.nInitialOuter,1);
    resOut = cell(growthCfg.nInitialOuter,1);
    for k = 1:growthCfg.nInitialOuter
        [Jout(k),resOut{k}] = evaluateNewOrbitCandidate(Xout(k,:),currentElements, ...
            fixed,growthWeights,growthScales,periodSec,growthCfg.altitudeRangeKm, ...
            semiMajorAxis,growthCfg.nNestedInnerLHS,newOrbitIndex);
    end

    for b = 1:growthCfg.nBOOuter
        gprOut = fitrgp(Xout,Jout,'KernelFunction','ardsquaredexponential', ...
            'Standardize',true,'BasisFunction','constant');
        bestSoFarOut = max(Jout);
        candU = lhsdesign(300,2);
        Xc = lbOuter+candU.*(ubOuter-lbOuter);
        [muC,sdC] = predict(gprOut,Xc);
        eiC = expectedImprovement(muC,sdC,bestSoFarOut);
        [~,ixBest] = max(eiC);
        x0 = Xc(ixBest,:);
        try
            xNextOut = fmincon(@(xx) negativeExpectedImprovement(gprOut,xx,bestSoFarOut), ...
                x0,[],[],[],[],lbOuter,ubOuter,[],optsFminGrowth);
        catch
            xNextOut = x0;
        end
        xNextOut = min(max(xNextOut,lbOuter),ubOuter);
        [Jnext,resNext] = evaluateNewOrbitCandidate(xNextOut,currentElements, ...
            fixed,growthWeights,growthScales,periodSec,growthCfg.altitudeRangeKm, ...
            semiMajorAxis,growthCfg.nNestedInnerLHS,newOrbitIndex);
        Xout(end+1,:) = xNextOut; %#ok<AGROW>
        Jout(end+1) = Jnext; %#ok<AGROW>
        resOut{end+1} = resNext; %#ok<AGROW>
    end

    % Same hard veto as the inner loop, applied here too: a new orbital
    % plane would otherwise always be adopted (it is normally the only
    % remaining growth action at this point), so this is the last line of
    % defense against a worst-case-revisit regression. If nothing in the
    % outer-BO pool survives the veto, there is no valid growth action
    % left at all and the loop stops.
    guardMaskOut = false(numel(resOut),1);
    for gk = 1:numel(resOut)
        guardMaskOut(gk) = passesWorstRevisitGuard(currentMetrics, ...
            resOut{gk}.combinedMetrics,growthCfg);
    end

    if ~any(guardMaskOut)
        stopReason = sprintf(['no candidate satellite or orbital plane could be found that ' ...
            'keeps the worst-case lake revisit within %.0f min of its current value'], ...
            growthCfg.worstRevisitRegressionToleranceMin);
        fprintf(['\n[Step %d] Every new-orbit candidate would regress the worst-case lake ' ...
            'revisit beyond the %.0f-min tolerance. Stopping.\n'], ...
            growthStep,growthCfg.worstRevisitRegressionToleranceMin);
        break;
    end

    JoutGuarded = Jout;
    JoutGuarded(~guardMaskOut) = -Inf;
    [~,bestOutIdx] = max(JoutGuarded);
    outerResult = resOut{bestOutIdx};

    % A new orbital plane is adopted once found (subject to the
    % worst-case-revisit veto above); DECIDEGROWTHSTEP is still used only
    % to report the achieved gains for the step table.
    [~,~,gainsOuter] = decideGrowthStep(currentMetrics,outerResult.combinedMetrics,growthCfg);
    reasonOuter = sprintf(['Accept (new orbital plane adopted; coverage %+.1f%%, ' ...
        'revisit %+.1f min, largest gap %+.2f hr, Wasserstein=%s)'], ...
        100*gainsOuter.coverageGain,gainsOuter.revisitGainMin,gainsOuter.gapReductionHr, ...
        formatMaybeNaN(outerResult.wasserstein));

    growthLog(end+1) = struct( ...
        'step',growthStep,'action',outerResult.action, ...
        'inclination',outerResult.inclination,'raan',outerResult.raan, ...
        'trueAnomalyDeg',outerResult.trueAnomalyDeg, ...
        'altitudeOffsetKm',outerResult.altitudeOffsetKm, ...
        'newSatMetrics',outerResult.newSatMetrics, ...
        'combinedMetrics',outerResult.combinedMetrics, ...
        'wasserstein',outerResult.wasserstein, ...
        'oldGapHr',gainsOuter.oldGapHr,'newGapHr',gainsOuter.newGapHr, ...
        'gapReductionHr',gainsOuter.gapReductionHr, ...
        'decision',reasonOuter); %#ok<AGROW>
    printGrowthStepRow(growthStep,outerResult,reasonOuter);

    currentElements = outerResult.combinedElements;
    currentMetrics = outerResult.combinedMetrics;
    activeOrbit = newOrbitIndex;
    orbitCount = orbitCount+1;

end % while growthStep < growthCfg.maxSteps

fprintf('\nGrowth loop finished: %s.\n',stopReason);
fprintf('Final constellation: %d orbit(s), %d satellite(s) total.\n', ...
    orbitCount,numel(currentElements));

%% ---------------------------------------------------------
% STEP-WISE OUTPUT TABLE
%% ---------------------------------------------------------
if ~isempty(growthLog)
    growthTable = growthLogToTable(growthLog);
    fprintf('\n================ CONSTELLATION GROWTH STEP TABLE ================\n');
    disp(growthTable);
end

%% ---------------------------------------------------------
% FINAL GROWN-CONSTELLATION SUMMARY
%% ---------------------------------------------------------
fprintf('\n========= FINAL GROWN CONSTELLATION METRICS =========\n');
printMetricsBlock(currentMetrics,[]);
fprintf('Six lakes covered: %d/6; worst lake revisit: %s min\n', ...
    currentMetrics.numLakesCovered,formatMaybeNaN(currentMetrics.worstRevisitMin));
finalGapHr = largestScheduleGapHours(currentMetrics.observationIntervals,periodSec);
fprintf('Largest daily pass-time gap: %.2f hr (baseline was %.2f hr)\n', ...
    finalGapHr,baselineGapHr);
fprintf('======================================================\n\n');

%% ---------------------------------------------------------
% SAVE GROWTH RESULTS
%% ---------------------------------------------------------
growthResults = struct();
growthResults.growthCfg = growthCfg;
growthResults.growthLog = growthLog;
growthResults.finalElements = currentElements;
growthResults.finalMetrics = currentMetrics;
growthResults.baselineParams = bestParams;
growthResults.baselineMetrics = bestMetrics;
growthResults.stopReason = stopReason;

growthSaveFile = sprintf('cubesat_constellation_growth_%s.mat', ...
    datestr(now,'yyyymmdd_HHMMSS'));
save(growthSaveFile,'growthResults');
fprintf('Saved constellation growth results to: %s\n',growthSaveFile);

if ~isempty(growthLog)
    try
        growthCsvFile = sprintf('cubesat_constellation_growth_%s.csv', ...
            datestr(now,'yyyymmdd_HHMMSS'));
        writetable(growthTable,growthCsvFile);
        fprintf('Saved step-wise growth table to: %s\n\n',growthCsvFile);
    catch ME
        warning('Could not write growth step CSV: %s',ME.message);
    end
end

%% ---------------------------------------------------------
% TEMPORAL PASS-TIME HISTOGRAM: BASELINE VS GROWN CONSTELLATION
%% ---------------------------------------------------------
% Illustrates how the added satellites/orbits fill gaps in the baseline's
% Nepal pass-time-of-day distribution — the effect the Wasserstein and
% largest-gap terms are designed to reward. Styled to match every other
% figure in this script (dark background, white text/grid) and colored to
% match the final constellation's static/animated maps below: NSATFINAL and
% SATCOLORSFINALGROWN are defined here (rather than lower down, where they
% were previously first introduced) so this figure, the static map, and the
% animated map all draw satellite P%02d/S%02d from the same LINES(...)
% colormap — a bar's color in this histogram now identifies the same
% satellite as that color's ground track/footprint in the maps that follow.
nSatFinal = numel(currentElements);
satColorsFinalGrown = lines(nSatFinal);

baselineObsHrs = mod(mean(bestMetrics.observationIntervals,2),periodSec)/3600;
histEdges = 0:1:24;   % 1-hour bins (finer than the previous 2-hour bins)

figHist = figure('Name','Nepal Pass-Time Distribution - Before vs After Growth', ...
    'Position',[150,150,1050,480]);
set(figHist,'Color',[0.1 0.1 0.1]);

% --- Baseline panel (single satellite): colored to match Section H's
% green baseline ground-track/marker color for the same cross-figure
% consistency. ---
axHistBase = subplot(1,2,1);
histogram(axHistBase,baselineObsHrs,histEdges,'FaceColor',[0.3 0.85 0.3],'EdgeColor','w');
formatHistAxes(axHistBase,sprintf('Baseline (1 orbit x 1 sat), n=%d', ...
    size(bestMetrics.observationIntervals,1)));
[baselineGapHrPanel,baseGapStartHr,baseGapEndHr] = largestScheduleGapHours( ...
    bestMetrics.observationIntervals,periodSec);
shadeLargestGap(axHistBase,baseGapStartHr,baseGapEndHr,baselineGapHrPanel);

% --- Grown-constellation panel: stacked per-satellite breakdown so every
% satellite's contribution to the daily schedule is visible on its own,
% using the SAME per-satellite colors as the final constellation's
% static/animated ground-track maps. ---
[~,finalGapStartHr,finalGapEndHr] = largestScheduleGapHours( ...
    currentMetrics.observationIntervals,periodSec);

axHistFinal = subplot(1,2,2);
binCenters = histEdges(1:end-1)+diff(histEdges)/2;
stackedCounts = zeros(numel(histEdges)-1,nSatFinal);
for si = 1:nSatFinal
    obsHrsSi = mod(mean(currentMetrics.satelliteMetrics(si).observationIntervals,2), ...
        periodSec)/3600;
    stackedCounts(:,si) = histcounts(obsHrsSi,histEdges).';
end
hBarFinal = bar(axHistFinal,binCenters,stackedCounts,'stacked','BarWidth',1);
for si = 1:nSatFinal
    hBarFinal(si).FaceColor = satColorsFinalGrown(si,:);
    hBarFinal(si).EdgeColor = 'w';
    hBarFinal(si).LineWidth = 0.5;
    hBarFinal(si).DisplayName = sprintf('P%02d S%02d', ...
        currentElements(si).orbitIndex,currentElements(si).satIndex);
end
legend(axHistFinal,'Location','eastoutside','TextColor','w', ...
    'Color',[0.1 0.1 0.1],'FontSize',7);
formatHistAxes(axHistFinal,sprintf('Grown constellation (%d orbit(s), %d sat(s)), n=%d', ...
    orbitCount,nSatFinal,size(currentMetrics.observationIntervals,1)));
shadeLargestGap(axHistFinal,finalGapStartHr,finalGapEndHr,finalGapHr);

sgtitle(figHist,'Effect of incremental growth on Nepal pass-time spread','Color','w');

%% ---------------------------------------------------------
% FINAL GROWN CONSTELLATION: 3-D VIEWER, STATIC MAP, ANIMATED MAP
%% ---------------------------------------------------------
% Rebuilds the grown constellation at the same high resolution Section H
% uses for the single-satellite baseline, and reproduces the same three
% views: a 3-D satellite-scenario viewer, a static full-world ground-track
% map with the swath union, and an animated ground-track build-up —
% generalized to every satellite in the final constellation instead of
% just one. Section H itself is untouched; this is a separate, additive
% rebuild from CURRENTELEMENTS only. NSATFINAL and SATCOLORSFINALGROWN are
% defined earlier, above the pass-time histogram, so that figure shares
% this section's per-satellite color coding.

scFinalGrown = satelliteScenario(startTime,stopTime,baseSampleTime);

satFinalGrown = cell(1,nSatFinal);
for si = 1:nSatFinal
    e = currentElements(si);
    satName = sprintf('GrownSat_P%02d_S%02d',e.orbitIndex,e.satIndex);
    satFinalGrown{si} = satellite(scFinalGrown,e.semiMajorAxis,eccentricity, ...
        e.inclination,e.raan,0,e.trueAnomaly,Name=satName);
    satFinalGrown{si}.MarkerColor = [0 1 0];
    satFinalGrown{si}.MarkerSize = 8;
    satFinalGrown{si}.Orbit.LineWidth = 2;
end

gsFinalGrown = cell(1,numel(groundStations));
for gi = 1:numel(groundStations)
    gsFinalGrown{gi} = groundStation(scFinalGrown,groundStations(gi).lat, ...
        groundStations(gi).lon,Name=groundStations(gi).name);
end

% GLOF-lake markers for the 3-D viewer, as in Section H (display only).
glacierMarkersFinalGrown = cell(1,numel(glacierNames));
for gi = 1:numel(glacierNames)
    glacierMarkersFinalGrown{gi} = groundStation(scFinalGrown, ...
        glacierLat(gi),glacierLon(gi),Name=glacierNames{gi});
    glacierMarkersFinalGrown{gi}.MarkerColor = [0 1 1];
    glacierMarkersFinalGrown{gi}.MarkerSize = 6;
end

accFinalGrown = cell(nSatFinal,numel(gsFinalGrown));
for si = 1:nSatFinal
    for gi = 1:numel(gsFinalGrown)
        accFinalGrown{si,gi} = access(satFinalGrown{si},gsFinalGrown{gi});
        accFinalGrown{si,gi}.LineColor = [0 1 0];
        accFinalGrown{si,gi}.LineWidth = 2;
    end
end

sensorFinalGrown = cell(1,nSatFinal);
for si = 1:nSatFinal
    sensorFinalGrown{si} = conicalSensor(satFinalGrown{si},'MaxViewAngle',sensorHalfAngle, ...
        'Name',sprintf('Grown Constellation Swath P%02d S%02d', ...
        currentElements(si).orbitIndex,currentElements(si).satIndex));
end

viewerFinalGrown = satelliteScenarioViewer(scFinalGrown,'Dimension','3D', ...
    'Name','Final Grown Constellation - Viewer');
for si = 1:nSatFinal
    fieldOfView(sensorFinalGrown{si});
end
playScenarioSafely(scFinalGrown,sensorFinalGrown,'Final Grown Constellation viewer');

% High-resolution geographic tracks for every satellite in the final
% constellation (mirrors Section H's BESTTRACKS construction).
finalGrownTracks = cell(1,nSatFinal);
for si = 1:nSatFinal
    [posF,~,~] = states(satFinalGrown{si},'CoordinateFrame','geographic');
    latRawF = posF(1,:);
    lonRawF = posF(2,:);
    altRawF = posF(3,:)/1000;
    validF = isfinite(latRawF) & isfinite(lonRawF) & isfinite(altRawF);
    latRawF = latRawF(validF);
    lonRawF = lonRawF(validF);
    altRawF = altRawF(validF);
    [latSegF,lonSegF,altSegF] = segmentGroundTrack(latRawF,lonRawF,altRawF,10);
    finalGrownTracks{si} = struct('lat',latRawF,'lon',lonRawF,'alt',altRawF, ...
        'latSeg',latSegF,'lonSeg',lonSegF,'altSeg',altSegF);
end

% Authoritative swath geometry for the final constellation: union all
% satellite swaths (mirrors Section H's UNIONBEST construction).
satSwathPolysFinalGrown = cell(1,nSatFinal);
maxOrbitIdxFinal = max([currentElements.orbitIndex]);
planeSwathPolysFinalGrown = cell(1,maxOrbitIdxFinal);
for oi = 1:maxOrbitIdxFinal
    planeSwathPolysFinalGrown{oi} = polyshape();
end
unionFinalGrown = polyshape();

for si = 1:nSatFinal
    satPolyF = computeSwathPolygon(finalGrownTracks{si}.latSeg, ...
        finalGrownTracks{si}.lonSeg,finalGrownTracks{si}.altSeg,fixed);
    satSwathPolysFinalGrown{si} = satPolyF;
    if ~isempty(satPolyF.Vertices)
        unionFinalGrown = union(unionFinalGrown,satPolyF);
        oiPlane = currentElements(si).orbitIndex;
        planeSwathPolysFinalGrown{oiPlane} = union(planeSwathPolysFinalGrown{oiPlane},satPolyF);
    end
end

if isempty(unionFinalGrown.Vertices)
    areaKm2FinalGrown = 0;
else
    areaKm2FinalGrown = area(unionFinalGrown)*kmPerDegLat*kmPerDegLonMean;
end

%% ---------------------------------------------------------
% FINAL CONSTELLATION - STATIC WORLD 2-D MAP
%% ---------------------------------------------------------
figure('Name','Final Grown Constellation - Full World Ground Tracks with Swath', ...
    'Position',[100,100,1200,700]);
set(gcf,'Color',[0.1 0.1 0.1]);
axFinalGrown = axes;
hold(axFinalGrown,'on');
set(axFinalGrown,'Color',[0 0 0]);
grid(axFinalGrown,'on');
axFinalGrown.GridColor = [0.3 0.3 0.3];
axFinalGrown.XColor = 'w';
axFinalGrown.YColor = 'w';
box(axFinalGrown,'on');
load coastlines;
plot(axFinalGrown,coastlon,coastlat,'Color',[0.6 0.6 0.6], ...
    'LineWidth',0.5,'DisplayName','Coastlines');
plot(axFinalGrown,nepalBorderLon,nepalBorderLat,'Color',[0.8 0.6 0.4], ...
    'LineWidth',1.5,'DisplayName','Nepal Border');

for si = 1:nSatFinal
    dispName = sprintf('Ground Track P%02d S%02d', ...
        currentElements(si).orbitIndex,currentElements(si).satIndex);
    plot(axFinalGrown,finalGrownTracks{si}.lonSeg,finalGrownTracks{si}.latSeg, ...
        'Color',satColorsFinalGrown(si,:),'LineWidth',1.0,'DisplayName',dispName);
end

if ~isempty(unionFinalGrown.Vertices)
    plot(axFinalGrown,unionFinalGrown,'FaceColor',[1 1 1],'FaceAlpha',0.35, ...
        'EdgeColor','none','DisplayName','Swath Inside Nepal');
end

for gi = 1:numel(groundStations)
    plot(axFinalGrown,groundStations(gi).lon,groundStations(gi).lat,'p', ...
        'MarkerSize',12,'MarkerFaceColor','yellow','MarkerEdgeColor','red', ...
        'DisplayName',groundStations(gi).name);
end

plot(axFinalGrown,glacierLon,glacierLat,'d','MarkerSize',8, ...
    'MarkerFaceColor','cyan','MarkerEdgeColor','k','LineStyle','none', ...
    'DisplayName','GLOF Lakes');
for gi = 1:numel(glacierLat)
    text(axFinalGrown,glacierLon(gi)+0.15,glacierLat(gi),glacierNames{gi}, ...
        'Color','w','FontSize',7,'Interpreter','none','HandleVisibility','off');
end
xlabel(axFinalGrown,'Longitude (°)','Color','w');
ylabel(axFinalGrown,'Latitude (°)','Color','w');
title(axFinalGrown,sprintf('Final Grown Constellation: %d orbit(s), %d satellite(s)', ...
    orbitCount,nSatFinal),'Color','w');
xlim(axFinalGrown,[-180 180]);
ylim(axFinalGrown,[-90 90]);
legend(axFinalGrown,'Location','southoutside','Orientation','horizontal', ...
    'NumColumns',4,'TextColor','w','Color',[0.1 0.1 0.1]);

fprintf('\nFinal grown constellation static map: %d orbit(s), %d satellite(s), %.2f km² union swath.\n', ...
    orbitCount,nSatFinal,areaKm2FinalGrown);

%% ---------------------------------------------------------
% FINAL CONSTELLATION - ANIMATED GROUND-TRACK BUILD-UP
%% ---------------------------------------------------------
% Generalizes Section H2's single-satellite growing-footprint animation to
% every satellite in the grown constellation: each satellite's marker and
% trailing ground track move together (one persistent color per
% satellite, matching the static map above), and each satellite's own
% accumulated Nepal footprint keeps growing, uninterrupted, for as long as
% it spends over the region during the day — visualizing how the
% constellation's daily coverage builds up. All satellites share one
% wall-clock frame counter because they were propagated over the identical
% START/STOP/SAMPLE time grid; each satellite's own NaN-segmented track
% array is indexed through a precomputed raw-sample-to-segmented-index map
% so the shared clock still draws every satellite's line correctly.

rawToSegIdxFinal = cell(1,nSatFinal);
inPassFlagRawFinal = cell(1,nSatFinal);
nRawLenFinal = zeros(1,nSatFinal);
nRawSamplesFinal = numel(finalGrownTracks{1}.lat);

for si = 1:nSatFinal
    nRawLenFinal(si) = numel(finalGrownTracks{si}.lon);

    nonNanPos = find(~isnan(finalGrownTracks{si}.latSeg));
    if isempty(nonNanPos)
        nonNanPos = 1;
    end
    if numel(nonNanPos) < nRawSamplesFinal
        nonNanPos = [nonNanPos,repmat(nonNanPos(end),1,nRawSamplesFinal-numel(nonNanPos))]; %#ok<AGROW>
    end
    rawToSegIdxFinal{si} = nonNanPos(1:nRawSamplesFinal);

    rClampVec = min(1:nRawSamplesFinal,nRawLenFinal(si));
    lonRawSi = finalGrownTracks{si}.lon(rClampVec);
    latRawSi = finalGrownTracks{si}.lat(rClampVec);
    inPassFlagRawFinal{si} = ...
        lonRawSi >= fixed.lonMin & lonRawSi <= fixed.lonMax & ...
        latRawSi >= fixed.latMin & latRawSi <= fixed.latMax & ...
        ~isnan(lonRawSi) & ~isnan(latRawSi);
end

anyInPassRawFinal = false(1,nRawSamplesFinal);
for si = 1:nSatFinal
    anyInPassRawFinal = anyInPassRawFinal | inPassFlagRawFinal{si};
end

strideOutsideFinal = 80;
strideInsideFinal = 5;
framePauseFinal = 0.01;

frameIdxListFinal = [];
rFrame = 1;
while rFrame <= nRawSamplesFinal
    frameIdxListFinal(end+1) = rFrame; %#ok<AGROW>
    if anyInPassRawFinal(rFrame)
        rFrame = rFrame+strideInsideFinal;
    else
        rFrame = rFrame+strideOutsideFinal;
    end
end
frameIdxListFinal = unique([frameIdxListFinal,nRawSamplesFinal]);
frameIdxListFinal = frameIdxListFinal(frameIdxListFinal >= 1 & frameIdxListFinal <= nRawSamplesFinal);
nFramesFinal = numel(frameIdxListFinal);

figAnimFinal = figure('Name','Final Grown Constellation - Animated Ground-Track Build-up', ...
    'Position',[200,100,1200,750],'Interruptible','on','BusyAction','queue');
set(figAnimFinal,'Color',[0.1 0.1 0.1]);
axAnimFinal = axes(figAnimFinal);
hold(axAnimFinal,'on');
set(axAnimFinal,'Color',[0 0 0]);
grid(axAnimFinal,'on');
axAnimFinal.GridColor = [0.3 0.3 0.3];
axAnimFinal.XColor = 'w';
axAnimFinal.YColor = 'w';
box(axAnimFinal,'on');
xlim(axAnimFinal,[-180 180]);
ylim(axAnimFinal,[-90 90]);
axAnimFinal.XLimMode = 'manual';
axAnimFinal.YLimMode = 'manual';
xlabel(axAnimFinal,'Longitude (°)','Color','w');
ylabel(axAnimFinal,'Latitude (°)','Color','w');

load coastlines;
plot(axAnimFinal,coastlon,coastlat,'Color',[0.6 0.6 0.6], ...
    'LineWidth',0.5,'DisplayName','Coastlines');
plot(axAnimFinal,nepalBorderLon,nepalBorderLat,'Color',[0.8 0.6 0.4], ...
    'LineWidth',1.5,'DisplayName','Nepal Border');

for gi = 1:numel(groundStations)
    plot(axAnimFinal,groundStations(gi).lon,groundStations(gi).lat,'p', ...
        'MarkerSize',12,'MarkerFaceColor','yellow','MarkerEdgeColor','red', ...
        'DisplayName',groundStations(gi).name);
end

plot(axAnimFinal,glacierLon,glacierLat,'d','MarkerSize',8, ...
    'MarkerFaceColor','cyan','MarkerEdgeColor','k','LineStyle','none', ...
    'DisplayName','GLOF Lakes');
for gi = 1:numel(glacierLat)
    text(axAnimFinal,glacierLon(gi)+0.15,glacierLat(gi),glacierNames{gi}, ...
        'Color','w','FontSize',7,'Interpreter','none','HandleVisibility','off');
end

titleAnimFinal = title(axAnimFinal,'Final Grown Constellation  |  t = 00:00:00','Color','w');

hTrackFinal = gobjects(1,nSatFinal);
hSatMarkerFinal = gobjects(1,nSatFinal);
for si = 1:nSatFinal
    dispName = sprintf('Track P%02d S%02d', ...
        currentElements(si).orbitIndex,currentElements(si).satIndex);
    hTrackFinal(si) = plot(axAnimFinal,NaN,NaN,'Color',satColorsFinalGrown(si,:), ...
        'LineWidth',1.0,'DisplayName',dispName);
    hSatMarkerFinal(si) = plot(axAnimFinal,NaN,NaN,'o','MarkerSize',6, ...
        'MarkerFaceColor',satColorsFinalGrown(si,:),'MarkerEdgeColor','w', ...
        'HandleVisibility','off');
end

hFootprintFinal = gobjects(1,nSatFinal);
accumFootprintFinal = cell(1,nSatFinal);
for si = 1:nSatFinal
    hFootprintFinal(si) = gobjects(1);
    accumFootprintFinal{si} = polyshape();
end

legend(axAnimFinal,'Location','southoutside','Orientation','horizontal', ...
    'NumColumns',4,'TextColor','w','Color',[0.1 0.1 0.1]);

warnStateFinalAnim = warning('off','MATLAB:polyshape:repairedBySimplify');

for f = 1:nFramesFinal

    rFrame = frameIdxListFinal(f);
    elapsedSecFinal = (rFrame-1)*baseSampleTime;
    hh = floor(elapsedSecFinal/3600);
    mm = floor(mod(elapsedSecFinal,3600)/60);
    ss = floor(mod(elapsedSecFinal,60));

    for si = 1:nSatFinal

        segIdxFinal = rawToSegIdxFinal{si}(rFrame);

        set(hTrackFinal(si), ...
            'XData',finalGrownTracks{si}.lonSeg(1:segIdxFinal), ...
            'YData',finalGrownTracks{si}.latSeg(1:segIdxFinal));

        rUseFinal = min(rFrame,nRawLenFinal(si));
        lonNow = finalGrownTracks{si}.lon(rUseFinal);
        latNow = finalGrownTracks{si}.lat(rUseFinal);
        if ~isnan(lonNow) && ~isnan(latNow)
            set(hSatMarkerFinal(si),'XData',lonNow,'YData',latNow);
        end

        if inPassFlagRawFinal{si}(rFrame)

            altNow = finalGrownTracks{si}.alt(rUseFinal);
            R_local_km = altNow*tanHalfAngle;
            rLat = R_local_km*fixed.degPerKmLat;
            rLon = R_local_km/(fixed.kmPerDegLat*cosd(latNow));
            circLon = lonNow+rLon*fixed.circleCos;
            circLat = latNow+rLat*fixed.circleSin;
            newCirc = polyshape(circLon,circLat,'Simplify',false);
            accumFootprintFinal{si} = union(accumFootprintFinal{si},newCirc);
            clippedFootprint = intersect(accumFootprintFinal{si},nepalPoly);

            if ~isempty(clippedFootprint.Vertices)
                if isgraphics(hFootprintFinal(si))
                    delete(hFootprintFinal(si));
                end
                hFootprintFinal(si) = plot(axAnimFinal,clippedFootprint, ...
                    'FaceColor',satColorsFinalGrown(si,:),'FaceAlpha',0.45, ...
                    'EdgeColor','none','DisplayName', ...
                    sprintf('Footprint P%02d S%02d',currentElements(si).orbitIndex, ...
                    currentElements(si).satIndex));
            end

        end

    end

    set(titleAnimFinal,'String',sprintf( ...
        'Final Grown Constellation  |  t = %02d:%02d:%02d',hh,mm,ss));

    drawnow;
    if framePauseFinal > 0
        pause(framePauseFinal);
    end

end

warning(warnStateFinalAnim);

fprintf('Final grown-constellation animation complete: %d orbit(s), %d satellite(s), %.2f km² union swath.\n\n', ...
    orbitCount,nSatFinal,areaKm2FinalGrown);

end % if growthCfg.enabled

%% =========================================================
% LOCAL HELPER FUNCTIONS
%% =========================================================

function elements = expandConstellationElements(params,cfg)
% EXPANDCONSTELLATIONELEMENTS
% Expand the three primary orbit variables into a constellation-ready list.
% Current configuration is 1 orbit x 1 satellite, so this is numerically
% identical to the original single-satellite model.
% For a uniform extension, RAAN and true anomaly offsets are spaced as:
%   ΔRAAN = 360/N_orbits
%   Δnu   = 360/N_sat_per_orbit

    nOrbit = cfg.numOrbits;
    nSat = cfg.satsPerOrbit;
    nTotal = nOrbit*nSat;
    elements = repmat(struct('inclination',NaN,'raan',NaN,'trueAnomaly',NaN, ...
        'orbitIndex',NaN,'satIndex',NaN),1,nTotal);

    if isempty(cfg.raanOffsetsDeg)
        raanOffsets = (0:nOrbit-1)*cfg.raanSpacingDeg;
    else
        raanOffsets = cfg.raanOffsetsDeg(:).';
        if numel(raanOffsets) ~= nOrbit
            error('constellation.raanOffsetsDeg must contain one value per orbit.');
        end
    end
    if isempty(cfg.trueAnomalyOffsetsDeg)
        nuOffsets = (0:nSat-1)*cfg.trueAnomalySpacingDeg;
    else
        nuOffsets = cfg.trueAnomalyOffsetsDeg(:).';
        if numel(nuOffsets) ~= nSat
            error('constellation.trueAnomalyOffsetsDeg must contain one value per satellite/plane.');
        end
    end
    if isempty(cfg.inclinationOffsetsDeg)
        incOffsets = zeros(1,nOrbit);
    else
        incOffsets = cfg.inclinationOffsetsDeg(:).';
        if numel(incOffsets) ~= nOrbit
            error('constellation.inclinationOffsetsDeg must contain one value per orbit.');
        end
    end

    k = 0;
    for o = 1:nOrbit
        for ss = 1:nSat
            k = k+1;
            elements(k).inclination = params(1)+incOffsets(o);
            elements(k).raan = mod(params(2)+raanOffsets(o),360);
            elements(k).trueAnomaly = mod(params(3)+nuOffsets(ss),360);
            elements(k).orbitIndex = o;
            elements(k).satIndex = ss;
        end
    end
end

function metrics = runSim(params,fixed)
% RUNSIM  Evaluate one candidate mission using the constellation data path.
% Steps:
%   1) Expand base orbit parameters into configured satellites.
%   2) Propagate all satellites on the optimization time grid.
%   3) Build each satellite swath polygon and the full constellation union.
%   4) Merge Nepal observation intervals across all satellites.
%   5) Merge access intervals across every satellite-ground-station pair.
%   6) Compute constellation-wide lake hits plus per-satellite/per-plane data.

    elements = expandConstellationElements(params,fixed.constellation);
    nSat = numel(elements);
    nGS = numel(fixed.groundStations);

    sc = satelliteScenario(fixed.startTime,fixed.stopTime,fixed.sampleTime);
    sat = cell(1,nSat);
    for si = 1:nSat
        e = elements(si);
        sat{si} = satellite(sc,fixed.semiMajorAxis,fixed.eccentricity, ...
            e.inclination,e.raan,0,e.trueAnomaly, ...
            Name=sprintf('EvalSat_P%02d_S%02d',e.orbitIndex,e.satIndex));
    end

    gs = cell(1,nGS);
    for gi = 1:nGS
        gs{gi} = groundStation(sc,fixed.groundStations(gi).lat, ...
            fixed.groundStations(gi).lon,Name=fixed.groundStations(gi).name);
    end

    tracks = cell(1,nSat);
    satSwathPolys = cell(1,nSat);
    planeSwathPolys = cell(1,fixed.constellation.numOrbits);
    for oi = 1:fixed.constellation.numOrbits
        planeSwathPolys{oi} = polyshape();
    end

    satObservationIntervals = cell(1,nSat);
    satelliteGlacierMetrics = cell(1,nSat);
    glacierHitsBySat = cell(1,nSat);

    unionSwath = polyshape();
    allObservationIntervals = zeros(0,2);

    for si = 1:nSat
        [pos,~,~] = states(sat{si},'CoordinateFrame','geographic');
        lat = pos(1,:); lon = pos(2,:); alt = pos(3,:)/1000;
        valid = isfinite(lat) & isfinite(lon) & isfinite(alt);
        lat = lat(valid); lon = lon(valid); alt = alt(valid);

        [latSeg,lonSeg,altSeg] = segmentGroundTrack(lat,lon,alt,10);
        tracks{si} = struct('lat',lat,'lon',lon,'alt',alt, ...
            'latSeg',latSeg,'lonSeg',lonSeg,'altSeg',altSeg);

        satPoly = computeSwathPolygon(latSeg,lonSeg,altSeg,fixed);
        satSwathPolys{si} = satPoly;
        if ~isempty(satPoly.Vertices)
            unionSwath = union(unionSwath,satPoly);
        end

        oi = elements(si).orbitIndex;
        planeSwathPolys{oi} = union(planeSwathPolys{oi},satPoly);

        obs = computeObservationIntervals(lat,lon,fixed.nepalBuffered,fixed.sampleTime);
        satObservationIntervals{si} = obs;
        allObservationIntervals = [allObservationIntervals;obs]; %#ok<AGROW>

        [glMetricsThisSat,glacierHit] = computeGlacierMetrics( ...
            lat,lon,alt,fixed.glacierNames,fixed.glacierLat, ...
            fixed.glacierLon,fixed.theta_half,fixed.sampleTime,fixed);
        satelliteGlacierMetrics{si} = glMetricsThisSat;
        glacierHitsBySat{si} = glacierHit;
    end

    if isempty(unionSwath.Vertices)
        areaKm2 = 0;
    else
        areaKm2 = area(unionSwath)*fixed.kmPerDegLat*fixed.kmPerDegLonMean;
    end

    mergedObservationIntervals = mergeIntervals(allObservationIntervals);

    % Build and merge every satellite-ground-station access interval.
    accessObjects = cell(nSat,nGS);
    allAccessIntervals = zeros(0,2);
    groundStationIntervals = cell(1,nGS);
    satelliteContactIntervals = cell(1,nSat);
    for gi = 1:nGS
        groundStationIntervals{gi} = zeros(0,2);
    end
    for si = 1:nSat
        satelliteContactIntervals{si} = zeros(0,2);
        for gi = 1:nGS
            acc = access(sat{si},gs{gi});
            accessObjects{si,gi} = acc;
            intervals = accessIntervalsToSeconds(acc,fixed.startTime);
            allAccessIntervals = [allAccessIntervals;intervals]; %#ok<AGROW>
            groundStationIntervals{gi} = [groundStationIntervals{gi};intervals]; %#ok<AGROW>
            satelliteContactIntervals{si} = [satelliteContactIntervals{si};intervals]; %#ok<AGROW>
        end
    end
    mergedAccessIntervals = mergeIntervals(allAccessIntervals);

    metrics = computeCoverageMetricsFromIntervals( ...
        mergedObservationIntervals,mergedAccessIntervals, ...
        fixed.nepalAreaKm2,areaKm2,fixed.startTime,fixed.sampleTime);

    % Constellation-wide lake metrics are logical OR across satellite hit masks.
    nSamples = numel(tracks{1}.lat);
    combinedHit = false(nSamples,numel(fixed.glacierLat));
    for si = 1:nSat
        combinedHit = combinedHit | glacierHitsBySat{si};
    end
    glacierMetrics = summarizeGlacierHits(combinedHit,fixed.glacierNames, ...
        fixed.sampleTime);

    coveredFlags = [glacierMetrics.covered];
    if all(coveredFlags)
        worstRevisitMin = max([glacierMetrics.meanRevisitMin]);
    else
        worstRevisitMin = NaN;
    end

    % Per-satellite metrics.
    satMetricsCell = cell(1,nSat);
    for si = 1:nSat
        obs = mergeIntervals(satObservationIntervals{si});
        aSat = 0;
        if ~isempty(satSwathPolys{si}.Vertices)
            aSat = area(satSwathPolys{si})*fixed.kmPerDegLat*fixed.kmPerDegLonMean;
        end
        cInt = mergeIntervals(satelliteContactIntervals{si});
        satMetricsCell{si} = summarizeEntityMetrics( ...
            obs,cInt,fixed.startTime,fixed.sampleTime, ...
            aSat,fixed.nepalAreaKm2,satelliteGlacierMetrics{si});
        satMetricsCell{si}.orbitIndex = elements(si).orbitIndex;
        satMetricsCell{si}.satelliteIndex = elements(si).satIndex;
    end
    satelliteMetrics = [satMetricsCell{:}];

    % Per-plane metrics: union geometry and merged time intervals of all sats in plane.
    planeMetricsCell = cell(1,fixed.constellation.numOrbits);
    for oi = 1:fixed.constellation.numOrbits
        memberIdx = find([elements.orbitIndex] == oi);
        obsPlane = zeros(0,2);
        contactPlane = zeros(0,2);
        glacierPlaneHit = false(nSamples,numel(fixed.glacierLat));
        for k = memberIdx
            obsPlane = [obsPlane;satObservationIntervals{k}]; %#ok<AGROW>
            contactPlane = [contactPlane;satelliteContactIntervals{k}]; %#ok<AGROW>
            glacierPlaneHit = glacierPlaneHit | glacierHitsBySat{k};
        end
        aPlane = 0;
        if ~isempty(planeSwathPolys{oi}.Vertices)
            aPlane = area(planeSwathPolys{oi})*fixed.kmPerDegLat*fixed.kmPerDegLonMean;
        end
        planeGlMetrics = summarizeGlacierHits(glacierPlaneHit,fixed.glacierNames,fixed.sampleTime);
        planeMetricsCell{oi} = summarizeEntityMetrics( ...
            mergeIntervals(obsPlane),mergeIntervals(contactPlane), ...
            fixed.startTime,fixed.sampleTime,aPlane,fixed.nepalAreaKm2,planeGlMetrics);
        planeMetricsCell{oi}.orbitIndex = oi;
        planeMetricsCell{oi}.satelliteCount = numel(memberIdx);
    end
    planeMetrics = [planeMetricsCell{:}];

    % Per-ground-station metrics.
    gsMetricsCell = cell(1,nGS);
    for gi = 1:nGS
        gsIntervals = mergeIntervals(groundStationIntervals{gi});
        gsMetricsCell{gi} = struct( ...
            'name',fixed.groundStations(gi).name, ...
            'totalContactMin',sum(gsIntervals(:,2)-gsIntervals(:,1),'all')/60, ...
            'numContacts',size(gsIntervals,1), ...
            'intervals',gsIntervals);
    end
    groundStationMetrics = [gsMetricsCell{:}];

    metrics.glacierMetrics = glacierMetrics;
    metrics.numLakesCovered = sum(coveredFlags);
    metrics.worstRevisitMin = worstRevisitMin;
    metrics.minPassesPerDayLakes = min([glacierMetrics.passesPerDay]);
    metrics.meanPassesPerDayLakes = mean([glacierMetrics.passesPerDay]);
    metrics.satelliteMetrics = satelliteMetrics;
    metrics.planeMetrics = planeMetrics;
    metrics.groundStationMetrics = groundStationMetrics;
    metrics.constellationElements = elements;
    metrics.accessIntervals = mergedAccessIntervals;
    metrics.observationIntervals = mergedObservationIntervals;
end

function [latSeg,lonSeg,altSeg] = segmentGroundTrack(lat,lon,alt,jumpThreshold)
% SEGMENTGROUNDTRACK
% Insert NaNs at longitude wrap discontinuities so plotting/polygon building
% never connects the end of one world-map segment to the next.
    jumpIdx = find(abs(diff(lon)) > jumpThreshold);
    if isempty(jumpIdx)
        latSeg = lat;
        lonSeg = lon;
        altSeg = alt;
        return;
    end

    nSeg = numel(jumpIdx)+1;
    starts = [1,jumpIdx+1];
    ends = [jumpIdx,numel(lat)];
    segLen = ends-starts+1;
    totalLen = sum(segLen)+nSeg-1;

    latSeg = nan(1,totalLen);
    lonSeg = nan(1,totalLen);
    altSeg = nan(1,totalLen);

    out = 1;
    for s = 1:nSeg
        idx = starts(s):ends(s);
        n = numel(idx);
        latSeg(out:out+n-1) = lat(idx);
        lonSeg(out:out+n-1) = lon(idx);
        altSeg(out:out+n-1) = alt(idx);
        out = out+n;
        if s < nSeg
            out = out+1;
        end
    end
end

function [glacierMetrics,hit] = computeGlacierMetrics(latRaw,lonRaw,altRaw, ...
    glacierNames,glacierLat,glacierLon,theta_half,sampleTime,fixed)
% COMPUTEGLACIERMETRICS
% Point-target test per sample and lake:
%   d^2 = ((lonLake-lonSat)/rLon)^2 + ((latLake-latSat)/rLat)^2 <= 1.
% Vectorization is across all lakes and all valid samples.

    nLakes = numel(glacierLat);
    nSamples = numel(latRaw);
    r = altRaw(:)*tan(theta_half);
    rLat = r*fixed.degPerKmLat;
    rLon = r./(fixed.kmPerDegLat*cosd(latRaw(:)));
    valid = isfinite(latRaw(:)) & isfinite(lonRaw(:)) & isfinite(altRaw(:)) & ...
        rLat > 0 & rLon > 0;

    hit = false(nSamples,nLakes);
    if any(valid)
        latV = latRaw(valid); latV = latV(:);
        lonV = lonRaw(valid); lonV = lonV(:);
        rLatV = rLat(valid); rLatV = rLatV(:);
        rLonV = rLon(valid); rLonV = rLonV(:);
        dLon = (glacierLon(:).' - lonV)./rLonV;
        dLat = (glacierLat(:).' - latV)./rLatV;
        hit(valid,:) = (dLon.^2 + dLat.^2) <= 1;
    end

    glacierMetrics = summarizeGlacierHits(hit,glacierNames,sampleTime);
end

function glacierMetrics = summarizeGlacierHits(hit,glacierNames,sampleTime)
% SUMMARIZEGLACIERHITS  Convert sample-level lake hits into pass/revisit/dwell metrics.
    nLakes = numel(glacierNames);
    nSamples = size(hit,1);
    timeVec = (0:sampleTime:(nSamples-1)*sampleTime).';
    glacierMetrics = repmat(struct('name','','passesPerDay',0, ...
        'meanRevisitMin',NaN,'maxRevisitMin',NaN,'totalDwellMin',0,'covered',false),1,nLakes);

    for g = 1:nLakes
        h = logical(hit(:,g)).';
        d = diff([false,h,false]);
        startIdx = find(d==1);
        endIdx = find(d==-1)-1;
        if isempty(startIdx)
            glacierMetrics(g).name = glacierNames{g};
            continue;
        end
        hitStartTimes = timeVec(startIdx);
        hitEndTimes = timeVec(endIdx);
        hitDurations = hitEndTimes-hitStartTimes;
        passesPerDay = numel(startIdx);
        glacierMetrics(g).name = glacierNames{g};
        glacierMetrics(g).passesPerDay = passesPerDay;
        glacierMetrics(g).totalDwellMin = sum(hitDurations)/60;
        glacierMetrics(g).covered = true;
        if passesPerDay > 1
            gaps = hitStartTimes(2:end)-hitEndTimes(1:end-1);
            glacierMetrics(g).meanRevisitMin = mean(gaps)/60;
            glacierMetrics(g).maxRevisitMin = max(gaps)/60;
        end
    end
end

function [glacierMetrics,hitBySat] = computeConstellationGlacierMetrics(tracks, ...
    glacierNames,glacierLat,glacierLon,theta_half,sampleTime,fixed)
% COMPUTECONSTELLATIONGLACIERMETRICS  Aggregate lake hits across all satellites.
    nSat = numel(tracks);
    hitBySat = cell(1,nSat);
    combined = false(numel(tracks{1}.lat),numel(glacierLat));
    for si = 1:nSat
        [~,hitBySat{si}] = computeGlacierMetrics( ...
            tracks{si}.lat,tracks{si}.lon,tracks{si}.alt,glacierNames, ...
            glacierLat,glacierLon,theta_half,sampleTime,fixed);
        combined = combined | hitBySat{si};
    end
    glacierMetrics = summarizeGlacierHits(combined,glacierNames,sampleTime);
end

function poly = computeSwathPolygon(latSeg,lonSeg,altSeg,fixed,referenceCirclePts)
% COMPUTESWATHPOLYGON  Union instantaneous footprint ellipses and clip to Nepal.
% Footprint radius on the surface:
%   r = h*tan(theta_half)
% Degree-space semi-axes:
%   rLat = r/111.32,
%   rLon = r/(111.32*cosd(latitude)).

    if nargin < 5 || isempty(referenceCirclePts)
        cCos = fixed.circleCos;
        cSin = fixed.circleSin;
    else
        theta = linspace(0,2*pi,referenceCirclePts+1);
        theta(end) = [];
        cCos = cos(theta);
        cSin = sin(theta);
    end

    valid = isfinite(latSeg) & isfinite(lonSeg) & isfinite(altSeg);
    latP = latSeg(valid); lonP = lonSeg(valid); altP = altSeg(valid);
    inRegion = lonP >= fixed.lonMin & lonP <= fixed.lonMax & ...
        latP >= fixed.latMin & latP <= fixed.latMax;
    latP = latP(inRegion); lonP = lonP(inRegion); altP = altP(inRegion);

    if isempty(lonP)
        poly = polyshape();
        return;
    end

    n = numel(lonP);
    polyArr = polyshape.empty(0,n);
    r = altP(:)*fixed.tanHalfAngle;
    rLat = r*fixed.degPerKmLat;
    rLon = r./(fixed.kmPerDegLat*cosd(latP(:)));

    warnState = warning('off','MATLAB:polyshape:repairedBySimplify');
    for i = 1:n
        polyArr(i) = polyshape( ...
            lonP(i)+rLon(i)*cCos, ...
            latP(i)+rLat(i)*cSin,'Simplify',false);
    end
    warning(warnState);
    poly = intersect(union(polyArr),fixed.nepalPoly);
end


function areaKm2 = computeConstellationSwathArea(tracks,fixed,referenceCirclePts)
% COMPUTECONSTELLATIONSWATHAREA  Recompute union coverage at a requested
% footprint polygon resolution for the final numerical sensitivity test.
    unionPoly = polyshape();
    for si = 1:numel(tracks)
        p = computeSwathPolygon(tracks{si}.latSeg,tracks{si}.lonSeg, ...
            tracks{si}.altSeg,fixed,referenceCirclePts);
        if ~isempty(p.Vertices)
            unionPoly = union(unionPoly,p);
        end
    end
    if isempty(unionPoly.Vertices)
        areaKm2 = 0;
    else
        areaKm2 = area(unionPoly)*fixed.kmPerDegLat*fixed.kmPerDegLonMean;
    end
end

function intervals = computeObservationIntervals(lat,lon,nepalBuffered,sampleTime)
% OBSERVATIONINTERVALS  Return [start,end] seconds when the buffered Nepal
% target is under the satellite nadir track.
    timeVec = (0:sampleTime:(numel(lat)-1)*sampleTime).';
    inside = inpolygon(lon(:),lat(:),nepalBuffered.Vertices(:,1),nepalBuffered.Vertices(:,2));
    d = diff([false;inside;false]);
    starts = find(d==1); ends = find(d==-1)-1;
    if isempty(starts)
        intervals = zeros(0,2);
    else
        intervals = [timeVec(starts),timeVec(ends)];
    end
end

function intervals = accessIntervalsToSeconds(acc,startTime)
% ACCESSINTERVALSTOSECONDS  Convert MATLAB satelliteScenario access windows to
% seconds from the mission start.
    ai = accessIntervals(acc);
    if isempty(ai)
        intervals = zeros(0,2);
        return;
    end
    refUTC = datetime(startTime,'TimeZone','UTC');
    t0 = posixtime(refUTC);
    starts = posixtime(ai.StartTime)-t0;
    ends = posixtime(ai.EndTime)-t0;
    intervals = [starts(:),ends(:)];
end

function merged = mergeIntervals(intervals)
% MERGEINTERVALS  Sort and union overlapping/touching time intervals.
    if isempty(intervals)
        merged = zeros(0,2);
        return;
    end
    intervals = sortrows(intervals,[1 2]);
    merged = zeros(size(intervals,1),2);
    k = 1;
    merged(1,:) = intervals(1,:);
    for i = 2:size(intervals,1)
        if intervals(i,1) <= merged(k,2)
            merged(k,2) = max(merged(k,2),intervals(i,2));
        else
            k = k+1;
            merged(k,:) = intervals(i,:);
        end
    end
    merged = merged(1:k,:);
end

function metrics = computeCoverageMetricsFromIntervals(observationInt,accessInt, ...
    nepalAreaKm2,areaKm2,~,~)
% COMPUTECOVERAGEMETRICSFROMINTERVALS  Compute constellation-wide temporal metrics
% from already merged observation and access intervals.
    %#ok<INUSD>
    if isempty(observationInt)
        overlapDur = 0;
        avgRevisitMin = NaN;
        maxRevisitMin = NaN;
    else
        overlapDurations = observationInt(:,2)-observationInt(:,1);
        overlapDur = sum(overlapDurations);
        if size(observationInt,1) > 1
            gaps = observationInt(2:end,1)-observationInt(1:end-1,2);
            avgRevisitMin = mean(gaps)/60;
            maxRevisitMin = max(gaps)/60;
        else
            avgRevisitMin = NaN;
            maxRevisitMin = NaN;
        end
    end

    if isempty(accessInt)
        sgsTotalMin = 0;
        sgsAvgMin = NaN;
        sgsCount = 0;
    else
        d = accessInt(:,2)-accessInt(:,1);
        sgsTotalMin = sum(d)/60;
        sgsAvgMin = mean(d)/60;
        sgsCount = size(accessInt,1);
    end

    latency = nan(size(observationInt,1),1);
    if ~isempty(accessInt) && ~isempty(observationInt)
        sgsStarts = accessInt(:,1);
        for i = 1:size(observationInt,1)
            k = find(sgsStarts >= observationInt(i,2),1,'first');
            if ~isempty(k)
                latency(i) = sgsStarts(k)-observationInt(i,2);
            end
        end
    end
    latency = latency(isfinite(latency));
    if isempty(latency)
        avgLatencyMin = NaN;
    else
        avgLatencyMin = mean(latency)/60;
    end

    metrics = struct( ...
        'passesPerDay',size(observationInt,1), ...
        'avgRevisitMin',avgRevisitMin, ...
        'maxRevisitMin',maxRevisitMin, ...
        'sgsTotalMin',sgsTotalMin, ...
        'sgsAvgMin',sgsAvgMin, ...
        'sgsCount',sgsCount, ...
        'overlapTotalMin',overlapDur/60, ...
        'overlapPct',100*overlapDur/(24*60), ...
        'coverageFraction',areaKm2/nepalAreaKm2, ...
        'avgLatencyMin',avgLatencyMin, ...
        'accumAreaKm2',areaKm2);
end

function out = summarizeEntityMetrics(obs,contact,startTime,sampleTime,areaKm2,nepalAreaKm2,glacierMetrics)
% SUMMARIZEENTITYMETRICS  Per-satellite/per-plane temporal and coverage metrics.
    out = computeCoverageMetricsFromIntervals(obs,contact,nepalAreaKm2,areaKm2,startTime,sampleTime);
    out.glacierMetrics = glacierMetrics;
end

function [metrics,details] = computeConstellationMetrics(tracks,accObjects,groundStations, ...
    startTime,sampleTime,nepalAreaKm2,areaKm2,nepalBuffered,~)
% COMPUTECONSTELLATIONMETRICS  High-resolution constellation-wide metrics.
    %#ok<INUSD>
    nSat = numel(tracks);
    nGS = numel(groundStations);
    allObs = zeros(0,2);
    for si = 1:nSat
        allObs = [allObs;computeObservationIntervals(tracks{si}.lat,tracks{si}.lon, ...
            nepalBuffered,sampleTime)]; %#ok<AGROW>
    end
    allObs = mergeIntervals(allObs);

    allAccess = zeros(0,2);
    perGS = cell(1,nGS);
    perSat = cell(1,nSat);
    for gi = 1:nGS, perGS{gi}=zeros(0,2); end
    for si = 1:nSat, perSat{si}=zeros(0,2); end
    for si = 1:nSat
        for gi = 1:nGS
            iv = accessIntervalsToSeconds(accObjects{si,gi},startTime);
            allAccess = [allAccess;iv]; %#ok<AGROW>
            perGS{gi} = [perGS{gi};iv]; %#ok<AGROW>
            perSat{si} = [perSat{si};iv]; %#ok<AGROW>
        end
    end
    allAccess = mergeIntervals(allAccess);
    metrics = computeCoverageMetricsFromIntervals(allObs,allAccess, ...
        nepalAreaKm2,areaKm2,startTime,sampleTime);

    details = struct('observationIntervals',allObs,'accessIntervals',allAccess, ...
        'groundStationIntervals',cell(1,nGS),'satelliteContactIntervals',cell(1,nSat));
    for gi = 1:nGS
        details.groundStationIntervals{gi} = mergeIntervals(perGS{gi});
    end
    for si = 1:nSat
        details.satelliteContactIntervals{si} = mergeIntervals(perSat{si});
    end
end

function [passShapes,passAreas,passIndexRanges] = computePerPassSwathShapesIdx( ...
    latSeg,lonSeg,altSeg,~,fixed)
% COMPUTEPERPASSSWATHSHAPESIDX  Build one clipped swath polygon per orbital segment.
    segStarts = [1,find(isnan(latSeg))+1];
    segEnds = [find(isnan(latSeg))-1,numel(latSeg)];
    passShapes = cell(1,0);
    passAreas = zeros(1,0);
    passIndexRanges = zeros(0,2);

    for s = 1:numel(segStarts)
        if segEnds(s) < segStarts(s), continue; end
        idx = segStarts(s):segEnds(s);
        passPoly = computeSwathPolygon(latSeg(idx),lonSeg(idx),altSeg(idx),fixed);
        if isempty(passPoly.Vertices), continue; end
        aS = area(passPoly)*fixed.kmPerDegLat*fixed.kmPerDegLonMean;
        if aS > 0
            passShapes{end+1} = passPoly; %#ok<AGROW>
            passAreas(end+1) = aS; %#ok<AGROW>
            passIndexRanges(end+1,:) = [segStarts(s),segEnds(s)]; %#ok<AGROW>
        end
    end
end

function score = computeScore(metrics,weights,scales)
% COMPUTESCORE
% Normalize each metric by its scale, apply its signed mission weight, and
% sum the terms. Missing revisit/latency values receive a large penalty.

    revisit = metrics.avgRevisitMin;
    if isnan(revisit)
        revisit = 3*scales.avgRevisitMin;
    end
    latency = metrics.avgLatencyMin;
    if isnan(latency)
        latency = 3*scales.avgLatencyMin;
    end
    glacierWorstRevisit = metrics.worstRevisitMin;
    if isnan(glacierWorstRevisit)
        glacierWorstRevisit = 3*scales.worstRevisitMin;
    end

    score = ...
        weights.coverageFraction*(metrics.coverageFraction/scales.coverageFraction) + ...
        weights.passesPerDay*(metrics.passesPerDay/scales.passesPerDay) + ...
        weights.sgsTotalMin*(metrics.sgsTotalMin/scales.sgsTotalMin) + ...
        weights.avgRevisitMin*(revisit/scales.avgRevisitMin) + ...
        weights.avgLatencyMin*(latency/scales.avgLatencyMin) + ...
        weights.numLakesCovered*(metrics.numLakesCovered/scales.numLakesCovered) + ...
        weights.worstRevisitMin*(glacierWorstRevisit/scales.worstRevisitMin) + ...
        weights.minPassesPerDayLakes*(metrics.minPassesPerDayLakes/scales.minPassesPerDayLakes);
end

function ei = expectedImprovement(mu,sigma,bestSoFar)
% EXPECTEDIMPROVEMENT  Vectorized EI formula used by Bayesian optimization.
    sigma = max(sigma,1e-9);
    delta = mu-bestSoFar;
    z = delta./sigma;
    ei = delta.*normcdf(z) + sigma.*normpdf(z);
    ei(sigma <= 1e-8) = 0;
end

function val = negativeExpectedImprovement(gprMdl,x,bestSoFar)
% NEGATIVEEXPECTEDIMPROVEMENT  Single predict() call per fmincon evaluation.
    [mu,sigma] = predict(gprMdl,x);
    val = -expectedImprovement(mu,sigma,bestSoFar);
end

function printMetricsBlock(m,perPassAreas)
% PRINTMETRICSBLOCK  Preserve the established console report format.
    fprintf('Revisit time with swath:\n');
    if ~isnan(m.avgRevisitMin)
        fprintf('  Average: %.2f min\n',m.avgRevisitMin);
        fprintf('  Maximum: %.2f min\n',m.maxRevisitMin);
    else
        fprintf('  N/A (less than 2 passes)\n');
    end
    fprintf('Number of passes per day within swath: %d\n',m.passesPerDay);
    fprintf('SGS contact:\n');
    if m.sgsCount > 0
        fprintf('  Total contact time: %.2f min\n',m.sgsTotalMin);
        fprintf('  Average per contact: %.2f min\n',m.sgsAvgMin);
        fprintf('  Number of contacts: %d\n',m.sgsCount);
    else
        fprintf('  No SGS contacts\n');
    end
    fprintf('Swath overlap time:\n');
    fprintf('  Total: %.2f min (%.2f%% of day)\n',m.overlapTotalMin,m.overlapPct);
    fprintf('Coverage fraction (accumulated area / Nepal area): %.4f\n',m.coverageFraction);
    fprintf('Data latency (average time from observation to next SGS contact):\n');
    if ~isnan(m.avgLatencyMin)
        fprintf('  Average: %.2f min\n',m.avgLatencyMin);
    else
        fprintf('  N/A\n');
    end

    if ~isempty(perPassAreas)
        for p = 1:numel(perPassAreas)
            fprintf('Pass %d coverage = %.2f km^2\n',p,perPassAreas(p));
        end
        totalAreaKm2 = sum(perPassAreas);
        termsStr = strjoin(arrayfun(@(x)sprintf('%.2f',x), ...
            perPassAreas,'UniformOutput',false),' + ');
        fprintf(['Accumulated coverage over 24 h (Σ per-pass areas, overlaps counted) ' ...
            '= %s = %.2f km^2\n'],termsStr,totalAreaKm2);
        fprintf('Unique union coverage = %.2f km^2  (overlap redundancy = %.2f km^2)\n', ...
            m.accumAreaKm2,totalAreaKm2-m.accumAreaKm2);
    else
        fprintf('Accumulated coverage over 24 h: %.2f km^2\n',m.accumAreaKm2);
    end
end

function s = formatHMS(tSec)
% FORMATHMS  Format seconds as H hr M min S sec.
    h = floor(tSec/3600);
    m = floor(mod(tSec,3600)/60);
    sec = floor(mod(tSec,60));
    if h > 0
        s = sprintf('%d hr %d min %d sec',h,m,sec);
    elseif m > 0
        s = sprintf('%d min %d sec',m,sec);
    else
        s = sprintf('%d sec',sec);
    end
end

function metrics = evaluateElementsList(elements,fixed)
% EVALUATEELEMENTSLIST  General constellation evaluator for an explicit,
% possibly heterogeneous list of satellites (each with its own
% inclination/RAAN/altitude/true-anomaly, grouped by ORBITINDEX). Mirrors
% RUNSIM's computation path exactly, but reads each satellite's own
% SEMIMAJORAXIS field instead of one global value, and infers the set of
% orbital planes from the elements themselves instead of
% FIXED.CONSTELLATION. Used only by the incremental constellation-growth
% loop (Section I); RUNSIM and the validated 1x1x1 baseline (Sections A-H)
% are left completely unmodified.

    nSat = numel(elements);
    nGS = numel(fixed.groundStations);
    orbitIdxList = unique([elements.orbitIndex]);
    maxOrbitIdx = max(orbitIdxList);

    sc = satelliteScenario(fixed.startTime,fixed.stopTime,fixed.sampleTime);
    sat = cell(1,nSat);
    for si = 1:nSat
        e = elements(si);
        sat{si} = satellite(sc,e.semiMajorAxis,fixed.eccentricity, ...
            e.inclination,e.raan,0,e.trueAnomaly, ...
            Name=sprintf('GrowSat_P%02d_S%02d',e.orbitIndex,e.satIndex));
    end

    gs = cell(1,nGS);
    for gi = 1:nGS
        gs{gi} = groundStation(sc,fixed.groundStations(gi).lat, ...
            fixed.groundStations(gi).lon,Name=fixed.groundStations(gi).name);
    end

    tracks = cell(1,nSat);
    satSwathPolys = cell(1,nSat);
    planeSwathPolys = cell(1,maxOrbitIdx);
    for oi = orbitIdxList
        planeSwathPolys{oi} = polyshape();
    end

    satObservationIntervals = cell(1,nSat);
    satelliteGlacierMetrics = cell(1,nSat);
    glacierHitsBySat = cell(1,nSat);

    unionSwath = polyshape();
    allObservationIntervals = zeros(0,2);

    for si = 1:nSat
        [pos,~,~] = states(sat{si},'CoordinateFrame','geographic');
        lat = pos(1,:); lon = pos(2,:); alt = pos(3,:)/1000;
        valid = isfinite(lat) & isfinite(lon) & isfinite(alt);
        lat = lat(valid); lon = lon(valid); alt = alt(valid);

        [latSeg,lonSeg,altSeg] = segmentGroundTrack(lat,lon,alt,10);
        tracks{si} = struct('lat',lat,'lon',lon,'alt',alt, ...
            'latSeg',latSeg,'lonSeg',lonSeg,'altSeg',altSeg);

        satPoly = computeSwathPolygon(latSeg,lonSeg,altSeg,fixed);
        satSwathPolys{si} = satPoly;
        if ~isempty(satPoly.Vertices)
            unionSwath = union(unionSwath,satPoly);
        end

        oi = elements(si).orbitIndex;
        planeSwathPolys{oi} = union(planeSwathPolys{oi},satPoly);

        obs = computeObservationIntervals(lat,lon,fixed.nepalBuffered,fixed.sampleTime);
        satObservationIntervals{si} = obs;
        allObservationIntervals = [allObservationIntervals;obs]; %#ok<AGROW>

        [glMetricsThisSat,glacierHit] = computeGlacierMetrics( ...
            lat,lon,alt,fixed.glacierNames,fixed.glacierLat, ...
            fixed.glacierLon,fixed.theta_half,fixed.sampleTime,fixed);
        satelliteGlacierMetrics{si} = glMetricsThisSat;
        glacierHitsBySat{si} = glacierHit;
    end

    if isempty(unionSwath.Vertices)
        areaKm2 = 0;
    else
        areaKm2 = area(unionSwath)*fixed.kmPerDegLat*fixed.kmPerDegLonMean;
    end

    mergedObservationIntervals = mergeIntervals(allObservationIntervals);

    allAccessIntervals = zeros(0,2);
    groundStationIntervals = cell(1,nGS);
    satelliteContactIntervals = cell(1,nSat);
    for gi = 1:nGS
        groundStationIntervals{gi} = zeros(0,2);
    end
    for si = 1:nSat
        satelliteContactIntervals{si} = zeros(0,2);
        for gi = 1:nGS
            acc = access(sat{si},gs{gi});
            intervalsAcc = accessIntervalsToSeconds(acc,fixed.startTime);
            allAccessIntervals = [allAccessIntervals;intervalsAcc]; %#ok<AGROW>
            groundStationIntervals{gi} = [groundStationIntervals{gi};intervalsAcc]; %#ok<AGROW>
            satelliteContactIntervals{si} = [satelliteContactIntervals{si};intervalsAcc]; %#ok<AGROW>
        end
    end
    mergedAccessIntervals = mergeIntervals(allAccessIntervals);

    metrics = computeCoverageMetricsFromIntervals( ...
        mergedObservationIntervals,mergedAccessIntervals, ...
        fixed.nepalAreaKm2,areaKm2,fixed.startTime,fixed.sampleTime);

    nSamples = numel(tracks{1}.lat);
    combinedHit = false(nSamples,numel(fixed.glacierLat));
    for si = 1:nSat
        combinedHit = combinedHit | glacierHitsBySat{si};
    end
    glacierMetrics = summarizeGlacierHits(combinedHit,fixed.glacierNames, ...
        fixed.sampleTime);

    coveredFlags = [glacierMetrics.covered];
    if all(coveredFlags)
        worstRevisitMin = max([glacierMetrics.meanRevisitMin]);
    else
        worstRevisitMin = NaN;
    end

    satMetricsCell = cell(1,nSat);
    for si = 1:nSat
        obs = mergeIntervals(satObservationIntervals{si});
        aSat = 0;
        if ~isempty(satSwathPolys{si}.Vertices)
            aSat = area(satSwathPolys{si})*fixed.kmPerDegLat*fixed.kmPerDegLonMean;
        end
        cInt = mergeIntervals(satelliteContactIntervals{si});
        satMetricsCell{si} = summarizeEntityMetrics( ...
            obs,cInt,fixed.startTime,fixed.sampleTime, ...
            aSat,fixed.nepalAreaKm2,satelliteGlacierMetrics{si});
        satMetricsCell{si}.orbitIndex = elements(si).orbitIndex;
        satMetricsCell{si}.satelliteIndex = elements(si).satIndex;
        satMetricsCell{si}.observationIntervals = obs;
    end
    satelliteMetrics = [satMetricsCell{:}];

    planeMetricsCell = cell(1,numel(orbitIdxList));
    for pIdx = 1:numel(orbitIdxList)
        oi = orbitIdxList(pIdx);
        memberIdx = find([elements.orbitIndex] == oi);
        obsPlane = zeros(0,2);
        contactPlane = zeros(0,2);
        glacierPlaneHit = false(nSamples,numel(fixed.glacierLat));
        for mIdx = memberIdx
            obsPlane = [obsPlane;satObservationIntervals{mIdx}]; %#ok<AGROW>
            contactPlane = [contactPlane;satelliteContactIntervals{mIdx}]; %#ok<AGROW>
            glacierPlaneHit = glacierPlaneHit | glacierHitsBySat{mIdx};
        end
        aPlane = 0;
        if ~isempty(planeSwathPolys{oi}.Vertices)
            aPlane = area(planeSwathPolys{oi})*fixed.kmPerDegLat*fixed.kmPerDegLonMean;
        end
        planeGlMetrics = summarizeGlacierHits(glacierPlaneHit,fixed.glacierNames,fixed.sampleTime);
        planeMetricsCell{pIdx} = summarizeEntityMetrics( ...
            mergeIntervals(obsPlane),mergeIntervals(contactPlane), ...
            fixed.startTime,fixed.sampleTime,aPlane,fixed.nepalAreaKm2,planeGlMetrics);
        planeMetricsCell{pIdx}.orbitIndex = oi;
        planeMetricsCell{pIdx}.satelliteCount = numel(memberIdx);
    end
    planeMetrics = [planeMetricsCell{:}];

    gsMetricsCell = cell(1,nGS);
    for gi = 1:nGS
        gsIntervals = mergeIntervals(groundStationIntervals{gi});
        gsMetricsCell{gi} = struct( ...
            'name',fixed.groundStations(gi).name, ...
            'totalContactMin',sum(gsIntervals(:,2)-gsIntervals(:,1),'all')/60, ...
            'numContacts',size(gsIntervals,1), ...
            'intervals',gsIntervals);
    end
    groundStationMetrics = [gsMetricsCell{:}];

    metrics.glacierMetrics = glacierMetrics;
    metrics.numLakesCovered = sum(coveredFlags);
    metrics.worstRevisitMin = worstRevisitMin;
    metrics.minPassesPerDayLakes = min([glacierMetrics.passesPerDay]);
    metrics.meanPassesPerDayLakes = mean([glacierMetrics.passesPerDay]);
    metrics.satelliteMetrics = satelliteMetrics;
    metrics.planeMetrics = planeMetrics;
    metrics.groundStationMetrics = groundStationMetrics;
    metrics.constellationElements = elements;
    metrics.accessIntervals = mergedAccessIntervals;
    metrics.observationIntervals = mergedObservationIntervals;
end

function [J,result] = evaluateSameOrbitCandidate(x,baseElements,orbitIdx,fixed, ...
    growthWeights,growthScales,periodSec,collisionOffsetKm,altitudeRangeKm,semiMajorAxis0)
% EVALUATESAMEORBITCANDIDATE
% x = [trueAnomalyDeg, altitudeOffsetKm]. Evaluates adding one satellite to
% orbital plane ORBITIDX, enforcing a minimum altitude separation
% (COLLISIONOFFSETKM) from every satellite already sharing that plane.
    trueAnomalyDeg = mod(x(1),360);
    altOffsetKm = x(2);

    members = find([baseElements.orbitIndex] == orbitIdx);
    existingOffsets = [baseElements(members).altitudeOffsetKm];
    altOffsetKm = nudgeToSatisfyCollision(altOffsetKm,existingOffsets, ...
        collisionOffsetKm,altitudeRangeKm);

    newSat = struct( ...
        'inclination',baseElements(members(1)).inclination, ...
        'raan',baseElements(members(1)).raan, ...
        'trueAnomaly',trueAnomalyDeg, ...
        'semiMajorAxis',semiMajorAxis0+altOffsetKm*1e3, ...
        'orbitIndex',orbitIdx, ...
        'satIndex',max([baseElements(members).satIndex])+1, ...
        'altitudeOffsetKm',altOffsetKm);

    combinedElements = [baseElements,newSat];
    combinedMetrics = evaluateElementsList(combinedElements,fixed);
    newSatMetrics = combinedMetrics.satelliteMetrics(end);

    existingObs = zeros(0,2);
    for k = 1:numel(baseElements)
        existingObs = [existingObs;combinedMetrics.satelliteMetrics(k).observationIntervals]; %#ok<AGROW>
    end
    existingObs = mergeIntervals(existingObs);
    wass = wassersteinPassTimeDistance(existingObs,newSatMetrics.observationIntervals,periodSec);

    oldGapHr = largestScheduleGapHours(existingObs,periodSec);
    newGapHr = largestScheduleGapHours(combinedMetrics.observationIntervals,periodSec);
    gapReductionHr = oldGapHr-newGapHr;

    J = computeGrowthObjective(combinedMetrics,wass,gapReductionHr,growthWeights,growthScales);

    result = struct( ...
        'action','AddSatellite', ...
        'orbitIndex',orbitIdx, ...
        'inclination',newSat.inclination, ...
        'raan',newSat.raan, ...
        'trueAnomalyDeg',trueAnomalyDeg, ...
        'altitudeOffsetKm',altOffsetKm, ...
        'combinedElements',combinedElements, ...
        'combinedMetrics',combinedMetrics, ...
        'newSatMetrics',newSatMetrics, ...
        'wasserstein',wass, ...
        'oldGapHr',oldGapHr,'newGapHr',newGapHr,'gapReductionHr',gapReductionHr, ...
        'J',J);
end

function [J,result] = evaluateNewOrbitCandidate(x,baseElements,fixed, ...
    growthWeights,growthScales,periodSec,altitudeRangeKm,semiMajorAxis0, ...
    nInnerSamples,newOrbitIndex)
% EVALUATENEWORBITCANDIDATE
% x = [inclinationDeg, raanDeg]. Places the first satellite of this new
% orbital plane with a small deterministic LHS sweep over
% [trueAnomaly, altitudeOffset] (no prior satellite in the plane to avoid
% colliding with, so a lightweight nested search is sufficient) and keeps
% the best-scoring placement found.
    inclDeg = min(x(1),179.9);  % SGP4/SDP4 rejects inclination == 180 deg exactly
    raanDeg = mod(x(2),360);

    rng(round(inclDeg*100)+round(raanDeg*100)+newOrbitIndex,'twister');
    lbIn = [0,-altitudeRangeKm];
    ubIn = [360, altitudeRangeKm];
    lhsU = lhsdesign(nInnerSamples,2);
    Xin = lbIn+lhsU.*(ubIn-lbIn);

    bestJin = -Inf;
    bestResIn = [];
    for k = 1:nInnerSamples
        nuK = Xin(k,1);
        altOffK = Xin(k,2);
        newSat = struct( ...
            'inclination',inclDeg, ...
            'raan',raanDeg, ...
            'trueAnomaly',nuK, ...
            'semiMajorAxis',semiMajorAxis0+altOffK*1e3, ...
            'orbitIndex',newOrbitIndex, ...
            'satIndex',1, ...
            'altitudeOffsetKm',altOffK);
        combinedElements = [baseElements,newSat];
        combinedMetrics = evaluateElementsList(combinedElements,fixed);
        newSatMetrics = combinedMetrics.satelliteMetrics(end);

        existingObs = zeros(0,2);
        for kk = 1:numel(baseElements)
            existingObs = [existingObs;combinedMetrics.satelliteMetrics(kk).observationIntervals]; %#ok<AGROW>
        end
        existingObs = mergeIntervals(existingObs);
        wass = wassersteinPassTimeDistance(existingObs,newSatMetrics.observationIntervals,periodSec);

        oldGapHr = largestScheduleGapHours(existingObs,periodSec);
        newGapHr = largestScheduleGapHours(combinedMetrics.observationIntervals,periodSec);
        gapReductionHr = oldGapHr-newGapHr;

        Jk = computeGrowthObjective(combinedMetrics,wass,gapReductionHr,growthWeights,growthScales);

        if Jk > bestJin
            bestJin = Jk;
            bestResIn = struct( ...
                'action','AddOrbit', ...
                'orbitIndex',newOrbitIndex, ...
                'inclination',inclDeg, ...
                'raan',raanDeg, ...
                'trueAnomalyDeg',nuK, ...
                'altitudeOffsetKm',altOffK, ...
                'combinedElements',combinedElements, ...
                'combinedMetrics',combinedMetrics, ...
                'newSatMetrics',newSatMetrics, ...
                'wasserstein',wass, ...
                'oldGapHr',oldGapHr,'newGapHr',newGapHr,'gapReductionHr',gapReductionHr, ...
                'J',Jk);
        end
    end
    J = bestJin;
    result = bestResIn;
end

function offsetKm = nudgeToSatisfyCollision(offsetKm,existingOffsetsKm,minSepKm,rangeKm)
% NUDGETOSATISFYCOLLISION
% Move a candidate altitude offset to the nearest feasible point that is
% at least MINSEPKM from every existing satellite sharing the plane,
% staying within [-RANGEKM, RANGEKM]. Candidate breakpoints are generated
% at +/- MINSEPKM around each existing offset and at the search-range
% edges; the closest feasible breakpoint to the original candidate is
% returned. This direct search (rather than iteratively nudging toward
% whichever satellite is currently nearest) avoids oscillating forever
% when two existing satellites are closer together than 2*MINSEPKM.
    if isempty(existingOffsetsKm)
        return;
    end
    tol = 1e-9;
    isFeasible = @(v) all(abs(v-existingOffsetsKm) >= minSepKm-tol);

    if isFeasible(offsetKm)
        return;
    end

    original = offsetKm;
    breakpoints = [-rangeKm,rangeKm, ...
        existingOffsetsKm(:).'+minSepKm, ...
        existingOffsetsKm(:).'-minSepKm];
    breakpoints = min(max(breakpoints,-rangeKm),rangeKm);

    feasibleMask = false(size(breakpoints));
    for i = 1:numel(breakpoints)
        feasibleMask(i) = isFeasible(breakpoints(i));
    end
    candidates = breakpoints(feasibleMask);

    if isempty(candidates)
        % No point in range clears every separation: keep the original
        % rather than moving to an arbitrarily chosen infeasible point.
        offsetKm = original;
        return;
    end

    [~,ix] = min(abs(candidates-original));
    offsetKm = candidates(ix);
end

function d = wassersteinPassTimeDistance(intervalsA,intervalsB,periodSec)
% WASSERSTEINPASSTIMEDISTANCE
% 1-D earth-mover (Wasserstein-1) distance between the pass-time-of-day
% distributions of two sets of observation intervals. Each interval
% contributes one sample: its midpoint time, wrapped into [0,periodSec).
% The distance is the area between the two empirical CDFs on their merged
% support, normalized by periodSec/2 so typical values fall roughly in
% [0,1] (higher = more temporally complementary).
    if isempty(intervalsA) || isempty(intervalsB)
        d = NaN;
        return;
    end
    tA = sort(mod(mean(intervalsA,2),periodSec));
    tB = sort(mod(mean(intervalsB,2),periodSec));
    grid = unique([0;tA(:);tB(:);periodSec]);
    nA = numel(tA);
    nB = numel(tB);
    FA = arrayfun(@(t) sum(tA<=t)/nA,grid);
    FB = arrayfun(@(t) sum(tB<=t)/nB,grid);
    d = sum(abs(FA(1:end-1)-FB(1:end-1)).*diff(grid));
    d = d/(periodSec/2);
end

function J = computeGrowthObjective(metrics,wassersteinDist,gapReductionHr,weights,scales)
% COMPUTEGROWTHOBJECTIVE
% Extends COMPUTESCORE with a temporal-complementarity bonus and a direct
% largest-gap-reduction bonus:
%   J = missionScore(metrics) + wW*(Dw/scaleW) + wG*(Dgap/scaleG)
% A NaN Wasserstein distance (fewer than one candidate pass) contributes no
% bonus rather than being penalized; likewise a NaN gap reduction. GAPRE-
% DUCTIONHR is the largest daily pass-time gap (hours) BEFORE the candidate
% minus the gap AFTER it, so a positive value (the gap shrank) is rewarded
% and a negative value (the gap grew) is penalized — a genuine, direct
% criterion targeting exactly the largest empty stretch, rather than the
% overall-distribution proxy the Wasserstein term provides.
    J = computeScore(metrics,weights,scales);
    if ~isnan(wassersteinDist) && isfield(weights,'wasserstein') && isfield(scales,'wasserstein')
        J = J + weights.wasserstein*(wassersteinDist/scales.wasserstein);
    end
    if ~isnan(gapReductionHr) && isfield(weights,'gap') && isfield(scales,'gap')
        J = J + weights.gap*(gapReductionHr/scales.gap);
    end
end

function [accept,reasonStr,gains] = decideGrowthStep(oldMetrics,newMetrics,cfg)
% DECIDEGROWTHSTEP
% Accept a candidate satellite addition if it clears any one of the
% user-defined improvement thresholds on coverage, revisit, latency, or
% the largest daily pass-time gap (a direct, wraparound-aware measure of
% the single biggest blank stretch in the Nepal pass schedule — see
% LARGESTSCHEDULEGAPHOURS — tracked here as its own criterion rather than
% relying solely on the softer Wasserstein distributional proxy).
    coverageGain = newMetrics.coverageFraction-oldMetrics.coverageFraction;

    oldRevisit = oldMetrics.avgRevisitMin;
    newRevisit = newMetrics.avgRevisitMin;
    if isnan(oldRevisit), oldRevisit = Inf; end
    if isnan(newRevisit), newRevisit = Inf; end
    revisitGainMin = oldRevisit-newRevisit;

    oldLatency = oldMetrics.avgLatencyMin;
    newLatency = newMetrics.avgLatencyMin;
    if isnan(oldLatency), oldLatency = Inf; end
    if isnan(newLatency), newLatency = Inf; end
    latencyGainMin = oldLatency-newLatency;

    oldGapHr = largestScheduleGapHours(oldMetrics.observationIntervals,cfg.periodSec);
    newGapHr = largestScheduleGapHours(newMetrics.observationIntervals,cfg.periodSec);
    gapReductionHr = oldGapHr-newGapHr;

    gains = struct('coverageGain',coverageGain,'revisitGainMin',revisitGainMin, ...
        'latencyGainMin',latencyGainMin, ...
        'passesGain',newMetrics.passesPerDay-oldMetrics.passesPerDay, ...
        'oldGapHr',oldGapHr,'newGapHr',newGapHr,'gapReductionHr',gapReductionHr);

    reasons = {};
    accept = false;
    if coverageGain >= cfg.coverageGainThreshold
        accept = true;
        reasons{end+1} = sprintf('coverage +%.1f%%',100*coverageGain);
    end
    if revisitGainMin >= cfg.revisitGainThresholdMin
        accept = true;
        reasons{end+1} = sprintf('revisit -%.1f min',revisitGainMin);
    end
    if latencyGainMin >= cfg.latencyGainThresholdMin
        accept = true;
        reasons{end+1} = sprintf('latency -%.1f min',latencyGainMin);
    end
    if gapReductionHr >= cfg.gapGainThresholdHr
        accept = true;
        reasons{end+1} = sprintf('largest gap -%.2f hr',gapReductionHr);
    end

    if accept
        reasonStr = ['Accept (',strjoin(reasons,', '),')'];
    else
        reasonStr = sprintf(['Reject (insufficient improvement: coverage %+.1f%%, ' ...
            'revisit %+.1f min, latency %+.1f min, largest gap %+.2f hr)'], ...
            100*coverageGain,revisitGainMin,latencyGainMin,gapReductionHr);
    end
end

function ok = passesWorstRevisitGuard(oldMetrics,newMetrics,cfg)
% PASSESWORSTREVISITGUARD
% Hard veto used ahead of DECIDEGROWTHSTEP for both the inner- and outer-BO
% candidate pools: blocks any candidate that would let the worst-case
% (per-lake) revisit time regress by more than
% CFG.WORSTREVISITREGRESSIONTOLERANCEMIN, or that would drop lake coverage
% back below 6/6 when it was previously complete. Going from "not all
% lakes covered" (NaN) to "all covered" (finite) is always an improvement,
% never a regression, and is always allowed through.
    oldW = oldMetrics.worstRevisitMin;
    newW = newMetrics.worstRevisitMin;

    if isnan(oldW) && isnan(newW)
        ok = true;                  % no information either way
    elseif isnan(oldW) && ~isnan(newW)
        ok = true;                  % newly achieved full lake coverage
    elseif ~isnan(oldW) && isnan(newW)
        ok = false;                 % lost full lake coverage
    else
        ok = (newW-oldW) <= cfg.worstRevisitRegressionToleranceMin;
    end
end

function printGrowthStepRow(stepNum,res,decisionStr)
% PRINTGROWTHSTEPROW  Console row in the step-wise report format described
% in the design brief: candidate parameters, new-satellite-only metrics,
% combined constellation metrics, the Wasserstein score, the largest
% pass-time gap before/after (so its shrinkage can be watched step by
% step), and the decision.
    if strcmp(res.action,'AddSatellite')
        actionStr = 'Add Satellite (same orbit)';
    else
        actionStr = 'Add New Orbit';
    end
    fprintf('--- Step %d: %s ---\n',stepNum,actionStr);
    fprintf('  Candidate: incl=%.2f deg, RAAN=%.2f deg, nu=%.2f deg, altitude offset=%+.1f km\n', ...
        res.inclination,res.raan,res.trueAnomalyDeg,res.altitudeOffsetKm);
    nsm = res.newSatMetrics;
    fprintf('  New-sat metrics (alone): revisit=%s min, covered area=%.1f km^2\n', ...
        formatMaybeNaN(nsm.avgRevisitMin),nsm.accumAreaKm2);
    cm = res.combinedMetrics;
    fprintf(['  Combined constellation: coverage=%.2f%%, revisit=%s min, passes/day=%d, ' ...
        'SGS contact=%.1f min, latency=%s min\n'], ...
        100*cm.coverageFraction,formatMaybeNaN(cm.avgRevisitMin), ...
        cm.passesPerDay,cm.sgsTotalMin,formatMaybeNaN(cm.avgLatencyMin));
    fprintf('  Wasserstein complementarity: %s\n',formatMaybeNaN(res.wasserstein));
    fprintf('  Largest pass-time gap: %.2f hr -> %.2f hr (%+.2f hr)\n', ...
        res.oldGapHr,res.newGapHr,res.gapReductionHr);
    fprintf('  Decision: %s\n\n',decisionStr);
end

function s = formatMaybeNaN(v)
% FORMATMAYBENAN  Render a scalar as "%.2f", or "N/A" when it is NaN.
    if isnan(v)
        s = 'N/A';
    else
        s = sprintf('%.2f',v);
    end
end

function T = growthLogToTable(growthLog)
% GROWTHLOGTOTABLE  Flatten the step-wise growth log into a plain table
% for console display and CSV export.
    n = numel(growthLog);
    step = zeros(n,1);
    action = cell(n,1);
    inclination = zeros(n,1);
    raan = zeros(n,1);
    trueAnomaly = zeros(n,1);
    altitudeOffsetKm = zeros(n,1);
    newSatRevisitMin = zeros(n,1);
    newSatAreaKm2 = zeros(n,1);
    coverageFraction = zeros(n,1);
    avgRevisitMin = zeros(n,1);
    passesPerDay = zeros(n,1);
    sgsTotalMin = zeros(n,1);
    avgLatencyMin = zeros(n,1);
    wasserstein = zeros(n,1);
    largestGapHr = zeros(n,1);
    gapReductionHr = zeros(n,1);
    decision = cell(n,1);

    for i = 1:n
        r = growthLog(i);
        step(i) = r.step;
        action{i} = r.action;
        inclination(i) = r.inclination;
        raan(i) = r.raan;
        trueAnomaly(i) = r.trueAnomalyDeg;
        altitudeOffsetKm(i) = r.altitudeOffsetKm;
        newSatRevisitMin(i) = r.newSatMetrics.avgRevisitMin;
        newSatAreaKm2(i) = r.newSatMetrics.accumAreaKm2;
        coverageFraction(i) = r.combinedMetrics.coverageFraction;
        avgRevisitMin(i) = r.combinedMetrics.avgRevisitMin;
        passesPerDay(i) = r.combinedMetrics.passesPerDay;
        sgsTotalMin(i) = r.combinedMetrics.sgsTotalMin;
        avgLatencyMin(i) = r.combinedMetrics.avgLatencyMin;
        wasserstein(i) = r.wasserstein;
        largestGapHr(i) = r.newGapHr;
        gapReductionHr(i) = r.gapReductionHr;
        decision{i} = r.decision;
    end

    T = table(step,action,inclination,raan,trueAnomaly,altitudeOffsetKm, ...
        newSatRevisitMin,newSatAreaKm2,coverageFraction,avgRevisitMin, ...
        passesPerDay,sgsTotalMin,avgLatencyMin,wasserstein,largestGapHr, ...
        gapReductionHr,decision);
end
function [gapHr,gapStartHr,gapEndHr] = largestScheduleGapHours(intervals,periodSec)
% LARGESTSCHEDULEGAPHOURS
% The single largest empty stretch, in hours, between the end of one pass
% and the start of the next in a day's merged Nepal pass schedule
% (INTERVALS: an Nx2 array of non-overlapping [start,end] second offsets,
% as returned by MERGEINTERVALS / stored in .observationIntervals). The
% day is treated as a 24-hour CIRCULAR clock, so the wraparound stretch
% from the last pass of the day back to the first pass of the next day is
% included — this is the piece the existing AVGREVISITMIN/MAXREVISITMIN
% metrics miss, since they only look at gaps between consecutive passes
% within the linear list and never check the wraparound gap, which is
% exactly why a blank period spanning midnight (or otherwise touching the
% edge of the 24-hour window) would not previously show up in any printed
% metric even though it is visible in the pass-time histogram.
% Returns the gap length in hours and the local-time-of-day bounds
% [gapStartHr,gapEndHr) of that gap; gapEndHr may exceed 24 to flag a
% wraparound gap (subtract 24 for its portion after midnight).
    periodHr = periodSec/3600;
    if isempty(intervals)
        gapHr = periodHr;      % no passes at all: the entire day is empty
        gapStartHr = 0;
        gapEndHr = periodHr;
        return;
    end
    iv = sortrows(intervals,1);
    starts = iv(:,1);
    ends = iv(:,2);
    n = size(iv,1);

    gapStartsSec = ends;
    gapEndsSec = [starts(2:n);starts(1)+periodSec];
    gapsSec = gapEndsSec-gapStartsSec;

    [gapSec,ix] = max(gapsSec);
    gapHr = max(gapSec,0)/3600;
    gapStartHr = gapStartsSec(ix)/3600;
    gapEndHr = gapEndsSec(ix)/3600;
end

function formatHistAxes(ax,titleStr)
% FORMATHISTAXES
% Apply the same dark, white-on-black styling used by every other figure
% in this script (Section H's static map, the final-constellation static
% map, and both animated ground-track builds) to a pass-time histogram
% axes, so this figure reads as part of the same report instead of a
% plain default-white MATLAB figure dropped in among dark ones.
    set(ax,'Color',[0 0 0]);
    hold(ax,'on');
    grid(ax,'on');
    ax.GridColor = [0.3 0.3 0.3];
    ax.XColor = 'w';
    ax.YColor = 'w';
    box(ax,'on');
    xlabel(ax,'Local time of day (hr)','Color','w');
    ylabel(ax,'Number of passes','Color','w');
    title(ax,titleStr,'Color','w');
    xlim(ax,[0 24]);
    xticks(ax,0:2:24);
end

function shadeLargestGap(ax,gapStartHr,gapEndHr,gapHr)
% SHADELARGESTGAP
% Shade the largest empty pass-time-of-day interval on a histogram axes
% (splitting the band in two if it wraps past midnight) and label it with
% its duration — a direct visualization of the same largest-gap metric
% now tracked as a genuine selection criterion through the constellation
% growth loop (see LARGESTSCHEDULEGAPHOURS, DECIDEGROWTHSTEP,
% COMPUTEGROWTHOBJECTIVE), rather than only inferred by eye from the bars.
    yl = ylim(ax);
    if yl(2) <= 0
        yl(2) = 1;
    end
    bandColor = [0.85 0.2 0.2];
    if gapEndHr <= 24
        patch(ax,[gapStartHr,gapEndHr,gapEndHr,gapStartHr],[yl(1),yl(1),yl(2),yl(2)], ...
            bandColor,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
    else
        patch(ax,[gapStartHr,24,24,gapStartHr],[yl(1),yl(1),yl(2),yl(2)], ...
            bandColor,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
        patch(ax,[0,gapEndHr-24,gapEndHr-24,0],[yl(1),yl(1),yl(2),yl(2)], ...
            bandColor,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
    end
    text(ax,0.5,yl(2)-0.05*(yl(2)-yl(1)),sprintf('Largest gap: %.2f hr',gapHr), ...
        'Color',[1 0.75 0.75],'FontSize',8,'FontWeight','bold', ...
        'VerticalAlignment','top','HandleVisibility','off');
    ylim(ax,yl);
end

function playScenarioSafely(sc,sensors,viewerLabel)
% PLAYSCENARIOSAFELY
% Start SATELLITESCENARIO animation playback with a short retry-with-
% backoff loop instead of a single unguarded PLAY call. The 3-D viewer's
% (browser-based) globe can time out waiting to finish initializing before
% PLAY is invoked — more likely the more satellites/sensors/ground
% stations there are to place — and when that happens the viewer can be
% left showing only the bare globe with every satellite, orbit, sensor
% footprint, and ground-station marker gone. Re-issuing FIELDOFVIEW for
% every sensor (which redraws that satellite's orbit/footprint) before
% each retry is cheap and gives the viewer another chance to finish
% initializing and recover the full scene; the retries use a growing
% pause so a slow-to-start viewer gets progressively more time. This
% wrapper never lets a viewer timeout propagate as an uncaught error and
% stop the rest of the script, and it never leaves the failure silent —
% MATLAB's dedicated satellite scenario viewer only supports one open
% scenario, so a failed PLAY does not affect any later processing.
    maxAttempts = 1;
    pauseSecondsByAttempt = [3,6,10];
    played = false;
    for attempt = 1:maxAttempts
        pause(pauseSecondsByAttempt(min(attempt,numel(pauseSecondsByAttempt))));
        for si = 1:numel(sensors)
            try
                fieldOfView(sensors{si});
            catch
                % A sensor may already have an active field-of-view display;
                % that is not a failure, so continue with the remaining ones.
            end
        end
        try
            play(sc);
            played = true;
            break;
        catch ME
            warning('%s: 3-D viewer playback attempt %d/%d timed out or failed (%s).', ...
                viewerLabel,attempt,maxAttempts,ME.message);
        end
    end
    if ~played
        warning(['%s: 3-D viewer playback could not be started automatically after %d ' ...
            'attempts. The viewer window is still open with every satellite, orbit, sensor ' ...
            'footprint, and ground station already placed (FIELDOFVIEW was re-issued on each ' ...
            'retry), and can be replayed manually from its own toolbar. This does not affect ' ...
            'any metric, map, or figure computed below, all of which are derived directly from ' ...
            'the propagated satellite states rather than from the 3-D viewer.'], ...
            viewerLabel,maxAttempts);
    end
end
