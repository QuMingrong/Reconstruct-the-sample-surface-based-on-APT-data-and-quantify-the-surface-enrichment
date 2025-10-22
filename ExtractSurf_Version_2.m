%% ExtractSurf_HEA_Layers_App.m
% HEA-only surface + per-cell ranked peeling (any layer on demand) + exports
% + persistent interactive app for N-layer surface & 1 nm enrichment with curves
% Author: (your name)   Date: (today)

clear; clc;

%% 1) Select files
[aptFile, aptPath] = uigetfile({'*.csv'}, 'Select APT point cloud file (apt.csv)');
if isequal(aptFile,0), error('Cancelled.'); end
[mzwFile, mzwPath] = uigetfile({'*.xlsx;*.csv'}, 'Select m/z window table (APT_mz_windows.xlsx)');
if isequal(mzwFile,0), error('Cancelled.'); end
outDefault = fullfile(aptPath, 'HEA_surface_output.xlsx');
[outFile, outPath] = uiputfile({'*.xlsx'}, 'Save surface Excel as...', outDefault);
if isequal(outFile,0), error('Cancelled.'); end
outXLSX = fullfile(outPath, outFile);

%% 2) Parameters (surface building)
prompt = {'XY grid step h (nm):', ...
          'Quantile q (0–1, 1=strict max):', ...
          'Min HEA points / XY cell (Nmin):', ...
          'Layer-1 definition (topmost/closest):'};
defAns = {'0.7','0.98','20','topmost'};
answ = inputdlg(prompt,'HEA-only Surface Parameters',1,defAns);
if isempty(answ), error('Cancelled.'); end
h    = str2double(answ{1});
q    = str2double(answ{2});
Nmin = str2double(answ{3});
layer_order = lower(strtrim(answ{4}));
if ~ismember(layer_order,["topmost","closest"])
    layer_order = "topmost";
end

% global HEA list
HEAset = ["Ir","Ru","Rh","Pt","Pd"];

%% 3) Load data + map m/z -> species
T = readtable(fullfile(aptPath, aptFile), 'VariableNamingRule','preserve', 'TextType','string');
cn = matlab.lang.makeUniqueStrings(strtrim(string(T.Properties.VariableNames)));
T.Properties.VariableNames = cn;
xCol = find(ismember(upper(cn), ["X","CX","COORDX"]), 1);
yCol = find(ismember(upper(cn), ["Y","CY","COORDY"]), 1);
zCol = find(ismember(upper(cn), ["Z","CZ","COORDZ"]), 1);
mzCol = find(contains(lower(cn),["mz","m/z","mass","tocharge"]), 1);
assert(~(isempty(xCol)||isempty(yCol)||isempty(zCol)||isempty(mzCol)), 'X/Y/Z/mz not detected.');
T.X = double(T{:,xCol});  T.Y = double(T{:,yCol});  T.Z = double(T{:,zCol});  T.mz = double(T{:,mzCol});

[~,~,extW] = fileparts(mzwFile);
if strcmpi(extW,'.csv')
    W = readtable(fullfile(mzwPath, mzwFile), 'VariableNamingRule','preserve', 'TextType','string');
else
    W = readtable(fullfile(mzwPath, mzwFile),'FileType','spreadsheet','VariableNamingRule','preserve','TextType','string');
end
W.Properties.VariableNames = matlab.lang.makeUniqueStrings(strtrim(string(W.Properties.VariableNames)));
assert(all(ismember(["Species","mz_min","mz_max"], string(W.Properties.VariableNames))), ...
       'Window table must have: Species, mz_min, mz_max');

label = assignSpeciesFromWindows(T.mz, W);
T.Species = categorical(label);

%% 4) Build HEA-only upper envelope (surface) + normals
Th = T(ismember(string(T.Species), HEAset), :);
assert(~isempty(Th), 'No HEA (Ir/Ru/Rh/Pt/Pd) points after filtering.');
XYZ = [Th.X, Th.Y, Th.Z];

xb = floor(min(XYZ(:,1))/h)*h : h : ceil(max(XYZ(:,1))/h)*h;
yb = floor(min(XYZ(:,2))/h)*h : h : ceil(max(XYZ(:,2))/h)*h;
nx = numel(xb)-1; ny = numel(yb)-1;

zq  = nan(nx,ny);  nxG = nan(nx,ny);  nyG = nan(nx,ny);  nzG = nan(nx,ny);
for ix = 1:nx
  for iy = 1:ny
    in = XYZ(:,1)>=xb(ix) & XYZ(:,1)<xb(ix+1) & XYZ(:,2)>=yb(iy) & XYZ(:,2)<yb(iy+1);
    if nnz(in) < Nmin, continue; end
    z_local   = XYZ(in,3);
    zq(ix,iy) = quantile(z_local, q);       % upper envelope per cell
    Xi = [XYZ(in,1), XYZ(in,2), ones(nnz(in),1)];
    ai = Xi \ z_local;                      % plane fit z ~ ax+by+c
    n  = [-ai(1), -ai(2), 1]; n = n/norm(n);
    nxG(ix,iy)=n(1); nyG(ix,iy)=n(2); nzG(ix,iy)=n(3);
  end
end
[Xg, Yg] = ndgrid(xb(1:end-1)+h/2, yb(1:end-1)+h/2);
Zg = zq;

%% 5) Signed-normal distance for ALL atoms; per-cell HEA ranking (NO fixed K)
FnZ = griddedInterpolant({xb(1:end-1)+h/2, yb(1:end-1)+h/2}, Zg,  'nearest','nearest');
Fnx= griddedInterpolant({xb(1:end-1)+h/2, yb(1:end-1)+h/2}, nxG, 'nearest','nearest');
Fny= griddedInterpolant({xb(1:end-1)+h/2, yb(1:end-1)+h/2}, nyG, 'nearest','nearest');
Fnz= griddedInterpolant({xb(1:end-1)+h/2, yb(1:end-1)+h/2}, nzG, 'nearest','nearest');

x=T.X; y=T.Y; z=T.Z;
zsurf=FnZ(x,y); nxv=Fnx(x,y); nyv=Fny(x,y); nzv=Fnz(x,y);
valid = isfinite(zsurf)&isfinite(nxv)&isfinite(nyv)&isfinite(nzv);

nvec = [nxv(valid), nyv(valid), nzv(valid)];
nvec = nvec ./ vecnorm(nvec,2,2);
sref = [x(valid), y(valid), zsurf(valid)];
d    = sum( ([x(valid),y(valid),z(valid)] - sref) .* nvec, 2 );

% per-atom grid index
nxg = size(Zg,1); nyg = size(Zg,2);
ix_all = max(1, min(nxg, floor((x - xb(1))/h)+1));
iy_all = max(1, min(nyg, floor((y - yb(1))/h)+1));
Ivalid = find(valid);

% ---- build HEA per-cell order list (rank by d) ----
[cellOrderByD, maxRankAvailable] = buildCellOrder(T, Ivalid, ix_all, iy_all, d, HEAset, layer_order);

% convenience masks for first 5 layers (so你不改其它段也能用 L{1..5})
K_pre = min(5, maxRankAvailable);
L = repmat({false(height(T),1)}, 1, K_pre);
for k = 1:K_pre
    L{k} = mask_from_rank(cellOrderByD, k, height(T));
end
L_top5 = false(height(T),1);
for k=1:K_pre, L_top5 = L_top5 | L{k}; end

% expose to base for the app
assignin('base','T',T);
assignin('base','outXLSX',outXLSX);
assignin('base','HEAset',HEAset);
assignin('base','layer_order',string(layer_order));
assignin('base','cellOrderByD',cellOrderByD);
assignin('base','maxRankAvailable',maxRankAvailable);
assignin('base','L',L);                  % first 1..min(5,available)
assignin('base','L_top5',L_top5);
assignin('base','Xg',Xg); assignin('base','Yg',Yg); assignin('base','Zg',Zg);

%% 6) Export default layers (1..5 or up to available) + Z/D stats
colHex = struct('Ir','#FF0000','Pt','#0000FF','Ru','#FF00FF','Pd','#00CCFF','Rh','#000000');

% SurfaceGrid
maskGrid = isfinite(Zg);
writetable(table(Xg(maskGrid),Yg(maskGrid),Zg(maskGrid), 'VariableNames',{'X','Y','Z'}), ...
           outXLSX,'Sheet','SurfaceGrid');

% per-layer exports (HEA only)
Kexp = K_pre;    % default export first up to 5 layers
Layer = (1:Kexp)'; Count=zeros(Kexp,1); MeanZ=nan(Kexp,1); MedianZ=nan(Kexp,1);
StdZ=nan(Kexp,1); MinZ=nan(Kexp,1); MaxZ=nan(Kexp,1);
d_full = nan(height(T),1); d_full(Ivalid)=d;
MeanD=nan(Kexp,1); MedianD=nan(Kexp,1); StdD=nan(Kexp,1);

for k=1:Kexp
    Mk = mask_from_rank(cellOrderByD, k, height(T));
    Tk = T(Mk & ismember(string(T.Species),HEAset), {'X','Y','Z','Species'});
    Count(k)=height(Tk);
    if Count(k)>0
        MeanZ(k)=mean(Tk.Z,'omitnan'); MedianZ(k)=median(Tk.Z,'omitnan');
        StdZ(k)=std(Tk.Z,0,'omitnan');  MinZ(k)=min(Tk.Z);  MaxZ(k)=max(Tk.Z);
    end
    dk = d_full(Mk & ismember(string(T.Species),HEAset));
    MeanD(k)=mean(dk,'omitnan'); MedianD(k)=median(dk,'omitnan'); StdD(k)=std(dk,0,'omitnan');

    % color hex column
    colorHexCol = strings(height(Tk),1);
    for i=1:height(Tk)
        s = char(string(Tk.Species(i)));
        if isfield(colHex,s), colorHexCol(i)=colHex.(s); else, colorHexCol(i)="#808080"; end
    end
    Tk.ColorHex = colorHexCol;

    writetable(Tk, outXLSX, 'Sheet', sprintf('SurfaceLayer%d',k));

    % species summary
    if height(Tk)>0
        S = groupsummary(Tk,'Species');
        S.Properties.VariableNames = {'Species','GroupCount'};
        S.Fraction = S.GroupCount / sum(S.GroupCount);
    else
        S = table(categorical([]),[],[], 'VariableNames',{'Species','GroupCount','Fraction'});
    end
    writetable(S, outXLSX, 'Sheet', sprintf('Layer%d_Summary',k));
end
writetable(table(Layer,Count,MeanZ,MedianZ,StdZ,MinZ,MaxZ), outXLSX, 'Sheet','LayerZ_Stats');
writetable(table(Layer,MeanD,MedianD,StdD), outXLSX, 'Sheet','LayerD_Stats');

%% 7) Plots — NO green mesh overlays in layer figures
% A) surface mesh only
figure('Color','w');
surf(Xg,Yg,Zg,'EdgeColor','none','FaceAlpha',0.95);
axis equal tight; view(3); camlight headlight; lighting gouraud;
xlabel('X (nm)'); ylabel('Y (nm)'); zlabel('Z (nm)');
title('HEA-only Surface (mesh only)'); colorbar;

% B) each of first Kexp layers — HEA colored scatter (no surface overlay)
for k=1:Kexp
    Mk = mask_from_rank(cellOrderByD, k, height(T));
    Tk = T(Mk & ismember(string(T.Species),HEAset), :);

    figure('Color','w'); hold on;
    legendNames = {};
    for s = HEAset
        in = Tk.Species == categorical(s);
        if any(in)
            c = pickcolor(colHex, s);
            scatter3(Tk.X(in),Tk.Y(in),Tk.Z(in),10,'filled', ...
                     'MarkerFaceColor',c,'MarkerEdgeColor','none');
            legendNames{end+1} = char(s); %#ok<AGROW>
        end
    end
    axis equal tight; grid on; view(3);
    xlabel('X (nm)'); ylabel('Y (nm)'); zlabel('Z (nm)');
    title(sprintf('Surface Layer %d (HEA only)', k));
    if ~isempty(legendNames), legend(legendNames,'Location','northeastoutside'); end
end

% C) (example) combined top N layers scatter (change Ntop if needed)
Ntop = min(2, maxRankAvailable);
Mtop = false(height(T),1);
for k=1:Ntop, Mtop = Mtop | mask_from_rank(cellOrderByD,k,height(T)); end
Ts = T(Mtop & ismember(string(T.Species),HEAset),:);
figure('Color','w'); hold on;
legendNames = {};
for s = HEAset
    in = Ts.Species == categorical(s);
    if any(in)
        c = pickcolor(colHex, s);
        scatter3(Ts.X(in),Ts.Y(in),Ts.Z(in),8,'filled', ...
                 'MarkerFaceColor',c,'MarkerEdgeColor','none');
        legendNames{end+1}=char(s);
    end
end
axis equal tight; grid on; view(3);
xlabel('X (nm)'); ylabel('Y (nm)'); zlabel('Z (nm)');
title(sprintf('Surface Layers 1-%d (HEA only, element-colored)', Ntop));
if ~isempty(legendNames), legend(legendNames,'Location','northeastoutside'); end

%% 8) Open interactive app (on-demand layers & 1 nm enrichment)
createSurfaceApp();

%% ================== Helper functions ==================
function label = assignSpeciesFromWindows(mz, W)
    sp = string(W.Species);
    mm = [double(W.mz_min), double(W.mz_max)];
    mz = mz(:);
    label = strings(size(mz));
    speciesList = unique(sp,'stable');
    for k = 1:numel(speciesList)
        s = speciesList(k); rows = find(sp==s);
        mask = false(size(mz));
        for r = 1:numel(rows)
            lo = mm(rows(r),1); hi = mm(rows(r),2);
            mask = mask | (mz>=lo & mz<=hi);
        end
        assignable = (label=="") & mask;
        label(assignable) = s;
    end
    label(label=="") = "Unknown";
end

function c = pickcolor(colHex, s)
    s = char(string(s));
    if isfield(colHex,s), c = hex2rgb(colHex.(s)); else, c=[.5 .5 .5]; end
end
function rgb = hex2rgb(hex)
    hex = char(hex); if hex(1) == '#', hex = hex(2:end); end
    rgb = [hex2dec(hex(1:2)) hex2dec(hex(3:4)) hex2dec(hex(5:6))]/255;
end

function [cellOrderByD, maxRankAvailable] = buildCellOrder(T, Ivalid, ix_all, iy_all, d, HEAset, layer_order)
    % 只在 HEA 中排名；为避免把第二层“误判掉”，不再强行 d>=0
    % 如需更保守，可把 keep = isHEA_valid 改成 keep = (d >= d_min) & isHEA_valid
    % 建议 d_min = -0.15 ~ -0.30（单位 nm）作为容差
    d_min = -Inf;   % ← 想开容差就改成 -0.2 之类
    isHEA_all   = ismember(string(T.Species), HEAset);
    isHEA_valid = isHEA_all(Ivalid);

    keep = isHEA_valid & (d >= d_min);   % ★ 不再用 d>=0 硬门槛

    % 有效原子所在格
    nx = max(ix_all); ny = max(iy_all);
    cellID_valid = sub2ind([nx,ny], ix_all(Ivalid), iy_all(Ivalid));

    % 只遍历“含 HEA 候选”的格
    uCells = unique(cellID_valid(keep));
    G_idx = []; G_d = []; G_cell = []; G_rank = [];

    for c = reshape(uCells,1,[])
        J = find(cellID_valid == c & keep);
        if isempty(J), continue; end

        % 排序方向：'topmost' = d 降序；'closest' = d 升序
        switch lower(layer_order)
            case 'closest', [~,ord] = sort(d(J), 'ascend');
            otherwise,      [~,ord] = sort(d(J), 'descend');  % 默认 topmost
        end
        J = J(ord);

        ranks = (1:numel(J))';
        G_idx  = [G_idx;  Ivalid(J)];            %#ok<AGROW>
        G_d    = [G_d;    d(J)];                 %#ok<AGROW>
        G_cell = [G_cell; repmat(c,numel(J),1)]; %#ok<AGROW>
        G_rank = [G_rank; ranks];                %#ok<AGROW>
    end

    cellOrderByD.idxGlobal = G_idx;   % 索引到 T 的全局行
    cellOrderByD.d         = G_d;
    cellOrderByD.cellID    = G_cell;
    cellOrderByD.rank      = G_rank;

    maxRankAvailable = 0;
    if ~isempty(G_rank), maxRankAvailable = max(G_rank); end
end


function M = mask_from_rank(cellOrderByD, k, N)
    M = false(N,1);
    if isempty(fieldnames(cellOrderByD)), return; end
    rows = (cellOrderByD.rank == k);
    if any(rows)
        M(cellOrderByD.idxGlobal(rows)) = true;
    end
end

function [r, b] = fcc_shell_params(d1, Smax, epsR)
    Smax = min(Smax,5);
    a = sqrt(2)*d1;
    rFCC = [a/sqrt(2), a, a*sqrt(3/2), a*sqrt(2), a*sqrt(5/2)];
    r = rFCC(1:Smax).';
    nextR = a*sqrt(3);
    b = zeros(Smax,1);
    if Smax>=2, b(1) = (r(1)+r(2))/2; else, b(1) = r(1)+0.5*(nextR-r(1)); end
    for s=2:Smax-1, b(s) = (r(s)+r(s+1))/2; end
    b(end) = (r(end)+nextR)/2;
    b = b + epsR;
end

function [tbl, pairCounts] = kde_pdf_pairs(T, centerMask, neighborsMask, neighList, r_min, r_max, bw, nGrid)
    neighList = string(neighList(:)');
    centers = find(centerMask);
    tbl = table(string.empty, [], [], [], ...
        'VariableNames', {'NeighborSpecies','r','pdf_kde','Npairs'});
    pairCounts = zeros(1, numel(neighList));
    if isempty(centers), return; end

    XYZ=[T.X,T.Y,T.Z]; nbrIdx=find(neighborsMask);
    if isempty(nbrIdx), return; end
    Mdl=createns(XYZ(nbrIdx,:), 'NSMethod','kdtree');

    D = cell(numel(neighList),1);
    idxCell = rangesearch(Mdl, XYZ(centers,:), r_max);
    for ii=1:numel(centers)
        c=centers(ii); glb=nbrIdx(idxCell{ii}); glb(glb==c)=[];
        if isempty(glb), continue; end
        dist=vecnorm(XYZ(glb,:)-XYZ(c,:),2,2);
        sp=string(T.Species(glb));
        sel=(dist>=r_min & dist<=r_max);
        dist=dist(sel); sp=sp(sel);
        for k=1:numel(neighList)
            if ~isempty(dist), D{k} = [D{k}; dist(sp==neighList(k))]; end %#ok<AGROW>
        end
    end

    r_grid = linspace(r_min, r_max, nGrid).';
    for k=1:numel(neighList)
        d = D{k}; pairCounts(k)=numel(d);
        if isempty(d)
            f=zeros(nGrid,1); Npairs=0;
        else
            try
                [f,~]=ksdensity(d, r_grid, 'Bandwidth', bw, ...
                                'Support', [r_min r_max], ...
                                'BoundaryCorrection','reflection');
            catch
                d_ext=[2*r_min-d; d; 2*r_max-d];
                [f,~]=ksdensity(d_ext, r_grid, 'Bandwidth', bw, ...
                                'Support', [r_min r_max]);
            end
            A=trapz(r_grid,f); if A>0, f=f/A; end
            Npairs=numel(d);
        end
        tmp=table(repmat(neighList(k),numel(r_grid),1), r_grid, f(:), ...
                  repmat(Npairs,numel(r_grid),1), ...
                  'VariableNames',{'NeighborSpecies','r','pdf_kde','Npairs'});
        tbl=[tbl; tmp]; %#ok<AGROW>
    end
end

function createSurfaceApp()
% Persistent interactive app: choose "surface = top N layers", center element,
% neighbor mode & baseline, radius R; compute within-R enrichment & 3 curves.
    persistent ui;
    if ~isempty(ui) && isvalid(ui.f), figure(ui.f); return; end

    % pull from base
    T            = evalin('base','T');
    outXLSX      = evalin('base','outXLSX');
    HEAset       = evalin('base','HEAset');
    layer_order  = evalin('base','layer_order');
    cellOrderByD = evalin('base','cellOrderByD');
    maxRankAvail = evalin('base','maxRankAvailable');

    % ui
    ui.f = uifigure('Name','Surface Layers & 1 nm Enrichment','Position',[100 100 980 620]);

    uilabel(ui.f,'Text','Surface = top N layers','Position',[20 580 180 22]);
    ui.Nsurf = uieditfield(ui.f,'numeric','Limits',[1 max(1,maxRankAvail)], ...
        'RoundFractionalValues','on','Value',min(3,maxRankAvail),'Position',[200 576 80 28]);

    uilabel(ui.f,'Text','Center element','Position',[20 540 120 22]);
    avail = intersect(HEAset, string(categories(removecats(T.Species))),'stable');
    if isempty(avail), avail = string(categories(removecats(T.Species))); end
    ui.center = uidropdown(ui.f,'Items',cellstr(avail),'Value',char(avail(1)),'Position',[140 536 140 28]);

    uilabel(ui.f,'Text','Neighbor mode','Position',[320 540 120 22]);
    ui.mode = uidropdown(ui.f,'Items',{'HEA (exclude center)','All species (exclude center)'}, ...
        'Value','HEA (exclude center)','Position',[440 536 200 28]);

    uilabel(ui.f,'Text','Baseline','Position',[660 540 80 22]);
    ui.base = uidropdown(ui.f,'Items',{'Global baseline','Surface baseline'}, ...
        'Value','Global baseline','Position',[740 536 200 28]);

    uilabel(ui.f,'Text','Radius r (nm)','Position',[320 580 100 22]);
    ui.R = uieditfield(ui.f,'numeric','Limits',[0.1 5],'Value',1.0,'Position',[420 576 80 28]);

    ui.btn = uibutton(ui.f,'Text','Compute 1 nm enrichment','Position',[740 576 200 32], ...
        'ButtonPushedFcn',@onCompute);

    % result table
    uilabel(ui.f,'Text','Enrichment within r (nm)','FontWeight','bold','Position',[20 510 220 20]);
    ui.tabRes = uitable(ui.f,'Position',[20 280 920 220], ...
        'ColumnName',{'Neighbor','PairsWithinR','P_withinR','pi_baseline','Enrichment E','Tendency'}, ...
        'ColumnEditable',false(1,6));

    % per-layer Z stats (first up to 10 for quick view)
    uilabel(ui.f,'Text','Per-layer Z stats (HEA only)','FontWeight','bold','Position',[20 250 220 20]);
    ui.tabZ = uitable(ui.f,'Position',[20 20 450 220], ...
        'ColumnName',{'Layer','Count','MeanZ','MedianZ','StdZ','MinZ','MaxZ'}, ...
        'ColumnEditable',false(1,7));

    ui.info = uitextarea(ui.f,'Position',[500 20 440 220],'Editable','off');

    % prefill Z stats for first min(5,maxRankAvail) layers
    Kshow = min(10, maxRankAvail);
    [Layer, Count, MeanZ, MedianZ, StdZ, MinZ, MaxZ] = layerStats(T, cellOrderByD, HEAset, Kshow);
    ui.tabZ.Data = table(Layer, Count, MeanZ, MedianZ, StdZ, MinZ, MaxZ);

    function onCompute(~,~)
        % read UI
        Ns = max(1, min(maxRankAvail, round(ui.Nsurf.Value)));
        centerA = string(ui.center.Value);
        R = ui.R.Value;
        modeStr = ui.mode.Value;
        baseStr = ui.base.Value;

        % surface mask = union layers 1..Ns
        surfaceMask = false(height(T),1);
        for k=1:Ns, surfaceMask = surfaceMask | mask_from_rank(cellOrderByD,k,height(T)); end

        % center & neighbors
        centerMask = surfaceMask & (string(T.Species) == centerA);

        if contains(modeStr,'HEA')
            universe  = HEAset;
            neighList = setdiff(universe, centerA, 'stable');
            neighborsMask = ismember(string(T.Species), neighList);   % [PATCH-3a: REPLACE LINES]

        else
            allSp = string(categories(removecats(T.Species)));
            neighList = setdiff(allSp, centerA, 'stable');
            neighborsMask = ismember(string(T.Species), neighList);   % [PATCH-3b: REPLACE LINES]

        end

        % table within R
        [ResTbl, metaTxt] = enrichment_within_radius(T, centerMask, neighborsMask, neighList, R, surfaceMask, baseStr);
        ui.tabRes.Data = ResTbl; ui.info.Value = metaTxt;
        try, writetable(ResTbl, outXLSX, 'Sheet', sprintf('SurfN%d_Enrich_%s', Ns, char(centerA))); catch, end

        if ~any(centerMask), return; end

        % ---- curves (r from 0) ----
        r_min = 0.22; r_max = 1.00; nGrid = 400; bw = 0.03;
        [tblKDE, ~] = kde_pdf_pairs(T, centerMask, neighborsMask, neighList, r_min, r_max, bw, nGrid);

        r  = unique(tblKDE.r);  nR = numel(r);  M = numel(neighList);
        F  = zeros(nR, M);
        for k = 1:M
            rows = (tblKDE.NeighborSpecies == neighList(k));
            F(:,k) = tblKDE.pdf_kde(rows);
        end

       % baseline (global vs surface)   % [PATCH-1: REPLACE BLOCK]
       if strcmp(baseStr,'Surface baseline')
           baseMask = neighborsMask & surfaceMask;
       else
           baseMask = neighborsMask;
       end
       spBase = string(T.Species(baseMask));

       % 只保留与曲线一致的物种集合（例如 HEA 模式下去掉中心 Ir）
       spBase = spBase(ismember(spBase, neighList));

       % 计数并在 neighList 上归一化 -> sum(pi_s) = 1
       counts = arrayfun(@(s) sum(spBase == s), neighList);
       if sum(counts)==0
           % 安全兜底：若基线集合为空，则用均匀先验
           pi_s = ones(numel(neighList),1) / numel(neighList);
       else
           pi_s = counts(:) / sum(counts);
       end

        g = F * pi_s; g(g<=0)=eps;
        P = (F .* pi_s.') ./ g;     % posterior
        E = P ./ pi_s.';            % enrichment

        % r->0 visual extension
        r_plot = [0; r]; F_plot = [zeros(1,M); F]; P_plot = [zeros(1,M); P]; E_plot = [zeros(1,M); E];

        % shell reference lines
        [rFCC, ~] = fcc_shell_params(0.25,5,0.04);

        % Plot 1: PDF
        f1 = figure('Color','w'); hold on;
        for i=1:numel(rFCC), xline(rFCC(i),'--','Color',[.7 .7 .7]); end
        xline(r_min,':','Color',[.6 .6 .6],'HandleVisibility','off');
        for k=1:M, plot(r_plot, F_plot(:,k), 'LineWidth', 1.8, 'DisplayName', char(neighList(k))); end
        grid on; box on; xlabel('r (nm)'); ylabel('Probability density (1/nm)');
        title(sprintf('Pair-distance PDF (KDE) — Center: %s', centerA));
        legend('Location','northeastoutside'); xlim([0 r_max]);

        % Plot 2: Posterior
        f2 = figure('Color','w'); hold on;
        for i=1:numel(rFCC), xline(rFCC(i),'--','Color',[.7 .7 .7]); end
        xline(r_min,':','Color',[.6 .6 .6],'HandleVisibility','off');
        for k=1:M, plot(r_plot, P_plot(:,k), 'LineWidth', 1.8, 'DisplayName', sprintf('%s (P)', char(neighList(k)))); end
        grid on; box on; xlabel('r (nm)'); ylabel('P(neighbor = s | r)'); ylim([0 1]);
        title(sprintf('Posterior around %s (%s baseline)', centerA, lower(baseStr)));
        legend('Location','northeastoutside'); xlim([0 r_max]);

        % Plot 3: Enrichment
        f3 = figure('Color','w'); hold on; yline(1,'k:');
        for i=1:numel(rFCC), xline(rFCC(i),'--','Color',[.7 .7 .7]); end
        xline(r_min,':','Color',[.6 .6 .6],'HandleVisibility','off');
        for k=1:M, plot(r_plot, E_plot(:,k), 'LineWidth', 1.8, 'DisplayName', sprintf('%s (E)', char(neighList(k)))); end
        grid on; box on; xlabel('r (nm)'); ylabel(sprintf('Enrichment E_s(r)=P/\\pi_s (%s)', lower(baseStr)));
        title(sprintf('Enrichment vs %s baseline — Center: %s', lower(baseStr), centerA));
        legend('Location','northeastoutside'); xlim([0 r_max]);

        % also refresh per-layer Z stats (first up to 10) for quick view
        [Layer, Count, MeanZ, MedianZ, StdZ, MinZ, MaxZ] = layerStats(T, cellOrderByD, HEAset, min(10,maxRankAvail));
        ui.tabZ.Data = table(Layer, Count, MeanZ, MedianZ, StdZ, MinZ, MaxZ);

        % optional: save three figures next to Excel
        %#ok<*NASGU>
        % print(f1, fullfile(fileparts(outXLSX), 'PDF.png'), '-dpng','-r300');  % if you want
        % print(f2, fullfile(fileparts(outXLSX), 'POST.png'),'-dpng','-r300');
        % print(f3, fullfile(fileparts(outXLSX), 'ENRICH.png'),'-dpng','-r300');
    end
end

function [Layer, Count, MeanZ, MedianZ, StdZ, MinZ, MaxZ] = layerStats(T, cellOrderByD, HEAset, Kshow)
    N = height(T);
    Layer = (1:Kshow)'; Count=zeros(Kshow,1); MeanZ=nan(Kshow,1); MedianZ=nan(Kshow,1);
    StdZ=nan(Kshow,1); MinZ=nan(Kshow,1); MaxZ=nan(Kshow,1);
    for k=1:Kshow
        Mk = mask_from_rank(cellOrderByD,k,N);
        Tk = T(Mk & ismember(string(T.Species),HEAset),:);
        Count(k)=height(Tk);
        if Count(k)>0
            MeanZ(k)=mean(Tk.Z,'omitnan'); MedianZ(k)=median(Tk.Z,'omitnan');
            StdZ(k)=std(Tk.Z,0,'omitnan');  MinZ(k)=min(Tk.Z);  MaxZ(k)=max(Tk.Z);
        end
    end
end

function [ResTbl, metaTxt] = enrichment_within_radius(T, centerMask, neighborsMask, neighList, R, surfaceMask, baseStr)
    neighList = string(neighList(:)');
    centers = find(centerMask);
    ResTbl = table();
    if isempty(centers)
        metaTxt = "No centers."; return;
    end

    XYZ = [T.X, T.Y, T.Z];
    nbrIdx = find(neighborsMask);
    if isempty(nbrIdx)
        metaTxt = "No neighbors under current mode."; return;
    end

    Mdl = createns(XYZ(nbrIdx,:), 'NSMethod','kdtree');
    K = numel(neighList);
    counts = zeros(1,K);
    idxCell = rangesearch(Mdl, XYZ(centers,:), R);
    for ii=1:numel(centers)
        c = centers(ii);
        glb = nbrIdx(idxCell{ii});
        glb(glb==c) = [];
        if isempty(glb), continue; end
        sp = string(T.Species(glb));
        for k=1:K
            counts(k) = counts(k) + sum(sp == neighList(k));
        end
    end
    totInR = sum(counts);
    if totInR==0
        metaTxt = sprintf('No neighbors found within %.3f nm.', R);
        ResTbl = table(neighList', zeros(K,1), zeros(K,1), zeros(K,1), zeros(K,1), strings(K,1), ...
            'VariableNames', {'Neighbor','PairsWithinR','P_withinR','pi_baseline','Enrichment','Tendency'});
        return;
    end
    P_within = counts / totInR;

    % baseline   % [PATCH-2: REPLACE BLOCK]
    if strcmp(baseStr,'Surface baseline')
        baseMask = neighborsMask & surfaceMask;
    else
        baseMask = neighborsMask;
    end
    spBase = string(T.Species(baseMask));
    
    % 仅统计 neighList 中的物种，并归一化
    spBase = spBase(ismember(spBase, neighList));
    counts = arrayfun(@(s) sum(spBase == s), neighList);
    if sum(counts)==0
        pi = ones(1,K)/K;           % 兜底：均匀基线
    else
        pi = counts / sum(counts);  % sum(pi)=1
    end
    
    E = P_within ./ (pi + eps);

    % tendency
    tol = 0.05; tend = strings(K,1);
    for k=1:K
        if E(k) > 1+tol, tend(k)="Aggregate";
        elseif E(k) < 1-tol, tend(k)="Repel";
        else, tend(k)="Neutral";
        end
    end

    ResTbl = table(neighList', counts', P_within', pi', E', tend, ...
        'VariableNames', {'Neighbor','PairsWithinR','P_withinR','pi_baseline','Enrichment','Tendency'});

    metaTxt = sprintf(['Centers: %d (in selected surface)\n' ...
                       'Radius: %.3f nm\n' ...
                       'Neighbors universe: %s\n' ...
                       'Baseline: %s\n' ...
                       'Total neighbor pairs within R: %d'], ...
                       numel(centers), R, strjoin(neighList,','), baseStr, totInR);
end
