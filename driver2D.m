
% P3 CutFEM for Laplace with two *unfitted circular* Dirichlet boundaries
% hot:  (0,0), R=1,   u=50
% cold: (2.5,0), R=0.5,u=0
% Outer boundary: do-nothing (Neumann).
% - TRUE-arc Nitsche on each circle
% - Volume masking (skip GPs inside circles)
% - STRONG clamping of all nodes strictly inside circles (prevents singular K)
%
% Requires:
%   - mesh5_gmsh.m
%   - P3 mesh (TRIANGLES10/LINES4), e.g. 'bigdisk.m'

clear; clc; 

%% ===== Controls =====
gmshFile   = 'circlep3r.m';
order      = 3;               % **P3**
kappa      = 1.0;
fconst     = 0.0;             % Laplace
gammaN     = 30;              % Nitsche penalty
tol        = 1e-12;

% Internal Dirichlet circles:
C(1).c = [0.0, 0.0];  C(1).R = 1.0;  C(1).g = 50.0;  % hot
C(2).c = [2.5, 0.0];  C(2).R = 0.5;  C(2).g =  0.0;  % cold

pinOneNode = false;   % not needed now that we clamp interiors

%% ===== Read mesh (P3) =====
[nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh(gmshFile, order, []); %#ok<ASGLU>
x = nodeInfo.velocity.x(:);  y = nodeInfo.velocity.y(:);
nodes = [x, y];  Nn = size(nodes,1);
tris10 = elemInfo.velElements;               % Nt×10
Nt = size(tris10,1);

%% ===== Identify nodes STRICTLY inside each circle → clamp later =====
inside1 = hypot(nodes(:,1)-C(1).c(1), nodes(:,2)-C(1).c(2)) < (C(1).R - 5e-13);
inside2 = hypot(nodes(:,1)-C(2).c(1), nodes(:,2)-C(2).c(2)) < (C(2).R - 5e-13);
clamp_ids  = [find(inside1); find(inside2)];
clamp_vals = [ C(1).g * ones(nnz(inside1),1);  C(2).g * ones(nnz(inside2),1) ];

% (circles disjoint ⇒ no overlap; if they could overlap, prefer nearer-center circle)

%% ===== P3 basis & quadrature =====
p3basis = @(xi,eta) p3IsoShape_local(xi,eta); % [N(10), dNxi(10), dNeta(10)]
[g_pt, g_wt] = triGaussPoints12_local();      % volume quad (12-pt)

%% ===== Assemble: volume (MASK inside circles) =====
K = sparse(Nn,Nn);  F = zeros(Nn,1);

for it = 1:Nt
    tri = tris10(it,:);  X = nodes(tri,:);
    Ke = zeros(10,10);   Fe = zeros(10,1);

    for kq = 1:numel(g_wt)
        xi  = g_pt(kq,1);  eta = g_pt(kq,2);
        [N, dNxi, dNeta] = p3basis(xi,eta);

        % physical Gauss point
        xq = (N.' * X); % 1x2

        % Heaviside mask: skip inside either circle
        if hypot(xq(1)-C(1).c(1), xq(2)-C(1).c(2)) <= (C(1).R - 1e-12), continue; end
        if hypot(xq(1)-C(2).c(1), xq(2)-C(2).c(2)) <= (C(2).R - 1e-12), continue; end

        J    = [dNxi.'; dNeta.'] * X;   % 2x2
        detJ = det(J);  if detJ <= 0, error('Nonpositive detJ in element %d.', it); end
        invJT = inv(J)';

        dNdx = [dNxi, dNeta] * invJT;   % 10x2

        Ke = Ke + kappa * (dNdx * dNdx.') * (detJ * g_wt(kq));
        Fe = Fe + (N * (fconst * detJ * g_wt(kq)));
    end

    K(tri,tri) = K(tri,tri) + Ke;
    F(tri)     = F(tri)     + Fe;
end

%% ===== Nitsche on both circular interfaces (TRUE arc integration) =====
for it = 1:Nt
    tri10 = tris10(it,:);  X10 = nodes(tri10,:);   % 10x2
    Xv    = X10(1:3,:);                             % vertex triangle

    for ic = 1:numel(C)
        cc = C(ic).c;  RR = C(ic).R;  gGamma = C(ic).g;

        % chord endpoints on the vertex triangle
        [ok, P, Q] = chord_in_triangle_linear(Xv, cc, RR, tol);
        if ~ok, continue; end

        % integrate along TRUE short arc P→Q
        thetaP = atan2(P(2)-cc(2), P(1)-cc(1));
        thetaQ = atan2(Q(2)-cc(2), Q(1)-cc(1));
        dth    = wrapToPi_local(thetaQ - thetaP); Ls = RR * abs(dth);
        if Ls < 1e-14, continue; end

        [t1d, w1d] = gauss_legendre_1D_local(3);
        for kq = 1:numel(w1d)
            th = thetaP + 0.5*(t1d(kq)+1)*dth;
            xp = cc + RR*[cos(th), sin(th)];
            nHat = (xp - cc) / RR;
            w    = w1d(kq) * (RR * abs(dth));      % ds

            % directional size via vertex triangle
            h_e = projected_he_vertices(Xv, nHat);

            % invert iso-P3 map to (xi,eta)
            [xi,eta,okN] = invert_isoP3_newton(xp, X10, p3basis);
            if ~okN, continue; end

            % basis + grads at xp
            [N, dNxi, dNeta] = p3basis(xi,eta);
            J    = [dNxi.'; dNeta.'] * X10;  invJT = inv(J)';
            dNdx = [dNxi, dNeta] * invJT;    % 10x2
            dn   = dNdx * nHat.';            % 10x1

            % symmetric Nitsche
            K(tri10,tri10) = K(tri10,tri10) ...
               + ( - (dn*N.' + N*dn.') + (gammaN/h_e)*(N*N.') ) * w * kappa;
            F(tri10) = F(tri10) + ( (gammaN/h_e)*gGamma*N - dn*gGamma ) * w * kappa;
        end
    end
end

%% ===== STRONG CLAMP: all nodes strictly inside the circles =====
% This removes disconnected interior DOFs (prevents singular K).
for idx = 1:numel(clamp_ids)
    i = clamp_ids(idx);  uD = clamp_vals(idx);
    F = F - K(:,i) * uD;
    K(:,i) = 0;  K(i,:) = 0;  K(i,i) = 1;  F(i) = uD;
end

%% ===== (Optional) pin one safe exterior node (not necessary now) =====
if pinOneNode
    r1 = hypot(x - C(1).c(1), y - C(1).c(2)) - C(1).R;
    r2 = hypot(x - C(2).c(1), y - C(2).c(2)) - C(2).R;
    mask = (r1 > 5e-3) & (r2 > 5e-3);  cand = find(mask);
    dmin = min(r1(mask), r2(mask)); [~,kmin] = min(dmin); iPin = cand(kmin);
    [PsiA, ~] = analytical_two_cylinders_local(nodes); uPin = PsiA(iPin);
    F = F - K(:,iPin)*uPin; K(:,iPin)=0; K(iPin,:)=0; K(iPin,iPin)=1; F(iPin)=uPin;
end

%% ===== Solve =====
u = K \ F;

%% ===== Analytical (overwrite inside with Dirichlet, clip [0,50]) =====
[PsiRaw, ~] = analytical_two_cylinders_local(nodes);
PsiShown = PsiRaw;
PsiShown(inside1) = C(1).g;
PsiShown(inside2) = C(2).g;
PsiShown = min(max(PsiShown,0),50);

% FE display clamp (visual only)
uShow = u;
uShow(inside1) = C(1).g;
uShow(inside2) = C(2).g;
uShow = min(max(uShow,0),50);

%% ===== Plots: surfaces + contours =====
figure('Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% FE surface
nexttile(1);
trisurf(tris10(:,1:3), nodes(:,1), nodes(:,2), uShow, 'EdgeColor','k'); view(2); colorbar;
axis equal tight off; title('CutFEM (P3) u (inside clamped)'); caxis([0 50]);
hold on; viscircle_local(C(1).c, C(1).R, 'Color',[0 0 0], 'LineStyle',':');
hold on; viscircle_local(C(2).c, C(2).R, 'Color',[0 0 0], 'LineStyle',':');

% Analytical surface
nexttile(2);
trisurf(tris10(:,1:3), nodes(:,1), nodes(:,2), PsiShown, 'EdgeColor','k'); view(2); colorbar;
axis equal tight off; title('Analytical \Psi (inside set to Dirichlet)'); caxis([0 50]);
hold on; viscircle_local(C(1).c, C(1).R, 'Color',[0 0 0], 'LineStyle',':');
hold on; viscircle_local(C(2).c, C(2).R, 'Color',[0 0 0], 'LineStyle',':');

% Contours
bb = [min(nodes(:,1)) max(nodes(:,1)) min(nodes(:,2)) max(nodes(:,2))];
nx = 350; ny = 350;
[xg, yg] = meshgrid(linspace(bb(1),bb(2),nx), linspace(bb(3),bb(4),ny));
ug   = griddata(nodes(:,1), nodes(:,2), uShow,   xg, yg, 'linear');
psig = griddata(nodes(:,1), nodes(:,2), PsiShown, xg, yg, 'linear');
levels = [10 20 30 40];

% nexttile(3);
figure; hold on;
contour(xg, yg, ug, levels, 'LineWidth', 1.2, 'ShowText','on');
axis equal tight; title('CutFEM u: contours 10/20/30/40');
hold on; viscircle_local(C(1).c, C(1).R, 'Color',[0 0 0], 'LineStyle','-');
hold on; viscircle_local(C(2).c, C(2).R, 'Color',[0 0 0], 'LineStyle','-');
R = 10;
th = linspace(0, 2*pi, 720);
plot(R*cos(th), R*sin(th), 'k-', 'LineWidth', 1.8);   % outline
axis equal;  % keep circle round
hold off;

% nexttile(4);
figure;
contour(xg, yg, psig, levels, 'LineWidth', 1.2, 'ShowText','on');
axis equal tight; title('Analytical \Psi: contours 10/20/30/40');
hold on; viscircle_local(C(1).c, C(1).R, 'Color',[0 0 0], 'LineStyle','-');
hold on; viscircle_local(C(2).c, C(2).R, 'Color',[0 0 0], 'LineStyle','-');
R = 10;
th = linspace(0, 2*pi, 720);
plot(R*cos(th), R*sin(th), 'k-', 'LineWidth', 1.8);   % outline
axis equal;  % keep circle round
hold off;
%% ===== Compact interface neighborhood viz (fixed) =====
icShow      = 1;        % which circle (1 = hot, 2 = cold)
ringDepth   = 1;        % k-ring by shared EDGE
maxCutShow  = 24;       % cap # of cut parents
tol         = 1e-12;
shrinkCut   = 0.90;     % shrink split polygons (visual only)
shrinkNon   = 0.98;     % slight shrink for non-cut context

cc = C(icShow).c;  RR = C(icShow).R;
V3 = tris10(:,1:3);                 % P1 vertex connectivity
Nt = size(V3,1);

% --- detect cut parents and chord endpoints
isCut  = false(Nt,1);
PQ     = zeros(Nt,4);
cutRank = [];   % [elemId, |dist(centroid to circle) - R|]
for it = 1:Nt
    Xv = nodes(V3(it,:),:);
    [ok,P,Q] = chord_in_triangle_linear(Xv, cc, RR, tol);
    if ok
        isCut(it) = true;
        PQ(it,:)  = [P Q];
        ctri = mean(Xv,1);
        cutRank = [cutRank; it, abs(norm(ctri-cc) - RR)]; %#ok<AGROW>
    end
end
if isempty(cutRank)
    warning('No cut triangles detected for circle %d.', icShow);
    return;
end
cutRank = sortrows(cutRank,2);
cutKeep = cutRank(1:min(maxCutShow, size(cutRank,1)),1);
cutMask = false(Nt,1); cutMask(cutKeep) = true;

% --- build EDGE adjacency (robust)
A = edge_adjacency_clean(V3);   % sparse Nt x Nt, 1 if share an edge

% --- k-ring grow (neighbors by edges)
showMask = cutMask;
front = cutMask;
for r = 1:ringDepth
    nbr = (A * front) > 0;   % logical column
    front = nbr & ~showMask;
    showMask = showMask | front;
end

% --- also include inside neighbors that touch kept cuts by an edge
centIn = false(Nt,1);
for it = 1:Nt
    Xv = nodes(V3(it,:),:);
    centIn(it) = (norm(mean(Xv,1) - cc) < RR - 1e-10);
end
insideNbr = ((A * cutMask) > 0) & centIn;
showMask = showMask | insideNbr;

elist = find(showMask);
cutInShow = intersect(find(isCut), elist);
nonCut    = setdiff(elist, cutInShow);

% split non-cuts into inside/outside by centroid
isIn_non  = false(size(nonCut));
for k = 1:numel(nonCut)
    it = nonCut(k);
    Xv = nodes(V3(it,:),:);
    isIn_non(k) = (norm(mean(Xv,1)-cc) < RR - 1e-10);
end
nonCutIn  = nonCut(isIn_non);
nonCutOut = nonCut(~isIn_non);

% --- CUT geometry (each sub-triangle its own color)
Vcut = []; Tcut = []; Fcol = [];
rng(123);
facePalette = hsv(max(1, 6*numel(cutInShow))) * 0.85;
facePtr = 1;

for it = cutInShow.'
    Xv = nodes(V3(it,:),:);
    P = PQ(it,1:2); Q = PQ(it,3:4);

    [polyA, polyB] = split_triangle_by_line_PQ_strict(Xv, P, Q);
    if isempty(polyA) || isempty(polyB), continue; end

    % inside/outside by centroid to circle center
    if norm(mean(polyA,1)-cc) < norm(mean(polyB,1)-cc)
        polyIn=polyA; polyOut=polyB;
    else
        polyIn=polyB; polyOut=polyA;
    end

    % shrink (visual gap)
    ctr = mean(Xv,1);
    polyIn  = ctr + shrinkCut*(polyIn  - ctr);
    polyOut = ctr + shrinkCut*(polyOut - ctr);

    % fan triangulate (returns (3m x 2) stacked)
    tIn  = fan_tris(polyIn);
    tOut = fan_tris(polyOut);

    if ~isempty(tOut)
        base = size(Vcut,1);
        Vcut = [Vcut; tOut]; %#ok<AGROW>
        nF = size(tOut,1)/3;
        Tloc = reshape((base+1):(base+3*nF), 3, []).';
        Tcut = [Tcut; Tloc]; %#ok<AGROW>
        Fcol = [Fcol; facePalette(facePtr:facePtr+nF-1,:)]; facePtr = facePtr + nF; %#ok<AGROW>
    end
    if ~isempty(tIn)
        base = size(Vcut,1);
        Vcut = [Vcut; tIn]; %#ok<AGROW>
        nF = size(tIn,1)/3;
        Tloc = reshape((base+1):(base+3*nF), 3, []).';
        Tcut = [Tcut; Tloc]; %#ok<AGROW>
        Fcol = [Fcol; facePalette(facePtr:facePtr+nF-1,:)]; facePtr = facePtr + nF; %#ok<AGROW>
    end
end

FVCDcut = zeros(size(Vcut,1),3);
for f = 1:size(Tcut,1)
    FVCDcut(Tcut(f,:),:) = repmat(Fcol(f,:), 3, 1);
end

% --- NON-CUT context (uniform color, dotted outlines, semi-transparent)
colOut = [0.90 0.94 1.00];
colIn  = [1.00 0.92 0.96];

Vout=[]; Tout=[];
for it = nonCutOut
    Xv = nodes(V3(it,:),:);
    ctr = mean(Xv,1); Xs = ctr + shrinkNon*(Xv - ctr);
    base = size(Vout,1);
    Vout = [Vout; Xs]; %#ok<AGROW>
    Tout = [Tout; base+(1:3)]; %#ok<AGROW>
end

Vin=[]; Tin=[];
for it = nonCutIn
    Xv = nodes(V3(it,:),:);
    ctr = mean(Xv,1); Xs = ctr + shrinkNon*(Xv - ctr);
    base = size(Vin,1);
    Vin = [Vin; Xs]; %#ok<AGROW>
    Tin = [Tin; base+(1:3)]; %#ok<AGROW>
end

% --- PLOT: draw non-cuts first (transparent), then cut pieces on top
figure('Color','w'); hold on; axis equal; view(2); box on;

if ~isempty(Tout)
    patch('Faces',Tout,'Vertices',Vout, ...
          'FaceColor',colOut,'EdgeColor',[0.1 0.2 0.5], ...
          'LineStyle',':','LineWidth',0.9,'FaceAlpha',0.45);
end
if ~isempty(Tin)
    patch('Faces',Tin,'Vertices',Vin, ...
          'FaceColor',colIn,'EdgeColor',[0.5 0.1 0.4], ...
          'LineStyle',':','LineWidth',0.9,'FaceAlpha',0.45);
end
if ~isempty(Tcut)
    patch('Faces',Tcut,'Vertices',Vcut, ...
          'FaceVertexCData',FVCDcut,'FaceColor','flat', ...
          'EdgeColor',[0.15 0.15 0.15],'LineWidth',0.7);
end

% circle outline
th = linspace(0,2*pi,400);
plot(cc(1)+RR*cos(th), cc(2)+RR*sin(th), 'k-', 'LineWidth', 1.2);

title(sprintf('Interface %d neighborhood (edge %d-ring + inside neighbors)', icShow, ringDepth));
xlabel('x'); ylabel('y'); axis tight;

%% ------------ helpers (viz only) ------------
function A = edge_adjacency_clean(V3)
    % robust edge adjacency: two elems adjacent iff they share a full edge
    Nt = size(V3,1);
    E1 = sort(V3(:,[1 2]),2);
    E2 = sort(V3(:,[2 3]),2);
    E3 = sort(V3(:,[3 1]),2);
    A = spalloc(Nt,Nt, 6*Nt);
    % accumulate per-edge band
    A = A + band_adj(E1) + band_adj(E2) + band_adj(E3);
    A = A - diag(diag(A));  % zero diagonal
end

function B = band_adj(Eband)
    % Eband: Nt x 2 sorted node ids for the same local edge across all elems
    Nt = size(Eband,1);
    [Euniq,~,ic] = unique(Eband,'rows','stable'); %#ok<ASGLU>
    buckets = accumarray(ic, (1:Nt).', [size(Euniq,1),1], @(ix){ix});
    ii=[]; jj=[];
    for b = 1:numel(buckets)
        els = buckets{b};
        if numel(els) >= 2
            % fully connect this small set (edge is shared by 2 elements in a manifold mesh)
            [p,q] = ndgrid(els,els);
            mask = p ~= q;
            ii = [ii; p(mask)]; %#ok<AGROW>
            jj = [jj; q(mask)]; %#ok<AGROW>
        end
    end
    B = sparse(ii,jj,1,Nt,Nt);
end

function [polyA, polyB] = split_triangle_by_line_PQ_strict(Xv, P, Q)
    % exact partition by infinite line through P->Q; returns two polygons
    v = Q - P; vz = [v 0];
    function s = sgn(X)
        s = cross(vz, [X-P 0]); s = s(3);
        if abs(s) < 1e-14, s = 0; end
    end
    V = Xv; E = [1 2; 2 3; 3 1];
    s = zeros(3,1);
    for i=1:3, s(i)=sgn(V(i,:)); end
    A=[]; B=[];
    for e=1:3
        i=E(e,1); j=E(e,2);
        Xi=V(i,:); Xj=V(j,:);
        si=s(i);   sj=s(j);
        if si>=0, A=[A; Xi]; else, B=[B; Xi]; end %#ok<AGROW>
        if (si>0 && sj<0) || (si<0 && sj>0)
            ti=abs(si); tj=abs(sj); t=ti/(ti+tj);
            Xint=Xi + t*(Xj - Xi);
            A=[A; Xint]; B=[B; Xint]; %#ok<AGROW>
        elseif si==0 && sj~=0
            % Xi on the line → include Xi in both to avoid tiny gaps
            A=[A; Xi]; B=[B; Xi]; %#ok<AGROW>
        end
    end
    A = unique_loop(A);  B = unique_loop(B);
    if poly_area(A) < 0, A = flipud(A); end
    if poly_area(B) < 0, B = flipud(B); end
    polyA=A; polyB=B;
end

function P2 = unique_loop(P)
    if isempty(P), P2=P; return; end
    keep = true(size(P,1),1);
    for i=2:size(P,1)
        if norm(P(i,:)-P(i-1,:)) < 1e-12, keep(i) = false; end
    end
    if norm(P(end,:)-P(1,:)) < 1e-12, keep(end) = false; end
    P2 = P(keep,:);
end

function A = poly_area(P)
    x=P(:,1); y=P(:,2); x2=[x; x(1)]; y2=[y; y(1)];
    A = 0.5*sum(x2(1:end-1).*y2(2:end) - x2(2:end).*y2(1:end-1));
end

function triStack = fan_tris(P)
    if size(P,1) < 3, triStack = zeros(0,2); return; end
    c = mean(P,1);
    ang = atan2(P(:,2)-c(2), P(:,1)-c(1));
    [~,ord] = sort(ang); P = P(ord,:);
    triStack = [];
    for i=2:size(P,1)-1
        triStack = [triStack; P(1,:); P(i,:); P(i+1,:)]; %#ok<AGROW>
    end
end

%% ------------- helpers (viz only) -----------------
function vals = eval_u_stack(triStack, X10, ue, p3basis)
    % triStack: (3m x 2) stacked vertices; return (3m x 1) values
    m = size(triStack,1);
    vals = zeros(m,1);
    for i = 1:m
        xp = triStack(i,:);
        [xi,eta,ok] = invert_isoP3_newton(xp, X10, p3basis);
        if ~ok
            % very rare; fall back to nearest of the 3 corners
            d = vecnorm((X10(1:3,:) - xp),2,2); [~,id] = min(d);
            Ni = zeros(10,1); Ni(id) = 1;
            vals(i) = (Ni.' * ue);
        else
            [N,~,~] = p3basis(xi,eta);
            vals(i) = (N.' * ue);
        end
    end
end

function [polyA, polyB] = split_triangle_by_line_PQ(Xv, P, Q)
% exact partition of triangle Xv by infinite line through P->Q
    v = Q - P;  vz = [v 0];
    function s = sgn(X)
        w = X - P; s = cross(vz, [w 0]); s = s(3);
    end
    V = Xv; E = [1 2; 2 3; 3 1];
    s = zeros(3,1);
    for i=1:3, s(i) = sgn(V(i,:)); end
    A = []; B = [];
    for e = 1:3
        i = E(e,1); j = E(e,2);
        Xi = V(i,:); Xj = V(j,:);
        si = s(i);   sj = s(j);
        if si >= 0, A = [A; Xi]; else, B = [B; Xi]; end
        if (si>0 && sj<0) || (si<0 && sj>0)
            % interpolate intersection with signed distances
            ti = abs(si); tj = abs(sj); t = ti/(ti+tj);
            Xint = Xi + t*(Xj - Xi);
            A = [A; Xint]; B = [B; Xint];
        end
    end
    A = unique_close_loop(A);  B = unique_close_loop(B);
    if poly_area(A) < 0, A = flipud(A); end
    if poly_area(B) < 0, B = flipud(B); end
    polyA = A; polyB = B;
end

function P2 = unique_close_loop(P)
    if isempty(P), P2 = P; return; end
    keep = true(size(P,1),1);
    for i = 2:size(P,1)
        if norm(P(i,:) - P(i-1,:)) < 1e-12, keep(i) = false; end
    end
    if norm(P(end,:) - P(1,:)) < 1e-12, keep(end) = false; end
    P2 = P(keep,:);
end



function triStack = fan_triangulate(P)
% returns (3m x 2) stacked triangle vertices (unshared)
    if size(P,1) < 3, triStack = zeros(0,2); return; end
    c = mean(P,1);
    ang = atan2(P(:,2)-c(2), P(:,1)-c(1));
    [~,ord] = sort(ang); P = P(ord,:);
    triStack = [];
    for i = 2:size(P,1)-1
        triStack = [triStack; P(1,:); P(i,:); P(i+1,:)]; %#ok<AGROW>
    end
end


%% ==================== Local helpers ====================

function [N, dNxi, dNeta] = p3IsoShape_local(xi, eta)
    % P3 (10-node) triangle (Gmsh-like ordering)
    L1 = 1 - xi - eta;  L2 = xi;  L3 = eta;
    N = zeros(10,1); dNxi=zeros(10,1); dNeta=zeros(10,1);
    N(1) = L1*(3*L1 - 1)*(3*L1 - 2)/2;  N(2) = L2*(3*L2 - 1)*(3*L2 - 2)/2;  N(3) = L3*(3*L3 - 1)*(3*L3 - 2)/2;
    function [dNx,dNe] = dCorner(L,dLx,dLe), f=(3*L-1).*(3*L-2)/2; df=(9*(2*L-1))/2; dNx=dLx.*f + L.*(df.*dLx); dNe=dLe.*f + L.*(df.*dLe); end
    [dNxi(1), dNeta(1)] = dCorner(L1,-1,-1); [dNxi(2), dNeta(2)] = dCorner(L2, 1, 0); [dNxi(3), dNeta(3)] = dCorner(L3, 0, 1);
    c = 9/2;
    N(4)=c*L1*L2*(3*L1-1); N(5)=c*L1*L2*(3*L2-1); N(6)=c*L2*L3*(3*L2-1);
    N(7)=c*L2*L3*(3*L3-1); N(8)=c*L3*L1*(3*L3-1); N(9)=c*L3*L1*(3*L1-1); N(10)=27*L1*L2*L3;
    function [dNx,dNe]=dProd3(s,X,Y,Z,dXx,dXe,dYx,dYe,dZx,dZe), dNx=s*( dXx*Y*Z + X*dYx*Z + X*Y*dZx ); dNe=s*( dXe*Y*Z + X*dYe*Z + X*Y*dZe ); end
    dL1x=-1; dL1e=-1; dL2x=1; dL2e=0; dL3x=0; dL3e=1;
    [dNxi(4), dNeta(4)]   = dProd3(c,L1,L2,(3*L1-1), dL1x,dL1e, dL2x,dL2e, 3*dL1x,3*dL1e);
    [dNxi(5), dNeta(5)]   = dProd3(c,L1,L2,(3*L2-1), dL1x,dL1e, dL2x,dL2e, 3*dL2x,3*dL2e);
    [dNxi(6), dNeta(6)]   = dProd3(c,L2,L3,(3*L2-1), dL2x,dL2e, dL3x,dL3e, 3*dL2x,3*dL2e);
    [dNxi(7), dNeta(7)]   = dProd3(c,L2,L3,(3*L3-1), dL2x,dL2e, dL3x,dL3e, 3*dL3x,3*dL3e);
    [dNxi(8), dNeta(8)]   = dProd3(c,L3,L1,(3*L3-1), dL3x,dL3e, dL1x,dL1e, 3*dL3x,3*dL3e);
    [dNxi(9), dNeta(9)]   = dProd3(c,L3,L1,(3*L1-1), dL3x,dL3e, dL1x,dL1e, 3*dL1x,3*dL1e);
    [dNxi(10), dNeta(10)] = dProd3(27,L1,L2,L3, dL1x,dL1e, dL2x,dL2e, dL3x,dL3e);
end

function [g_pt, g_wt] = triGaussPoints12_local()
    g_pt = [
      0.24928674517091, 0.24928674517091;
      0.24928674517091, 0.50142650965818;
      0.50142650965818, 0.24928674517091;
      0.06308901449150, 0.06308901449150;
      0.06308901449150, 0.87382197101700;
      0.87382197101700, 0.06308901449150;
      0.31035245103378,  0.63650249912140;
      0.63650249912140,  0.05314504984482;
      0.05314504984482,  0.31035245103378;
      0.63650249912140,  0.31035245103378;
      0.31035245103378,  0.05314504984482;
      0.05314504984482,  0.63650249912140 ];
    g_wt = [
      0.05839313786319; 0.05839313786319; 0.05839313786319;
      0.02542245318510; 0.02542245318510; 0.02542245318510;
      0.04142553780919; 0.04142553780919; 0.04142553780919;
      0.04142553780919; 0.04142553780919; 0.04142553780919 ];
end

function [x,w] = gauss_legendre_1D_local(n)
    switch n
      case 2, x=[-1/sqrt(3); 1/sqrt(3)]; w=[1;1];
      case 3, x=[-sqrt(3/5); 0; sqrt(3/5)]; w=[5/9; 8/9; 5/9];
      case 4, x=[-0.8611363116; -0.3399810436; 0.3399810436; 0.8611363116]; w=[0.3478548451; 0.6521451549; 0.6521451549; 0.3478548451];
      case 5, x=[-0.9061798459; -0.5384693101; 0; 0.5384693101; 0.9061798459]; w=[0.2369268851; 0.4786286705; 0.5688888889; 0.4786286705; 0.2369268851];
      otherwise, error('n must be 2..5');
    end
end

function ang = wrapToPi_local(ang), ang = mod(ang + pi, 2*pi) - pi; end

function [xi,eta,ok] = invert_isoP3_newton(xp, X10, p3basis_fun)
    xi=1/3; eta=1/3; ok=false;
    for it=1:30
        [N,dNxi,dNeta]=p3basis_fun(xi,eta);
        xhat = (N.'*X10); r = xp - xhat;
        if norm(r) < 1e-12, ok=true; return; end
        J = [dNxi.'; dNeta.'] * X10; if rcond(J) < 1e-12, return; end
        step = (J \ r.').'; xi=xi+step(1); eta=eta+step(2);
        if xi < -0.35, xi=-0.35; end; if eta < -0.35, eta=-0.35; end
        if xi+eta > 1.35, t=(xi+eta-1.35); xi=xi-t/2; eta=eta-t/2; end
    end
    [N,~,~]=p3basis_fun(xi,eta);
    if norm(xp-(N.'*X10)) < 1e-8, ok=true; end
end

function [ok, P, Q] = chord_in_triangle_linear(Xv, cxy, R, tol)
% Intersect linear (vertex) triangle with circle; return chord endpoints P,Q.
    pts = [];
    E = [1 2; 2 3; 3 1];
    for e = 1:3
        i = E(e,1); j = E(e,2);
        P0 = Xv(i,:); P1 = Xv(j,:);
        tt = seg_circle_t_local(P0,P1,cxy,R,tol);
        for k = 1:numel(tt)
            t = tt(k); pts = [pts; P0 + t*(P1 - P0)]; %#ok<AGROW>
        end
    end
    for k=1:3
        if abs(norm(Xv(k,:) - cxy) - R) < 1e-10
            pts = [pts; Xv(k,:)]; %#ok<AGROW>
        end
    end
    if isempty(pts), ok=false; P=[0 0]; Q=[0 0]; return; end
    pts = unique(round(pts,12),'rows');
    if size(pts,1) < 2, ok=false; P=[0 0]; Q=[0 0]; return; end
    if size(pts,1) > 2
        [ii,jj] = farthest_pair_local(pts); pts = pts([ii jj],:);
    end
    P = pts(1,:); Q = pts(2,:); ok = norm(P-Q) > 1e-14;
end

function [iBest,jBest] = farthest_pair_local(P)
    m = size(P,1); iBest=1; jBest=2; dBest=-inf;
    for i = 1:m-1
        Pi = P(i,:);
        for j = i+1:m
            d = (Pi(1)-P(j,1))^2 + (Pi(2)-P(j,2))^2;
            if d > dBest, dBest=d; iBest=i; jBest=j; end
        end
    end
end

function tt = seg_circle_t_local(P0, P1, cxy, R, tol)
    v  = P1 - P0; a=dot(v,v); b=2*dot(v, P0 - cxy); c=dot(P0 - cxy, P0 - cxy) - R^2;
    disc=b*b-4*a*c; tt=[];
    if disc < -1e-14, return; end
    disc = max(disc, 0);
    t1 = (-b - sqrt(disc)) / (2*a);
    t2 = (-b + sqrt(disc)) / (2*a);
    if t1>-tol && t1<1+tol, tt(end+1)=min(max(t1,0),1); end %#ok<AGROW>
    if t2>-tol && t2<1+tol, tt(end+1)=min(max(t2,0),1); end %#ok<AGROW>
    tt=unique(tt);
end

function he = projected_he_vertices(Xv, nHat)
% h_e = 2|T| / sum_k |n·n_k| |e_k| using the vertex triangle
    A = 0.5*abs(det([Xv(2,:)-Xv(1,:); Xv(3,:)-Xv(1,:)]));
    E=[1 2; 2 3; 3 1]; denom=0;
    for k=1:3
        P=Xv(E(k,1),:); Q=Xv(E(k,2),:); evec=Q-P; len=norm(evec);
        if len < 1e-14, continue; end
        tEdge=evec/len; nEdge=[tEdge(2), -tEdge(1)];
        denom = denom + abs(dot(nHat,nEdge)) * len;
    end
    he = 2*A / max(denom,1e-14);
end

function [Psi, R] = analytical_two_cylinders_local(nodes)
% Greenberg/bipolar closed form (x1=3, x2=2; hot=50, cold=0)
    x1=3; x2=2;
    a = (x1*x2 + 1 + sqrt((x1*x1 -1)*(x2*x2 - 1)))/(x1 + x2);
    R = (x1*x2 - 1 - sqrt((x1*x1 -1)*(x2*x2 - 1)))/(x1 - x2);
    x = nodes(:,1); y = nodes(:,2);
    ratio = ((x - a).^2 + y.^2) ./ ((a.*x - 1).^2 + (a^2).*(y.^2));
    Psi = 50/log(R) * (log(R) - 0.5*log(ratio));
end

function viscircle_local(cxy, R, varargin)
    th = linspace(0,2*pi,400);
    plot(cxy(1)+R*cos(th), cxy(2)+R*sin(th), varargin{:});
end