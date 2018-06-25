function ubs_dots
% Plot undolded band structure

%% Init. parameters
KPATH = [0 1/2 0; ...
        0 0 0; ...
        0 0 1/2]; % k-point path
FOLDS = [1 2 3]; % multiplicity in the corresponding direction used when constructing the super-cell
ERANGE = [-0.5 0.5]; % energy range for plot (Ry)
finpt = 'fold2Bloch.out'; % input file name
Ef = 0.0; % Fermi energy (Ry)
ha2ev = 2*13.605698066; % Ha -> eV conversion factor
pwr = 1/1; % power for result plotting
         % 1 - linear scale, 1/2 - sqrt, etc.
         % 0 - folded bands (needs wth = 0)
msz = 10; % marker size for plot
lwdth = 0.5; % plot line width
PLTSZ = [300 300 4*150 2*200]; % plot size
wth = 0.00; % threshold weight
clrmp = jet;    % flipud(gray)
                % flipud(pink)
                % flipud(hot)
                % flipud(autumn)
                % cool
                % flipud(bone)
                % flipud(jet)
                % jet
G = [0.33333 0.0000000 0.0000000;
   0.000000  0.16667 0.0000000;
   0.000000  0.000000  0.11111]; % Reciprocal latt. vect. [Bohr^-1] from *.out

%% INITIALIZATION
[KEIG, EIG, W] = readinput(finpt); % read input data from file
% EIG - energy eigenvalues
% KEIG - k-list for eigenvalues
% W - list of characters

%% MAIN
L = [];
ENE = [];
WGHT = [];
for i=1 : 3
    G(i,:)=G(i,:)*FOLDS(i); % rescale reciprocal lattice vectors 
end                         % from supercell to primitive cell
dl = 0; % cumulative length of the path
KPATH = coordTransform(KPATH,G);
KEIG = coordTransform(KEIG,G);
for ikp = 1 : size(KPATH,1)-1
    for j = 1 : length(EIG)
        if EIG(j) > ERANGE(1) && EIG(j) < ERANGE(2) && W(j) >= wth
            dist = dp2l( KEIG(j,:) , KPATH(ikp,:) , KPATH(ikp+1,:) );
            B = KPATH(ikp,:) - KPATH(ikp+1,:);
            dk = sqrt(dot(B,B));
            if dist < eps % k-point is on the path
                A = KPATH(ikp,:) - KEIG(j,:);
                x = dot(A,B)/dk;
                if x > 0  &&  x-dk < eps % k-point is within the path range
                    L = [L; x+dl]; % append k-point coordinate along the path
                    ENE = [ENE; EIG(j)]; % append energy list
                    WGHT = [WGHT; W(j)];
                end
            end
        end
    end
    dl = dl + dk;
end

%% Plot results
hFig = figure(1);

% Fig 1(a)
subplot(1,2,1);
set(hFig, 'Position', PLTSZ)
colormap(clrmp);
WGHTRS = rescale(WGHT,pwr);
scatter(L,(ENE-Ef)*ha2ev, WGHTRS*msz, WGHTRS,'LineWidth',lwdth);
axis([0 max(L) min((ENE-Ef)*ha2ev) max((ENE-Ef)*ha2ev)])
xlabel('Wave vector (Y - G - X)')
ylabel('Energy (eV)')
box on

% Fig 1(b)
subplot(1,2,2);
DAT = linspace(0,1,10);
DATX = ones(size(DAT));
DATRS = rescale(DAT,pwr);
scatter(DATX,DAT, DATRS*msz, DATRS,'LineWidth',lwdth);
ylabel('Spectral weight')

% -------------------------------------------------------------------------
function W = coordTransform(V,G)
% transform vector V(:,3) in G(3,3) coord. system -> W(:,3) in Cartesian coordinates
% G vector elements are in columns!
W = zeros(size(V));
G = G'; % transform G
for i = 1:length(V)
    W(i,:) = G(1,:)*V(i,1) + G(2,:)*V(i,2) + G(3,:)*V(i,3);
end;
% -------------------------------------------------------------------------
function WRESCL = rescale(W,pwr)
% rescale weights W = [0...1] using a power functio W^pwr
WRESCL=W.^(pwr); % rescale if needed to enhance
if abs( max(WRESCL)-min(WRESCL) ) > eps
    WRESCL = (WRESCL-min(WRESCL))/(max(WRESCL)-min(WRESCL)) + eps;
end; % need eps to make plot "heapy"
% -------------------------------------------------------------------------
function [KEIG, EIG, W] = readinput(filename)
% read input data
DATA = importdata(filename);
KEIG = DATA(:,1:3);
EIG = DATA(:,4);
W = DATA(:,5);
% -------------------------------------------------------------------------
function RES = dp2l(X0,X1,X2) % distance from point {X0} to line {X1}-{X2}
% see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
denom = X2 - X1;
denomabs = sqrt(dot(denom,denom));
if denomabs < eps
    display(X1); display(X2);
    error('X1 = X2');
end;
numer = cross( X0-X1 , X0-X2 );
numerabs = sqrt(dot(numer,numer));
RES = numerabs/denomabs;
% -------------------------------------------------------------------------
