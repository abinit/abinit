%EffMass   Effective Mass of semiconductors
%    EffMass(filename, Band, nDataSet) analyse the eigenenergies of datasets in output 
%    file of ABINIT to get the curvature of band structure at Conduction Band Minimum 
%    and Valence Band Maximum. According to the formula,
%          ( 1/Mass_eff )ij = ( 1/h_bar^2 ) * 2nd_Derivative( CBM/VBM eigenenergies ) 
%                              with respect to ki to kj.
%    We can easily acquire the three main components of effective mass tensors,
%          Mass_eff_xx,  Mass_eff_yy, Mass_eff_zz.
%
%    [EffectiveMass, K, eigenE] = EffMass(filename, Band, nDataSet)
%    Varible Band indicates which band requires to be analysed.
%    Varible nDataSet indicates the number of dataset contained in 'filename' output file.
%    Because the data is derived from band structure calculation following a SCF total
%    energy calculation, the nDataSet refers to dataset 2 to (1+nDataSet).
%
%    Written by Mike Zhang, 8/11/2002, copyright reserved.
%    Used by your own warranty, Mike Zhang(mz24cn@hotmail.com)
%    Copyright (C) 2002-2020 ABINIT group (Mike Zhang)
%    This file is distributed under the terms of the
%    GNU General Public License, see ~abinit/COPYING
%    or http://www.gnu.org/copyleft/gpl.txt .
%    For the initials of contributors, see ~abinit/doc/developers/contributors.txt .


function [EffectiveMass, K, eigenE] = EffMass(filename, Band, nDataSet)

% Currently we presumably nBand and nKpt of every dataset are same.
nBand = mexData([1 1], filename, '== DATASET  2 ==', 'newkpt: treating', '', ''); 
nKpt = mexData([1 1], filename, '== DATASET  2 ==', 'Eigenvalues (   eV  ) for nkpt=', '', ''); 
% Find the HOMO(Highest Occupied Molecule Orbit) and LUMO(Lowest Unoccupied Molecule Orbit)
occ = mexData([1 nBand], filename, '-outvars: echo values of preprocessed input variables --------', 'occ', '', '');
i = 1;
while occ(i) ~= 0
    i = i + 1;
end
HOMO = i - 1;
LUMO = i;

nGroup = fix(nKpt / 2);  % Output nGroup groups Effective Mass data

eigenE = mexData([nKpt nBand nDataSet], filename, '== DATASET  2 ==', '(reduced coord)', 'Eigenvalues (   eV  ) for nkpt=', '');
% Search CBM/VBM
indexExtreE(nDataSet) = 0;
for i = 1:nDataSet
    indexExtreE(i) = 1;
    for n = 1:nKpt
        if eigenE(n, LUMO, i) <= eigenE(indexExtreE(i), LUMO, i)
            indexExtreE(i) = n;
        end
    end
end
eigenEha = eigenE / 27.2113961;        % eigenEha is in atom unit (Bohr)

k = mexData([nKpt 3 nDataSet], filename, '== DATASET  2 ==', 'kpt=', 'Eigenvalues (   eV  ) for nkpt=', '');
% Now alter k from reduced coordinates to cartesian coordinates in atom unit
unitK1 = mexData([1 3], filename, '== DATASET  2 ==', 'G(1)=', '', '');
unitK2 = mexData([1 3], filename, '== DATASET  2 ==', 'G(2)=', '', '');
unitK3 = mexData([1 3], filename, '== DATASET  2 ==', 'G(3)=', '', '');
K(nKpt, 3, nDataSet) = 0;
for i = 1:3
    K(:, i, :) = (k(:, 1, :)*unitK1(i) + k(:, 2, :)*unitK2(i) + k(:, 3, :)*unitK3(i)) * 2 * pi;
end
scalarK(nKpt, nDataSet) = 0;
for i = 1:nDataSet
    scalarK(:, i) = sqrt((K(:, 1, i)-K(1, 1, i)).^2 + (K(:, 2, i)-K(1, 2, i)).^2 + (K(:, 3, i)-K(1, 3, i)).^2);
end

% Now solve the three main components of effective mass tensor
Curvature(nGroup, 3, nDataSet) = 0;
% Enter main loop of nDataSet datasets
for i = 1:nDataSet
    for j = 1:nGroup
        nMin = indexExtreE(i) - j;
        if nMin < 1
            nMin = 1;
        end
        nMax = indexExtreE(i) + j;
        if nMax > nKpt
            nMax = nKpt;
        end
        if nMax-nMin < 3
            if nMin == 1
                nMax = 3;
            else
                nMin = nMax -3;
            end
        end
        Curvature(j, :, i) = polyfit(scalarK(nMin:nMax, i), eigenEha(nMin:nMax, Band, i), 2);        % K(nMin:nMax, 1, i)
    end
end
EffectiveMass = 0.5 ./ Curvature(:, 1, :);  % 2nd_Derivative( CBM/VBM eigenenergies ) = 2 * Curvature
EffectiveMass
