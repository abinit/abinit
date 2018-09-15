---
authors: RC
---

# List of 230 3D symmetry space groups

This page presents the list of the 230 3D symmetry space groups with
characteristics, and also, if applicable, derivatives of the space group. 

This list of symmetry groups is part of the on-line help of the ABINIT code.
Besides *classical* crystallographic information, it provides the values of
certain symmetry-related variables used in the input file of the ABINIT code.

The table entries are as follows:

SPGROUP
:   *SPace GROUP number* as found in the International Tables of Crystallography, 1983 edition [[cite:Hahn1983]].

SPGAXOR
:   *SPace Group AXis ORientation*, the orientation of the unit cell axis. The allowed values are: 

    * Trigonal groups:
        * 1 represents the hexagonal axes
        * 2 represents the rhombohedral axes
    * Orthorhombic groups:
        * 1 abc -> abc
        * 2 abc -> cab
        * 3 abc -> bca
        * 4 abc -> acb
        * 5 abc -> bac
        * 6 abc -> cba
    * Monoclinic: 3 or 9 possibilities depending on the space group

SPGORIG
:   *SPace Group ORIGin*, the position of the origin in the unit cell. 
    The allowed values are 1 and 2. they correspond to the actual choices in the 
    International Tables of Crystallography, 1986 edition. 

BRVLTT
:   *BRaVais LaTTice*, the type of Bravais lattice. The allowed values are:

    * 1 = Primitive with no associated translations
    * 2 = Inner centered with (a/2 + b/2 + c/2) associated translation
    * 3 = Face centered with (a/2 + b/2; b/2 + c/2; c/2 + a/2) associated translations
    * 4 = C - centered with (a/2 + b/2) associated translation
    * 5 = A - centered with (b/2 + c/2) associated translation
    * 6 = B - centered with (c/2 + a/2) associated translation
    * 7 = Rhombohedral lattice

INTERNATIONAL
: the INTERNATIONAL notation of the space group

SCHOENFLIES
: the equivalent SCHOENFLIES notation of the space group

MULTIPLICITY (x y z)
: the maximum multiplicity, after applying Bravais lattice translations,  corresponding to a general position (x y z)


## List of 230 3D symmetry space groups

This section presents the list of the 230 3D symmetry space groups with
characteristics, and also, if applicable, derivatives of the space group. 

<TABLE BORDER = "2" cellspacing = "2" cellpadding = "5">
<tr>
<td width = 10%>SPGROUP </td>
<td width = 10%> SPGAXOR </td>
<td width = 10%> SPGORIG </td>
<td width = 10%> BRVLTT </td>
<td width = 22%> INTERNATIONAL </td>
<td width = 20%> SCHOENFLIES </td>
<td width = 10%> ORDER </td>
</tr>

<tr><td>  1  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P1  </td><td>  C1<SUP>1</SUP>  </td><td>  1  </td></tr>
<tr><td>  2  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-1  </td><td>  Ci<SUP>1</SUP>  </td><td>  2  </td></tr>
<tr><td>  3  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P2_b = P121  </td><td>  C2<SUP>1</SUP>  </td><td>  2  </td></tr>
<tr><td>  3  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  P2_c = P112  </td><td>  C2<SUP>1</SUP>  </td><td>  2  </td></tr>
<tr><td>  3  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2_a = P211  </td><td>  C2<SUP>1</SUP>  </td><td>  2  </td></tr>
<tr><td>  4  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P21_b = P1211  </td><td>  C2<SUP>2</SUP>  </td><td>  2  </td></tr>
<tr><td>  4  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  P21_c = P1121  </td><td>  C2<SUP>2</SUP>  </td><td>  2  </td></tr>
<tr><td>  4  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P21_a = P2111  </td><td>  C2<SUP>2</SUP>  </td><td>  2  </td></tr>
<tr><td>  5  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  C2_b1 = C121  </td><td>  C2<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  5  </td><td>  5  </td><td>  1  </td><td>  5  </td><td>  C2_b2 = A121  </td><td>  C2<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  5  </td><td>  6  </td><td>  1  </td><td>  2  </td><td>  C2_b3 = I121  </td><td>  C2<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  5  </td><td>  7  </td><td>  1  </td><td>  5  </td><td>  C2_c1 = A112  </td><td>  C2<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  5  </td><td>  8  </td><td>  1  </td><td>  6  </td><td>  C2_c2 = B112 = B2  </td><td>  C2<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  5  </td><td>  9  </td><td>  1  </td><td>  2  </td><td>  C2_c3 = I112  </td><td>  C2<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  5  </td><td>  2  </td><td>  1  </td><td>  6  </td><td>  C2_a1 = B211  </td><td>  C2<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  5  </td><td>  3  </td><td>  1  </td><td>  4  </td><td>  C2_a2 = C211  </td><td>  C2<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  5  </td><td>  4  </td><td>  1  </td><td>  2  </td><td>  C2_a3 = I211  </td><td>  C2<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  6  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pm_b = P1m1  </td><td>  Cs<SUP>1</SUP>  </td><td>  2  </td></tr>
<tr><td>  6  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pm_c = P11m  </td><td>  Cs<SUP>1</SUP>  </td><td>  2  </td></tr>
<tr><td>  6  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pm_a = Pm11  </td><td>  Cs<SUP>1</SUP>  </td><td>  2  </td></tr>
<tr><td>  7  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pc_b1 = P1c1  </td><td>  Cs<SUP>2</SUP>  </td><td>  2  </td></tr>
<tr><td>  7  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pc_b2 = P1n1  </td><td>  Cs<SUP>2</SUP>  </td><td>  2  </td></tr>
<tr><td>  7  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  Pc_b3 = P1a1  </td><td>  Cs<SUP>2</SUP>  </td><td>  2  </td></tr>
<tr><td>  7  </td><td>  7  </td><td>  1  </td><td>  1  </td><td>  Pc_c1 = P11a  </td><td>  Cs<SUP>2</SUP>  </td><td>  2  </td></tr>
<tr><td>  7  </td><td>  8  </td><td>  1  </td><td>  1  </td><td>  Pc_c2 = P11n  </td><td>  Cs<SUP>2</SUP>  </td><td>  2  </td></tr>
<tr><td>  7  </td><td>  9  </td><td>  1  </td><td>  1  </td><td>  Pc_c3 = P11b = Pb  </td><td>  Cs<SUP>2</SUP>  </td><td>  2  </td></tr>
<tr><td>  7  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pc_a1 = Pb11  </td><td>  Cs<SUP>2</SUP>  </td><td>  2  </td></tr>
<tr><td>  7  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pc_a2 = Pn11  </td><td>  Cs<SUP>2</SUP>  </td><td>  2  </td></tr>
<tr><td>  7  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pc_a3 = Pc11  </td><td>  Cs<SUP>2</SUP>  </td><td>  2  </td></tr>
<tr><td>  8  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  Cm_b1 = C1m1  </td><td>  Cs<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  8  </td><td>  5  </td><td>  1  </td><td>  5  </td><td>  Cm_b2 = A1m1  </td><td>  Cs<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  8  </td><td>  6  </td><td>  1  </td><td>  2  </td><td>  Cm_b3 = I1m1  </td><td>  Cs<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  8  </td><td>  7  </td><td>  1  </td><td>  5  </td><td>  Cm_c1 = A11m  </td><td>  Cs<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  8  </td><td>  8  </td><td>  1  </td><td>  6  </td><td>  Cm_c2 = B11m = Bm  </td><td>  Cs<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  8  </td><td>  9  </td><td>  1  </td><td>  2  </td><td>  Cm_c3 = I11m  </td><td>  Cs<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  8  </td><td>  2  </td><td>  1  </td><td>  6  </td><td>  Cm_a1 = Bm11  </td><td>  Cs<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  8  </td><td>  3  </td><td>  1  </td><td>  4  </td><td>  Cm_a2 = Cm11  </td><td>  Cs<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  8  </td><td>  4  </td><td>  1  </td><td>  2  </td><td>  Cm_a3 = Im11  </td><td>  Cs<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  9  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  Cc_b1 = C1c1  </td><td>  Cs<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  9  </td><td>  5  </td><td>  1  </td><td>  5  </td><td>  Cc_b2 = A1n1  </td><td>  Cs<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  9  </td><td>  6  </td><td>  1  </td><td>  2  </td><td>  Cc_b3 = I1a1  </td><td>  Cs<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  9  </td><td>  7  </td><td>  1  </td><td>  5  </td><td>  Cc_c1 = A11a  </td><td>  Cs<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  9  </td><td>  8  </td><td>  1  </td><td>  6  </td><td>  Cc_c2 = B11n  </td><td>  Cs<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  9  </td><td>  9  </td><td>  1  </td><td>  2  </td><td>  Cc_c3 = I11b  </td><td>  Cs<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  9  </td><td>  2  </td><td>  1  </td><td>  6  </td><td>  Cc_a1 = Bb11  </td><td>  Cs<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  9  </td><td>  3  </td><td>  1  </td><td>  4  </td><td>  Cc_a2 = Cn11  </td><td>  Cs<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  9  </td><td>  4  </td><td>  1  </td><td>  2  </td><td>  Cc_a3 = Ic11  </td><td>  Cs<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  10  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P2/m_b = P12/m1  </td><td>  C2h<SUP>1</SUP>  </td><td>  4  </td></tr>
<tr><td>  10  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  P2/m_c = P112/m  </td><td>  C2h<SUP>1</SUP>  </td><td>  4  </td></tr>
<tr><td>  10  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2/m_a = P2/m11  </td><td>  C2h<SUP>1</SUP>  </td><td>  4  </td></tr>
<tr><td>  11  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>/m_b = P12<SUB>1</SUB>/m1  </td><td>  C2h<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  11  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>/m_c = P112<SUB>1</SUB>/m  </td><td>  C2h<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  11  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>/m_a = P2<SUB>1</SUB>/m11  </td><td>  C2h<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  12  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  C2/m_b1 = C12/m1  </td><td>  C2h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  12  </td><td>  5  </td><td>  1  </td><td>  5  </td><td>  C2/m_b2 = A12/m1  </td><td>  C2h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  12  </td><td>  6  </td><td>  1  </td><td>  2  </td><td>  C2/m_b3 = I12/m1  </td><td>  C2h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  12  </td><td>  7  </td><td>  1  </td><td>  5  </td><td>  C2/m_c1 = A112/m  </td><td>  C2h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  12  </td><td>  8  </td><td>  1  </td><td>  6  </td><td>  C2/m_c2 = B112/m = B2/m  </td><td>  C2h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  12  </td><td>  9  </td><td>  1  </td><td>  2  </td><td>  C2/m_c3 = I112/m  </td><td>  C2h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  12  </td><td>  2  </td><td>  1  </td><td>  6  </td><td>  C2/m_a1 = B2/m11  </td><td>  C2h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  12  </td><td>  3  </td><td>  1  </td><td>  4  </td><td>  C2/m_a2 = C2/m11  </td><td>  C2h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  12  </td><td>  4  </td><td>  1  </td><td>  2  </td><td>  C2/m_a3 = I2/m11  </td><td>  C2h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  13  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P2/c_b1 = P12/c1  </td><td>  C2h<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  13  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  P2/c_b2 = P12/n1  </td><td>  C2h<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  13  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  P2/c_b3 = P12/a1  </td><td>  C2h<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  13  </td><td>  7  </td><td>  1  </td><td>  1  </td><td>  P2/c_c1 = P112/a  </td><td>  C2h<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  13  </td><td>  8  </td><td>  1  </td><td>  1  </td><td>  P2/c_c2 = P112/n  </td><td>  C2h<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  13  </td><td>  9  </td><td>  1  </td><td>  1  </td><td>  P2/c_c3 = P112/b = P2/b  </td><td>  C2h<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  13  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2/c_a1 = P2/b11  </td><td>  C2h<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  13  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  P2/c_a2 = P2/n11  </td><td>  C2h<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  13  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  P2/c_a3 = P2/c11  </td><td>  C2h<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  14  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>/c_b1 = P12<SUB>1</SUB>/c1  </td><td>  C2h<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  14  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>/c_b2 = P12<SUB>1</SUB>/n1  </td><td>  C2h<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  14  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>/c_b3 = P12<SUB>1</SUB>/a1  </td><td>  C2h<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  14  </td><td>  7  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>/c_c1 = P112<SUB>1</SUB>/a  </td><td>  C2h<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  14  </td><td>  8  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>/c_c2 = P112<SUB>1</SUB>/n  </td><td>  C2h<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  14  </td><td>  9  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>/c_c3 = P112<SUB>1</SUB>/b = P2<SUB>1</SUB>/b  </td><td>  C2h<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  14  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>/c_a1 = P2<SUB>1</SUB>/b11  </td><td>  C2h<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  14  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>/c_a2 = P2<SUB>1</SUB>/n11  </td><td>  C2h<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  14  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>/c_a3 = P2<SUB>1</SUB>/c11  </td><td>  C2h<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  15  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  C2/c_b1 = C12/c1  </td><td>  C2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  15  </td><td>  5  </td><td>  1  </td><td>  5  </td><td>  C2/c_b2 = A12/n1  </td><td>  C2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  15  </td><td>  6  </td><td>  1  </td><td>  2  </td><td>  C2/c_b3 = I12/a1  </td><td>  C2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  15  </td><td>  7  </td><td>  1  </td><td>  4  </td><td>  C2/c_c1 = A112/a  </td><td>  C2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  15  </td><td>  8  </td><td>  1  </td><td>  6  </td><td>  C2/c_c2 = B112/n  </td><td>  C2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  15  </td><td>  9  </td><td>  1  </td><td>  2  </td><td>  C2/c_c3 = I112/b  </td><td>  C2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  15  </td><td>  2  </td><td>  1  </td><td>  6  </td><td>  C2/c_a1 = B2/b11  </td><td>  C2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  15  </td><td>  3  </td><td>  1  </td><td>  4  </td><td>  C2/c_a2 = C2/n11  </td><td>  C2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  15  </td><td>  4  </td><td>  1  </td><td>  2  </td><td>  C2/c_a3 = I2/c11  </td><td>  C2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  16  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P222  </td><td>  D2<SUP>1</SUP>  </td><td>  4  </td></tr>
<tr><td>  17  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P2221  </td><td>  D2<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  17  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>22  </td><td>  D2<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  17  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  P22<SUB>1</SUB>2  </td><td>  D2<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  18  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>2<SUB>1</SUB>2  </td><td>  D2<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  18  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P22<SUB>1</SUB>2<SUB>1</SUB>  </td><td>  D2<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  18  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>22<SUB>1</SUB>  </td><td>  D2<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  19  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>2<SUB>1</SUB>2<SUB>1</SUB>  </td><td>  D2<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  20  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  C222<SUB>1</SUB>  </td><td>  D2<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  20  </td><td>  2  </td><td>  1  </td><td>  5  </td><td>  A2<SUB>1</SUB>22  </td><td>  D2<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  20  </td><td>  3  </td><td>  1  </td><td>  6  </td><td>  B22<SUB>1</SUB>2  </td><td>  D2<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  21  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  C222  </td><td>  D2<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  21  </td><td>  2  </td><td>  1  </td><td>  5  </td><td>  A222  </td><td>  D2<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  21  </td><td>  3  </td><td>  1  </td><td>  6  </td><td>  B222  </td><td>  D2<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  22  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  F222  </td><td>  D2<SUP>7</SUP>  </td><td>  16  </td></tr>
<tr><td>  23  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I222  </td><td>  D2<SUP>8</SUP>  </td><td>  8  </td></tr>
<tr><td>  24  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I2<SUB>1</SUB>2<SUB>1</SUB>2<SUB>1</SUB>  </td><td>  D2<SUP>9</SUP>  </td><td>  8  </td></tr>
<tr><td>  25  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pmm2  </td><td>  C2v<SUP>1</SUP>  </td><td>  4  </td></tr>
<tr><td>  25  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2mm  </td><td>  C2v<SUP>1</SUP>  </td><td>  4  </td></tr>
<tr><td>  25  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pm2m  </td><td>  C2v<SUP>1</SUP>  </td><td>  4  </td></tr>
<tr><td>  26  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pmc2<SUB>1</SUB>  </td><td>  C2v<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  26  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pcm2<SUB>1</SUB>  </td><td>  C2v<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  26  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>ma  </td><td>  C2v<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  26  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>am  </td><td>  C2v<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  26  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pb2<SUB>1</SUB>m  </td><td>  C2v<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  26  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pm2<SUB>1</SUB>b  </td><td>  C2v<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  27  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pcc2  </td><td>  C2v<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  27  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2aa  </td><td>  C2v<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  27  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pb2b  </td><td>  C2v<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  28  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pma2  </td><td>  C2v<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  28  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pbm2  </td><td>  C2v<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  28  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2mb  </td><td>  C2v<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  28  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  P2cm  </td><td>  C2v<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  28  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pc2m  </td><td>  C2v<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  28  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pm2a  </td><td>  C2v<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  29  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pca2<SUB>1</SUB>  </td><td>  C2v<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  29  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pbc2<SUB>1</SUB>  </td><td>  C2v<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  29  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>ab  </td><td>  C2v<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  29  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>ca  </td><td>  C2v<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  29  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pc2<SUB>1</SUB>b  </td><td>  C2v<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  29  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pb2<SUB>1</SUB>a  </td><td>  C2v<SUP>5</SUP>  </td><td>  4  </td></tr>
<tr><td>  30  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pnc2  </td><td>  C2v<SUP>6</SUP>  </td><td>  4  </td></tr>
<tr><td>  30  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pcn2  </td><td>  C2v<SUP>6</SUP>  </td><td>  4  </td></tr>
<tr><td>  30  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2na  </td><td>  C2v<SUP>6</SUP>  </td><td>  4  </td></tr>
<tr><td>  30  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  P2an  </td><td>  C2v<SUP>6</SUP>  </td><td>  4  </td></tr>
<tr><td>  30  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pb2n  </td><td>  C2v<SUP>6</SUP>  </td><td>  4  </td></tr>
<tr><td>  30  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pn2b  </td><td>  C2v<SUP>6</SUP>  </td><td>  4  </td></tr>
<tr><td>  31  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pmn2<SUB>1</SUB>  </td><td>  C2v<SUP>7</SUP>  </td><td>  4  </td></tr>
<tr><td>  31  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pnm2<SUB>1</SUB>  </td><td>  C2v<SUP>7</SUP>  </td><td>  4  </td></tr>
<tr><td>  31  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>mn  </td><td>  C2v<SUP>7</SUP>  </td><td>  4  </td></tr>
<tr><td>  31  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>nm  </td><td>  C2v<SUP>7</SUP>  </td><td>  4  </td></tr>
<tr><td>  31  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pn2<SUB>1</SUB>m  </td><td>  C2v<SUP>7</SUP>  </td><td>  4  </td></tr>
<tr><td>  31  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pm2<SUB>1</SUB>n  </td><td>  C2v<SUP>7</SUP>  </td><td>  4  </td></tr>
<tr><td>  32  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pba2  </td><td>  C2v<SUP>8</SUP>  </td><td>  4  </td></tr>
<tr><td>  32  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2cb  </td><td>  C2v<SUP>8</SUP>  </td><td>  4  </td></tr>
<tr><td>  32  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pc2a  </td><td>  C2v<SUP>8</SUP>  </td><td>  4  </td></tr>
<tr><td>  33  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pna2<SUB>1</SUB>  </td><td>  C2v<SUP>9</SUP>  </td><td>  4  </td></tr>
<tr><td>  33  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pbn2<SUB>1</SUB>  </td><td>  C2v<SUP>9</SUP>  </td><td>  4  </td></tr>
<tr><td>  33  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>nb  </td><td>  C2v<SUP>9</SUP>  </td><td>  4  </td></tr>
<tr><td>  33  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>cn  </td><td>  C2v<SUP>9</SUP>  </td><td>  4  </td></tr>
<tr><td>  33  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pc2<SUB>1</SUB>n  </td><td>  C2v<SUP>9</SUP>  </td><td>  4  </td></tr>
<tr><td>  33  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pn2<SUB>1</SUB>a  </td><td>  C2v<SUP>9</SUP>  </td><td>  4  </td></tr>
<tr><td>  34  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pnn2  </td><td>  C2v<SUP>10</SUP>  </td><td>  4  </td></tr>
<tr><td>  34  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  P2nn  </td><td>  C2v<SUP>10</SUP>  </td><td>  4  </td></tr>
<tr><td>  34  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pn2n  </td><td>  C2v<SUP>10</SUP>  </td><td>  4  </td></tr>
<tr><td>  35  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  Cmm2  </td><td>  C2v<SUP>11</SUP>  </td><td>  8  </td></tr>
<tr><td>  35  </td><td>  2  </td><td>  1  </td><td>  5  </td><td>  A2mm  </td><td>  C2v<SUP>11</SUP>  </td><td>  8  </td></tr>
<tr><td>  35  </td><td>  3  </td><td>  1  </td><td>  6  </td><td>  Bm2m  </td><td>  C2v<SUP>11</SUP>  </td><td>  8  </td></tr>
<tr><td>  36  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  Cmc2<SUB>1</SUB>  </td><td>  C2v<SUP>12</SUP>  </td><td>  8  </td></tr>
<tr><td>  36  </td><td>  5  </td><td>  1  </td><td>  4  </td><td>  Ccm2<SUB>1</SUB>  </td><td>  C2v<SUP>12</SUP>  </td><td>  8  </td></tr>
<tr><td>  36  </td><td>  2  </td><td>  1  </td><td>  5  </td><td>  A2<SUB>1</SUB>ma  </td><td>  C2v<SUP>12</SUP>  </td><td>  8  </td></tr>
<tr><td>  36  </td><td>  6  </td><td>  1  </td><td>  5  </td><td>  A2<SUB>1</SUB>am  </td><td>  C2v<SUP>12</SUP>  </td><td>  8  </td></tr>
<tr><td>  36  </td><td>  3  </td><td>  1  </td><td>  6  </td><td>  Bb2<SUB>1</SUB>m  </td><td>  C2v<SUP>12</SUP>  </td><td>  8  </td></tr>
<tr><td>  36  </td><td>  4  </td><td>  1  </td><td>  6  </td><td>  Bm2<SUB>1</SUB>b  </td><td>  C2v<SUP>12</SUP>  </td><td>  8  </td></tr>
<tr><td>  37  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  Ccc2  </td><td>  C2v<SUP>13</SUP>  </td><td>  8  </td></tr>
<tr><td>  37  </td><td>  2  </td><td>  1  </td><td>  5  </td><td>  A2aa  </td><td>  C2v<SUP>13</SUP>  </td><td>  8  </td></tr>
<tr><td>  37  </td><td>  3  </td><td>  1  </td><td>  6  </td><td>  Bb2b  </td><td>  C2v<SUP>13</SUP>  </td><td>  8  </td></tr>
<tr><td>  38  </td><td>  1  </td><td>  1  </td><td>  5  </td><td>  Amm2  </td><td>  C2v<SUP>14</SUP>  </td><td>  8  </td></tr>
<tr><td>  38  </td><td>  5  </td><td>  1  </td><td>  6  </td><td>  Bmm2  </td><td>  C2v<SUP>14</SUP>  </td><td>  8  </td></tr>
<tr><td>  38  </td><td>  2  </td><td>  1  </td><td>  6  </td><td>  B2mm  </td><td>  C2v<SUP>14</SUP>  </td><td>  8  </td></tr>
<tr><td>  38  </td><td>  6  </td><td>  1  </td><td>  4  </td><td>  C2mm  </td><td>  C2v<SUP>14</SUP>  </td><td>  8  </td></tr>
<tr><td>  38  </td><td>  3  </td><td>  1  </td><td>  4  </td><td>  Cm2m  </td><td>  C2v<SUP>14</SUP>  </td><td>  8  </td></tr>
<tr><td>  38  </td><td>  4  </td><td>  1  </td><td>  5  </td><td>  Am2m  </td><td>  C2v<SUP>14</SUP>  </td><td>  8  </td></tr>
<tr><td>  39  </td><td>  1  </td><td>  1  </td><td>  5  </td><td>  Abm2  </td><td>  C2v<SUP>15</SUP>  </td><td>  8  </td></tr>
<tr><td>  39  </td><td>  5  </td><td>  1  </td><td>  6  </td><td>  Bma2  </td><td>  C2v<SUP>15</SUP>  </td><td>  8  </td></tr>
<tr><td>  39  </td><td>  2  </td><td>  1  </td><td>  6  </td><td>  B2cm  </td><td>  C2v<SUP>15</SUP>  </td><td>  8  </td></tr>
<tr><td>  39  </td><td>  6  </td><td>  1  </td><td>  4  </td><td>  C2mb  </td><td>  C2v<SUP>15</SUP>  </td><td>  8  </td></tr>
<tr><td>  39  </td><td>  3  </td><td>  1  </td><td>  4  </td><td>  Cm2a  </td><td>  C2v<SUP>15</SUP>  </td><td>  8  </td></tr>
<tr><td>  39  </td><td>  4  </td><td>  1  </td><td>  5  </td><td>  Ac2m  </td><td>  C2v<SUP>15</SUP>  </td><td>  8  </td></tr>
<tr><td>  40  </td><td>  1  </td><td>  1  </td><td>  5  </td><td>  Ama2  </td><td>  C2v<SUP>16</SUP>  </td><td>  8  </td></tr>
<tr><td>  40  </td><td>  5  </td><td>  1  </td><td>  6  </td><td>  Bbm2  </td><td>  C2v<SUP>16</SUP>  </td><td>  8  </td></tr>
<tr><td>  40  </td><td>  2  </td><td>  1  </td><td>  6  </td><td>  B2mb  </td><td>  C2v<SUP>16</SUP>  </td><td>  8  </td></tr>
<tr><td>  40  </td><td>  6  </td><td>  1  </td><td>  4  </td><td>  C2cm  </td><td>  C2v<SUP>16</SUP>  </td><td>  8  </td></tr>
<tr><td>  40  </td><td>  3  </td><td>  1  </td><td>  4  </td><td>  Cc2m  </td><td>  C2v<SUP>16</SUP>  </td><td>  8  </td></tr>
<tr><td>  40  </td><td>  4  </td><td>  1  </td><td>  5  </td><td>  Am2a  </td><td>  C2v<SUP>16</SUP>  </td><td>  8  </td></tr>
<tr><td>  41  </td><td>  1  </td><td>  1  </td><td>  5  </td><td>  Aba2  </td><td>  C2v<SUP>17</SUP>  </td><td>  8  </td></tr>
<tr><td>  41  </td><td>  5  </td><td>  1  </td><td>  6  </td><td>  Bba2  </td><td>  C2v<SUP>17</SUP>  </td><td>  8  </td></tr>
<tr><td>  41  </td><td>  2  </td><td>  1  </td><td>  6  </td><td>  B2cb  </td><td>  C2v<SUP>17</SUP>  </td><td>  8  </td></tr>
<tr><td>  41  </td><td>  6  </td><td>  1  </td><td>  4  </td><td>  C2cb  </td><td>  C2v<SUP>17</SUP>  </td><td>  8  </td></tr>
<tr><td>  41  </td><td>  3  </td><td>  1  </td><td>  4  </td><td>  Cc2a  </td><td>  C2v<SUP>17</SUP>  </td><td>  8  </td></tr>
<tr><td>  41  </td><td>  4  </td><td>  1  </td><td>  5  </td><td>  Ac2a  </td><td>  C2v<SUP>17</SUP>  </td><td>  8  </td></tr>
<tr><td>  42  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  Fmm2  </td><td>  C2v<SUP>18</SUP>  </td><td>  16  </td></tr>
<tr><td>  42  </td><td>  2  </td><td>  1  </td><td>  3  </td><td>  F2mm  </td><td>  C2v<SUP>18</SUP>  </td><td>  16  </td></tr>
<tr><td>  42  </td><td>  3  </td><td>  1  </td><td>  3  </td><td>  Fm2m  </td><td>  C2v<SUP>18</SUP>  </td><td>  16  </td></tr>
<tr><td>  43  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  Fdd2  </td><td>  C2v<SUP>19</SUP>  </td><td>  16  </td></tr>
<tr><td>  43  </td><td>  2  </td><td>  1  </td><td>  3  </td><td>  F2dd  </td><td>  C2v<SUP>19</SUP>  </td><td>  16  </td></tr>
<tr><td>  43  </td><td>  3  </td><td>  1  </td><td>  3  </td><td>  Fd2d  </td><td>  C2v<SUP>19</SUP>  </td><td>  16  </td></tr>
<tr><td>  44  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  Imm2  </td><td>  C2v<SUP>20</SUP>  </td><td>  8  </td></tr>
<tr><td>  44  </td><td>  2  </td><td>  1  </td><td>  2  </td><td>  I2mm  </td><td>  C2v<SUP>20</SUP>  </td><td>  8  </td></tr>
<tr><td>  44  </td><td>  3  </td><td>  1  </td><td>  2  </td><td>  Im2m  </td><td>  C2v<SUP>20</SUP>  </td><td>  8  </td></tr>
<tr><td>  45  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  Iba2  </td><td>  C2v<SUP>21</SUP>  </td><td>  8  </td></tr>
<tr><td>  45  </td><td>  2  </td><td>  1  </td><td>  2  </td><td>  I2cb  </td><td>  C2v<SUP>21</SUP>  </td><td>  8  </td></tr>
<tr><td>  45  </td><td>  3  </td><td>  1  </td><td>  2  </td><td>  Ic2a  </td><td>  C2v<SUP>21</SUP>  </td><td>  8  </td></tr>
<tr><td>  46  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  Ima2  </td><td>  C2v<SUP>22</SUP>  </td><td>  8  </td></tr>
<tr><td>  46  </td><td>  5  </td><td>  1  </td><td>  2  </td><td>  Ibm2  </td><td>  C2v<SUP>22</SUP>  </td><td>  8  </td></tr>
<tr><td>  46  </td><td>  2  </td><td>  1  </td><td>  2  </td><td>  I2mb  </td><td>  C2v<SUP>22</SUP>  </td><td>  8  </td></tr>
<tr><td>  46  </td><td>  6  </td><td>  1  </td><td>  2  </td><td>  I2cm  </td><td>  C2v<SUP>22</SUP>  </td><td>  8  </td></tr>
<tr><td>  46  </td><td>  3  </td><td>  1  </td><td>  2  </td><td>  Ic2m  </td><td>  C2v<SUP>22</SUP>  </td><td>  8  </td></tr>
<tr><td>  46  </td><td>  4  </td><td>  1  </td><td>  2  </td><td>  Im2a  </td><td>  C2v<SUP>22</SUP>  </td><td>  8  </td></tr>
<tr><td>  47  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pmmm  </td><td>  D2h<SUP>1</SUP>  </td><td>  8  </td></tr>
<tr><td>  48  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pnnn_1  </td><td>  D2h<SUP>2</SUP>  </td><td>  8  </td></tr>
<tr><td>  48  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  Pnnn_2  </td><td>  D2h<SUP>2</SUP>  </td><td>  8  </td></tr>
<tr><td>  49  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pccm  </td><td>  D2h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  49  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pmaa  </td><td>  D2h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  49  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pbmb  </td><td>  D2h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  50  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pban_1  </td><td>  D2h<SUP>4</SUP>  </td><td>  8  </td></tr>
<tr><td>  50  </td><td>  5  </td><td>  2  </td><td>  1  </td><td>  Pban_2  </td><td>  D2h<SUP>4</SUP>  </td><td>  8  </td></tr>
<tr><td>  50  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pncb_1  </td><td>  D2h<SUP>4</SUP>  </td><td>  8  </td></tr>
<tr><td>  50  </td><td>  6  </td><td>  2  </td><td>  1  </td><td>  Pncb_2  </td><td>  D2h<SUP>4</SUP>  </td><td>  8  </td></tr>
<tr><td>  50  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pcna_1  </td><td>  D2h<SUP>4</SUP>  </td><td>  8  </td></tr>
<tr><td>  50  </td><td>  4  </td><td>  2  </td><td>  1  </td><td>  Pcna_2  </td><td>  D2h<SUP>4</SUP>  </td><td>  8  </td></tr>
<tr><td>  51  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pmma  </td><td>  D2h<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  51  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pmmb  </td><td>  D2h<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  51  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pbmm  </td><td>  D2h<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  51  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  Pcmm  </td><td>  D2h<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  51  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pmcm  </td><td>  D2h<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  51  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pmam  </td><td>  D2h<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  52  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pnna  </td><td>  D2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  52  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pnnb  </td><td>  D2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  52  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pbnn  </td><td>  D2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  52  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  Pcnn  </td><td>  D2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  52  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pncn  </td><td>  D2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  52  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pnan  </td><td>  D2h<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  53  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pmna  </td><td>  D2h<SUP>7</SUP>  </td><td>  8  </td></tr>
<tr><td>  53  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pnmb  </td><td>  D2h<SUP>7</SUP>  </td><td>  8  </td></tr>
<tr><td>  53  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pbmn  </td><td>  D2h<SUP>7</SUP>  </td><td>  8  </td></tr>
<tr><td>  53  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  Pcnm  </td><td>  D2h<SUP>7</SUP>  </td><td>  8  </td></tr>
<tr><td>  53  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pncm  </td><td>  D2h<SUP>7</SUP>  </td><td>  8  </td></tr>
<tr><td>  53  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pman  </td><td>  D2h<SUP>7</SUP>  </td><td>  8  </td></tr>
<tr><td>  54  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pcca  </td><td>  D2h<SUP>8</SUP>  </td><td>  8  </td></tr>
<tr><td>  54  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pccb  </td><td>  D2h<SUP>8</SUP>  </td><td>  8  </td></tr>
<tr><td>  54  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pbaa  </td><td>  D2h<SUP>8</SUP>  </td><td>  8  </td></tr>
<tr><td>  54  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  Pcaa  </td><td>  D2h<SUP>8</SUP>  </td><td>  8  </td></tr>
<tr><td>  54  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pbcb  </td><td>  D2h<SUP>8</SUP>  </td><td>  8  </td></tr>
<tr><td>  54  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pbab  </td><td>  D2h<SUP>8</SUP>  </td><td>  8  </td></tr>
<tr><td>  55  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pbam  </td><td>  D2h<SUP>9</SUP>  </td><td>  8  </td></tr>
<tr><td>  55  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pmcb  </td><td>  D2h<SUP>9</SUP>  </td><td>  8  </td></tr>
<tr><td>  55  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pcma  </td><td>  D2h<SUP>9</SUP>  </td><td>  8  </td></tr>
<tr><td>  56  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pccn  </td><td>  D2h<SUP>10</SUP>  </td><td>  8  </td></tr>
<tr><td>  56  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pnaa  </td><td>  D2h<SUP>10</SUP>  </td><td>  8  </td></tr>
<tr><td>  56  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pbnb  </td><td>  D2h<SUP>10</SUP>  </td><td>  8  </td></tr>
<tr><td>  57  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pbcm  </td><td>  D2h<SUP>11</SUP>  </td><td>  8  </td></tr>
<tr><td>  57  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pcam  </td><td>  D2h<SUP>11</SUP>  </td><td>  8  </td></tr>
<tr><td>  57  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pmca  </td><td>  D2h<SUP>11</SUP>  </td><td>  8  </td></tr>
<tr><td>  57  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  Pmab  </td><td>  D2h<SUP>11</SUP>  </td><td>  8  </td></tr>
<tr><td>  57  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pbma  </td><td>  D2h<SUP>11</SUP>  </td><td>  8  </td></tr>
<tr><td>  57  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pcmb  </td><td>  D2h<SUP>11</SUP>  </td><td>  8  </td></tr>
<tr><td>  58  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pnnm  </td><td>  D2h<SUP>12</SUP>  </td><td>  8  </td></tr>
<tr><td>  58  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pmnn  </td><td>  D2h<SUP>12</SUP>  </td><td>  8  </td></tr>
<tr><td>  58  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pnmn  </td><td>  D2h<SUP>12</SUP>  </td><td>  8  </td></tr>
<tr><td>  59  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pmmn_1  </td><td>  D2h<SUP>13</SUP>  </td><td>  8  </td></tr>
<tr><td>  59  </td><td>  5  </td><td>  2  </td><td>  1  </td><td>  Pmmn_2  </td><td>  D2h<SUP>13</SUP>  </td><td>  8  </td></tr>
<tr><td>  59  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pnmm_1  </td><td>  D2h<SUP>13</SUP>  </td><td>  8  </td></tr>
<tr><td>  59  </td><td>  6  </td><td>  2  </td><td>  1  </td><td>  Pnmm_2  </td><td>  D2h<SUP>13</SUP>  </td><td>  8  </td></tr>
<tr><td>  59  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pmnm_1  </td><td>  D2h<SUP>13</SUP>  </td><td>  8  </td></tr>
<tr><td>  59  </td><td>  4  </td><td>  2  </td><td>  1  </td><td>  Pmnm_2  </td><td>  D2h<SUP>13</SUP>  </td><td>  8  </td></tr>
<tr><td>  60  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pbcn  </td><td>  D2h<SUP>14</SUP>  </td><td>  8  </td></tr>
<tr><td>  60  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pcan  </td><td>  D2h<SUP>14</SUP>  </td><td>  8  </td></tr>
<tr><td>  60  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pnca  </td><td>  D2h<SUP>14</SUP>  </td><td>  8  </td></tr>
<tr><td>  60  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  Pnab  </td><td>  D2h<SUP>14</SUP>  </td><td>  8  </td></tr>
<tr><td>  60  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pbna  </td><td>  D2h<SUP>14</SUP>  </td><td>  8  </td></tr>
<tr><td>  60  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pcnb  </td><td>  D2h<SUP>14</SUP>  </td><td>  8  </td></tr>
<tr><td>  61  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pbca  </td><td>  D2h<SUP>15</SUP>  </td><td>  8  </td></tr>
<tr><td>  61  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pcab  </td><td>  D2h<SUP>15</SUP>  </td><td>  8  </td></tr>
<tr><td>  62  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pnma  </td><td>  D2h<SUP>16</SUP>  </td><td>  8  </td></tr>
<tr><td>  62  </td><td>  5  </td><td>  1  </td><td>  1  </td><td>  Pmnb  </td><td>  D2h<SUP>16</SUP>  </td><td>  8  </td></tr>
<tr><td>  62  </td><td>  2  </td><td>  1  </td><td>  1  </td><td>  Pbnm  </td><td>  D2h<SUP>16</SUP>  </td><td>  8  </td></tr>
<tr><td>  62  </td><td>  6  </td><td>  1  </td><td>  1  </td><td>  Pcmn  </td><td>  D2h<SUP>16</SUP>  </td><td>  8  </td></tr>
<tr><td>  62  </td><td>  3  </td><td>  1  </td><td>  1  </td><td>  Pmcn  </td><td>  D2h<SUP>16</SUP>  </td><td>  8  </td></tr>
<tr><td>  62  </td><td>  4  </td><td>  1  </td><td>  1  </td><td>  Pnam  </td><td>  D2h<SUP>16</SUP>  </td><td>  8  </td></tr>
<tr><td>  63  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  Cmcm  </td><td>  D2h<SUP>17</SUP>  </td><td>  16  </td></tr>
<tr><td>  63  </td><td>  5  </td><td>  1  </td><td>  4  </td><td>  Ccmm  </td><td>  D2h<SUP>17</SUP>  </td><td>  16  </td></tr>
<tr><td>  63  </td><td>  2  </td><td>  1  </td><td>  5  </td><td>  Amma  </td><td>  D2h<SUP>17</SUP>  </td><td>  16  </td></tr>
<tr><td>  63  </td><td>  6  </td><td>  1  </td><td>  5  </td><td>  Amam  </td><td>  D2h<SUP>17</SUP>  </td><td>  16  </td></tr>
<tr><td>  63  </td><td>  3  </td><td>  1  </td><td>  6  </td><td>  Bbmm  </td><td>  D2h<SUP>17</SUP>  </td><td>  16  </td></tr>
<tr><td>  63  </td><td>  4  </td><td>  1  </td><td>  6  </td><td>  Bmmb  </td><td>  D2h<SUP>17</SUP>  </td><td>  16  </td></tr>
<tr><td>  64  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  Cmca  </td><td>  D2h<SUP>18</SUP>  </td><td>  16  </td></tr>
<tr><td>  64  </td><td>  5  </td><td>  1  </td><td>  4  </td><td>  Ccmb  </td><td>  D2h<SUP>18</SUP>  </td><td>  16  </td></tr>
<tr><td>  64  </td><td>  2  </td><td>  1  </td><td>  5  </td><td>  Abma  </td><td>  D2h<SUP>18</SUP>  </td><td>  16  </td></tr>
<tr><td>  64  </td><td>  6  </td><td>  1  </td><td>  5  </td><td>  Acam  </td><td>  D2h<SUP>18</SUP>  </td><td>  16  </td></tr>
<tr><td>  64  </td><td>  3  </td><td>  1  </td><td>  6  </td><td>  Bbcm  </td><td>  D2h<SUP>18</SUP>  </td><td>  16  </td></tr>
<tr><td>  64  </td><td>  4  </td><td>  1  </td><td>  6  </td><td>  Bmab  </td><td>  D2h<SUP>18</SUP>  </td><td>  16  </td></tr>
<tr><td>  65  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  Cmmm  </td><td>  D2h<SUP>19</SUP>  </td><td>  16  </td></tr>
<tr><td>  65  </td><td>  2  </td><td>  1  </td><td>  5  </td><td>  Ammm  </td><td>  D2h<SUP>19</SUP>  </td><td>  16  </td></tr>
<tr><td>  65  </td><td>  3  </td><td>  1  </td><td>  6  </td><td>  Bmmm  </td><td>  D2h<SUP>19</SUP>  </td><td>  16  </td></tr>
<tr><td>  66  </td><td>  1  </td><td>  1  </td><td>  6  </td><td>  Cccm  </td><td>  D2h<SUP>20</SUP>  </td><td>  16  </td></tr>
<tr><td>  66  </td><td>  2  </td><td>  1  </td><td>  5  </td><td>  Amaa  </td><td>  D2h<SUP>20</SUP>  </td><td>  16  </td></tr>
<tr><td>  66  </td><td>  3  </td><td>  1  </td><td>  6  </td><td>  Bbmb  </td><td>  D2h<SUP>20</SUP>  </td><td>  16  </td></tr>
<tr><td>  67  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  Cmma  </td><td>  D2h<SUP>21</SUP>  </td><td>  16  </td></tr>
<tr><td>  67  </td><td>  5  </td><td>  1  </td><td>  4  </td><td>  Cmmb  </td><td>  D2h<SUP>21</SUP>  </td><td>  16  </td></tr>
<tr><td>  67  </td><td>  2  </td><td>  1  </td><td>  5  </td><td>  Abmm  </td><td>  D2h<SUP>21</SUP>  </td><td>  16  </td></tr>
<tr><td>  67  </td><td>  6  </td><td>  1  </td><td>  5  </td><td>  Acmm  </td><td>  D2h<SUP>21</SUP>  </td><td>  16  </td></tr>
<tr><td>  67  </td><td>  3  </td><td>  1  </td><td>  6  </td><td>  Bmcm  </td><td>  D2h<SUP>21</SUP>  </td><td>  16  </td></tr>
<tr><td>  67  </td><td>  4  </td><td>  1  </td><td>  6  </td><td>  Bmam  </td><td>  D2h<SUP>21</SUP>  </td><td>  16  </td></tr>
<tr><td>  68  </td><td>  1  </td><td>  1  </td><td>  4  </td><td>  Ccca_1  </td><td>  D2h<SUP>22</SUP>  </td><td>  16  </td></tr>
<tr><td>  68  </td><td>  1  </td><td>  2  </td><td>  4  </td><td>  Ccca_2  </td><td>  D2h<SUP>22</SUP>  </td><td>  16  </td></tr>
<tr><td>  68  </td><td>  5  </td><td>  1  </td><td>  4  </td><td>  Cccb_1  </td><td>  D2h<SUP>22</SUP>  </td><td>  16  </td></tr>
<tr><td>  68  </td><td>  5  </td><td>  2  </td><td>  4  </td><td>  Cccb_2  </td><td>  D2h<SUP>22</SUP>  </td><td>  16  </td></tr>
<tr><td>  68  </td><td>  2  </td><td>  1  </td><td>  5  </td><td>  Abaa_1  </td><td>  D2h<SUP>22</SUP>  </td><td>  16  </td></tr>
<tr><td>  68  </td><td>  2  </td><td>  2  </td><td>  5  </td><td>  Abaa_2  </td><td>  D2h<SUP>22</SUP>  </td><td>  16  </td></tr>
<tr><td>  68  </td><td>  6  </td><td>  1  </td><td>  5  </td><td>  Acaa_1  </td><td>  D2h<SUP>22</SUP>  </td><td>  16  </td></tr>
<tr><td>  68  </td><td>  6  </td><td>  2  </td><td>  5  </td><td>  Acaa_2  </td><td>  D2h<SUP>22</SUP>  </td><td>  16  </td></tr>
<tr><td>  68  </td><td>  3  </td><td>  1  </td><td>  6  </td><td>  Bbcb_1  </td><td>  D2h<SUP>22</SUP>  </td><td>  16  </td></tr>
<tr><td>  68  </td><td>  3  </td><td>  2  </td><td>  6  </td><td>  Bbcb_2  </td><td>  D2h<SUP>22</SUP>  </td><td>  16  </td></tr>
<tr><td>  68  </td><td>  4  </td><td>  1  </td><td>  6  </td><td>  Bbab_1  </td><td>  D2h<SUP>22</SUP>  </td><td>  16  </td></tr>
<tr><td>  68  </td><td>  4  </td><td>  2  </td><td>  6  </td><td>  Bbab_2  </td><td>  D2h<SUP>22</SUP>  </td><td>  16  </td></tr>
<tr><td>  69  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  Fmmm  </td><td>  D2h<SUP>23</SUP>  </td><td>  32  </td></tr>
<tr><td>  70  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  Fddd_1  </td><td>  D2h<SUP>24</SUP>  </td><td>  32  </td></tr>
<tr><td>  70  </td><td>  1  </td><td>  2  </td><td>  3  </td><td>  Fddd_2  </td><td>  D2h<SUP>24</SUP>  </td><td>  32  </td></tr>
<tr><td>  71  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  Immm  </td><td>  D2h<SUP>25</SUP>  </td><td>  16  </td></tr>
<tr><td>  72  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  Ibam  </td><td>  D2h<SUP>26</SUP>  </td><td>  16  </td></tr>
<tr><td>  72  </td><td>  2  </td><td>  1  </td><td>  2  </td><td>  Imcb  </td><td>  D2h<SUP>26</SUP>  </td><td>  16  </td></tr>
<tr><td>  72  </td><td>  3  </td><td>  1  </td><td>  2  </td><td>  Icma  </td><td>  D2h<SUP>26</SUP>  </td><td>  16  </td></tr>
<tr><td>  73  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  Ibca  </td><td>  D2h<SUP>27</SUP>  </td><td>  16  </td></tr>
<tr><td>  73  </td><td>  2  </td><td>  1  </td><td>  2  </td><td>  Icab  </td><td>  D2h<SUP>27</SUP>  </td><td>  16  </td></tr>
<tr><td>  74  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  Imma  </td><td>  D2h<SUP>28</SUP>  </td><td>  16  </td></tr>
<tr><td>  74  </td><td>  5  </td><td>  1  </td><td>  2  </td><td>  Immb  </td><td>  D2h<SUP>28</SUP>  </td><td>  16  </td></tr>
<tr><td>  74  </td><td>  2  </td><td>  1  </td><td>  2  </td><td>  Ibmm  </td><td>  D2h<SUP>28</SUP>  </td><td>  16  </td></tr>
<tr><td>  74  </td><td>  6  </td><td>  1  </td><td>  2  </td><td>  Icmm  </td><td>  D2h<SUP>28</SUP>  </td><td>  16  </td></tr>
<tr><td>  74  </td><td>  3  </td><td>  1  </td><td>  2  </td><td>  Imcm  </td><td>  D2h<SUP>28</SUP>  </td><td>  16  </td></tr>
<tr><td>  74  </td><td>  4  </td><td>  1  </td><td>  2  </td><td>  Imam  </td><td>  D2h<SUP>28</SUP>  </td><td>  16  </td></tr>
<tr><td>  75  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4  </td><td>  C4<SUP>1</SUP>  </td><td>  4  </td></tr>
<tr><td>  76  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>1</SUB>  </td><td>  C4<SUP>2</SUP>  </td><td>  4  </td></tr>
<tr><td>  77  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>  </td><td>  C4<SUP>3</SUP>  </td><td>  4  </td></tr>
<tr><td>  78  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>3</SUB>  </td><td>  C4<SUP>4</SUP>  </td><td>  4  </td></tr>
<tr><td>  79  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4  </td><td>  C4<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  80  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4<SUB>1</SUB>  </td><td>  C4<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  81  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-4  </td><td>  S4<SUP>1</SUP>  </td><td>  4  </td></tr>
<tr><td>  82  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I-4  </td><td>  S4<SUP>2</SUP>  </td><td>  8  </td></tr>
<tr><td>  83  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4/m  </td><td>  C4h<SUP>1</SUP>  </td><td>  8  </td></tr>
<tr><td>  84  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>/m  </td><td>  C4h<SUP>2</SUP>  </td><td>  8  </td></tr>
<tr><td>  85  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4/n_1  </td><td>  C4h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  85  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  P4/n_2  </td><td>  C4h<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  86  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>/n_1  </td><td>  C4h<SUP>4</SUP>  </td><td>  8  </td></tr>
<tr><td>  86  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  P4<SUB>2</SUB>/n_2  </td><td>  C4h<SUP>4</SUP>  </td><td>  8  </td></tr>
<tr><td>  87  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4/m  </td><td>  C4h<SUP>5</SUP>  </td><td>  16  </td></tr>
<tr><td>  88  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4<SUB>1</SUB>/a_1  </td><td>  C4h<SUP>6</SUP>  </td><td>  16  </td></tr>
<tr><td>  88  </td><td>  1  </td><td>  2  </td><td>  2  </td><td>  I4<SUB>1</SUB>/a_2  </td><td>  C4h<SUP>6</SUP>  </td><td>  16  </td></tr>
<tr><td>  89  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>2  </td><td>  D4<SUP>1</SUP>  </td><td>  8  </td></tr>
<tr><td>  90  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P42<SUB>1</SUB>2  </td><td>  D4<SUP>2</SUP>  </td><td>  8  </td></tr>
<tr><td>  91  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>1</SUB>22  </td><td>  D4<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  92  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>1</SUB>2<SUB>1</SUB>2  </td><td>  D4<SUP>4</SUP>  </td><td>  8  </td></tr>
<tr><td>  93  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4222  </td><td>  D4<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  94  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P422<SUB>1</SUB>2  </td><td>  D4<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  95  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>3</SUB>22  </td><td>  D4<SUP>7</SUP>  </td><td>  8  </td></tr>
<tr><td>  96  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>3</SUB>2<SUB>1</SUB>2  </td><td>  D4<SUP>8</SUP>  </td><td>  8  </td></tr>
<tr><td>  97  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I422  </td><td>  D4<SUP>9</SUP>  </td><td>  16  </td></tr>
<tr><td>  98  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4<SUB>1</SUB>22  </td><td>  D4<SUP>10</SUP>  </td><td>  16  </td></tr>
<tr><td>  99  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4mm  </td><td>  C4v<SUP>1</SUP>  </td><td>  8  </td></tr>
<tr><td>  100  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4bm  </td><td>  C4v<SUP>2</SUP>  </td><td>  8  </td></tr>
<tr><td>  101  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P42cm  </td><td>  C4v<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  102  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P42nm  </td><td>  C4v<SUP>4</SUP>  </td><td>  8  </td></tr>
<tr><td>  103  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4cc  </td><td>  C4v<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  104  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4nc  </td><td>  C4v<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  105  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>mc  </td><td>  C4v<SUP>7</SUP>  </td><td>  8  </td></tr>
<tr><td>  106  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>bc  </td><td>  C4v<SUP>8</SUP>  </td><td>  8  </td></tr>
<tr><td>  107  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4mm  </td><td>  C4v<SUP>9</SUP>  </td><td>  16  </td></tr>
<tr><td>  108  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4cm  </td><td>  C4v<SUP>10</SUP>  </td><td>  16  </td></tr>
<tr><td>  109  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4<SUB>1</SUB>md  </td><td>  C4v<SUP>11</SUP>  </td><td>  16  </td></tr>
<tr><td>  110  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4<SUB>1</SUB>cd  </td><td>  C4v<SUP>12</SUP>  </td><td>  16  </td></tr>
<tr><td>  111  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-42m  </td><td>  D2d<SUP>1</SUP>  </td><td>  8  </td></tr>
<tr><td>  112  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-42c  </td><td>  D2d<SUP>2</SUP>  </td><td>  8  </td></tr>
<tr><td>  113  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-42<SUB>1</SUB>m  </td><td>  D2d<SUP>3</SUP>  </td><td>  8  </td></tr>
<tr><td>  114  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-42<SUB>1</SUB>c  </td><td>  D2d<SUP>4</SUP>  </td><td>  8  </td></tr>
<tr><td>  115  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-4m2  </td><td>  D2d<SUP>5</SUP>  </td><td>  8  </td></tr>
<tr><td>  116  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-4c2  </td><td>  D2d<SUP>6</SUP>  </td><td>  8  </td></tr>
<tr><td>  117  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-4b2  </td><td>  D2d<SUP>7</SUP>  </td><td>  8  </td></tr>
<tr><td>  118  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-4n2  </td><td>  D2d<SUP>8</SUP>  </td><td>  8  </td></tr>
<tr><td>  119  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I-4m2  </td><td>  D2d<SUP>9</SUP>  </td><td>  16  </td></tr>
<tr><td>  120  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I-4c2  </td><td>  D2d<SUP>10</SUP>  </td><td>  16  </td></tr>
<tr><td>  121  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I-42m  </td><td>  D2d<SUP>11</SUP>  </td><td>  16  </td></tr>
<tr><td>  122  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I-42d  </td><td>  D2d<SUP>12</SUP>  </td><td>  16  </td></tr>
<tr><td>  123  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4/mmm  </td><td>  D4h<SUP>1</SUP>  </td><td>  16  </td></tr>
<tr><td>  124  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4/mcc  </td><td>  D4h<SUP>2</SUP>  </td><td>  16  </td></tr>
<tr><td>  125  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4/nbm_1  </td><td>  D4h<SUP>3</SUP>  </td><td>  16  </td></tr>
<tr><td>  125  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  P4/nbm_2  </td><td>  D4h<SUP>3</SUP>  </td><td>  16  </td></tr>
<tr><td>  126  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4/nnc_1  </td><td>  D4h<SUP>4</SUP>  </td><td>  16  </td></tr>
<tr><td>  126  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  P4/nnc_2  </td><td>  D4h<SUP>4</SUP>  </td><td>  16  </td></tr>
<tr><td>  127  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4/mbm  </td><td>  D4h<SUP>5</SUP>  </td><td>  16  </td></tr>
<tr><td>  128  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4/mnc  </td><td>  D4h<SUP>6</SUP>  </td><td>  16  </td></tr>
<tr><td>  129  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4/nmm_1  </td><td>  D4h<SUP>7</SUP>  </td><td>  16  </td></tr>
<tr><td>  129  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  P4/nmm_2  </td><td>  D4h<SUP>7</SUP>  </td><td>  16  </td></tr>
<tr><td>  130  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4/ncc_1  </td><td>  D4h<SUP>8</SUP>  </td><td>  16  </td></tr>
<tr><td>  130  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  P4/ncc_2  </td><td>  D4h<SUP>8</SUP>  </td><td>  16  </td></tr>
<tr><td>  131  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>/mmc  </td><td>  D4h<SUP>9</SUP>  </td><td>  16  </td></tr>
<tr><td>  132  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>/mcm  </td><td>  D4h<SUP>10</SUP>  </td><td>  16  </td></tr>
<tr><td>  133  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>/nbc_1  </td><td>  D4h<SUP>11</SUP>  </td><td>  16  </td></tr>
<tr><td>  133  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  P4<SUB>2</SUB>/nbc_2  </td><td>  D4h<SUP>11</SUP>  </td><td>  16  </td></tr>
<tr><td>  134  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>/nnm_1  </td><td>  D4h<SUP>12</SUP>  </td><td>  16  </td></tr>
<tr><td>  134  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  P4<SUB>2</SUB>/nnm_2  </td><td>  D4h<SUP>12</SUP>  </td><td>  16  </td></tr>
<tr><td>  135  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>/mbc  </td><td>  D4h<SUP>13</SUP>  </td><td>  16  </td></tr>
<tr><td>  136  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>/mnm  </td><td>  D4h<SUP>14</SUP>  </td><td>  16  </td></tr>
<tr><td>  137  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>/nmc_1  </td><td>  D4h<SUP>15</SUP>  </td><td>  16  </td></tr>
<tr><td>  137  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  P4<SUB>2</SUB>/nmc_2  </td><td>  D4h<SUP>15</SUP>  </td><td>  16  </td></tr>
<tr><td>  138  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>/ncm_1  </td><td>  D4h<SUP>16</SUP>  </td><td>  16  </td></tr>
<tr><td>  138  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  P4<SUB>2</SUB>/ncm_2  </td><td>  D4h<SUP>16</SUP>  </td><td>  16  </td></tr>
<tr><td>  139  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4/mmm  </td><td>  D4h<SUP>17</SUP>  </td><td>  32  </td></tr>
<tr><td>  140  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4/mcm  </td><td>  D4h<SUP>18</SUP>  </td><td>  32  </td></tr>
<tr><td>  141  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4<SUB>1</SUB>/amd_1  </td><td>  D4h<SUP>19</SUP>  </td><td>  32  </td></tr>
<tr><td>  141  </td><td>  1  </td><td>  2  </td><td>  2  </td><td>  I4<SUB>1</SUB>/amd_2  </td><td>  D4h<SUP>19</SUP>  </td><td>  32  </td></tr>
<tr><td>  142  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4<SUB>1</SUB>/acd_1  </td><td>  D4h<SUP>20</SUP>  </td><td>  32  </td></tr>
<tr><td>  142  </td><td>  1  </td><td>  2  </td><td>  2  </td><td>  I4<SUB>1</SUB>/acd_2  </td><td>  D4h<SUP>20</SUP>  </td><td>  32  </td></tr>
<tr><td>  143  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P3  </td><td>  C3<SUP>1</SUP>  </td><td>  3  </td></tr>
<tr><td>  144  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P3<SUB>1</SUB>  </td><td>  C3<SUP>2</SUP>  </td><td>  3  </td></tr>
<tr><td>  145  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P3<SUB>2</SUB>  </td><td>  C3<SUP>3</SUP>  </td><td>  3  </td></tr>
<tr><td>  146  </td><td>  1  </td><td>  1  </td><td>  7  </td><td>  R3_H  </td><td>  C3<SUP>4</SUP>  </td><td>  9  </td></tr>
<tr><td>  146  </td><td>  2  </td><td>  1  </td><td>  7  </td><td>  R3_R  </td><td>  C3<SUP>4</SUP>  </td><td>  3  </td></tr>
<tr><td>  147  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-3  </td><td>  C3i<SUP>1</SUP>  </td><td>  6  </td></tr>
<tr><td>  148  </td><td>  1  </td><td>  1  </td><td>  7  </td><td>  R-3_H  </td><td>  C3i<SUP>2</SUP>  </td><td> 18  </td></tr>
<tr><td>  148  </td><td>  2  </td><td>  1  </td><td>  7  </td><td>  R-3_R  </td><td>  C3i<SUP>2</SUP>  </td><td>  6  </td></tr>
<tr><td>  149  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P312  </td><td>  D3<SUP>1</SUP>  </td><td>  6  </td></tr>
<tr><td>  150  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P321  </td><td>  D3<SUP>2</SUP>  </td><td>  6  </td></tr>
<tr><td>  151  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P3<SUB>1</SUB>12  </td><td>  D3<SUP>3</SUP>  </td><td>  6  </td></tr>
<tr><td>  152  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P3<SUB>1</SUB>21  </td><td>  D3<SUP>4</SUP>  </td><td>  6  </td></tr>
<tr><td>  153  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P3<SUB>2</SUB>12  </td><td>  D3<SUP>5</SUP>  </td><td>  6  </td></tr>
<tr><td>  154  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P3<SUB>2</SUB>21  </td><td>  D3<SUP>6</SUP>  </td><td>  6  </td></tr>
<tr><td>  155  </td><td>  1  </td><td>  1  </td><td>  7  </td><td>  R32_H  </td><td>  D3<SUP>7</SUP>  </td><td>  18  </td></tr>
<tr><td>  155  </td><td>  2  </td><td>  1  </td><td>  7  </td><td>  R32_R  </td><td>  D3<SUP>7</SUP>  </td><td>  6  </td></tr>
<tr><td>  156  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P3m1  </td><td>  C3v<SUP>1</SUP>  </td><td>  6  </td></tr>
<tr><td>  157  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P31m  </td><td>  C3v<SUP>2</SUP>  </td><td>  6  </td></tr>
<tr><td>  158  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P3c1  </td><td>  C3v<SUP>3</SUP>  </td><td>  6  </td></tr>
<tr><td>  159  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P31c  </td><td>  C3v<SUP>4</SUP>  </td><td>  6  </td></tr>
<tr><td>  160  </td><td>  1  </td><td>  1  </td><td>  7  </td><td>  R3m_H  </td><td>  C3v<SUP>5</SUP>  </td><td>  18  </td></tr>
<tr><td>  160  </td><td>  2  </td><td>  1  </td><td>  7  </td><td>  R3m_R  </td><td>  C3v<SUP>5</SUP>  </td><td>  6  </td></tr>
<tr><td>  161  </td><td>  1  </td><td>  1  </td><td>  7  </td><td>  R3c_H  </td><td>  C3v<SUP>6</SUP>  </td><td>  18  </td></tr>
<tr><td>  161  </td><td>  2  </td><td>  1  </td><td>  7  </td><td>  R3c_R  </td><td>  C3v<SUP>6</SUP>  </td><td>  6  </td></tr>
<tr><td>  162  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-31m  </td><td>  D3d<SUP>1</SUP>  </td><td>  12  </td></tr>
<tr><td>  163  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-31c  </td><td>  D3d<SUP>2</SUP>  </td><td>  12  </td></tr>
<tr><td>  164  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-3m1  </td><td>  D3d<SUP>3</SUP>  </td><td>  12  </td></tr>
<tr><td>  165  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-3c1  </td><td>  D3d<SUP>4</SUP>  </td><td>  12  </td></tr>
<tr><td>  166  </td><td>  1  </td><td>  1  </td><td>  7  </td><td>  R-3m_H  </td><td>  D3d<SUP>5</SUP>  </td><td>  36  </td></tr>
<tr><td>  166  </td><td>  2  </td><td>  1  </td><td>  7  </td><td>  R-3m_R  </td><td>  D3d<SUP>5</SUP>  </td><td>  12  </td></tr>
<tr><td>  167  </td><td>  1  </td><td>  1  </td><td>  7  </td><td>  R-3c_H  </td><td>  D3d<SUP>6</SUP>  </td><td>  36  </td></tr>
<tr><td>  167  </td><td>  2  </td><td>  1  </td><td>  7  </td><td>  R-3c_R  </td><td>  D3d<SUP>6</SUP>  </td><td>  12  </td></tr>
<tr><td>  168  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6  </td><td>  C6<SUP>1</SUP>  </td><td>  6  </td></tr>
<tr><td>  169  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>1</SUB>  </td><td>  C6<SUP>2</SUP>  </td><td>  6  </td></tr>
<tr><td>  170  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>5</SUB>  </td><td>  C6<SUP>3</SUP>  </td><td>  6  </td></tr>
<tr><td>  171  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>2</SUB>  </td><td>  C6<SUP>4</SUP>  </td><td>  6  </td></tr>
<tr><td>  172  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>4</SUB>  </td><td>  C6<SUP>5</SUP>  </td><td>  6  </td></tr>
<tr><td>  173  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>3</SUB>  </td><td>  C6<SUP>6</SUP>  </td><td>  6  </td></tr>
<tr><td>  174  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-6  </td><td>  C3h<SUP>1</SUP>  </td><td>  6  </td></tr>
<tr><td>  175  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6/m  </td><td>  C6h<SUP>1</SUP>  </td><td>  12  </td></tr>
<tr><td>  176  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>3</SUB>/m  </td><td>  C6h<SUP>2</SUP>  </td><td>  12  </td></tr>
<tr><td>  177  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P622  </td><td>  D6<SUP>1</SUP>  </td><td>  12  </td></tr>
<tr><td>  178  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>1</SUB>22  </td><td>  D6<SUP>2</SUP>  </td><td>  12  </td></tr>
<tr><td>  179  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>5</SUB>22  </td><td>  D6<SUP>3</SUP>  </td><td>  12  </td></tr>
<tr><td>  180  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>2</SUB>22  </td><td>  D6<SUP>4</SUP>  </td><td>  12  </td></tr>
<tr><td>  181  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>4</SUB>22  </td><td>  D6<SUP>5</SUP>  </td><td>  12  </td></tr>
<tr><td>  182  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>3</SUB>22  </td><td>  D6<SUP>6</SUP>  </td><td>  12  </td></tr>
<tr><td>  183  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6mm  </td><td>  C6v<SUP>1</SUP>  </td><td>  12  </td></tr>
<tr><td>  184  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6cc  </td><td>  C6v<SUP>2</SUP>  </td><td>  12  </td></tr>
<tr><td>  185  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>3</SUB>cm  </td><td>  C6v<SUP>3</SUP>  </td><td>  12  </td></tr>
<tr><td>  186  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>3</SUB>mc  </td><td>  C6v<SUP>4</SUP>  </td><td>  12  </td></tr>
<tr><td>  187  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-6m2  </td><td>  D3h<SUP>1</SUP>  </td><td>  12  </td></tr>
<tr><td>  188  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-6c2  </td><td>  D3h<SUP>2</SUP>  </td><td>  12  </td></tr>
<tr><td>  189  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-62m  </td><td>  D3h<SUP>3</SUP>  </td><td>  12  </td></tr>
<tr><td>  190  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-62c  </td><td>  D3h<SUP>4</SUP>  </td><td>  12  </td></tr>
<tr><td>  191  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6/mmm  </td><td>  D6h<SUP>1</SUP>  </td><td>  24  </td></tr>
<tr><td>  192  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6/mcc  </td><td>  D6h<SUP>2</SUP>  </td><td>  24  </td></tr>
<tr><td>  193  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>3</SUB>/mcm  </td><td>  D6h<SUP>3</SUP>  </td><td>  24  </td></tr>
<tr><td>  194  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P6<SUB>3</SUB>/mmc  </td><td>  D6h<SUP>4</SUP>  </td><td>  24  </td></tr>
<tr><td>  195  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P23  </td><td>  T<SUP>1</SUP>  </td><td>  12  </td></tr>
<tr><td>  196  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  F23  </td><td>  T<SUP>2</SUP>  </td><td>  48  </td></tr>
<tr><td>  197  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I23  </td><td>  T<SUP>3</SUP>  </td><td>  24  </td></tr>
<tr><td>  198  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P2<SUB>1</SUB>3  </td><td>  T<SUP>4</SUP>  </td><td>  12  </td></tr>
<tr><td>  199  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I2<SUB>1</SUB>3  </td><td>  T<SUP>5</SUP>  </td><td>  24  </td></tr>
<tr><td>  200  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pm-3  </td><td>  Th<SUP>1</SUP>  </td><td>  24  </td></tr>
<tr><td>  201  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pn-3_1  </td><td>  Th<SUP>2</SUP>  </td><td>  24  </td></tr>
<tr><td>  201  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  Pn-3_2  </td><td>  Th<SUP>2</SUP>  </td><td>  24  </td></tr>
<tr><td>  202  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  Fm-3  </td><td>  Th<SUP>3</SUP>  </td><td>  96  </td></tr>
<tr><td>  203  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  Fd-3_1  </td><td>  Th<SUP>4</SUP>  </td><td>  96  </td></tr>
<tr><td>  203  </td><td>  1  </td><td>  2  </td><td>  3  </td><td>  Fd-3_2  </td><td>  Th<SUP>4</SUP>  </td><td>  96  </td></tr>
<tr><td>  204  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  Im-3  </td><td>  Th<SUP>5</SUP>  </td><td>  48  </td></tr>
<tr><td>  205  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pa-3  </td><td>  Th<SUP>6</SUP>  </td><td>  24  </td></tr>
<tr><td>  206  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  Ia-3  </td><td>  Th<SUP>7</SUP>  </td><td>  48  </td></tr>
<tr><td>  207  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P432  </td><td>  O<SUP>1</SUP>  </td><td>  24  </td></tr>
<tr><td>  208  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>2</SUB>32  </td><td>  O<SUP>2</SUP>  </td><td>  24  </td></tr>
<tr><td>  209  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  F432  </td><td>  O<SUP>3</SUP>  </td><td>  96  </td></tr>
<tr><td>  210  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  F4<SUB>1</SUB>32  </td><td>  O<SUP>4</SUP>  </td><td>  96  </td></tr>
<tr><td>  211  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I432  </td><td>  O<SUP>5</SUP>  </td><td>  48  </td></tr>
<tr><td>  212  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>3</SUB>32  </td><td>  O<SUP>6</SUP>  </td><td>  24  </td></tr>
<tr><td>  213  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P4<SUB>1</SUB>32  </td><td>  O<SUP>7</SUP>  </td><td>  24  </td></tr>
<tr><td>  214  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I4<SUB>1</SUB>32  </td><td>  O<SUP>8</SUP>  </td><td>  48  </td></tr>
<tr><td>  215  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-43m  </td><td>  Td<SUP>1</SUP>  </td><td>  24  </td></tr>
<tr><td>  216  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  F-43m  </td><td>  Td<SUP>2</SUP>  </td><td>  96  </td></tr>
<tr><td>  217  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I-43m  </td><td>  Td<SUP>3</SUP>  </td><td>  48  </td></tr>
<tr><td>  218  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  P-43n  </td><td>  Td<SUP>4</SUP>  </td><td>  24  </td></tr>
<tr><td>  219  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  F-43c  </td><td>  Td<SUP>5</SUP>  </td><td>  96  </td></tr>
<tr><td>  220  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  I-43d  </td><td>  Td<SUP>6</SUP>  </td><td>  48  </td></tr>
<tr><td>  221  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pm-3m  </td><td>  Oh<SUP>1</SUP>  </td><td>  48  </td></tr>
<tr><td>  222  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pn-3n_1  </td><td>  Oh<SUP>2</SUP>  </td><td>  48  </td></tr>
<tr><td>  222  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  Pn-3n_2  </td><td>  Oh<SUP>2</SUP>  </td><td>  48  </td></tr>
<tr><td>  223  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pm-3n  </td><td>  Oh<SUP>3</SUP>  </td><td>  48  </td></tr>
<tr><td>  224  </td><td>  1  </td><td>  1  </td><td>  1  </td><td>  Pn-3m_1  </td><td>  Oh<SUP>4</SUP>  </td><td>  48  </td></tr>
<tr><td>  224  </td><td>  1  </td><td>  2  </td><td>  1  </td><td>  Pn-3m_2  </td><td>  Oh<SUP>4</SUP>  </td><td>  48  </td></tr>
<tr><td>  225  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  Fm-3m  </td><td>  Oh<SUP>5</SUP>  </td><td>  192  </td></tr>
<tr><td>  226  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  Fm-3c  </td><td>  Oh<SUP>6</SUP>  </td><td>  192  </td></tr>
<tr><td>  227  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  Fd-3m_1  </td><td>  Oh<SUP>7</SUP>  </td><td>  192  </td></tr>
<tr><td>  227  </td><td>  1  </td><td>  2  </td><td>  3  </td><td>  Fd-3m_2  </td><td>  Oh<SUP>7</SUP>  </td><td>  192  </td></tr>
<tr><td>  228  </td><td>  1  </td><td>  1  </td><td>  3  </td><td>  Fd-3c_1  </td><td>  Oh<SUP>8</SUP>  </td><td>  192  </td></tr>
<tr><td>  228  </td><td>  1  </td><td>  2  </td><td>  3  </td><td>  Fd-3c_2  </td><td>  Oh<SUP>8</SUP>  </td><td>  192  </td></tr>
<tr><td>  229  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  Im-3m  </td><td>  Oh<SUP>9</SUP>  </td><td>  96  </td></tr>
<tr><td>  230  </td><td>  1  </td><td>  1  </td><td>  2  </td><td>  Ia-3d  </td><td>  Oh<SUP>10</SUP>  </td><td>  96  </td></tr>
</TABLE>
</body></html>
