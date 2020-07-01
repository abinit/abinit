#!/usr/bin/env python

import sys
import re
import os.path
import numpy as N
import sgmllib
from pprint import pprint

PTG_NUMS=[1, 2, 3, 6, 10, 16, 25, 47, 75, 81, 83, 89, 99, 111, 
        123, 143, 147, 149, 156, 162, 168, 174, 175, 177, 183, 189, 191, 195, 200, 207, 215, 221]

ptg_names =  [
["C1"  , "1"],   
["Ci"  , "-1"],
["C2"  , "2"],
["Cs"  , "m"],
["C2h" , "2/m"],	
["D2"  , "222"],	
["C2v" , "mm2"],
["D2h" , "mmm"],
["C4"  , "4"],
["S4"  , "-4"],
["C4h" , "4/m"],	
["D4"  , "422"],
["C4v" , "4mm"],
["D2d" , "-42m"],
["D4h" , "4/mmm"],
["C3"  , "3"],
["C3i" , "-3"],
["D3"  , "32"],
["C3v" , "3m"],
["D3d" , "-3m"],
["C6"  , "6"],
["C3h" , "-6"],
["C6h" , "6/m"], 
["D6"  , "622"],
["C6v" , "6mm"],	
["D3h" , "-62m"],	
["D6h" , "6/mmm"],
["T"   , "23"],
["Th"  , "m-3"],
["O"   , "432"],
["Td"  , "-43m"],
["Oh"  , "m-3m"],
]

_E3D = N.identity(3,  N.int)
_03D = N.zeros((3,3), N.int)

def dotc(v1, v2):
    """Scalar product between two complex vectors"""
    return N.dot(N.conjugate(v1), v2)

def rflat(iterables):
    """Iterator over all elements of a nested iterable. It's recursive!"""
    for item in iterables:
        if not hasattr(item, "__iter__"): 
            yield item
        else: # iterable object.
            for it in rflat(item): yield it

def npymat2Fmat(npymat):
    """Return a string with the F90 declaration of a numpy matrix."""
                                                                    
    shape = str(npymat.shape)
    shape = shape.replace("(","(/").replace(")","/)")

    fmat = npymat.T # From C to F ordering.
    vect = ""
    for idx, ele in enumerate(fmat.flat): 
        #print ele
        tk = str(ele)
        tk = tk.replace("+-","-")  # Have to Fix weird output of str when imag part is negative (+-)
        tk = tk.replace("++","+")  # Fortran compilers will likely complain!
        tk = tk.replace("(","")
        tk = tk.replace(")","")
        vect += tk
        if idx != (fmat.size-1): vect += ", " 

    vect = vect.replace("j","*j") # For complex numbers, j has to be defined in F source.
    vect = "(/" + vect + "/)"
    return " RESHAPE( %(vect)s, %(shape)s )" % locals()

#########################################################################################
class RotationException(Exception):
    pass


class Rotation(object):
    """Object describing a pure rotation (proper, improper, mirror symmetry"""

    def __init__(self, rotation, order=None, trcoords=None, versor=None):

        self.rotation = N.matrix(rotation, N.int)

        self.order = order
        self.trcoords = trcoords

        if versor is None:
            self.versor = [0,0,0]
        else:
            self.versor = versor

        dd = {0:"x",1:"y",2:"z"} 
        self.my_trcoords = ""
        for ridx, row in enumerate(self.rotation):
            for cidx, el in enumerate(row.flat):
                ax = dd[cidx] 
                if   el == -1: self.my_trcoords += "-" + ax
                elif el ==  0: pass
                elif el == +1: self.my_trcoords += "+" + ax
                else: 
                    raise RotationException("wrong element value" + str(el))
            if ridx < 2: self.my_trcoords += ", "


    # Might subclass N.matrix though.
    def __eq__(self, other):
        return N.allclose(self.rotation, other.rotation)
    def __neq__(self, other):
        return not self == other 

    # Implement the unary arithmetic operations (+, -)
    def __pos__(self): return self
    def __neg__(self): return Rotation(-self.rotation) 

    def __mul__(self, other):
        return Rotation(self.rotation * other.rotation)

    def __pow__(self, intexp, modulo=1):
       if intexp ==  0: return Rotation(_E3D)
       if intexp  >  0: return Rotation(self.rotation**intexp)
       if intexp == -1: return self.invert()
       if intexp  <  0: return self.__pow__(-intexp).invert()

    def _get_det(self):
        """Return the determinant of a symmetry matrix mat[3,3]. It must be +-1"""
        mat = self.rotation
        det =   mat[0,0]* ( mat[1,1]*mat[2,2] - mat[1,2]*mat[2,1] )\
              - mat[0,1]* ( mat[1,0]*mat[2,2] - mat[1,2]*mat[2,0] )\
              + mat[0,2]* ( mat[1,0]*mat[2,1] - mat[1,1]*mat[2,0] )

        if abs(det) != 1: 
            raise RotationException("abs(det) must be 1 while it is " + str(abs(det)))
        else:
            return det

    det = property(_get_det, doc="The determinant of the rotation")

    def _get_trace(self):
        return self.rotation.trace()[0,0]
    trace = property(_get_trace, doc="The trace of the rotation")

    def _isproper(self):
        return bool(self.det+1)
    isproper = property(_isproper, doc="True if proper rotation")

    def invert(self): 
        """
        Invert an orthogonal 3x3 matrix of INTEGER elements.
        Note use of integer arithmetic. Raise RotationException if not invertible.
        """

        det = self.det
        mm = self.rotation
        inv= N.matrix(N.zeros((3,3), N.int))

        inv[0,0] = mm[1,1] * mm[2,2] - mm[1,2] * mm[2,1]
        inv[0,1] = mm[0,2] * mm[2,1] - mm[0,1] * mm[2,2] 
        inv[0,2] = mm[0,1] * mm[1,2] - mm[0,2] * mm[1,1]

        inv[1,0] = mm[1,2] * mm[2,0] - mm[1,0] * mm[2,2]
        inv[1,1] = mm[0,0] * mm[2,2] - mm[0,2] * mm[2,0]
        inv[1,2] = mm[0,2] * mm[1,0] - mm[0,0] * mm[1,2] 

        inv[2,0] = mm[1,0] * mm[2,1] - mm[1,1] * mm[2,0]
        inv[2,1] = mm[0,1] * mm[2,0] - mm[0,0] * mm[2,1] 
        inv[2,2] = mm[0,0] * mm[1,1] - mm[0,1] * mm[1,0]

        # Make sure matrix is not singular
        if det != 0: 
            return Rotation(inv/det)
        else: 
            raise RotationException("Attempting to invert singular matrix")

    def _rottype(self):
        """
        Receive a 3x3 orthogonal matrix and reports its type:
         1 Identity
         2 Inversion
         3 Proper rotation of an angle <> 180 degrees
         4 Proper rotation of 180 degrees
         5 Mirror symmetry
         6 Improper rotation
        """
        rot = self.rotation # Just an alias. 

        # Treat identity and inversion first
        #identity  = Rotation(_E3D)

        if self.isE: return 1
        if self.isI: return 2
    
        if self.isproper: # Proper rotation
            t = 3 # try angle != 180
            #det180 = get_sym_det(rot+_E3D)
            if (self + identity).det == 0: t = 4 # 180 rotation
        else: # Mirror symmetry or Improper rotation
            t = 6
            #detmirror = get_sym_det(rot-_E3D)
            if (self - identity).det == 0: t = 5 # Mirror symmetry if an eigenvalue is 1

        return t

    def _isE(self):
        return N.allclose(self.rotation, _E3D)
    isE = property(_isE, doc="True if it is the identity")

    def _isI(self):
        return N.allclose(self.rotation, -_E3D)
    isI = property(_isI, doc="True if it is the inversion")

    def get_versor(self):
        raise NotImplementedError

        if self.isE or self.isI:
            versor = [0, 0, 0]
        return versor

    def _get_order(self):

        order = None
        root_invers = 0
        for ior in range(1,7):
            rn = self**ior
            if rn.isE:
                order = ior
                break
            if rn.isI: root_invers = ior
        if order is None: 
            raise RotationException("symmetry is not a root of unit!")
        return (order, root_invers)

    info = property(_get_order, doc="Order and root of unit")

    def _get_name(self):

        order, root_invers = self.info
        name = ""
        if self.det == -1: name = "-"

        name += str(order) # FIXME this one doesn't work yet.
        if root_invers != 0: 
            name += "-"
        else:
            name += "+"
        return name

    name = property(_get_name, doc="")

    def __str__(self):
        string  = "Rotation: " + str(self.order) + ", versor: " + str(self.versor) + ", " + str(self.trcoords) + "\n"
        string += str(self.rotation) 
        return string

    def toFortran(self, varname):
        """Return a string with the F90 declaration of the symmetry."""

        fmat = self.rotation.T # From C to F indexing
        string = "RESHAPE(("
        vect = str([e for e in fmat.flat])
        vect = vect.replace("[","(/")
        vect = vect.replace("]","/)")
        return " %(varname)s = RESHAPE( %(vect)s ,(/3,3/) )" % locals()

#########################################################################################
class IrreducibleRepr(object):
    """Class defining an irreducible representation"""  

    def __init__(self, name, dim, matrices):
        self.name = name
        self.dim  = dim
        self.matrices = matrices
        #self.matrices = [ N.matrix(mat) for mat in matrices ]

    def __str__(self):
        string = " Irred repr: " + self.name + " dimension= " + str(self.dim) + "\n"
        for mat in self.matrices: string +=  str(mat) + "\n" 
        return string

    def traces(self):
        return N.array( [ mat.trace()[0,0] for mat in self.matrices ] )

def mk_classes(rotations):
    """
    Find the classes of a group. Return a list containing the indeces
    of the operations in each class sorted in ascending order.
                                                                       
    A class is defined as the set of distinct items obtained by 
    considering for each element, S, of the group all its conjugate 
    X^-1 S X where X is one of the elements of the group.
    """

    class_ids = list()
    seen = [ 0 for i in range(len(rotations)) ]
    nclass = -1 
                                                                        
    for idx, m in enumerate(rotations):
        if seen[idx]: continue
        seen[idx] = 1

        nclass += 1
        class_ids.append([])
        class_ids[nclass].append(idx)

        for x in rotations:
            new =  x.invert() * m * x
            idx_found = -1
            for idx_search, search in enumerate(rotations):
                if search == new:
                    idx_found = idx_search
                    break
            if (idx_found == -1): sys.stderr.write("idx_found == -1")
            if not seen[idx_found]: 
                seen[idx_found] = 1
                class_ids[nclass].append(idx_found)
                                                                        
    # Now sort the indeces.
    sort_class_ids = list()
    for ids in class_ids:
        ids.sort()
        sort_class_ids.append(ids)
    return sort_class_ids


#########################################################################################
class PointGroupException(Exception):
    pass

class NotAGroup(PointGroupException):
    pass

class RotationNotFound(PointGroupException):
    pass

class PointGroup(list):
    """A PointGroup is a list of Rotations and it has irreducible representations""" 

    def __init__(self, rotations, name=None, irreprs=None):

        class_ids = mk_classes(rotations)

        # Always reorder rotations and irreprs according to class indeces.
        ord_rotations = [ None for ii in range(len(rotations)) ]
        idx = -1
        for ord_idx in rflat(class_ids): 
            idx += 1
            ord_rotations[idx] = rotations[ord_idx]

        ord_irreprs  = list() 
        for irr in irreprs:
            ord_matrices = [ None for ii in range(len(irr.matrices)) ]
            idx = -1
            for ord_idx in rflat(class_ids):
               idx += 1
               ord_matrices[idx] = irr.matrices[ord_idx]

            ord_irreprs.append( IrreducibleRepr(irr.name, irr.dim, ord_matrices) )

        list.__init__(self)
        for orot in ord_rotations: self.append(orot)

        if __debug__ and False:
            for rot in self: 
                print "info", rot.det, rot.info
                print "name",  rot.name
                print "order", rot.order

        self.class_ids = mk_classes(ord_rotations)
        self.nclass = len(self.class_ids)

        # Create name of each class.
        #self.class_names = [ "None" for ii in range(self.nclass) ]
        first_rot_ids = [ self.class_ids[ii][0] for ii in range(self.nclass) ]
        self.class_names = [ self[ii].name for ii in first_rot_ids ]

        self.nsym = len(self)
        self.name = str(name)

        self.irreprs = ord_irreprs
        #for ii in self.irreprs: print ii

        self.nirrepr = len(self.irreprs)

    def find(self, rot):
        """Return the index of rot."""
        try:
            return self.index(rot)
        except ValueError:
            raise RotationNotFound(rot)

    def findE(self):
        """Return the index of the identity."""
        try:
            return self.index(Rotation(_E3D))
        except RotationNotFound:
            raise

    def find_inverse(self, rot):
        """Return the index of the inverse of rot."""
        E = Rotation(_E3D)
        for s in self:
            if s * rot == E: return s

        sys.stderr.write("Warning: Inverse not found!!\n")
        raise RotationNotFound(rot)

    def isgroup(self):
        try:
            self.findE() 
            for rot in self: self.find_inverse(rot)
            return True
        except RotationNotFound:
            sys.stderr.write("Not a group! Identity or inverse are missing")
            return False

    def mk_mtable(self):
        """Check if it is a group, then build the multiplication table"""

        # Check if r1 * r2 is in group and build the multiplication table.
        mtable = dict()
        for idx1, r1 in enumerate(self):
            for idx2, r2 in enumerate(self): 

                try:
                    ij = (idx1, idx2)
                    mtable[ij] = self.index(r1 * r2)
                except RotationNotFound:
                    sys.stderr.write("Not a group. Not close wrt *")
                    raise
        return mtable

    def show_mtable(self):
        """Print out multiplication table."""
                                                                                     
        mtable = self.mk_mtable() 

        print 4*" " + (2*" ").join([str(i) for i in xrange(self.nsym)]) + "\n"
        for i in xrange(self.nsym):
            lij = [(i, j) for j in xrange(self.nsym)]
            print str(i) + (2*" ").join([str(mtable[ij]) for ij in lij]) + "\n"

    def show_character_table(self):
        vlen = 10
                                                
        print 100*"*"
        print ("Point Group" + self.name)
        cln = ""
        for clname in self.class_names:
            cln += str(clname).center(vlen)
        print "Class" + cln 
                                                
        mult = "Mult" 
        for cls in self.class_ids:
            mult += str(len(cls)).center(vlen)
        print mult
                                                
        for irrepr in self.irreprs:
            #print "irrepr ", irrepr
            row = irrepr.name.ljust(5)
            for icls in range(self.nclass): 
               sym_id = self.class_ids[icls][0] 
               mat = irrepr.matrices[sym_id]
               char = mat.trace()[0,0]
               row += str(char).center(vlen)
            print row
                                                
        print 100*"*"
        print 100*"*"

    def check(self):
        if not self.isgroup(): raise NotAGroup

        class_ids = mk_classes(self)
        #print class_ids
        check = -1
        for idx in rflat(class_ids):
            check = check +1
            if check!= idx: raise PointGroupException("Symmetries are not ordered by classes")

        mtable = self.mk_mtable()

        err = 0.0
        for idx1 in range(len(self)):
            for idx2 in range(len(self)):
                ij = (idx1, idx2)
                idx_prod = mtable[ij] 

                for irr in self.irreprs:
                    mat_prod = irr.matrices[idx1] * irr.matrices[idx2]
                    my_err = (mat_prod - irr.matrices[idx_prod]).max()
                    err = max(err, abs(my_err))

        print "Error in Group Representation", err

        character_of = dict()
        for irr in self.irreprs:
            traces = irr.traces()
            #character = [ traces[ii] 
            chr = list()
            for clids in self.class_ids:
                idx = clids[0]
                chr.append(traces[idx])
            #character_of[irr.name] = N.array(chr)
            character_of[irr.name] = traces
            #irr.name 

        err_otrace = 0.0
        for k1, v1 in character_of.iteritems():
            for k2, v2 in character_of.iteritems():
                my_err = dotc(v1, v2) / self.nsym
                if k2 == k1: my_err -= 1.0 
                err_otrace = max(err_otrace, abs(my_err))
        print "Error in orthogonality relation of traces ", err

    def dump_Fortran_sub(self, fh):

        subname = "ptg_"  + self.name.split()[0].strip()

        fh.write("!{\src2tex{textfont=tt}}\n") 
        fh.write("!!****f* ABINIT/%s\n" % subname) 
        fh.write("!!\n") 
        fh.write("!! NAME\n") 
        fh.write("!! %s\n" % subname) 
        fh.write("!!\n") 
        fh.write("!! FUNCTION\n") 
        fh.write("!!\n") 
        fh.write("!! COPYRIGHT\n") 
        fh.write("!! Copyright (C) 2010-2020 ABINIT group (MG)\n")
        fh.write("!! This file is distributed under the terms of the\n")
        fh.write("!! GNU General Public License, see ~abinit/COPYING\n")
        fh.write("!! or http://www.gnu.org/copyleft/gpl.txt .\n")
        fh.write("!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .\n")
        fh.write("!!\n")
        fh.write("!! INPUTS\n")
        fh.write("!!\n")
        fh.write("!! OUTPUT\n")
        fh.write("!!\n")
        fh.write("!! PARENTS\n")
        fh.write("!!\n")
        fh.write("!! CHILDREN\n")
        fh.write("!!\n")
        fh.write("!! SOURCE\n")
        fh.write("!!\n")
        fh.write("!" + 80*"*" + "\n")
        fh.write("! This include file has been automatically generated by the script " + os.path.basename(__file__) +"\n")
        fh.write("! Do not edit! Change the script source instead.\n") 
        fh.write("!" + 80*"*" + "\n")
        fh.write("\n")
        fh.write("! Point group name " + self.name + "\n")
        fh.write("\n")

        fh.write("#if defined HAVE_CONFIG_H\n#include \"config.h\"\n#endif\n\n")
        fh.write("#include \"abi_common.h\"\n\n")
        fh.write(" subroutine %s (nsym,nclass,sym,class_ids,class_names,Irr)\n" % subname )
        # Disable optimization if ifc is used.
        fh.write(" !DEC$ NOOPTIMIZE")
        fh.write(" use defs_basis\n")
        fh.write(" use m_profiling_abi\n")
        fh.write(" use m_defs_ptgroups,  only : irrep_t \n")
        fh.write(" implicit none\n")
        fh.write("!Arguments ------------------------------------\n")
        fh.write(" integer,intent(out) :: nclass,nsym \n")
        #fh.write(" character(len=5),intent(in) :: ptg_name \n")
        fh.write(" !arrays\n")
        fh.write(" integer,allocatable,intent(out) :: sym(:,:,:), class_ids(:,:)\n")
        fh.write(" character(len=5),allocatable,intent(out) :: class_names(:)\n")
        fh.write(" type(irrep_t),allocatable,intent(out) :: Irr(:)\n")
        #fh.write(" integer,pointer :: sym(:,:,:), class_ids(:,:)\n")
        #fh.write(" character(len=5),pointer :: class_names(:)\n")
        #fh.write(" type(irrep_t),pointer :: Irr(:)\n")
        fh.write(" !Local variables-------------------------------\n")
        fh.write(" complex(dpc) :: j=(0.0_dp,1.0_dp) \n")
        #fh.write(" character(len=500) :: msg \n")
        fh.write(" ! " + 80*"*" + "\n")


        # Write list of Symmetries first.
        fh.write("! List of symmetries packed in classes\n")
        nsym = self.nsym
        fh.write(" nsym = %s\n" % nsym)
        fh.write(" ABI_MALLOC(sym, (3,3,nsym))\n")
        isym=0
        for sym in self:
            isym += 1
            varname = "sym(:,:,%s)" % isym
            fh.write(sym.toFortran(varname)+"\n")
        fh.write("\n")

        # Write classes and their names.
        fh.write("! Number of classes and corresponding indeces\n")
        nclass = self.nclass
        fh.write(" nclass = %s\n" % nclass)
        fh.write(" ABI_MALLOC(class_ids, (2,nclass))\n")

        for iclas in range(self.nclass): 
            first = self.class_ids[iclas][0]  +1  # From C to Fortran
            last  = self.class_ids[iclas][-1] +1
            idx = iclas + 1
            fh.write(" class_ids(1,%(idx)s) = %(first)s\n" % locals())
            fh.write(" class_ids(2,%(idx)s) = %(last)s\n" % locals())

        # Write the name of each class to be reported in the table.
        fh.write("\n")
        fh.write("ABI_MALLOC(class_names,(%(nclass)s))\n" % locals())
        idx = 0
        for iclas in range(nclass):
            idx += 1
            name = self.class_names[iclas]
            if name: name = name.lstrip().rstrip()
            fh.write(" class_names(%(idx)s) = \"%(name)s\" \n" % locals() )

        # Write the irreducible representations.
        fh.write("\n")
        fh.write("! List of irreducible representations.\n")
        #fh.write(" nirrepr = %s\n" % self.nirrepr)
        fh.write(" ABI_MALLOC(Irr, (%(nclass)s))\n" % locals())
        idx = 0
        for irrepr in self.irreprs:
            idx += 1
            dim = irrepr.dim
            name = irrepr.name.lstrip().rstrip()
            fh.write(" Irr(%(idx)s)%%name = \"%(name)s\"\n" % locals())
            fh.write(" Irr(%(idx)s)%%dim = %(dim)s\n"  % locals())
            fh.write(" Irr(%(idx)s)%%nsym = %(nsym)s\n"  % locals())
            fh.write(" ABI_MALLOC(Irr(%(idx)s)%%mat, (%(dim)s,%(dim)s,%(nsym)s))\n" % locals())
            irp = 0
            for mat in irrepr.matrices:
                fmat = npymat2Fmat(mat)
                irp += 1
                fh.write(" Irr(%(idx)s)%%mat(:,:,%(irp)s) = %(fmat)s\n" % locals())
            fh.write("\n")

        fh.write(" RETURN\n ")
        fh.write(" if (.FALSE.) write(std_out,*) j\n")

        fh.write(" end subroutine %s \n" % subname)
        fh.write("!!***\n")

        return None

    def to_dict(self):
        d = {}
        #subname = "ptg_"  + self.name.split()[0].strip()
        #fh.write("! Point group name " + self.name + "\n")

        # List of symmetries packed in classes
        d["rotations"] = [o.rotation.tolist() for o in self]

        # Write classes and their names.
        #fh.write("! Number of classes and corresponding indeces\n")
        d["nclass"] = self.nclass

        #fh.write(" allocate(class_ids(2,nclass))\n")
        class_range = []
        for iclas in range(self.nclass): 
            first = self.class_ids[iclas][0]  #+1  # From C to Fortran
            last  = self.class_ids[iclas][-1] #+1
            #fh.write(" class_ids(1,%(idx)s) = %(first)s\n" % locals())
            #fh.write(" class_ids(2,%(idx)s) = %(last)s\n" % locals())
            class_range.append((first, last+1))

        assert last+1 == self.nsym
        d["class_range"] = class_range


        # Write the name of each class to be reported in the table.
        #fh.write("\n")
        #fh.write(" allocate(class_names(%(nclass)s))\n" % locals())
        #for iclas in range(self.nclass):
        #    name = self.class_names[iclas]
        #    if name: name = name.lstrip().rstrip()
        #    fh.write(" class_names(%(idx)s) = \"%(name)s\" \n" % locals() )

        d["class_names"] = [name.lstrip().rstrip() for name in self.class_names]

        # Write the irreducible representations.
        #fh.write("! List of irreducible representations.\n")
        #fh.write(" nirrepr = %s\n" % self.nirrepr)
        irreps = {}
        for irrepr in self.irreprs:
            name = irrepr.name.lstrip().rstrip()
            dim = irrepr.dim
            #fh.write(" Irr(%(idx)s)%%name = \"%(name)s\"\n" % locals())
            #fh.write(" Irr(%(idx)s)%%dim = %(dim)s\n"  % locals())
            #fh.write(" Irr(%(idx)s)%%nsym = %(nsym)s\n"  % locals())
            #fh.write(" allocate(Irr(%(idx)s)%%mat(%(dim)s,%(dim)s,%(nsym)s))\n" % locals())
            #for mat in irrepr.matrices:
                #fmat = npymat2Fmat(mat)
                #fh.write(" Irr(%(idx)s)%%mat(:,:,%(irp)s) = %(fmat)s\n" % locals())

            assert name not in irreps
            irreps[name] = {"dim": dim, "matrices": [mat.tolist() for mat in irrepr.matrices]}

        d["irreps"] = irreps
        return d


#########################################################################################
class PtGroupParserException(Exception):
    pass

class PtGroupParser(object):

    def __init__(self):
        self.text = list()

    def parse_sym_field(self):
        """Parse the section with the symmetry operations. Returns the list of symmetries found"""

        symmetries = list()

        #2 : -x,-y,z         => 2+ [[ 0 0 1 ]]
        re_symheader = re.compile("(\d+)\s?:\s+(.+)=>\s+([-]?[\d\w][+-]?)(.*)")

        for line in self.sym_field:
            m = re_symheader.match(line)
            if m: 
                #print "match ", line
                idx      = m.group(1) 
                trcoords = m.group(2).rstrip()
                order    = m.group(3)
                strversor   = m.group(4)
                if not strversor: strversor = "[ 0 0 0 ]"
                strversor = strversor.lstrip().rstrip()
                strversor = strversor.replace("[ ","[")
                strversor = strversor.replace(" ] ","]")
                #print "idx", idx, "trcoords", trcoords , "order", order, "versor", strversor, "\n"
                strmat = ""
                row = 0
            else:
                if row in (0,1): strmat += line + ","
                if row == 2: strmat += line 
                row += 1

            if row == 3: 
                #print strmat
                #matrix = strmat.replace(".",".,")
                matrix = strmat.replace(".",",")
                matrix = eval(matrix)
                versor = strversor.replace(" ",",")
                versor = eval(versor)
                #print matrix, versor

                # Istanciate the symmetry and append it to symmetries.
                rot = Rotation(matrix, order, trcoords, versor)
                symmetries.append(rot)

        return symmetries

    def parse_irrepr_field(self):

        irrow = list()
        for name, dim, text in zip(self.irrepr_names, self.irrepr_dims, self.irrepr_fields):

            # Normalize input such that each string is of the form "1 : matrix"
            newtext = list()
            sbuf = None
            for line in text:
                if ":" in line:
                    if sbuf: newtext.append(sbuf)
                    sbuf = line + " "
                else:
                    sbuf += line
            newtext.append(sbuf)

            matrices = [ None for x in range(self.nsym) ]

            for line in newtext: 
                #print name, line, self.nsym
                tokens = line.split(":")
                idx = int(tokens[0]) -1
                val = tokens[1]
                if dim != 1: # Convert string to python list. #FIXME this is not safe!
                    #print line
                    #print "before eval: ",val
                    val = val.replace("] ","], ")
                    if "j" not in val:
                        val = val.replace(". ","., ") 
                    else:
                        #print "before eval: ",val
                        val = val.replace("j ","j, ") # Dirty but it works 
                        #print "after eval : ", val

                val = eval(val)
                #print idx, val
                matrices[idx] = N.matrix(val)
                #print matrices[idx]

            irepr = IrreducibleRepr(name, dim, matrices)
            #print irepr
            irrow.append(irepr)

        return irrow

    def parse(self, fname):

        # Read data from file and remove blank lines
        lines = open(fname,"r").readlines()
        self.text = [ line for line in lines if line.lstrip() ]

        # Extract the group name from the first string.
        title = self.text.pop(0)
        substr = "Symmetry operations for the point group"
        if substr not in title: 
            raise PtGroupParser("Read wrong title " + title)

        self.pgroup_name = title.replace(substr, "",1).rstrip()

        # Separate the symmetry field from the one with the irred. representations. 

        # 1) Consume text until we reach the files with the irreducible repr.
        self.sym_field = []
        substr = "Irreducible representations for the point group" 
        line = self.text.pop(0)
        while substr not in line:
           if line.lstrip(): self.sym_field.append(line.rstrip())
           line = self.text.pop(0)

        # Extract the group name from the string
        #check = line.replace(substr, "",1).rstrip()
        #if check != self.pgroup_name :
        #    raise ValueError, "Wrong format" + self.pgroup_name + "!=" + check

        symmetries = self.parse_sym_field()
        self.nsym = len(symmetries)

        self.irrepr_fields= list()
        self.irrepr_names = list()
        self.irrepr_dims  = list()
        self.nirrepr = 0

        # Now extract the Irreducible representations. 
        # Each representation is signaled by a line of the form
        #Irrep A ( dimension  1 )
        re_irrep = re.compile("\s*Irrep(.*)\(\s+dimension\s+(\d)\s+\)")

        #line = self.text(0)
        #nlines = len(self.text)
        for line in self.text: 
           m = re_irrep.match(line)
           #print line
           if m:
               if self.nirrepr: self.irrepr_fields.append(newirr)
               self.nirrepr += 1
               irr_name =  m.group(1)
               irr_dim  =  int(m.group(2))
               #print line, irr_name, irr_dim
               self.irrepr_names.append(irr_name)
               self.irrepr_dims.append(irr_dim)
               newirr = list()
           else:
               newirr.append(line.rstrip())

        self.irrepr_fields.append(newirr) # Save the last irred repr.

        irreprs = self.parse_irrepr_field()
        #for irr in irreprs: print irr
        #print "remaining",self.text
        xx = PointGroup(symmetries, self.pgroup_name, irreprs)
        #xx.show_table()
        #xx.dump_Fortran_inc(sys.stdout)

        return xx


#########################################################################################
class HTMLStripper(sgmllib.SGMLParser):
    """Simple HTML parser that removes all HTML tags"""

    def __init__(self):
        sgmllib.SGMLParser.__init__(self)

    def handle_data(self, data):
        self.plain_text.append(data)
                                            
    def parse(self, input):
        """Remove HTML tags from input. Return string.""" 
        self.plain_text = list()
        for item in input: self.feed(item) 
        return "".join(self.plain_text)

def download_ptgroup_tables(ptgnames=None, stripHTML=False):
    """Download point group tables from the Bilbao server.
       Return the list of filenames containing the tables
    """

    import telnetlib
    HOST = "www.cryst.ehu.es"

    if stripHTML: 
       html_parser = HTMLStripper()

    fnames = list()

    ii=-1
    for pg in PTG_NUMS:
        ii += 1
        doublename = ptg_names[ii]

        if ptgnames and doublename[0] not in ptgnames:        
            continue
        
        print "Downloading table for point group: " + str(pg) + str(doublename) + "..."
        tn = telnetlib.Telnet(HOST, 80)
        command = "GET /cgi-bin/rep/programs/sam/point.py?sg=" + str(pg) + "&what=irreps\n"

        tn.write(command)
        table = tn.read_all()

        fname = "ptgroup_" + doublename[0] + ".html" 
        if stripHTML: 
            table = html_parser.parse(table)
            print table
            fname = "ptgroup_" + doublename[0] + ".txt" 

        # Write table on file.
        fh = open(fname,"w")
        fh.write(table)
        fh.close()
        fnames.append(fname)

    return fnames

def download_klgroup_table(ita_spgnum, basis, kcoords, label=None):
    """Download point group tables from the Bilbao server.
       Return the list of filenames containing the tables
    """

    import telnetlib
    HOST = "www.cryst.ehu.es"

    # To obtain the representations for a given k-vector and given space group in text form you have to give the command 
    # GET /cgi-bin/cryst/text/nph-repr?g=[gn]&b=[p|c|a]&x=[number]&y=[number]&z=[number]&l=[label]
    # where 
    # [gn] is the group number in ITA, 
    # b corresponds to the basis in which the k-vector coordinates will be given. 
    # The choices are: p - primitive, c - centered dual, a - adjusted coefficients. 
    # x,y,z are the values for the three coordinates of the k-vector, 
    # l is the label for the k-vector. 
    # NOTE: By now the program uses the default choice for the group setting when there is more than one conventional setting.

    #print "Downloading table for point group: " + str(pg) + str(doublename) + "..."
    tn = telnetlib.Telnet(HOST, 80)

    if label is None: label="None"
    if basis not in ("pca"): raise ValueError
    k1 = kcoords[0]
    k2 = kcoords[1]
    k3 = kcoords[2]

    command = \
    "GET /cgi-bin/cryst/text/nph-repr?g=%(ita_spgnum)s&b=%(basis)s&x=%(k1)s&y=%(k2)s&z=%(k3)s&l=%(label)s\n"  % locals()
    print "Executing command " + command

    tn.write(command)

    # Read table in a single string.
    k_table = tn.read_all()

    # Write table on file.
    fname = "k_test"
    fh = open(fname,"w")
    fh.write(k_table)
    fh.close()

    return k_table.splitlines()

WIDTH = 45

def parse_klgroup_table(text_lines, fh=None):

    WIDTH = 45
    if fh is None:
        write = sys.stdout.write
    else:
        write = fh.write

    # Remove blank lines from input text.
    text = [ line for line in text_lines if line.lstrip() ]

    line = text[0]
    while "Number of elements" not in line: line = text.pop(0)
    nsym = int(line.split(":")[1])
    #print "nsym= ", nsym

    line = text.pop(0)
    while "The k-vector coordinates relative to the standard dual basis are" not in line:
        line = text.pop(0)

    kdual = line.split(":")[1].split()
    strk = ""
    for kc in kdual: strk += str(kc) + " "
    write(strk.ljust(WIDTH) + " # kdual\n")

    line = text.pop(0)
    search_str = "The little group of the k-vector has the following "
    while search_str not in line:
        line = text.pop(0)

    nsym_ltgk = int(line.replace(search_str, ""))
    write(str(nsym_ltgk).ljust(WIDTH) + " # Symmetries of the little group of k\n")
    text.pop(0) # Remove next lines 

    # Read the rows of symmetries of the little group (C-ordering).
    nfields = nsym_ltgk/5 + 1

    for ifield in range(nfields):
        sym_ids = [ int(idx) for idx in text.pop(0).split()] # Remove line with symmetry indices.

        # coy this field in sym_mats because we are going to read the matrices by colums
        sym_mats = list()
        for il in xrange(3): sym_mats.append(text.pop(0).lstrip())
        
        # Extract the columns
        nmats = 4
        if ifield == (nfields-1) and (nsym_ltgk % 4 > 0):
          nmats = nsym_ltgk % 4 # Last field might have nmats < 4

        # Print the operations of the little group.
        for ii in xrange(nmats):
            cols = cut_cols(sym_mats, ncols=4)
            print cols
            eval_tnons = "" # Fractional translations (e.g. 3/4) have to be evaluated and converted to float
            for tk in cols[3].split(): eval_tnons += str(eval("1.*"+tk)) + " "
            print eval_tnons
            write("".join([ c for c in cols[:3]]) + "\n") 
            write(eval_tnons + "\n") 

    #####################################
    # Parse the section with the irreps #
    #####################################

    # Move to the beginning of the section.
    search_str = "The little group of the k-vector has " #6  allowed irreps.
    line = text.pop(0)
    while search_str not in line: line = text.pop(0)
    nirreps = int(line.replace(search_str, "").split()[0])
    line = text.pop(0) # Remove comment: The matrices, corresponding to all of the little group elements are :

    write(str(nirreps).ljust(WIDTH) + " # Number of irreps of the little group.\n")

    # Extract the fields with the Irreducible representations.
    # Each irrep starts with an header of the form: Irrep (M)(5) ,  dimension 2
    re_irrep_header = re.compile("\s*Irrep \((\w+)\)\((\d+)\) ,  dimension (\d+)\s*")

    irrep_name = list()
    irrep_dim = list()
    irrep_num = list()
    irrep_field = list()

    txt_buf = None
    str_end = "The full-group irreps for the generators of the space group are :"

    for line in text:
        #print line
        if str_end in line: break
        m = re_irrep_header.match(line)
        if m: 
            if txt_buf: irrep_field.append(txt_buf)
            irrep_name.append   (m.group(1))
            irrep_num.append(int(m.group(2)))
            irrep_dim.append(int(m.group(3)))
            txt_buf = list()
        else:
            txt_buf.append(line)
    irrep_field.append(txt_buf) # Add last field

    assert nirreps == irrep_num[-1]

    for irp in range(nirreps):
        st = str(irrep_num[irp]) + " " + str(irrep_dim[irp]) + " " + str(irrep_name[irp]) 
        write(st.ljust(WIDTH) + " # irrep_index, irrep_dim, irrep_name\n")
        #for line in irrep_field[irp]: print line
        irreps_k = parse_irrep_k(nsym_ltgk, irrep_dim[irp], irrep_field[irp])
        for isym in range(nsym_ltgk):
            mat_str = "".join(c for c in irreps_k[isym]).lstrip()
            #mat_str = mat_str.replace(" ","")
            #mat_str = mat_str.replace(")",") ")
            write(str(isym+1).ljust(5) + mat_str + "\n")

def parse_irrep_k(nsym, dim, text):
    """This function parses text in the form
                 1                             2                
    (1.000,  0.0) (0.000,  0.0)   (1.000,120.0) (0.000,  0.0)   
    (0.000,  0.0) (1.000,  0.0)   (0.000,  0.0) (1.000,240.0)   
    """

    irreps_k = [ None for isym in range(nsym) ]

    while len(text)>0:
        sym_ids = [ int(idx) for idx in text.pop(0).split()] # Read the indices of the symmetries 
        nsym_in_row = len(sym_ids)
        #print sym_ids

        rows = list()  # Get the irreps of this symmetries
        for irow in range(dim):
            row = text.pop(0)
            #print row
            rows.append(row)

        for isym in sym_ids: #range(nsym_in_row):
            cols = cut_cols(rows, separator=")", ncols=dim)
            #print " isym =", isym, "cols= " ,cols
            irreps_k[isym-1] = cols

    return irreps_k

def cut_cols(rows, separator=None, ncols=1):
    columns = list()

    for icol in xrange(ncols):
        col = ""
        for il in xrange(len(rows)):
            line = rows[il]
            #print line
            tokens =  line.split(separator, 1)
            head = tokens[0] 
            if len(tokens) > 1: 
                tail= tokens[1] 
            else:
                tail = ""  # separator might not be in line.
            rows[il] = tail
            if not separator:
                col += head + " "
            else:
                col += head + separator
        columns.append(col)

    return columns


def write_ltgroup_file(fname, ita_spgnum, basis, kpoints, labels):
    """Write file with the irreps of the little group that can be read by abinit."""

    nkpts = len(kpoints)

    fh = open(fname,"w")
    write = fh.write

    # Header (two rows)
    fvers = 1
    write("# Little group file generated by " + os.path.basename(__file__) + "\n") 
    write(str(fvers).ljust(WIDTH)      + " # Version\n")
    write(str(ita_spgnum).ljust(WIDTH) + " # ITA space group\n")
    write(str(basis).ljust(WIDTH)      + " # Basis\n")
    write(str(nkpts).ljust(WIDTH)      + " # Number of k-points in database\n")

    # Write the complete list of k-points first
    for kpt, kname in zip(kpoints, labels):
        str_kname = ""
        for kc in kpt: str_kname += str(kc) + " "
        str_kname += kname 
        write(str_kname.ljust(WIDTH) + "\n")

    for kpt, kname in zip(kpoints, labels):
        k_table = download_klgroup_table(ita_spgnum, basis, kpt, kname)
        parse_klgroup_table(k_table, fh)

    fh.close()
    return None

###############################################################################
###############################################################################
if __name__ == "__main__":

    #if True: 
    if False: 
        # Little group tables
        ita_spgnum = 227
        basis = "p"
        kcoords = [0.5, 0.0, 0.0]
        label = "M"

        #ktable_text = download_klgroup_table(ita_spgnum, basis, kcoords, label)
        #parse_klgroup_table(ktable_text) 
        kpoints = [kcoords]
        labels = [label]

        lgroup_fname = "lgroup_" + str(ita_spgnum)        

        write_ltgroup_file(lgroup_fname, ita_spgnum, basis, kpoints, labels)
        sys.exit(2)

    if False:
        fnames = download_ptgroup_tables(["C2v"], stripHTML=True)
        sys.exit(2)

    p = PtGroupParser() 

    #fcheck  = "pgroup_T"
    #fcheck  = "pgroup_O"
    fcheck  = ["pgroup_C1", "pgroup_Ci"]
    #fcheck = None
    #for ii, pg in enumerate(PTG_NUMS):
    #  doublename = ptg_names[ii]
    ##                                    
    #  fname = "ptgroup_" + doublename[0] + ".txt"
    #  fname_new = "ptg_" + doublename[0] + ".txt"
    #  print "renaming " + fname + "into" + fname_new
    #  os.rename(fname, fname_new)

    #print "done"
    #sys.exit(0)

    dirname = "./ptgroup_data/"

    #ttt = open("tmp","w")

    all_irreps = {}

    for ii, pg in enumerate(PTG_NUMS):
      doublename = ptg_names[ii]

      fname = dirname + "ptg_" + doublename[0] + ".txt"
      #if fcheck and fname not in fcheck: continue
      #print "analysing " + str(doublename) 

      pgroup = p.parse(fname)
      #pgroup.check()
      #pgroup.show_character_table()
      #print pgroup._mk_classes()

      #if True or doublename[0] == "C2":
      key = doublename[0]
      assert key not in all_irreps
      d = pgroup.to_dict()
      #pprint(d)
      all_irreps[key] = d

      #import json 
      #s = json.dumps(d)
      #print(s)

      # >>>> Write Fortran routines <<<<
      #finc_name = "ptg_" + doublename[0] + ".F90"
      #print "Writing Fortran  file: " + finc_name
      #fh = open(finc_name,"w")
      #fh = sys.stdout
      #pgroup.dump_Fortran_sub(fh)
      fh.close 

      #ttt.write(" CASE ('" + doublename[1] + "')\n")
      #subname = "ptg_"  + doublename[0].strip()
      #ttt.write("   call %s (ptg_name,nsym,nclass,sym,class_ids,class_names,Irr)\n" % subname )

    pprint(all_irreps)
    
