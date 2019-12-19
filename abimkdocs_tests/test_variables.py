# coding: utf-8
"""Tests abiref.bib file."""
from __future__ import division, print_function, unicode_literals, absolute_import

from .tools import patch_syspath, AbimkdocsTest
patch_syspath()

import os
import json

from collections import OrderedDict
from pprint import pprint
from abimkdocs.variables import get_codevars, ValueWithUnit, ValueWithConditions, MultipleValue, Range


class VariablesTest(AbimkdocsTest):

    def test_variables(self):
        codevars = get_codevars()
        assert codevars is get_codevars()
        assert "abinit" in codevars
        abinit_vars = codevars["abinit"]
        assert "ecut" in abinit_vars
        ecut = abinit_vars["ecut"]
        assert ecut == ecut
        assert ecut != None
        assert ecut.name == "ecut"
        assert ecut.vartype == "real"
        assert not ecut.isarray
        assert ecut != abinit_vars["pawecutdg"]
        assert ecut in {ecut: "foo"}
        assert ecut.executable == "abinit"
        assert ecut.wikilink == "[[abinit:ecut]]"
        assert not ecut.depends_on_dimension("ntypat")
        assert not ecut.is_internal
        assert abinit_vars["mpw"].is_internal
        acell = abinit_vars["acell"]
        assert acell.isarray
        assert not acell.depends_on_dimension("ntypat")
        xred = abinit_vars["xred"]
        assert xred.isarray and xred.varset == "basic" and "[[EVOLVING]]" in xred.characteristics
        assert xred.depends_on_dimension("natom")
        assert xred.depends_on_dimension("natrd")
        assert len(acell.dimensions) == 1 and acell.dimensions[0] == 3
        abinit_asr = abinit_vars["asr"]
        anaddb_asr = codevars["anaddb"]["asr"]
        assert abinit_asr.executable == "abinit"
        assert abinit_asr.name == "asr"
        assert anaddb_asr.executable == "anaddb"
        assert anaddb_asr.name == "asr"
        assert anaddb_asr.wikilink == "[[anaddb:asr]]"
        assert anaddb_asr != abinit_asr
        symsigma_parents = abinit_vars["symsigma"].get_parent_names()
        assert len(symsigma_parents) > 0 and "optdriver" in symsigma_parents

        # The text of this variable contaings greek symbols in HTML.
        var = abinit_vars["cd_frqim_method"]
        repr(var); str(var)

        # Test "tricky" variables.
        fxcartfactor = codevars["abinit"]["fxcartfactor"]
        str(fxcartfactor)
        assert fxcartfactor.to_abimarkdown()
        d = fxcartfactor.topic2relevances
        assert fxcartfactor.topic2relevances is d and len(d) == 2
        assert "expert" in d["TransPath"] and "expert" in d['GeoOpt']
        assert isinstance(fxcartfactor.defaultval, ValueWithUnit)
        assert fxcartfactor.defaultval.units == "(Bohr^2)/Hartree"
        assert fxcartfactor.defaultval.value == 1

        iomode = codevars["abinit"]["iomode"]
        str(iomode)
        assert iomode.to_abimarkdown()
        assert len(iomode.characteristics) == 1 and iomode.characteristics[0] == "[[DEVELOP]]"
        assert isinstance(iomode.defaultval, ValueWithConditions)
        # TODO: waiting for new version dict-based.
        #assert iomode.defaultval["defaultval"] == 0
        #assert iomode.defaultval['[[MPI_IO]] and [[paral_kgb]]==1'] == 1

        istwfk = codevars["abinit"]["istwfk"]
        assert istwfk.to_abimarkdown()
        assert istwfk.characteristics is None and istwfk.vartype == "integer" and istwfk.varset == "dev"
        assert istwfk.dimensions == ["[[nkpt]]"]
        assert isinstance(istwfk.defaultval, MultipleValue)
        #assert istwfk.defaultval.number is None
        #assert istwfk.defaultval.value == 0
        assert istwfk.requires is None
        assert "nkpt" in istwfk.get_parent_names()

        jdtset = codevars["abinit"]["jdtset"]
        assert isinstance(jdtset.defaultval, Range)
        #assert jdtset.defaultval.start == 1
        #assert jdtset.defaultval.stop == "[[ndtset]]"

        # Abipy test
        # Database methods.
        database = abinit_vars
        assert database.apropos("ecut")
        #assert len(database.json_dumps_varnames())

        print("vargeo section:\n", database.vars_with_varset("vargeo"))
        for section in database.my_varset_list:
            assert len(database.vars_with_varset(section))

        for charact in database.my_characteristics:
            #print("character:", charact)
            assert len(database.vars_with_char(charact))

        name2varset = database.name2varset
        assert name2varset["ecut"] == "basic" and name2varset["ionmov"] == "rlx"

        print("d:", database.group_by_varset("ecut"), "hello")
        assert database.group_by_varset("ecut") ==  {'basic': ['ecut']}

        #abinit_help("ecut", info=True)
        # Should not raise
        #abinit_help("foobar", info=True)
        #for codename, d in codevars.items():
        #    d.validate_vars()

        for var in codevars.iter_allvars():
            try:
                assert repr(var)
                assert str(var)
                assert var.to_string(verbose=2)
                assert var._repr_html_()
                assert var.info
                assert var.to_abimarkdown()
                # topic --> list of tribes
                assert len(var.topic2relevances)
                var.validate()
            except Exception as exc:
                print("Error in %s:\n%s" % (var.abivarname, str(exc)))
                raise

    def test_variables_in_tests(self):
        """
        Find variables that are not tested by comparing the database of variables
        with the input file of the test.
        Build dictionary with list of untested variables for the different codes.
        Finally compare the new dictionary with the reference one and fail if they don't match.
        """
        # Build database with all input variables indexed by code name.
        from abimkdocs.variables import get_codevars
        codevars = get_codevars()

        # Build AbinitTestSuite object.
        from doc import tests as tmod
        tests = tmod.abitests.select_tests(suite_args=[], regenerate=True, flat_list=True)

        # Build conter for the different codes, keys are the varnames from the database.
        from collections import Counter
        count_code = {}
        for code, d in codevars.items():
            count_code[code] = Counter({k: 0 for k in d})

        ierr = 0
        doc_vnames = codevars["anaddb"].get_all_vnames(with_internal=False)
        anaddb_f90vnames = self.get_anaddb_varnames_from_f90()
        #print("anaddb_vames:", anaddb_f90vnames)
        diff = anaddb_f90vnames - doc_vnames
        if diff:
            ierr += 1
            print("\nThe following variables are found in anaddb F90 code but not in variables_anaddb.py")
            pprint(diff)

        diff = doc_vnames - anaddb_f90vnames
        if diff:
            ierr += 1
            print("\nThe following variables are found in variables_anaddb.py but not in anaddb F90 code.")
            pprint(diff)

        doc_vnames = codevars["abinit"].get_all_vnames(with_internal=False)
        abinit_f90vnames = self.get_abinit_varnames_from_f90()
        #print("abinit_f90vnames:", abinit_f90vnames)

        diff = abinit_f90vnames - doc_vnames
        if diff:
            ierr += 1
            print("\nThe following variables are found in abinit F90 code but not in variables_abinit.py")
            pprint(diff)

        diff = doc_vnames - abinit_f90vnames
        if diff:
            ierr += 1
            print("\nThe following variables are found in variables_abinit.py but not in abinit F90 code")
            pprint(diff)

        # TODO: should parse chkvars and
        black_list = set([
            "atompaw", "cut3d", "multibinit", "fftprof", "conducti", "mrgscr", "tdep",
            "mrgddb", "mrggkk", "mrgdv", "band2eps", "ujdet", "fold2Bloch", "macroave", "testtransposer",
        ])
        for test in tests:
            if test.executable in black_list: continue
            vnset = test.get_varname_set()
            count_code[test.executable].update(vnset)

        untested = OrderedDict()
        for code, count in count_code.items():
            untested[code] = []
            # Add it if var is not tested and not internal.
            for vname, c in count.items():
                if c == 0 and not codevars[code][vname].is_internal:
                    #print(code, vname)
                    untested[code].append(vname)

        for code in sorted(untested.keys()):
            untested[code] = sorted(untested[code])
            if untested[code]:
                print("\nList of untested variables for code:", code)
                pprint(untested[code])

        #ref_json_path = os.path.join(os.path.dirname(__file__), "untested_variables.json")
        #update_ref = False
        #if update_ref:
        #    with open(ref_json_path, "wt") as fh:
        #        json.dump(untested, fh, indent=4)
        #else:
        #    with open(ref_json_path, "rt") as fh:
        #        ref_untested = json.load(fh)

        #    self.assertDictEqual(untested, ref_untested,
        #        msg="Detected mismatch between reference file and new list of untested variables.")

        assert ierr == 0
