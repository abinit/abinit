# coding: utf-8
"""Tests abiref.bib file."""
from __future__ import division, print_function, unicode_literals, absolute_import

from .tools import patch_syspath, AbimkdocsTest
patch_syspath()

import os
import json

from collections import OrderedDict
from abimkdocs.variables import get_variables_code, ValueWithUnit, ValueWithConditions, MultipleValue, Range

class VariablesTest(AbimkdocsTest):

    def test_variables(self):
        vars_code = get_variables_code()
        assert vars_code is get_variables_code()
        assert "abinit" in vars_code
        abinit_vars = vars_code["abinit"]
        assert "ecut" in abinit_vars
        ecut = abinit_vars["ecut"]
        assert ecut == ecut
        assert ecut != None
        assert ecut.name == "ecut"
        assert ecut.vartype == "real"
        assert ecut != abinit_vars["pawecutdg"]
        assert ecut in {ecut: "foo"}
        assert ecut.executable == "abinit"
        assert ecut.mdlink == "[[abinit:ecut]]"
        assert not ecut.is_internal
        assert abinit_vars["mpw"].is_internal
        abinit_asr = abinit_vars["asr"]
        anaddb_asr = vars_code["anaddb"]["asr"]
        assert abinit_asr.executable == "abinit"
        assert abinit_asr.name == "asr"
        assert anaddb_asr.executable == "anaddb"
        assert anaddb_asr.name == "asr"
        assert anaddb_asr.mdlink == "[[anaddb:asr]]"
        assert anaddb_asr != abinit_asr
        symsigma_parents = abinit_vars["symsigma"].get_parent_names()
        assert len(symsigma_parents) > 0 and "optdriver" in symsigma_parents

        # Test "tricky" variables.
        fxcartfactor = vars_code["abinit"]["fxcartfactor"]
        str(fxcartfactor)
        assert fxcartfactor.to_markdown()
        d = fxcartfactor.topic_tribes
        assert fxcartfactor.topic_tribes is d and len(d) == 2
        assert "expert" in d["TransPath"] and "expert" in d['GeoOpt']
        assert isinstance(fxcartfactor.defaultval, ValueWithUnit)
        assert fxcartfactor.defaultval.units == "(Bohr^2)/Hartree"
        assert fxcartfactor.defaultval.value == 1

        iomode = vars_code["abinit"]["iomode"]
        str(iomode)
        assert iomode.to_markdown()
        assert len(iomode.characteristics) == 1 and iomode.characteristics[0] == "[[DEVELOP]]"
        assert isinstance(iomode.defaultval, ValueWithConditions)
        # TODO: waiting for new version dict-based.
        #assert iomode.defaultval["defaultval"] == 0
        #assert iomode.defaultval['[[MPI_IO]] and [[paral_kgb]]==1'] == 1

        istwfk = vars_code["abinit"]["istwfk"]
        assert istwfk.to_markdown()
        assert istwfk.characteristics is None and istwfk.vartype == "integer" and istwfk.varset == "dev"
        assert istwfk.dimensions == ["[[nkpt]]"]
        assert isinstance(istwfk.defaultval, MultipleValue)
        #assert istwfk.defaultval.number is None
        #assert istwfk.defaultval.value == 0
        assert istwfk.requires is None
        assert "nkpt" in istwfk.get_parent_names()

        jdtset = vars_code["abinit"]["jdtset"]
        assert isinstance(jdtset.defaultval, Range)
        #assert jdtset.defaultval.start == 1
        #assert jdtset.defaultval.stop == "[[ndtset]]"

        errors = []
        for var in vars_code.iter_allvars():
            try:
                assert repr(var)
                assert str(var)
                #assert var.to_string(verbose=2)
                assert var.to_markdown()
                # topic --> list of tribes
                assert len(var.topic_tribes)
                var.validate()
            except Exception as exc:
                errors.append(str(exc))

        if errors:
            raise ValueError("\n".join(errors))

    def test_variables_in_tests(self):
        """
        Find variables that are not tested by comparing the database of variables
        with the input file of the test.
        Build dictionary with list of untested variables for the different codes.
        Finally compare the new dictionary with the reference one and fail if they don't match.
        """
        # Build database with all input variables indexed by code name.
        from abimkdocs.variables import get_variables_code
        variables_code = get_variables_code()

        # Build AbinitTestSuite object.
        from doc import tests as tmod
        tests = tmod.abitests.select_tests(suite_args=[], regenerate=True, flat_list=True)

        # Build conter for the different codes, keys are the varnames from the database.
        from collections import Counter
        count_code = {}
        for code, d in variables_code.items():
            count_code[code] = Counter({k: 0 for k in d})

        black_list = set([
            "atompaw", "cut3d", "multibinit", "fftprof", "conducti", "mrgscr",
            "mrgddb", "mrggkk", "mrgdv", "band2eps", "ujdet", "fold2Bloch", "macroave",
        ])
        for test in tests:
            if test.executable in black_list: continue
            #print(test)
            vnset = test.get_varname_set()
            #print(vnset)
            count_code[test.executable].update(vnset)

        untested = OrderedDict()
        for code, count in count_code.items():
            untested[code] = []
            # Add it if var is not tested and not internal.
            for vname, c in count.items():
                if c == 0 and not variables_code[code][vname].is_internal:
                    print(code, vname)
                    untested[code].append(vname)

        for code in sorted(untested.keys()):
            untested[code] = sorted(untested[code])
            print(code, untested[code])
        #assert 0

        ref_json_path = os.path.join(os.path.dirname(__file__), "untested_variables.json")
        update_ref = False
        if update_ref:
            with open(ref_json_path, "wt") as fh:
                json.dump(untested, fh, indent=4)
        else:
            with open(ref_json_path, "rt") as fh:
                ref_untested = json.load(fh)

            self.assertDictEqual(untested, ref_untested,
                msg="Detected mismatch between reference file and new list of untested variables.")
