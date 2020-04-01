TODO list for the documentation of ABINIT
=========================================

The first version of this TODO was written by Yann.
Several points have been adressed in the mkdocs version but there
are still points worth discussing

* General remarks:

    * In prevision of the Autotools support and packaging, the documentation
    should be reorganized in relevant sections => needs discussion.
    * Copyright and history information might be inserted as comments.
    * To be discussed: which files should exist in one format only, and which
    ones should be multi-format?
    * Dependencies should be stored somewhere, to known what to update and when.

* aim_help.html:

    * Is it up-to-date?

* aimhelp.tex,aimhelp.ps:

    * Do not match aim_help.html => which is the right one?

* anaddb_help.html:

    * Is it up-to-date?

* band2eps_help:

	* Can barely be called documentation (and executable) => restart from scratch.

* check_list:

	* To whom is it truly useful? 
        MG: Xavier.
        @Xavier: Could you move the file somewhere>

* conducti_manual.tex:

	* Is it complete and accurate?
	* Is it up-to-date?


* contributors:

	* Format better the file.
	* Remove copyright notice.

* cut3d_help:

	* Is it complete and accurate?
	* Is it up-to-date?

* cut3d_help.html:

	* Does not match cut3d_help => which is the right one?
	* Remove links at bottom.
	* Remove copyright notice.

* ddbs.upgrade:

	* Name: what about ddbs-upgrade-howto?
	* Reformat for markdown.

* elphon_manual.ps:

	* Violation of the GPL: where is the source? Should be relicensed?

* format_KSS:

	* Is it complete and accurate?
	* Is it up-to-date?
	* Reformat for markdown.

* gwmanual.txt:

	* Is it complete and accurate?
	* Is it up-to-date?
	* Reformat for markdown.

* known_problems.x.y.z:

	* Should be generated automatically from a database.
	* Should be presented in reverse-chronological order.

* make_help, make_dev_help, make_targz_help, makefile_macros_help:

	* Are they up-to-date?
	* Will be changed by the Autotools support.

* mrgddb_help.html:

	* Barely formatted.
	* Prepare text version for markdown.
	* Remove links at bottom.
	* Remove copyright notice.

* new_user_guide.html:

	* Is it complete and accurate?
	* Is it up-to-date?

* other_codes:

	* Probably out-of-date.
	* Set-up a database?

* paral_use:

	* Is it complete and accurate?
	* Is it up-to-date?

* piezoelectric.txt:

	* Only notes for now, will evolve a lot.
	* Reformat for markdown.

* planning:

	* Is it complete and accurate?
	* Is it up-to-date?
	* Reformat for markdown.
	* Remove copyright notice.

* problem_report:

	* Is it up-to-date?
	* Reformat for markdown.
	* Remove copyright notice.

* respfn_help.html:

	* Is it complete and accurate?
	* Is it up-to-date?

* spacegrouphelpfile.html, spgrdefinition.html, spgrdescription.html, spgrhead.html, spgrcopyright.html:

	* Should be rewritten. DONE by MG

* tuning:

	* Is it complete and accurate?
	* Is it up-to-date?

* welcome.html:

	* Is it complete and accurate?
	* Is it up-to-date?
	* Prepare text version for markdown (may take time).
	* Remove links at bottom.
	* Remove copyright notice.

* Features/features*:

	* File extensions vary => choose only one.
	* Barely formatted, almost plain-text.
	* Reformat for markdown.
	* Remove copyright notice.

* Images/*:

	* Only images for the tutorial => nothing to do.

* Installation_notes/install*.htm*:

	* File extensions vary => choose only one.
	* The source is not even valid HTML !!!
	* Improve layout.
	* Prepare text version for markdown (may take time).
	* Remove links at bottom.
	* Remove copyright notice.

* Macroave_Docs/README:

	* Ready for markdown.
	* Remove copyright notice.

* Macroave_Docs/macroave.(tex,ps,pdf):

	* Are they complete and accurate?
	* Are they up-to-date?

* Macroave_Docs/macroave.toc:

	* Temporary file => should be removed.

* Miscellaneous/*:

	* Are they useful?
	* Are they complete and accurate?
	* Are they up-to-date?
	* Reformatting: discuss types and priorities.

* Notes_for_coding/*:

	* Are they complete and accurate?
	* Are they up-to-date?
	* Should be reformatted to allow multi-format presentation.

* Presentation/*.tex:

	* Are they complete and accurate?
	* Are they up-to-date?

* Presentation/presentation.pdf.gz:

	* Small file anyway => decompress it.

* Psp_Infos/*.info:

	* All files should be reformatted in order to produce both HTML
	  and PDF versions of the information.

* Release_notes/release_notes*

	* File extensions vary => choose only one.
	* Plain-text version required by GNU coding standards.
	* Barely formatted, almost plain-text.
	* Reformat for markdown.
	* Remove copyright notice.

* Theory/*:

	* Are they complete and accurate?
	* Are they up-to-date?
	* Find more explicit names.

* Tutorial/*.html:

	* welcome.html: Rename it to something like summary.html?
	* The source is not even valid HTML !!! => fix it.
	* Upgrade to XHTML.
	* Should be part of a PDF manual too.
	* Remove links at top and bottom.
	* Remove copyright notice.


* README_conda with ${CONDA_ENV}
