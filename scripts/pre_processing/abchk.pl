#!/usr/bin/perl

# this program reads through an abinit input file, and checks to see
# whether each variable encountered can be matched to a variable in the
# abinit input list. The purpose is to avoid input errors of the type
# "rfstr" instead of the correct "rfstrs"; since abinit has defaults for
# many variables, if it reads an unknown name it seems to just ignore the
# input and continue on with the default.

# usage: abchk <input file name>
# output: variables in the input file that didn't match the glossary
# program returns the number of unmatched variables, which makes it
# also useful in shell scripts

# By: Josef Zwanziger, Dept. of Chemistry and Institute for Research
# in Materials, Dalhousie University, Halifax, NS B3H 4J3
# jzwanzig@dal.ca

# have fun!!


# here are many of the abinit variables. You can add more.
@abinit_variables = ("acell","angdeg","amu","algalch",
                     "brvltt","bdberry","berryopt","boxcenter","boxcutmin",
                     "cmlfile","chkexit","chkrprim","cpus","cpum","cpuh",
                     "delayperm","dilatmx","dtion","dsikpt","diecut",
		     "diegap","dielam","dielng","diemac","diemix","dosdeltae",
		     "d3e_pert1_atpol","d3e_pert1_dir","d3e_pert1_elfd","d3e_pert1_phon",
		     "d3e_pert2_atpol","d3e_pert2_dir","d3e_pert2_elfd","d3e_pert2_phon",
		     "d3e_pert3_atpol","d3e_pert3_dir","d3e_pert3_elfd","d3e_pert3_phon",
                     "ecut","ecutsm","efield","enunit",
                     "friction","fband","spinmagntarget",
		     "getddk","getden","getkss","getocc","getscr",
		     "getwfk","getwfq","get1den","get1wf","get1wfden",
                     "genafm","getcell","getxcart","getxred",
                     "iscf","ixc","iatcon","iatfix","iatfixx","iatfixy",
                     "iatfixz","ionmov","irdddk","irdkss","irdscr","irdwfk",
		     "irdwfq","ird1wf","iatsph","iprcel",
                     "jdtset",
                     "kpt","kptnrm","kptopt","kssform","kberry","kptbounds",
		     "kptrlatt","kptrlen",
		     "localrdwf",
                     "mdftemp","mditemp","mdwall","mffmem","mkmem",
		     "mkqmem","mk1mem","mixalch",
                     "natfix","natfixx","natfixy","natfixz","natcon",
                     "nconeq","ntime",
                     "natom","nband","ndtset","ngkpt","nkpt","nshiftk",
		     "natsph","nberry","nbdbuf","ndivk","ngfft","nline",
		     "npsp","nqpt","nspinor","ntypalch",
                     "nspol","nstep","nsym","ntypat","natrd","nobj",
                     "occopt","objatt","objbat","objaax","objbax","objan",
                     "objbn","objarf","objbrf","objaro","objbro","obfatr",
                     "objbtr","optcell","occ","optdriver",
		     "prtcml","prtden","prtdos","prteig","prtfsurf","prtgeo",
		     "prtkpt","prtpot","prtstm","prtvha","prtvhxc",
		     "prtvol","prtvxc","prtwf","prt1dm","prepanl","prtbbb",
		     "pspso",
		     "qpt","qptnrm",
                     "rprim","restartxf","rfasr","rfatpol","rfdir","rfelfd",
		     "rfphon","rfstrs","rfuser","ratsph",
                     "shiftk","symrel","spgaxor","spgorig","spgroup",
                     "spgroupma","signperm","strfact","strprecon",
                     "strtarget","sciss","so_typat","spinat","stmbias",
		     "symafm",
                     "tnons","toldfe","toldff","tolvrs","tolwfr","typat",
                     "tolmxf","td_maxene","td_mexcit","timopt",
		     "tphysel","tsmear",
                     "udtset",
                     "vel","vis","vacuum","vacwidth",
                     "wtk","wtatcon",
                     "xangst","xcart","xred",
                     "znucl");

#now set number of unknown variables to 0.
$unknown_variables = 0;
while (<>) { # loop over lines in the input file...
      @fields = split; # split each line on white space...
      if ( @fields[0] =~ /^[a-z]/ ) { # check for a match on lowercase characters. These are variables
         while( @fields[0] =~ /[0-9]\b/ ) { chop(@fields[0]); } # if the trailing character is a number, 
                                                                # strip it off. These are the variables in multiple data sets
	 $known = 0;
	 foreach $abinit_variable (@abinit_variables) { # if the variable matches a glossary entry EXACTLY, it's known.
		 if (@fields[0] =~ /\b$abinit_variable\b/) { $known = 1; }
	}
	if(!$known) { #if the variable is not known, tell us about it. 
		    print "do not know variable @fields[0]\n";
		    ++$unknown_variables;
		    };
	}
}
exit($unknown_variables); # finally, tell us the total number of unknown variables
