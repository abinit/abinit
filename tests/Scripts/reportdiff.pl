sub reportdiff {
	local($fn1,$fn2,$fn3) = @_;
#
# purpose : examine $fn1, using the reference values in $fn2, and echo result in $fn3
#
	$rc = open(FILEIN,"$fn1");
	if ($rc eq '') {
	 print "Unable to open input file $fn1\n";
	 }
	$rc = open(FILEREF,"$fn2");
	if ($rc eq '') {
	 print "Unable to open input file $fn2\n";
	 }
	$rc = open (FILEOUT,">>$fn3");
	if ($rc eq '') {
	 print "Unable to open input file $fn2\n";
	 }
#       Read both files, one line at a time
	$lnm1 = 0;
	$lnm2 = 0;
	while (1) {
	  $lin1 = <FILEIN>  ;   # read next line from file 1
	  $lin2 = <FILEREF> ;   # read next line from file 2
#
# ignore comment in report.in
	  while ( $lin2 =~ '^#' ) {
		$lin2 = <FILEREF> ;
	  }
#
	  last if($lin1 eq '' && $lin2 eq ''); # end of both files
	  @field1 = split(' ',$lin1);       # parse line into array
	  @field2 = split(' ',$lin2);
	  if ($field1[0] ne 'Summary') {
	    print FILEOUT "Error : all lines of file to be analysed must start with 'Summary' " ;
	    print FILEOUT "Analysed file  $fn1, line $lnm1, gives :" ;
	    print FILEOUT "$lin1" ;
	    print FILEOUT "$field1[0]" ;
	    print FILEOUT "exit" ;
	    exit ;
	    }
	  if ($field1[1] ne $field2[0]) {
	    print FILEOUT "Error : indices of cases do not match in analysed and reference files " ;
	    print FILEOUT "Analysed file  $fn1, line $lnm1, gives $field1[1] ;" ;
	    print FILEOUT "Reference file $fn2, line $lnm2, gives $field2[0] ;" ;
	    print FILEOUT "exit" ;
	    exit ;
	    }
	  if ($field1[3] eq "fatal") {
	    print FILEOUT "$field1[1]            failed" ;
	    }
	  elsif ($field1[3] eq "no") {
	    print FILEOUT "$field1[1] succeeded" ;
	    }
	  elsif ($field1[3] eq "different") {
	    $tolnlines = $field2[2] ;
	    $tolabs = $field2[4] ;
            $tolrel = $field2[6] ;
	    if ( $field1[9] > $tolabs*1.5 ) {
              print FILEOUT "$field1[1]            failed (too large absolute error : $field1[9] , accepted $tolabs ) " ;
              }
            elsif ( $field1[13] > $tolrel*1.5 ) {
              print FILEOUT "$field1[1]            failed (too large relative error : $field1[13] , accepted $tolrel ) " ;
              }
            elsif ( $field1[5] > $tolnlines ) {
              print FILEOUT "$field1[1]            failed (too many erroneous lines : $field1[5] , accepted $tolnlines ) " ;
              }
            elsif ( $field1[9] > $tolabs ) {
              print FILEOUT "$field1[1] within 1.5 of tolerance (absolute error : $field1[9] , accepted $tolabs ) " ;
              }
            elsif ( $field1[13] > $tolrel ) {
              print FILEOUT "$field1[1] within 1.5 of tolerance (relative error : $field1[13] , accepted $tolrel ) " ;
              }
	    else {
	      print FILEOUT "$field1[1] passed" ;
	      }
	    }
	  }
	close (FILEIN);
	close (FILEREF);
	close (FILEOUT);
	return;                 # go read next line from configuration file
	}
