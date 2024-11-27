"""
Implement a class used to analyze some data from ABINIT .abo file.
Can be used to count datasets, extract number of iterations...
"""
from __future__ import print_function, division, unicode_literals
from math import ceil, floor

# ---------------------------------------------------------------------

class AboFileAnalysis(object):
    """Main object containing data from abo file (for furher analysis)."""

    def __init__(self,file_name,option):
        """
        Arguments:
        file_name: name of the file to analyze
        option (string): which type of information do we extract from the file (possible values: "iterations")
        """
        self.file_type = "abo"
        self.file_name = file_name
        self.option = option

        if option != "":
            self.dtsets = self.extract(option=option)

#   -------------

    def extract(self,option):
        """
        Extract data from the abo file, dataset per dataset
        Argument:
        option (string): which type of information do we extract from the file
             possible values: "iterations" = extract the number of iterations of all cycles
        """

        dataset_list = []
        abo_lines = open(self.file_name, "rt").readlines()

        # Loop over file lines (loop over datasets)
        inDatasetMode = False
        for i, line in enumerate(abo_lines):

            # Open dataset mode
            if line.startswith("== DATASET"):
                if inDatasetMode:
                    dataset_list.append(current_dataset)
                    del(current_dataset)
                else:
                    inDatasetMode = True
                current_dataset = AboDataset(int(line.split()[2]))
            # Close dataset mode
            elif line.startswith("== END DATASET(S)"):
                inDatasetMode = False
                dataset_list.append(current_dataset)
                del(current_dataset)
  
            # Read data from current dataset
            else:

                # Read optdriver
                if "meta: {optdriver:" in line:
                    current_dataset.optdriver = int(line.split()[2].split(',')[0])

                if "iterations" in option:
                    # Read MD iteration number
                    # Look for:
                    #    At Broyd/MD step X, gradients are converged
                    if "At Broyd/MD step" in line and "converged" in line:
                        current_dataset.MD_niter = int(line.split()[3].split(',')[0])
                    # Look for:
                    #    ntime= X was not enough Broyd/MD steps to converge gradients
                    if "ntime" in line and "was not enough Broyd/MD steps" in line:
                        current_dataset.MD_niter = int(line.split()[1])
	
                    # Read SCF iteration number
                    # Look for:
                    #    At SCF step X, etot is converged
                    #    At SCF step X, forces are converged
                    #    At SCF step X, forces are sufficiently converged
                    #    At SCF step X,        nres2 [...] =>converged
                    #    At SCF step X, max residual [...] =>converged
                    #    At SCF step X, max grdnorm  [...] =>converged
                    if "At SCF step" in line and "converged" in line:
                        current_dataset.SCF_niter.append(int(line.split()[3].split(',')[0]))
                    # Look for:
                    #    At SCF step X, the difference between
                    #    is converged :  diff(etot_el-etot_pos)=
                    if "At SCF step" in line and "converged" in abo_lines[i+1]:
                        current_dataset.SCF_niter.append(int(line.split()[3].split(',')[0]))
                    # Look for:
                    #    nstep= X was not enough SCF cycles to converge;
                    #    nstep= X was not enough non-SCF iterations to converge;
                    if "nstep=" in line and "was not enough" in line:
                        current_dataset.SCF_niter.append(int(line.split()[1]))

        # Debug
        # print([[j.number,j.optdriver,j.MD_niter,j.SCF_niter] for j in dataset_list])            

        return dataset_list

#   -------------

    def compare_with(self,other_abo_file,option,percent_allowed_small=0,percent_allowed_large=0):
        """
        Compare the current abo file with another one
        Compare only specific parts specified by argument option (string)
        Arguments:
          option: what do we compare (possible values: "iterations")
          percent_allowed_small: percentage allowed for a change in the number of iterations
                                 for small numbers of iterations (n_iter<=8)
          percent_allowed_large: percentage allowed for a change in the number of iterations
                                 for large numbers of iterations (n_iter>8)
        """

        status = "succeeded"
        err_msg = "" ; err_msg_short = ""
        tol_small = float(percent_allowed_small)/100.
        tol_large = float(percent_allowed_large)/100.

        if status == "succeeded":
            if other_abo_file is None:
                status = "failed"
                raise ValueError("BUG: no abo file provided for the diff!")
          
        if status == "succeeded":
            if len(self.dtsets) != len(other_abo_file.dtsets):
                status = "failed"
                raise ValueError("ERROR: the two abo files have different dataset numbers!")

        if status == "succeeded":
            if "iterations" in option and "iterations" in self.option:                

                for i, dtset1 in enumerate(self.dtsets):
                    dtset2 = other_abo_file.dtsets[i]
                    jdt = dtset1.number

                    if dtset1.MD_niter is not None and dtset2.MD_niter is not None:
                        tol = tol_small if dtset1.MD_niter<=8 else tol_large
                        if dtset2.MD_niter > ceil(dtset1.MD_niter*(1.+tol)) or dtset2.MD_niter < floor(dtset1.MD_niter*(1.-tol)):
                            status = "failed"
                            err_msg = err_msg+'\n' if err_msg != "" else ""
                            err_msg += "Dataset %d, # of MD/relax iterations differs by more than %d%%!" % (jdt,int(tol*100))
                            err_msg_short += "(dtset %d, MD/relax cycle)" % (jdt)

                    if len(dtset1.SCF_niter)>0 and len(dtset2.SCF_niter)>0:
                        ncycle = len(dtset1.SCF_niter)	
                        for it, niter1 in enumerate(dtset1.SCF_niter): 
                            niter2 = dtset2.SCF_niter[it]
                            tol = tol_small if niter1<=8 else tol_large
                            if niter2 > ceil(niter1*(1.+tol)) or niter2 < floor(niter1*(1.-tol)):
                                status = "failed"
                                err_msg = err_msg+'\n' if err_msg != "" else ""
                                if ncycle == 1:
                                    err_msg += "Dataset %d, # of [non-]SCF iterations differs by more than %d%%!" % (jdt,int(tol*100))
                                    err_msg_short += "(dtset %d, SCF_iter)" % (i)
                                else:
                                    err_msg += "Dataset %d, MD/relax cycle %d, # of [non-]SCF iterations differs by more than %d%%!" % (jdt,it+1,int(tol*100))
                                    err_msg_short += "(dtset %d, MD/relax cycle %d, SCF_iter)" % (jdt,it+1)

        return status,err_msg,err_msg_short

# ---------------------------------------------------------------------

class AboDataset(object):
    """Object storing data extracted from ABINIT abo file for ONE dataset."""

    def __init__(self,number):
        self.number = number
        self.optddriver = 0
        self.MD_niter = None
        self.SCF_niter = [] # This is a list because several SCF can occur in a dataset
