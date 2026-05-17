"""
Implement a class used to analyze some data from ABINIT .abo file.
Can be used to count datasets, extract number of iterations...

This module provides the AboFileAnalysis class which can be used to count datasets,
extract iteration counts (SCF and MD), and compare different output files.
"""
from math import ceil, floor


class AboFileAnalysis:
    """
    Main object containing data from abo file (for further analysis).

    This class parses an ABINIT output file, identifies datasets, and extracts
    relevant information like SCF and MD iteration counts. It also provides
    functionality to compare iteration counts between two files.

    Attributes:
        file_type (str): Type of file (always "abo").
        file_name (str): Path to the file being analyzed.
        option (str): Extraction option string passed during initialization.
        dtsets (list[AboDataset]): List of AboDataset objects containing the extracted information.
    """

    def __init__(self, file_name, option):
        """
        Initialize the AboFileAnalysis object with a file and an extraction option.

        Args:
            file_name (str): Path to the ABINIT output file to analyze.
            option (str): Type of information to automatically extract.
                Pass an empty string to skip automatic extraction.
        """
        self.file_type = "abo"
        self.file_name = file_name
        self.option = option

        if option != "":
            self.dtsets = self.extract(option=option)

    def extract(self, option):
        """
        Extract data from the abo file, dataset per dataset.

        Args:
            option (str): Type of information to extract.
                Possible values: "iterations" = extract iterations of all cycles.

        Returns:
            list: List of AboDataset objects containing the extracted information.
        """
        dataset_list = []
        abo_lines = open(self.file_name).readlines()

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
                    current_dataset.optdriver = int(line.split()[2].split(",")[0])

                if "iterations" in option:
                    # Read MD iteration number
                    # Look for:
                    #    At Broyd/MD step X, gradients are converged
                    if "At Broyd/MD step" in line and "converged" in line:
                        current_dataset.MD_niter = int(line.split()[3].split(",")[0])
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
                        current_dataset.SCF_niter.append(int(line.split()[3].split(",")[0]))
                    # Look for:
                    #    At SCF step X, the difference between
                    #    is converged :  diff(etot_el-etot_pos)=
                    if "At SCF step" in line and "converged" in abo_lines[i+1]:
                        current_dataset.SCF_niter.append(int(line.split()[3].split(",")[0]))
                    # Look for:
                    #    nstep= X was not enough SCF cycles to converge;
                    #    nstep= X was not enough non-SCF iterations to converge;
                    if "nstep=" in line and "was not enough" in line:
                        current_dataset.SCF_niter.append(int(line.split()[1]))

        # Debug
        # print([[j.number,j.optdriver,j.MD_niter,j.SCF_niter] for j in dataset_list])

        return dataset_list

    def compare_with(self,other_abo_file,option,percent_allowed_small=0,percent_allowed_large=0):
        """
        Compare the current abo file with another one.
        Compare only specific parts specified by argument option.

        This method is primarily used in the test suite to ensure that changes in the code
        don't significantly affect the number of iterations required for convergence.

        Args:
            other_abo_file (AboFileAnalysis): The other AboFileAnalysis instance to compare with.
            option (str): What to compare. Only "iterations" is currently supported.
            percent_allowed_small (int): Percentage tolerance for small iteration counts (n_iter <= 8).
                Defaults to 0.
            percent_allowed_large (int): Percentage tolerance for large iteration counts (n_iter > 8).
                Defaults to 0.

        Returns:
            tuple: (status, err_msg, err_msg_short) where:
                status (str): "succeeded" if counts are within tolerance, "failed" otherwise.
                err_msg (str): Detailed multiline error message listing discrepancies.
                err_msg_short (str): Succinct summary of where discrepancies occurred.

        Raises:
            ValueError: If other_abo_file is None or if files have different numbers of datasets.
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
                print ("2 lengths of dtsets = ", len(self.dtsets), len(other_abo_file.dtsets))
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
                            err_msg = err_msg+"\n" if err_msg != "" else ""
                            err_msg += "Dataset %d, # of MD/relax iterations differs by more than %d%%!" % (jdt,int(tol*100))
                            err_msg_short += "(dtset %d, MD/relax cycle)" % (jdt)

                    if len(dtset1.SCF_niter)>0 and len(dtset2.SCF_niter)>0:
                        ncycle = len(dtset1.SCF_niter)
                        for it, niter1 in enumerate(dtset1.SCF_niter):
                            niter2 = dtset2.SCF_niter[it]
                            tol = tol_small if niter1<=8 else tol_large
                            if niter2 > ceil(niter1*(1.+tol)) or niter2 < floor(niter1*(1.-tol)):
                                status = "failed"
                                err_msg = err_msg+"\n" if err_msg != "" else ""
                                if ncycle == 1:
                                    err_msg += "Dataset %d, # of [non-]SCF iterations differs by more than %d%%!" % (jdt,int(tol*100))
                                    err_msg_short += "(dtset %d, SCF_iter)" % (i)
                                else:
                                    err_msg += "Dataset %d, MD/relax cycle %d, # of [non-]SCF iterations differs by more than %d%%!" % (jdt,it+1,int(tol*100))
                                    err_msg_short += "(dtset %d, MD/relax cycle %d, SCF_iter)" % (jdt,it+1)

        return status,err_msg,err_msg_short

class AboDataset:
    """
    Object storing data extracted from ABINIT abo file for ONE dataset.

    Attributes:
        number (int): Dataset number.
        optdriver (int): Optimization driver value.
        MD_niter (int): Number of MD/relax iterations.
        SCF_niter (list): List of SCF iteration counts.
    """

    def __init__(self, number):
        """
        Initialize an AboDataset instance.

        Args:
            number (int): The dataset number.
        """
        self.number = number
        self.optddriver = 0
        self.MD_niter = None
        self.SCF_niter = [] # This is a list because several SCF can occur in a dataset
