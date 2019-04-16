-- ================================================================================================ --
--                                         Basic parameters                                         --
-- ================================================================================================ --
-- Template of configuration file for oneview module

-- Name of the binary file to analyze
binary         = "/home/bsataric/ABinit/ABinit/abinit/tmp/src/98_main/abinit"

-- List of external libraries to analyze
external_libraries = {
-- "lib.so", "lib.so"
}
-- Path to the dataset directory used by the application
-- Use an empty string or remove the declaration if no 
-- dataset is required by the application
dataset        = ""

-- How dataset is handled. "link" (default) to create a symbolic link
-- from the experiment directory to the specified dataset, "copy" to
-- duplicate the dataset directory into the experiment directory
dataset_handler= "link"

-- Command used to run the binary. 
-- + <binary> will be replaced by the path to the binary file to run 
-- If some parameters are used, specify them in the string.
--   example: "<binary> -n4"
-- Use an empty string or remove the declaration if no 
-- specific command is required to run the application
run_command    = "<binary> < t03_x.files"

-- Script to use with job scheduler. 
-- If your application must be run on a cluster using a job scheduler (ex. SLURM), fill
-- this field with the path to a script which will be run using
-- the command detailed in 'batch_command'. The script must have been modified to replace
-- the application executable and its arguments with keyword <run_command>.
-- This configuration file must be filled with options to run your application.
-- The number of processes to use can be referred as <number_processes>.
-- The number of nodes to use can be referred as <number_nodes>.
-- The number of tasks per node to use can be referred as <number_tasks_nodes>.
-- It can be set to an empty string or removed if it not used.
-- batch_script    = ""

-- Command to use to submit job to the job scheduler. 
-- If your application must be run on a cluster using a job scheduler (ex. SLURM), fill
-- this field with the submission command. The script passed to the job scheduler
-- must be referred as <batch_script>.
-- Example using sbatch command (SLURM): batch_command = "sbatch <batch_script>"
-- batch_command   = ""

-- Number of MPI processes.
-- The variable can be referred as <number_processes> in field 'mpi_command' and
-- in the batch script.
number_processes = 1

-- Number of nodes to uses in an MPI application.
-- The variable can be referred as <number_nodes> in field 'mpi_command' and
-- in the batch script.
-- number_nodes = 1

-- Number of tasks per node to uses in an MPI application.
-- The variable can be referred as <number_tasks_nodes> in field 'mpi_command' and
-- in the batch script.
-- number_tasks_nodes = 1

-- Part of the command used to run the binary. 
-- It will be added at the beginning of the final command. It is used to specify
-- a launcher for the application, such as mpirun for MPI applications
--   example: "mpirun -n <number_processes>"
-- Use an empty string or remove the declaration if no 
-- specific command prefix is required to run the application
mpi_command    = "mpirun -n <number_processes> -host csl01"

-- Define the corresponding OpenMP variable to set the maximal number of threads.
-- If the application uses OpenMP, the value must be greater than 1.
-- If the application does not use OpenMP, the value can be equal to 1, nil
-- or the field can be removed from the configuration file.
-- omp_num_threads= 1

-- Directory where the binary must be run
-- + <dataset> will be replaced by the dataset directory located into the
--   experiment directory.
-- Use an empty string or remove the declaration if no 
-- specific directory is required to run the application
run_directory  = "/home/bsataric/ABinit/ABinit/abinit/tests/tutorial/Input/Work5A"

-- ================================================================================================ --
--                            Filter used to select loops to analyze                                --
-- ================================================================================================ --
-- !! Uncomment the filter you want to use !!
-- !! If no filter is specified, all innermost/single loops are used !!

-- This filter uses the first <value> loops, ordered by coverage
--filter = {
--   type = "number",
--   value = 1,
--}

-- This filter uses all loops whose coverage is greater than <value> (in percentage)
--filter = {
--   type = "coverage",
--   value = 1,
--}

-- This filter uses all loops while the cumulated coverage is lower than <value> (in percentage)
--filter = {
--   type = "cumulated_coverage",
--   value = 1,
--}

-- This filter uses all loops. In this case, the filter table can also be nil
--filter = {
--   type = "all",
--}

-- ================================================================================================ --
--                            Specify when the profiling should start                               --
-- ================================================================================================ --
-- !! Uncomment the table you want to use !!
-- !! If no table is specified, the default table has 'p' as unit and 30 as value !!
-- !! Report ONE always analyzes all loops !!

-- If the profiling should begin when the application is started
--profile_start = {
--   unit = "none",   -- Specify that no delay is needed
--   value = 0,         -- Useless value
--}

-- If the profiling should begin after a given time (in second)
--profile_start = {
--   unit = "s",      -- 's' for 'seconds'
--   value = 0,         -- delay in seconds
--}

-- If the profiling should begin after a given percentage of the application time.
-- A first run of the application is automatically performed to time the application.
--profile_start = {
--   unit = "p",      -- 'p' for 'percentage'
--   value = 0,         -- delay in percentage of the total application time
--}

-- ================================================================================================ --
--                                      Additional parameters                                       --
-- ================================================================================================ --
-- Frequencies to use in some dynamic analysis
-- Use an empty table if no frequency must be used.
-- !! Need permissions to change cpufreq files !!
-- !! Do not specify frequencies if MAQAO can not modify these files !!
-- frequencies    = {3300000,1200000}

-- Table describing additional hardware counters analyzed.
-- Each entry is a list of hardware counters analyzed during a run.
-- Hardware counters must be separated with a comma ','.
-- Either hardware counter codes or names can be used.
-- If a set of hardware counters needs sudo permissions, set sudo_mode at true
-- !! Check that all hardware counters in a same list can be analyzed in a single run !!
--additional_hwc = {
--   {names="<names_1>", sudo_mode = <boolean>},
--   {names="<names_2>"},
--}

-- Exclude some areas (loops or blocks) from the experiment
-- type defines if the area is a loop ("loop") or a basic block ("block")
-- id is the MAQAO internal identifier
-- module represents the binary file containing the area to exclude.
-- The main binary is referred as "binary" (or a nil / empty string).
-- Another source must be referred by its name and must be declared in the external_libraries table.
--excluded_areas = {
--   {type = "loop", id = <MAQAO loop identifier>, module = ""},
--}

-- DECAN analyzes several variants in a single run.
-- If decan_multi_variant is set to false, only one variant will be analyzed per run.
-- The experiment will be longer.
-- decan_multi_variant = false

-- Specify the location of the source code.
-- ONE-View looks for source code in directory specified in debug data. However when the
-- application is not compiled on the same machine than the run, source code is not located
-- in the same directory. This variable can be used to specify the source code of the application
-- source_code_location = ""

-- Specify maximal number of path a loop can have.
-- If a loop has more paths, it will not be analyzed by ONEVIEW. To analyze loops regardless of
-- the number of paths, the value must be 0. Default value is 4.
maximal_path_number = 4

-- Specify if sudo mode can be used during experiments.
-- Sudo mode is needed to analyze some hardware counters in report THREE or change the frequency
-- in reports TWO and THREE.
is_sudo_available = false

-- Specify additional options for LPROF.
-- If not empty, everything specified in this field is passed to LPROF through the command line
lprof_params = ""

-- Specify additional options for VPROF.
-- If not empty, everything specified in this field is passed to VPROF through the command line
vprof_params = ""

-- Specify additional options for DECAN.
-- If not empty, everything specified in this field is passed to DECAN through the command line
decan_params = ""

-- Specify additional options for CQA.
-- If not empty, everything specified in this field is passed to CQA through the CQA context
cqa_params = {}

-- If true, it specify that the binary should not be copied into the experiment directory but keep in place
-- keep_binary_location = true

-- Specify runs parameters to use when the scalability report must be generated.
-- Each entry in the table describes an experiment to run with its number of processes, threads, 
-- and the path to a binary when a specific binary must be used for the couple (nb_threads, nb_processes).
-- If nb_threads or nb_processes is not set, they are considered as equal to 1. If binary is not set, it is 
-- the main binary defined in the experiment parameters.
scalability_params = {
 {nb_processes = 1, nb_threads = 2, binary = nil},
 {nb_processes = 1, nb_threads = 4, binary = nil},

}

