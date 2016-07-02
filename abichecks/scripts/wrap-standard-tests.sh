#
# Wrapper for the standard tests of abinit
#

# Init
my_name="wrap-standard-tests"

# The leading './' is essential!
my_cnffile="./abichecks.env"

# Check arguments
if test "${#}" -lt "1"; then
  echo "Usage: ${my_name} machine_name [ start_test [ stop_test ] ]"
  echo ""
  exit 0
fi

# Set-up environment
if test -s "${my_cnffile}"; then
 . ${my_cnffile}
else
 echo "${my_name}: ${my_cnffile} not found - aborting now"
 exit 1
fi

# Save log if machine_name is chkinabi
machine_name="${1}"
shift
if test "${machine_name}" = "chkinabi"; then
 my_logfile='>& tmp-chkinabi.log'
else
 my_logfile=''
fi

${PERL} ${abinit_inpdir}/scripts/run-standard-tests.pl \
  ${machine_name} $@ ${my_logfile}
exit ${?}
