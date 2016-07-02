!! this file is aimed at beeing included into functions
!! having some dummy functions as argument

interface
  function dp_dum_v2dp(nv1,nv2,arg)
    use defs_basis
    use defs_datatypes
    integer, intent(in) :: nv1,nv2
    real(dp),intent(in) :: arg(nv1,nv2)
    real(dp)::dp_dum_v2dp
  end function dp_dum_v2dp
end interface


interface
  function v2dp_dum_v2dp(nv1,nv2,arg)
    use defs_basis
    use defs_datatypes
    integer, intent(in) :: nv1,nv2
    real(dp),intent(in) :: arg(nv1,nv2)
    real(dp)            :: v2dp_dum_v2dp(nv1,nv2)
  end function v2dp_dum_v2dp
end interface


interface
  subroutine sub_dum_dp_v2dp_v2dp(nv1,nv2,inarg1, inarg2, ioarg3)
    use defs_basis
    use defs_datatypes
    integer, intent(in)    :: nv1,nv2
    real(dp),intent(in)    :: inarg1
    real(dp),intent(inout) :: inarg2(nv1,nv2)
     real(dp),intent(inout):: ioarg3(nv1,nv2)
  end subroutine sub_dum_dp_v2dp_v2dp
end interface

!! end of dummy_functions.inc

