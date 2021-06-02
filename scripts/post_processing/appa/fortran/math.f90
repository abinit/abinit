!------------------------------------!
!---------AVERAGE ROUTINE------------!
!------------------------------------!

subroutine average(data,size_data,average_value)

! input values:  
  integer :: size_data
  real(8),dimension(size_data),intent(in) :: data
 !f2py optional , depend(data) :: size_data=len(data)
  real(8),intent(out) :: average_value
  average_value  = 0.0

  do i=1,size_data
    average_value  = average_value + data(i)
  end do

  average_value = average_value / size_data
  
end subroutine average


!------------------------------------!
!------STANDAR DEVIATION ROUTINE-----!
!------------------------------------!
subroutine standard_deviation(data,average_input,size_data,deviation)
!  input  values:    
  integer  ::  size_data
  real(8), dimension(size_data),intent(in)  ::  data
  real(8), intent(in), optional  ::  average_input  
  !f2py  optional  ,  depend(data)  ::  size_data=len(data)
  real(8),intent(out)  ::  deviation
  real(8) :: average_data = 0

  
  if (average_input/=0) then
    average_data = average_input
  else
    average_data = 0
    call  average(data,size_data,average_data)
  end if

  do  i=1,size_data
    deviation  =  deviation  +  (data(i)  -  average_data)**2
  end  do
    
  deviation  =  deviation/size_data
  deviation  =  deviation**0.5

end subroutine  standard_deviation


!------------------------------------!
!------MEAN SQUARED DISPLACEMENT-----!
!------------------------------------!
subroutine mean_square_displacement(position,indexAtom1,msd,size2,size3,nba1)
  integer :: size2,size3,nba1,nbtime
  !input  values:    
  !f2py  optional  ,  depend(position)    ::  size2=len(position[:,1,1])
  !f2py  optional  ,  depend(position)    ::  size3=len(position[1,:,1])
  !f2py  optional  ,  depend(indexAtom1)  ::  nba1=len(indexAtom1)
  integer, intent(in)  ::  indexAtom1(nba1)
  real, intent(in)  ::  position(size2,size3,3)
  !output 
  real(8),intent(out)  ::  msd((size2-1)/2)

  !local variable
  real(8), dimension(3)  ::  p1=(/1,1,1/),p2=(/1,1,1/)
  real(8) r
  integer i,k,tau

  nbtime = size2

  msd(:) = 0
  do tau=1,int((nbtime-1)/2) !loop on step
    do i=1,nbtime-tau-1 !loop on step
      do k=1,size(indexAtom1)
        p1 = position(i,indexAtom1(k),:)
        p2 = position(i+tau,indexAtom1(k),:) - p1(:)
        r = p2(1)**2+p2(2)**2+p2(3)**2
        msd(tau)= msd(tau) + r
      end do
    end do
    msd(tau)= msd(tau) / (nbtime-tau-1)
  end do
  msd(:) = msd(:)/(size(indexAtom1))
end subroutine mean_square_displacement
