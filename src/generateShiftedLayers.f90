! Read in vasp input file and generate new file with nLayers vertically stacked copies of original unit cell with random displacement in plane
program shiftedLayers
implicit none

integer, parameter :: nLayers=4
! random shift in plane 'R' or stacking on sites 'S'
character(1), parameter :: shiftType = 'S'

character(20) :: fileIn, fileOut
character(50) :: string

integer, dimension(2) :: nAtoms
integer :: i, j, k, n, iLayer, totalAtoms
integer, allocatable :: seed(:)

real :: posA
real, dimension(3) :: a,b,c, newPosition
real, dimension(3,nLayers) :: shiftVector 
real, allocatable :: atomPositions(:,:)
real, dimension(3) :: T2Stacking, T5Stacking, T6Stacking, T0Stacking, T7Stacking

write(fileIn,*) 'POSCAR'
write(fileOut,*) 'POSCAR_layers.vasp'

open(30,file=fileIn,status='old',recl=320)
open(35,file=fileOut,status='unknown',recl=320)

! read first two lines which are useless
read(30,'(a)') string
write(35,'(a)') string
read(30,'(a)') string
write(35,'(a)') string

! lattice constants
read(30,*) a
write(35,*) a
read(30,*) b
write(35,*) b
read(30,*) c
write(35,*) nLayers*c

! names of atoms
read(30,'(a)') string
write(35,'(a)') string

! numbers of atoms
read(30,*) nAtoms
write(35,*) nAtoms*nLayers

! type of coordinates (cartesian or fractional)
read(30,*) string
write(35,*) string
if(string(1:1) == "C") then
	write(*,*) "ERROR!! This program only works with fractional coordinates!!"
	stop
endif

! read atomPositions
totalAtoms=nAtoms(1)+nAtoms(2)
write(*,*) "total number of atoms",totalAtoms
allocate(atomPositions(totalAtoms,3))
do i=1,totalAtoms
	read(30,*) atomPositions(i,:)
end do
! z fractional coordinates get scaled to take into account the layers
atomPositions(:,3)=atomPositions(:,3)/real(nLayers)

! vectors for stacking according to paper by Nakanishi et al.
T0Stacking(:) = 0.0

T2Stacking(:) = atomPositions(12,:) - atomPositions(1,:) 
T2Stacking(:) = T2Stacking - nint(T2Stacking)

T5Stacking(:) = atomPositions(10,:) - atomPositions(1,:) 
T5Stacking(:) = T5Stacking - nint(T5Stacking)

T6Stacking(:) = atomPositions(8,:) - atomPositions(1,:) 
T6Stacking(:) = T6Stacking - nint(T6Stacking)

T7Stacking(:) = atomPositions(11,:) - atomPositions(1,:) 
T7Stacking(:) = T7Stacking - nint(T7Stacking)

shiftVector=0.0
if(shiftType == 'R') then
	! initialize random number generator
	call random_seed(size = n)
	allocate(seed(n))
	seed(:)=780781
	call random_seed(put=seed)
	deallocate(seed)
	call random_number(shiftVector)
	shiftVector=shiftVector/10.0
	! first layer = no change
	shiftVector(:,1) = 0.0
else if(shiftType == 'S') then
	! edit here to customize how layers are built
	! first layer = no change
	shiftVector(:,1) = 0.0
	! second layer
	shiftVector(:,2) = shiftVector(:,1) + T0Stacking(:)
	! third layer
	!shiftVector(:,3) = shiftVector(:,2) + T0Stacking(:)
	shiftVector(:,3) = shiftVector(:,2) + T2Stacking(:)
	! fourth layer
	shiftVector(:,4) = shiftVector(:,3) + T0Stacking(:)
	! fifth layer
	!shiftVector(:,5) = shiftVector(:,4) + T6Stacking(:)
	! sixth layer
	!shiftVector(:,6) = shiftVector(:,5) + T0Stacking(:)
	
else
	write(*,*) "ERROR!! Unknown stacking"
	stop
end if
! adjust z coordinates to shift layers up
!shiftVector(3,2:)=(/ (i,i=1,nLayers-1) /)
!shiftVector(3,2:)=shiftVector(3,2:)/real(nLayers)
shiftVector(3,2)=4.0/20.0
shiftVector(3,3)=0.5
shiftVector(3,4)=0.5+4.0/20.0

do i=1,nLayers
	write(*,*) "Layer ", i, "Shift = ",shiftVector(:,i)
end do

! loop over atom types
do i=1,2
	do iLayer=1,nLayers
		do j=1,nAtoms(i)
			k=j+(i-1)*nAtoms(1)
			newPosition(:)=atomPositions(k,:)+shiftVector(:,iLayer)
			newPosition(:)=newPosition(:)-floor(newPosition(:))
			write(35,*) newPosition(:)
		end do
	end do
end do

deallocate(atomPositions)
close(30)
close(35)
end program shiftedLayers
