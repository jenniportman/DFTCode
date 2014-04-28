program shiftedLayers
implicit none

character(20) :: fileIn, fileOut
character(50) :: string

integer, dimension(2) :: nAtoms
integer :: nLayers=3, totalAtoms
integer :: i, j, k, n
integer, allocatable :: seed(:)

real :: posA
real, dimension(3) :: a,b,c, newPosition
real, dimension(3) :: shiftVector1, shiftVector2
real, allocatable :: atomPositions(:,:)

write(fileIn,*) 'POSCAR_shifted.vasp'
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

! initialize random number generator
call random_seed(size = n)
allocate(seed(n))
seed(:)=124576
call random_seed(put=seed)
deallocate(seed)

! generate shift vectors for layers 2,3
call random_number(shiftVector1)
shiftVector1=shiftVector1/2.0
shiftVector1(3)=1.0/real(nLayers)
call random_number(shiftVector2)
shiftVector2=shiftVector2/2.0
shiftVector2(3)=2.0/real(nLayers)
write(*,*) "Shift of 2nd layer", shiftVector1
write(*,*) "Shift of 3nd layer", shiftVector2

! loop over atom types
do i=1,2
	! first layer = no change
	do j=1,nAtoms(i)
		k=j+(i-1)*nAtoms(1)
		write(35,*) atomPositions(k,:)
	end do

	! second layer = random shift of whole layer
	do j=1,nAtoms(i)
		k=j+(i-1)*nAtoms(1)
		newPosition(:)=atomPositions(k,:)+shiftVector1(:)
		newPosition(:)=newPosition(:)-nint(newPosition(:))
		write(35,*) newPosition(:)
	end do

	! third layer = random shift of whole layer
	do j=1,nAtoms(i)
		k=j+(i-1)*nAtoms(1)
		newPosition(:)=atomPositions(k,:)+shiftVector2(:)
		newPosition(:)=newPosition(:)-nint(newPosition(:))
		write(35,*) newPosition(:)
	end do

end do

deallocate(atomPositions)
close(30)
close(35)

end program shiftedLayers
