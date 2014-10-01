program sumfiles
use ifport
implicit none
integer, parameter :: nFiles=13, nPoints=3500 
character(80) :: inputFile, outputFile, outputFileA, outputFileB, outputFileC

integer :: i, j,k, nColumns
character(30), dimension(10) :: string
character(80) :: command
real, dimension(:,:), allocatable :: totalData
real, dimension(:), allocatable :: dataLine
real, dimension(nPoints) :: energy, densityA, densityB, densityC
real :: sumLine

write(outputFile,*) 'dosp_sum.dat'
open(60,file=outputFile,recl=360)
write(outputFileA,*) 'dosp_sum_A.dat'
open(61,file=outputFileA,recl=360)
write(outputFileB,*) 'dosp_sum_B.dat'
open(62,file=outputFileB,recl=360)
write(outputFileC,*) 'dosp_sum_C.dat'
open(63,file=outputFileC,recl=360)

command="sed -n '2p' dosp.001.dat | wc -w > col.dat"
i = system(command)
open(88,file="col.dat")
read(88,*) nColumns
close(88)
command="rm col.dat"
i = system(command)

allocate(dataLine(nColumns-1))
allocate(totalData(nColumns-1,nPoints))


totalData = 0.0
densityA = 0.0
densityB = 0.0
densityC = 0.0
do i=1,nFiles
	if( i < 10 ) then 
		write(inputFile,'(a,i1,a)') 'dosp.00',i,'.dat'
	else
		write(inputFile,'(a,i2,a)') 'dosp.0',i,'.dat'
	end if
	open(20,file=inputFile,status='old')
	read(20,*) string
	if( i == 1) write(60,*) string

	do j=1,nPoints
		read(20,*) energy(j), dataLine
		! sum densities of all atoms for different orbitals
		totalData(:,j) = totalData(:,j) + dataLine(:)

		! sum over orbitals for atoms of type A, B, C
		if(nColumns == 10) then
			if(i == 1) then
				densityA(j) = densityA(j) + sum(dataLine(:))
			else if( i > 1 .and. i < 8) then
				densityB(j) = densityB(j) + sum(dataLine(:))
			else if( i >= 8 .and. i < 14) then
				densityC(j) = densityC(j) + sum(dataLine(:))
			endif
		else
			sumLine = 0
			do k=1,10
				sumLine=sumLine+dataLine(1+4*(k-1))
			end do
			if(i == 1) then
				densityA(j) = densityA(j) + sumLine
			else if( i > 1 .and. i < 8) then
				densityB(j) = densityB(j) + sumLine
			else if( i >= 8 .and. i < 14) then
				densityC(j) = densityC(j) + sumLine
			endif
		endif
	end do
	close(20)
end do

if(nColumns == 10) then
	do i=1,nPoints
		write(60,*) energy(i), totalData(:,i)
	end do
else
	do i=1,nPoints
		write(60,*) energy(i), totalData(1,i), totalData(5,i), totalData(9,i), totalData(13,i), totalData(17,i), totalData(21,i), totalData(25,i), totalData(29,i), totalData(33,i)
	end do
endif

do i=1,nPoints
	write(61,*) energy(i), densityA(i)
	write(62,*) energy(i), densityB(i)
	write(63,*) energy(i), densityB(i)
end do

close(60)
close(61)
close(62)
close(63)
deallocate(dataLine)
deallocate(totalData)

end program sumfiles
