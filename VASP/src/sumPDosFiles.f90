! sums partial dos files into one big file for plotting
program sumfiles
implicit none
integer, parameter :: nFiles=13, nPoints=3500, nColumns=10
character(80) :: inputFile, outputFile

integer :: i, j
character(30), dimension(nColumns) :: string
real, dimension(nColumns-1,nPoints) :: totalData
real, dimension(nColumns-1) :: dataLine
real, dimension(nPoints) :: energy

write(outputFile,*) 'dosp_sum.dat'
open(60,file=outputFile,recl=360)

totalData = 0.0
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
		totalData(:,j) = totalData(:,j) + dataLine(:)
	end do
	close(20)
end do

do i=1,nPoints
	write(60,*) energy(i), totalData(:,i)
end do

close(60)

end program sumfiles
