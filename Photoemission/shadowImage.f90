program shadowImage
implicit none
real, parameter :: dx=1.0e0, dy=1.0e0, dz=1.0e-6
integer, parameter :: nSlices=800
integer :: nStep, nProc

integer :: i,stat,j,k
character(80) :: fileIn,fileOut
integer, dimension(nSlices) :: zHistogram
real :: x, y, z

write(*,*) "Insert time step number"
read(*,*) nStep
write(*,*) "Insert total number of processors"
read(*,*) nProc

write(*,*) "This will generate shadow images with dx", dx, "dy", dy

write(fileOut,'(i,a)') nStep,"-shadowImage.txt"
open(unit=41,file=trim(fileOut),status="unknown")

! Save z positions from all the separate files for the electrons in the x slice
zHistogram=0
k=0
do i=1,nProc
	stat=0
	if (i < 10) then
		write(fileIn,'(i,a,i1,a)') nStep,"-x-",i,".txt"
	else
		write(fileIn,'(i,a,i2,a)') nStep,"-x-",i,".txt"
	end if
	write(*,*) "reading file",trim(fileIn)
	open(unit=40,file=trim(fileIn),status="old")
	do while(stat==0)
		k=k+1
		read(40,*,iostat=stat) x,y,z
		if(abs(x) < dx .and. abs(y) < dy .and. stat==0) then
			j=int(z/dz)+1
			if(j > nSlices) then
				write(*,*) "skipping, z=",z
			else
				zHistogram(j)=zHistogram(j)+1
			end if
		end if
	end do
	k=k-1
	close(40)
end do

! Calculate shadow image slices
do i=1,nSlices
	write(41,*) i*dz,zHistogram(i)
end do

close(41)

write(*,*) "rows read",k
write(*,*) "total number z",sum(zHistogram(:))

end program shadowImage
