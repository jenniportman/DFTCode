program radialProfile
implicit none
real, parameter :: sigmaXBin=10.0e-6 !m
integer :: nProc, nBin
integer, dimension(4) :: nTimeSteps=(/ 1000,2000,4000,8000 /), nPart
integer :: i,j,stat, iBin
character(80) :: fileIn,fileOut
real :: x, xMin, xMax
real, allocatable :: radialData(:,:)

write(*,*) "Insert total number of processors"
read(*,*) nProc

write(fileOut,*) "radial_profile.txt"
open(unit=41,file=trim(fileOut),status="unknown",recl=600)

xMin=1.0
xMax=0.0

do i=1,nProc
	stat=0
	if (i < 10) then
		write(fileIn,'(i,a,i1,a)') nTimeSteps(4),"-x-",i,".txt"
	else
		write(fileIn,'(i,a,i2,a)') nTimeSteps(4),"-x-",i,".txt"
	end if
    !write(*,*) fileIn
	open(unit=40,file=trim(fileIn),status="old")
	do while(stat==0)
		read(40,*,iostat=stat) x
        if(stat == 0) then
            !write(*,*) x
            if(x > xMax) xMax=x
            if(x < xMin) xMin=x
        end if
	end do
	close(40)
end do

nBin=int((xMax-xMin)/sigmaXBin+1)
write(*,*) "N bins used",nBin
allocate(radialData(nBin,4))
radialData(:,:)=0
nPart(:)=0
do i=1,4
    do j=1,nProc
        stat=0
        if (j < 10) then
            write(fileIn,'(i,a,i1,a)') nTimeSteps(i),"-x-",j,".txt"
        else
            write(fileIn,'(i,a,i2,a)') nTimeSteps(i),"-x-",j,".txt"
        end if
        !write(*,*) fileIn
        open(unit=40,file=trim(fileIn),status="old")
        do while(stat==0)
            read(40,*,iostat=stat) x
            if(stat == 0) then
                nPart(i)=nPart(i)+1
                iBin=int((x-xMin)*nBin/(xMax-xMin))+1
                if(iBin > nBin) iBin=nBin
                radialData(iBin,i)=radialData(iBin,i)+1
            end if
        end do
        close(40)
    end do
end do

do i=1,nBin
	write(41,*) xMin+(xMax-xMin)/real(nBin)*i,radialData(i,:)/nPart(:)
end do

close(41)
end program radialProfile
