program emittance
implicit none
character(10) :: stringField
integer :: i,j,k,stat, nProc, nMacroParticle, nSteps
character(80) :: fileIn,fileIn2,fileOut,fileOut2, fileOut3
real, dimension(3) :: meanX,meanX2,meanP,meanP2,meanXP,rmsEmittance,xVect,pVect
real, dimension(3) :: sigmaX,sigmaP,sigmaXP
real, allocatable :: time(:)
real :: x1, x2

write(*,*) "Insert total number of processors"
read(*,*) nProc
write(*,*) "Insert number of electrons in each macroparticle"
read(*,*) nMacroParticle
write(*,*) "Insert total number of time steps"
read(*,*) nSteps

allocate(time(nSteps))

write(fileOut,*) "emittance_data_norm.txt"
open(unit=41,file=trim(fileOut),status="unknown",recl=400)
write(fileOut2,*) "emittance_data.txt"
open(unit=42,file=trim(fileOut2),status="unknown",recl=400)
write(fileOut3,*) "check.txt"
open(unit=44,file=trim(fileOut3),status="unknown",recl=600)

write(fileIn2,'(a)') "screen.txt"
open(unit=43,file=trim(fileIn2),status="old")
do i=1,nSteps
    read(43,*) x1, x2, time(i)
end do
close(43)

do i=10,nSteps,10
    meanX=0.0
    meanX2=0.0
    meanP=0.0
    meanP2=0.0
    meanXP=0.0
	j=0

    do k=1,nProc
        stat=0
        ! Generate file name and open it
        if (k < 10) then
            write(fileIn,'(i,a,i1,a)') i,"-x-",k,".txt"
        else
            write(fileIn,'(i,a,i2,a)') i,"-x-",k,".txt"
        end if
        open(unit=40,file=trim(fileIn),status="old")

        do while(stat==0)
            j=j+1
            read(40,*,iostat=stat) xVect,pVect
            if(stat == 0) then
                meanX(:)=meanX(:)+xVect(:)
                meanX2(:)=meanX2(:)+xVect(:)**2.0
                meanP(:)=meanP(:)+pVect(:)
                meanP2(:)=meanP2(:)+pVect(:)**2.0
                meanXP(:)=meanXP(:)+xVect(:)*pVect(:)
            end if
        end do
        close(40)
        j=j-1
    end do

    write(*,*) "Nel", j
    meanP(:)=meanP(:)/nMacroParticle
    meanP2(:)=meanP2(:)/nMacroParticle**2.d0
    meanXP(:)=meanXP(:)/nMacroParticle
    sigmaX(:)=(real(j)*meanX2(:)-meanX(:)*meanX(:))/(real(j)*real(j))
    sigmaP(:)=(real(j)*meanP2(:)-meanP(:)*meanP(:))/(real(j)*real(j))
    sigmaXP(:)=(real(j)*meanXP(:)-meanX(:)*meanP(:))/(real(j)*real(j))
	rmsEmittance(:)=sqrt(sigmaX(:)*sigmaP(:)-(sigmaXP)**2.d0)
    if(rmsEmittance(3) .ne. rmsEmittance(3)) then
        write(*,*) "Problem!", rmsEmittance(3)
    end if
	write(41,*) time(i),j*nMacroParticle,rmsEmittance
	write(42,*) time(i),j*nMacroParticle,rmsEmittance(:)/meanP(3)*real(j)
	write(44,*) time(i),j,sigmaX(:),sigmaP(:),sigmaXP(:)
	!write(44,*) time(i),j,meanX(:)/real(j),sqrt(sigmaX(:)),meanP(:)/real(j),sqrt(sigmaP(:))
end do

close(41)
close(42)
close(44)

end program emittance
