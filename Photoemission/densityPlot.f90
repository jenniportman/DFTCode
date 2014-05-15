program densityCalc
implicit none
real, parameter :: dx=1.0e0, dy=1.0e-4, dz=1.0e-6, nXSigma=3.0, nZSigma=4.0, nPzSigma=3.0
integer, parameter :: nSlices=800, nGrid=200

character(80) :: fileIn,fileIn2,fileOut,fileOut2,fileOut3,fileOut4, fileOut5

integer :: i,stat,j,k
integer :: nStep, nProc, nPart
integer :: ix, iz, ipz, nxzMax, npzMax, densityCount 
integer, dimension(nSlices) :: zHistogram
integer, dimension(nGrid+1, nGrid+1) :: densityArrayXZ, densityArrayPZ

real, dimension(nGrid+1, nGrid+1) :: vectorDx, vectorDz
real, dimension(nGrid+1, nGrid+1) :: densityArrayXZSmooth
real, dimension(15) :: inputVector
real, dimension(6) :: particleVector
real :: zMin, zMax, xMin, xMax, pzMin, pzMax, dzGrid, dxGrid, dpzGrid 
real :: sigmaPx, sigmaPz, meanPz, pxMaxV, pzMaxV
real :: densitySum, xzMaxSmooth

write(*,*) "Insert time step number"
read(*,*) nStep
write(*,*) "Insert total number of processors"
read(*,*) nProc

! Read in pulse parameters from file
write(fileIn2,'(a)') "screen.txt"
open(unit=39,file=trim(fileIn2),status="old")
do i=1,nStep-1
    read(39,*)
end do
read(39,*) inputVector
close(39)
nPart=int(inputVector(2))
meanPz=inputVector(12)
sigmaPz=inputVector(15)
sigmaPx=inputVector(13)

! Calculate grid ranges and spacing
!zMin=0.0
zMin=inputVector(6)-nZSigma*inputVector(9)
zMax=inputVector(6)+nZSigma*inputVector(9)
dzGrid=(zMax-zMin)/real(nGrid)

xMin=-nXSigma*inputVector(7)
xMax=nXSigma*inputVector(7)
dxGrid=(xMax-xMin)/real(nGrid)

pzMin=inputVector(12)-nPzSigma*inputVector(15)
pzMax=inputVector(12)+nPzSigma*inputVector(15)
dpzGrid=(pzMax-pzMin)/real(nGrid)

write(*,*) "grid ranges"
write(*,*) "z from", zMin, "to", zMax
write(*,*) "x from", xMin, "to", xMax
write(*,*) "pz from", pzMin, "to", pzMax

! All the output files
write(fileOut,'(i,a)') nStep,"-shadowImage.txt"
open(unit=41,file=trim(fileOut),status="unknown", recl=500)
write(fileOut2,'(i,a)') nStep,"-xz_2.txt"
open(unit=42,file=trim(fileOut2),status="unknown", recl=500)
write(fileOut3,'(i,a)') nStep,"-zpz_2.txt"
open(unit=43,file=trim(fileOut3),status="unknown", recl=500)
write(fileOut4,'(i,a)') nStep,"-vector.txt"
open(unit=44,file=trim(fileOut4),status="unknown", recl=500)
write(fileOut5,'(i,a)') nStep,"-xz_smooth.txt"
open(unit=45,file=trim(fileOut5),status="unknown", recl=500)

! read in all the separate files and calculate stuff
zHistogram=0
densityArrayXZ=0
densityArrayPZ=0
vectorDx=0
vectorDz=0
k=0
do i=1,nProc

    ! Generate file name and open it
	if (i < 10) then
		write(fileIn,'(i,a,i1,a)') nStep,"-x-",i,".txt"
	else
		write(fileIn,'(i,a,i2,a)') nStep,"-x-",i,".txt"
	end if
	write(*,*) "reading file",trim(fileIn)
	open(unit=40,file=trim(fileIn),status="old")

    ! Start reading file and bin data
	stat=0
	do while(stat==0)
		k=k+1
		read(40,*,iostat=stat) particleVector
        if(stat == 0) then
            ! shadow image stuff
            if(abs(particleVector(1)) < dx .and. abs(particleVector(2)) < dy) then
                j=int(particleVector(3)/dz)+1
                if(j > nSlices) then
                    !write(*,*) "skipping, z=",particleVector(3)
                else
                    zHistogram(j)=zHistogram(j)+1
                end if
			end if 
            
            ! density plot stuff
            if(particleVector(3) < zMax .and. particleVector(3) > zMin) then
                iz=int((particleVector(3)-zMin)/dzGrid)+1
                if(particleVector(1) < xMax .and. particleVector(1) > xMin) then
                    ix=int((particleVector(1)-xMin)/dxGrid)+1
                    densityArrayXZ(iz,ix)=densityArrayXZ(iz,ix)+1
                    vectorDx(iz,ix)=vectorDx(iz,ix)+particleVector(4)
                    vectorDz(iz,ix)=vectorDz(iz,ix)+particleVector(6)
                end if
                if(particleVector(6) < pzMax .and. particleVector(6) > pzMin) then
                    ipz=int((particleVector(6)-pzMin)/dpzGrid)+1
                    densityArrayPZ(iz,ipz)=densityArrayPZ(iz,ipz)+1
                end if
            end if

		end if
	end do
	k=k-1
	close(40)
end do

! Check that correct number of particles was read
!if (k /= nPart) then
    !write(*,*) "Error, number particles read doesnt match screen.txt data!"
    !stop
!end if

! Calculate shadow image slices and write to file
do i=1,nSlices
    write(41,*) i*dz,zHistogram(i)
end do

! Smooth arrays for nicer plots
do i=1,nGrid+1
    do j=1,nGrid+1
        densitySum=2*densityArrayXZ(i,j)
        densityCount=2
        if(i > 1) then
            densitySum=densitySum+densityArrayXZ(i-1,j)
            densityCount=densityCount+1
            if(j > 1) then
                densitySum=densitySum+densityArrayXZ(i-1,j-1)
                densityCount=densityCount+1
            end if
            if(j < nGrid+1) then
                densitySum=densitySum+densityArrayXZ(i-1,j+1)
                densityCount=densityCount+1
            end if
        end if
        if(i < nGrid+1) then
            densitySum=densitySum+densityArrayXZ(i+1,j)
            densityCount=densityCount+1
            if(j > 1) then
                densitySum=densitySum+densityArrayXZ(i+1,j-1)
                densityCount=densityCount+1
            end if
            if(j < nGrid+1) then
                densitySum=densitySum+densityArrayXZ(i+1,j+1)
                densityCount=densityCount+1
            end if
        end if
        if(j > 1) then
            densitySum=densitySum+densityArrayXZ(i,j-1)
            densityCount=densityCount+1
        end if
        if(j < nGrid+1) then
            densitySum=densitySum+densityArrayXZ(i,j+1)
            densityCount=densityCount+1
        end if
        densityArrayXZSmooth(i,j)=densitySum/real(densityCount)
    end do
end do

! Calculate density arrays and write to file
nxzMax=maxval(densityArrayXZ)
npzMax=maxval(densityArrayPZ)
pxMaxV=maxval(vectorDx)
pzMaxV=maxval(vectorDz)
xzMaxSmooth=maxval(densityArrayXZSmooth)

do i=1,nGrid+1
    do j=1,nGrid+1
        write(42,*) zMin+(i-1)*dzGrid, xMin+(j-1)*dxGrid, real(densityArrayXZ(i,j))/real(nxzMax)
        write(43,*) zMin+(i-1)*dzGrid, pzMin+(j-1)*dpzGrid, real(densityArrayPZ(i,j))/real(npzMax)
        if(vectorDz(i,j) /= 0 .or. vectorDx(i,j) /= 0) then 
            write(44,*) zMin+(i-1)*dzGrid, xMin+(j-1)*dxGrid, (vectorDz(i,j)-meanPz)/pzMaxV*dzGrid, vectorDx(i,j)/pxMaxV*dxGrid
        else
            write(44,*) zMin+(i-1)*dzGrid, xMin+(j-1)*dxGrid, (vectorDz(i,j))/pzMaxV*dzGrid, vectorDx(i,j)/pxMaxV*dxGrid
        end if
        write(45,*) zMin+(i-1)*dzGrid, xMin+(j-1)*dxGrid, densityArrayXZSmooth(i,j)/xzMaxSmooth
    end do
    write(42,*)
    write(43,*)
    write(44,*)
    write(45,*)
end do

close(41)
close(42)
close(43)
close(44)
close(45)

write(*,*) "rows read",k
write(*,*) "total number z",sum(zHistogram(:))

end program densityCalc
