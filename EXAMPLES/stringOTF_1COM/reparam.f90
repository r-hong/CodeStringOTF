       program reparam

!-------------------------------
!  reparametrization of a string 
!  nimage input number of images
!  nimage output number of images
!  nd coll vars dimension

!  c(nimage,nd) input str
!  cr(nimage,nd) output str
!-------------------------------

        implicit none
        integer, parameter:: nimage=3
        integer, parameter:: nimage2=3 !number of points in the string 
        integer, parameter:: nd=3  !! number of CVs (columns in file initialPcurve.dat) if 1COM (x,y,z) as CV ---> nd=3 
        integer:: i,j,k            !! if nd=6 (e.g. 2 COMs) then also CHANGE FORMAT at 101 
        real*8:: c(nimage,nd)
        real*8:: cr(nimage2,nd)
        real*8:: cd(nd)
        real*8:: scd,d
        real*8:: LL(nimage),s(nimage2)
        real*8:: L
        character*15::strin,strout

       strin="images_norep.ls"
       strout="images_rep.ls"

       open(10,file=strin,status="old")
       do i = 1,nimage
        read(10,*) (c(i,k),k=1,nd)
       enddo
       close(10)

       L=0.d0

       do i=1,nimage
        LL(i)=0.d0
       enddo

       do i=2,nimage

        scd = 0.d0
        do k = 1,nd
         cd(k)=c(i,k)-c(i-1,k)
         scd = scd + cd(k)*cd(k)
        enddo

         write(*,*) i,dsqrt(scd),"before"

          LL(i)=LL(i-1)+dsqrt(scd)

       ! write(*,*) i,LL(i),"LL 32"

       enddo

       L=LL(nimage)

       do k=1,nimage2

         s(k)=(k-1)*L/((nimage2)-1.d0)

!        write(*,*) k,s(k),L/((nimage2)-1.d0),"s 64"

       enddo

       do k = 1,nd
       cr(1,k)=c(1,k)
       cr(nimage2,k)=c(nimage,k)
       enddo


       do i=2,nimage2-1
 
        do j=2,nimage


         if(LL(j).ge.s(i)) then

           scd = 0.d0
           do k = 1,nd
            cd(k)=c(j,k)-c(j-1,k)
            scd = scd + cd(k)**2
           enddo
           d=dsqrt(scd)
!          write(*,*) "LOOP", i,j,LL(j),s(i),d

           do k = 1,nd
           cr(i,k)=(s(i)-LL(j-1))*cd(k)/d+c(j-1,k)
           enddo

           goto 10

         endif

        enddo

10      continue
       enddo

! check done repa 
       do i=2,nimage2
        scd = 0.d0
        do k = 1,nd
         cd(k)=cr(i,k)-cr(i-1,k)
         scd = scd + cd(k)**2
        enddo
         write(*,*) i,dsqrt(scd),"after"
       enddo

      open(10,file=strout,status="unknown")
       do i = 1,nimage2
        write(10,101) (cr(i,k),k=1,nd) !! to write single files for all process
       enddo
       close(10)

101   format(3f10.4)  ! if nd=6 then change this to format(6f10.4)


       end


