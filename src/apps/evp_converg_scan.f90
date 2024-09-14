program evp_converg_scan

implicit none
INTEGER:: nrchop,IS,mread,akread

call system('make swipe')
call system('rm ./converg/*')

DO mread = 1,10
    DO akread = 1,10
        OPEN(11,FILE='almostdegen_read.input',STATUS='UNKNOWN')
        WRITE(11,'(i4)') mread
        WRITE(11,'(i4)') akread
        WRITE(11,'(i4)') 1
        WRITE(11,'(i4)') 1
        WRITE(11,'(f10.3)') 1
        WRITE(11,*) 'End'
        CLOSE(11)

        DO nrchop = 100,400,20
            OPEN(12,FILE='./converg/nrchop.input',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IS)
            WRITE(12,*) nrchop
            CLOSE(12)
            call system('mpirun -np 8 ./bin/evp_converg_exec ')
        ENDDO
    ENDDO
ENDDO

end program evp_converg_scan