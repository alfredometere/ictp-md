!-------------------------------------------------------------------------------
! GENSYS v0.1 - XYZ File generator utility
! (C) Alfredo Metere, 2012-Aug-09. All Rights Reserved
! 
! This small utility is released under the GPLv3 License
! Info at: http://www.gnu.org/licenses/quick-guide-gplv3.html
!
! Description:
! This small utility can generate XYZ files describing homoatomic systems
! inside a box of arbitrary dimensions.
! The user has the freedom of choosing the atomic type (it is important for
! its atomic radius and therefore the system density) and the filling mode.
!
! Enjoy!
! 
! Credits: Axel Kohlmeyer
!
!-------------------------------------------------------------------------------
program gxyz

    implicit none

    ! Constants
    integer      (4), parameter             :: bufsize    = 256
    integer      (4), parameter             :: n_elements = 120 ! atoms.dat size
    character(len=4), parameter             :: f_ext      = ".xyz"

    ! Input from the command line variables
    character(len=bufsize)                  :: argv, appname
    real               (8)                  :: m,n,o ! Box Size
    character      (len=3)                  :: atomsym ! Atomic Symbol
    integer            (4)                  :: fmode ! Filling mode
    character(len=bufsize)                  :: outfile ! Output file name

    ! Storage variables for elements table file
    logical                                 :: db_status
    real         (8), dimension(n_elements) :: a_mass, a_radius
    integer      (4), dimension(n_elements) :: a_number
    character(len=3), dimension(n_elements) :: a_name

    ! Filling modes implemented for now (2012-08-10):
    ! 1 - Linear order
   
    ! Variables for internal usage
    integer      (4)                        :: argc ! Holder for iargc()
    integer      (4)                        :: atom_id ! ID of the chosen atom
    integer      (4)                        :: i ! Generic counter

    real         (8)                        :: rpart
    integer      (8)                        :: natoms
    real         (8), allocatable           :: x(:), y(:), z(:)

    ! Global Variables
    ! Not always a good idea, but safe in this case
    common appname, a_radius, a_mass, a_number, a_name

    ! Greeter    
    call welcome()
    
    ! Assignment of default values

    call getarg(0,appname)
    argc = iargc()

    atom_id = 0
    fmode = 1
    outfile = "default" // f_ext

    ! Command line processing conditional block
    select case (argc)
        case(4)
            call getarg(1,argv); read(argv,*) m
            call getarg(2,argv); read(argv,*) n
            call getarg(3,argv); read(argv,*) o
            call getarg(4,argv); read(argv,*) atomsym

            print *, "Assuming default for output file name: ", trim(outfile)
            print *, "Assuming default for fill mode", fmode
        
        case(5)
            call getarg(1,argv); read(argv,*) m
            call getarg(2,argv); read(argv,*) n
            call getarg(3,argv); read(argv,*) o
            call getarg(4,argv); read(argv,*) atomsym
            call getarg(5,argv); read(argv,*) outfile

            print *, "Assuming default for fill mode", fmode

        case(6)
            call getarg(1,argv); read(argv,*) m
            call getarg(2,argv); read(argv,*) n
            call getarg(3,argv); read(argv,*) o
            call getarg(4,argv); read(argv,*) atomsym
            call getarg(5,argv); read(argv,*) outfile
            call getarg(6,argv); read(argv,*) fmode

        case default
            call except(0) ! Help Exception
    end select

    ! Checks if the filling mode is valid and proceeds to locate the
!             particles accordingly
    select case(fmode)
        case (1)
            print *, "Linear filling mode selected"                                                   
        case default
            call except(4)
    end select

    ! Preliminary check for atoms.dat and loading it into memory
    write (*,'(A30)',advance='no') "Searching for [atoms.dat] ..."
    inquire(file="atoms.dat", exist=db_status)

    if (db_status .eqv. .false.) call except(2)

    print *, "Found!"

    write (*,'(A33)',advance='no') "Loading [atom.dat] in memory ..."
    
    open(1,file='atoms.dat',status='old')

    ! This loop transfers the values of atoms.dat in memory .AND.
    ! at the same time annotates the atom id specified in the command line
    do i = 1, n_elements
        read(1,*) a_number(i), a_name(i), a_radius(i), a_mass(i)

        ! Grabs the id for the atom
        if (atomsym .eq. a_name(i)) atom_id = a_number(i)
    end do

    close(1)
    print *, "Done!"
    print *

    ! Atom not found triggers exception
    if (atom_id == 0) call except(3)

    call genpos(m,n,o,atom_id,fmode)

! Subroutines are here!!! RUUUN!!! ---------------------------------------------
    contains

! WELCOME SUBROUTINE -----------------------------------------------------------
        subroutine welcome()
            ! Welcome message
            print *, "---------------------------------------------------------"
            print *, "- Welcome to GENSYS v.0.1                               -"
            print *, "- Utility for MD initial homoatomic system generator    -"
            print *, "- Developed as an utility for testing ICTP-MD           -"
            print *, "- (C) Alfredo Metere. All rights reserved.              -"
            print *, "-                                                       -"
            print *, "- This code is released under the GPLv3 license         -"
            print *, "- For more info about the license:                      -"
            print *, "- http://www.gnu.org/licenses/quick-guide-gplv3.html    -"
            print *, "-                                                       -"
            print *, "- Big thanks to Axel Kohlmeyer for his support!         -"
            print *, "---------------------------------------------------------"
            print *, "Enjoy!"
            print *
            print *
            return
        end subroutine welcome

! GENPOS SUBROUTINE ------------------------------------------------------------
        subroutine genpos(m,n,o,atom_id,fmode)

            implicit none

            ! Input Parameters
            integer(4), intent(in)                  :: fmode,atom_id
            real(8), intent(in)                     :: m,n,o

            ! Output Parameters


            ! COMMON VARS
            character,        dimension(bufsize)    :: appname
            real(8),          dimension(n_elements) :: a_mass, a_radius
            integer(4),       dimension(n_elements) :: a_number
            character(len=3), dimension(n_elements) :: a_name

            ! Internal variables
            integer(8)                              :: i,j,k ! Counters
            real(8)                                 :: dpart ! atomic diameter
            real(8)                                 :: radius ! atomic radius
            integer(8)                              :: sx,sy,sz ! Subcells sizes
            real(8)                                 :: nx,ny,nz ! 
            real(8)                                 :: fsize ! Expectd file size

            integer(8)                              :: natoms
            real(8), allocatable                    :: x(:),y(:),z(:)            

            ! IMPORTANT: The order of the variables determines the memory
!                        alignment. If you change this order, something weird
!                        will happen. (Tony Montana)
            common appname, a_radius, a_mass, a_number, a_name

            ! Take the excess integer part of the atomic diameter
            dpart = (a_radius(atom_id) * 2) + 0.5
            radius = dpart / 2.0

            write (*,'(A26,I40)') "Atom ID:",atom_id
            write (*,'(A26,A40)') "Atom name:",trim(a_name(atom_id))
            write (*,'(A26,F40.4)') "Atomic occupancy:",dpart
            ! Checks that the chosen box size is not too small for the
!             chosen element.
            if (m <= dpart .or. n <= dpart .or. o <= dpart) call except(1)

            ! Determine the maximum number of particles the box can contain
            
            sx = m / dpart
            sy = n / dpart
            sz = o / dpart

            natoms = sx*sy*sz
            write(*,'(A26)', advance='no') "Total mass of the system:"
            write(*,'(F40.4,A7)') a_mass(atom_id) * natoms," a.m.u."
            write(*,'(A26,I40)') "Number of atoms:", natoms

            ! Calculates expected file size
            write(*,'(A26)',advance='no') "Expected file size:"
            fsize = natoms * 24 + natoms * 3

            ! Humanize the file size readability
            select case (int(fsize))
                case (1025:1024**2-1)
                    fsize = fsize / 1024
                    write(*,'(F40.2,A3)') fsize,"KB"
                case (1024**2:1024**3-1)
                    fsize = fsize / 1024**2
                    write(*,'(F40.2,A3)') fsize,"MB"
                case (1024**3:)
                    fsize = fsize / 1024**3
                    write(*,'(F40.2,A3)') fsize,"GB"
                case default
                    write(*,'(F40.2,A3)') fsize,"B "
            end select
            
            allocate(x(sx),y(sy),z(sz))
            open (666, file = 'full_rd.out', access = 'sequential',&
                  &status = 'replace', form = 'formatted')
            do i = 1, sx
                do j = 1, sy
                    do k = 1, sz
                       x(i) = (i-1) + radius * i
                       y(j) = (j-1) + radius * j
                       z(k) = (k-1) + radius * k
                       write(666,'(F6.4,A1,F6.4,A1,F6.4)')&
                             & x(i)," ", y(j)," ", z(k)                  
                    end do
                end do                
            end do

            deallocate(x,y,z)
            return
        end subroutine genpos

! GENFILE SUBROUTINE -----------------------------------------------------------
        subroutine genfile(natoms,atomname,x,y,z,outfile)
            implicit none
            
            integer(4) :: natoms
            real(8), dimension(natoms) :: x,y,z
            character(len=bufsize) :: outfile
            character(len=3) :: atomname

            print *,"Hi, I am genfile subroutine! Nice to meet you"
            print *,natoms,atomname,x,y,z,outfile

            !write(*,'(F5.4,F5.4,F5.4)'), x(i), y(j), z(k)
            return
        end subroutine genfile

! EXCEPT SUBROUTINE ------------------------------------------------------------
        subroutine except(state)

            implicit none

            character(len=bufsize) :: app
            
            integer(4), intent(in) :: state            
            
            common app
        
            select case (state)

                case (0)
                    print *, "Usage"
                    print *
                    print *, "Synopsis: "
                    print *
                    print *, trim(app)," m n o esym (ofile fmode)"
                    print *
                    print *, "Options: "
                    print *
                    print *, "m n o - Real values defining the length of the"
                    print *, "------- box sides."
                    print *
                    print *, "esym  - Symbol of the preferred element"
                    print *, "------- Putting 'Z10' will create an archetypal"
                    print *, "------- system of Van der Waals radius 1.00 Å and"
                    print *, "------- mass 100 a.m.u."
                    print *, "------- Putting 'Z14' will create an archetypal"
                    print *, "------- system of Van der Waals radius 1.40 Å and"
                    print *, "------- mass 100 a.m.u."
                    print *
                    print *, "ofile - Optional: Output file name"
                    print *, "------- You should not specify any extension to"
                    print *, "------- the file name."
                    print *, "------- Default value is 'default.xyz'"
                    print *
                    print *, "fmode - Optional: Filling mode"
                    print *, "------- Choose 1 to linearly displace the atoms"
                    print *, "------- Choose 2 to randomly displace the atoms"
                    print *, "------- Default value is '1'"

                    print *
                    print *, "Example: "
                    print *, trim(app), " 26.4 32.345 12.4 Ar structfile 0"
                    print *
                    print *, "It will generate a system where argon atoms are"
                    print *, "randomly displaced in a 26.4 x 32.345 x 12.4 box"
                    print *, "and it will save this structure in a file called"
                    print *, "'structfile.xyz'"
                    stop

                case (1)
                    print *, "Exception",state
                    print *, "Box size unsufficiently big for the atomic radius"
                    stop

                case (2)
                    print *, "Exception",state
                    print *, "Cannot find [atoms.dat] in the working folder."
                    stop

                case (3)
                    print *, "Exception",state
                    print *, "Cannot find the element in the database. Retry"
                    stop

                case (4)
                    print *, "Exception",state
                    print *, "Filling mode is not implemented yet. Retry"
                    stop

                case default
                    print *, "Exception",state
                    print *, "Unhandled exception! Something is wrong with your"
                    print *, "hardware. This program is too perfect to fail! :D"
                    print *, "Just joking."
                    print *, "Double-check your compiler optimization settings,"
                    print *, "because the optimization might be too aggressive."
                    stop
            end select

            stop
            return    
        end subroutine except

end program gxyz
