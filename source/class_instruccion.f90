! ==========================================================================
module class_instruccion
! ==========================================================================
    ! Modulo dedicado a opciones de seteo del problema: case parameters
    use Aux, only : FindStringInFile, get_file_unit
    ! ----
    implicit none
    ! ----
    private
    ! ----
    public :: instruccion
    public :: read_from_configfile
    ! ----

    ! ----------------------------------------------------------------------
    type instruccion
        ! ----
        ! variables generales
        character(len=120) :: tipo  ! nombre de la instruccion a realizar
        integer :: opcion_archivo ! opcion si se da un archivo (1) o una lista de archivos (2)
        character(len=120) :: nombre_archivo ! input file
        integer :: num_mallas ! cantidad de mallas a procesar
        character(len=120), allocatable :: archs_mallas(:) ! nombre(s) de archivo(s) la(s) malla(s) a procesar
        ! ----
        ! variables metodo intersectar
        integer :: num_pasos_intsec
        logical :: period_intsec
        ! variables metodo simplificar
        ! (no hay)
        ! variables metodo equilibrar
        integer :: num_pasos_vibra ! numero de iteraciones para resolver el equilibrio
        integer, allocatable :: lista_iters_vibra(:)
        real(8), allocatable :: lista_drs_vibra(:)
        real(8) :: fza_ref ! referencia de fuerza para calcular desplazamientos
        real(8) :: fza_tol ! tolerancia de fuerza para calcular si un nodo esta en equilibrio
        integer :: opcionF ! opcion para el F macro (=1 es un solo F, =2 es una lista de valores de F)
        character(len=120) :: arch_Fs ! solo se usa en el caso de haber dado opcionF=2 (lista de deformaciones)
        integer :: num_Fs
        real(8), allocatable :: lista_Fs(:,:,:) ! shape=(2,2,num_Fs)
        ! variables metodo traccion
        ! (uso varias del metodo simplificar)
        real(8) :: dtime, dot_F11, dot_F22, F11_fin
        integer :: opcion_save
        real(8) :: dF11_save
        integer :: num_F11_saves
        real(8), allocatable :: lista_F11_tosave(:) ! shape=(num_F11_saves)
        logical, allocatable :: lista_F11_ifsaved(:) ! shape=(num_F11_saves)
        logical, allocatable :: completed_F11_saves(:) ! shape=(num_mallas)
        ! variables metodo uniaxial
        ! (usa la mayoria de los metodos previos)
        real(8) :: ten22
        ! ----
        contains
        ! ----
    end type
    ! ----------------------------------------------------------------------

! ==========================================================================
contains
! ==========================================================================

    ! ----------------------------------------------------------------------
    ! HUB
    function read_from_configfile(self, cfid, label) result(istat)
        ! funcion para leer los parametros de una instruccion cualquiera a partir de configfile
        implicit none
        !
        type(instruccion), intent(inout) :: self
        integer, intent(in) :: cfid
        character(len=3), intent(in) :: label
        !
        integer :: istat
        integer :: istat2
        !
        istat2 = FindStringInFile(label, cfid, .True.)
        read(cfid, *) self%tipo
        select case (trim(self%tipo))
            case ("Intersectar")
                istat = read_from_configfile_intersectar(self, cfid)
            case ("Simplificar")
                istat = read_from_configfile_simplificar(self, cfid)
            case ("Equilibrar")
                istat = read_from_configfile_equilibrar(self, cfid)
            case ("Traccion")
                istat = read_from_configfile_traccion(self, cfid)
            case ("Uniaxial")
                istat = read_from_configfile_uniaxial(self, cfid)
            case default
                write(*,*) "Error, tipo de instruccion desconocido"
                stop
        end select
        !
        istat = 0 ! todo ok
        !
    end function read_from_configfile
    ! ----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! LOCAL
    subroutine get_lista_mallas(self)
        ! ----
        implicit none
        ! ----
        type(instruccion), intent(inout) :: self
        ! ----
        integer :: fid
        integer :: i
        ! ----

        if (self%opcion_archivo==1) then ! si opcion=1, entonces tengo una sola malla
            self%num_mallas = 1
            allocate( self%archs_mallas(self%num_mallas) )
            self%archs_mallas(1) = self%nombre_archivo
        else if (self%opcion_archivo==2) then ! si opcion_archivo=2, entonces se lee una lista de mallas
            fid = get_file_unit()
            open(unit=fid, file=trim(self%nombre_archivo), status="old")
            ! -
                read(fid, *) self%num_mallas
                allocate( self%archs_mallas(self%num_mallas) )
                do i=1,self%num_mallas
                    read(fid, *) self%archs_mallas(i)
                end do
            ! -
            close(unit=fid)
        else
            write(*,*) "Error, opcion debe ser 1 o 2 y es: ", self%opcion_archivo
        end if

        ! ----
    end subroutine get_lista_mallas
    ! ----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    function read_from_configfile_intersectar(self, cfid) result(istat)
        ! leer los parametros propios del metodo intersectar
        ! ya se leyo el tipo de instruccion, restan los demas parametros
        implicit none
        !
        type(instruccion), intent(inout) :: self
        integer, intent(in) :: cfid
        integer :: mfid
        integer :: istat
        integer :: i
        !
        read(cfid, *) self%opcion_archivo, self%nombre_archivo
        call get_lista_mallas(self)
        read(cfid, *) self%num_pasos_intsec
        read(cfid, *) self%period_intsec
        !
        istat = 0 ! todo ok
        !
    end function read_from_configfile_intersectar
    ! ----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    function read_from_configfile_simplificar(self, cfid) result(istat)
        ! leer los parametros propios del metodo simplificar
        ! ya se leyo el tipo de instruccion, restan los demas parametros
        implicit none
        !
        type(instruccion), intent(inout) :: self
        integer, intent(in) :: cfid
        integer :: mfid
        integer :: istat
        integer :: i
        !
        read(cfid, *) self%opcion_archivo, self%nombre_archivo
        call get_lista_mallas(self)
        !
        istat = 0 ! todo ok
        !
    end function read_from_configfile_simplificar
    ! ----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! LOCAL
    subroutine get_lista_Fs(self)
        ! ----
        implicit none
        ! ----
        type(instruccion), intent(inout) :: self
        ! ----
        integer :: fid
        integer :: i
        ! ----

        if (self%opcionF==1) then ! si opcion=1, entonces tengo un solo F dado en lugar de self%arch_Fs
            self%num_Fs = 1
            allocate( self%lista_Fs(2,2,self%num_Fs) )
            read(self%arch_Fs, *) self%lista_Fs(:,:,1) ! leo los valores de F de la strin self%arch_Fs
        else if (self%opcionF==2) then ! si opcion_archivo=2, entonces se lee una lista de Fs del archivo dado en line
            fid = get_file_unit()
            open(unit=fid, file=trim(self%arch_Fs), status="old")
            ! -
                read(fid, *) self%num_Fs
                allocate( self%lista_Fs(2,2,self%num_Fs) )
                do i=1,self%num_Fs
                    read(fid, *) self%lista_Fs(:,:,i)
                end do
            ! -
            close(unit=fid)
        else
            write(*,*) "Error, opcion debe ser 1 o 2 y es: ", self%opcionF
        end if

        ! ----
    end subroutine get_lista_Fs
    ! ----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    function read_from_configfile_equilibrar(self, cfid) result(istat)
        ! leer los parametros propios del metodo equilibrar
        ! ya se leyo el tipo de instruccion, restan los demas parametros
        implicit none
        !
        type(instruccion), intent(inout) :: self
        integer, intent(in) :: cfid
        integer :: Ffid
        integer :: istat
        integer :: i
        !
        ! leo las mallas a procesar
        read(cfid, *) self%opcion_archivo, self%nombre_archivo
        call get_lista_mallas(self)
        ! leo parametros de solver
        read(cfid, *) self%num_pasos_vibra
        if (self%num_pasos_vibra > 0) then ! si es =0 entonces no hay pasos de equilibrio, se calcula el afin y nada mas
            allocate( self%lista_iters_vibra(self%num_pasos_vibra) )
            allocate( self%lista_drs_vibra(self%num_pasos_vibra) )
        end if
        read(cfid, *) self%lista_iters_vibra
        read(cfid, *) self%lista_drs_vibra
        read(cfid, *) self%fza_ref, self%fza_tol
        ! leo deformaciones macro a equilibrar
        read(cfid, *) self%opcionF
        read(cfid, *) self%arch_Fs
        call get_lista_Fs(self)
        !
        istat = 0 ! todo ok
        !
    end function read_from_configfile_equilibrar
    ! ----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! LOCAL
    subroutine assemble_saves(self)
        ! ---
        implicit none
        ! ---
        type(instruccion), intent(inout) :: self
        ! ---
        integer :: i
        ! ---

        if (self%opcion_save==1) then ! preparo una lista de saves
            self%num_F11_saves = int( (self%F11_fin-1.d0) / self%dF11_save ) + 1
            allocate( self%lista_F11_tosave(self%num_F11_saves) )
            allocate( self%lista_F11_ifsaved(self%num_F11_saves) )
            do i=1,self%num_F11_saves
                self%lista_F11_tosave(i) = 1.d0 + dfloat(i-1)*self%dF11_save
                self%lista_F11_ifsaved(i) = .false.
            end do
            allocate ( self%completed_F11_saves(self%num_mallas) )
            self%completed_F11_saves = .false.
        else
            self%num_F11_saves = 0
            allocate ( self%completed_F11_saves(self%num_mallas) )
            self%completed_F11_saves = .true.
        end if

        ! ---
    end subroutine
    ! ----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    function read_from_configfile_traccion(self, cfid) result(istat)
        ! leer los parametros propios del metodo traccion
        ! ya se leyo el tipo de instruccion, restan los demas parametros
        implicit none
        !
        type(instruccion), intent(inout) :: self
        integer, intent(in) :: cfid
        integer :: Ffid
        integer :: istat
        integer :: i
        !
        ! leo las mallas a procesar
        read(cfid, *) self%opcion_archivo, self%nombre_archivo
        call get_lista_mallas(self)
        ! leo parametros de solver
        read(cfid, *) self%num_pasos_vibra
        if (self%num_pasos_vibra > 0) then ! si es =0 entonces no hay pasos de equilibrio, se calcula el afin y nada mas
            allocate( self%lista_iters_vibra(self%num_pasos_vibra) )
            allocate( self%lista_drs_vibra(self%num_pasos_vibra) )
        end if
        read(cfid, *) self%lista_iters_vibra
        read(cfid, *) self%lista_drs_vibra
        read(cfid, *) self%fza_ref, self%fza_tol
        ! leo parametros temporales de la traccion
        read(cfid, *) self%dtime, self%dot_F11, self%dot_F22, self%F11_fin
        read(cfid, *) self%opcion_save, self%dF11_save
        call assemble_saves(self)
        !
        istat = 0 ! todo ok
        !
    end function read_from_configfile_traccion
    ! ----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    function read_from_configfile_uniaxial(self, cfid) result(istat)
        ! leer los parametros propios del metodo uniaxial
        ! ya se leyo el tipo de instruccion, restan los demas parametros
        implicit none
        !
        type(instruccion), intent(inout) :: self
        integer, intent(in) :: cfid
        integer :: Ffid
        integer :: istat
        integer :: i
        !
        ! leo las mallas a procesar
        read(cfid, *) self%opcion_archivo, self%nombre_archivo
        call get_lista_mallas(self)
        ! leo parametros de solver
        read(cfid, *) self%num_pasos_vibra
        if (self%num_pasos_vibra > 0) then ! si es =0 entonces no hay pasos de equilibrio, se calcula el afin y nada mas
            allocate( self%lista_iters_vibra(self%num_pasos_vibra) )
            allocate( self%lista_drs_vibra(self%num_pasos_vibra) )
        end if
        read(cfid, *) self%lista_iters_vibra
        read(cfid, *) self%lista_drs_vibra
        read(cfid, *) self%fza_ref, self%fza_tol
        ! leo parametros temporales de la traccion
        read(cfid, *) self%dtime, self%dot_F11, self%ten22, self%F11_fin
        read(cfid, *) self%opcion_save, self%dF11_save
        call assemble_saves(self)
        !
        istat = 0 ! todo ok
        !
    end function read_from_configfile_uniaxial
    ! ----------------------------------------------------------------------

! ==========================================================================
end module class_instruccion
! ==========================================================================
