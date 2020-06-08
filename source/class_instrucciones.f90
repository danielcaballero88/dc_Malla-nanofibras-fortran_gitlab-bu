! ==========
module class_instrucciones
    ! ==========
    use class_instruccion
    use Aux, only : FindStringInFile, get_file_unit
    implicit none
    !private
    ! ==========

    ! ----------
    type, public :: instrucciones
        ! -----
        private
        character(len=120) :: nomarch
        integer :: fid
        logical :: arch_abierto
        integer :: num
        integer, allocatable :: ilabels(:)
        type(instruccion), allocatable :: lista(:)
        ! -----
    contains
        ! -----
        procedure :: init => init
        procedure :: leer => read_instrucciones
        procedure :: imprimir => imprimir
        ! -----
    end type instrucciones
    ! ----------

    ! ==========
    contains
    ! ==========

    ! ----------
    subroutine init(this, nomarch)
        ! -----
        implicit none
        class(instrucciones), intent(inout) :: this
        character(len=120), intent(in) :: nomarch
        ! -----
        ! ---
        this%nomarch = nomarch
        this%fid = get_file_unit()
        open(file=trim(this%nomarch), unit=this%fid, status='old')
        this%arch_abierto = .true.
        ! ---
        ! -----
    end subroutine init
    ! ----------

    ! ----------
    subroutine read_labels(this)
        ! -----
        implicit none
        class(instrucciones), intent(inout) :: this
        integer :: istat
        ! -----
        ! ---
        ! Busco la etiqueta necesaria
        istat = FindStringInFile("* Numero de acciones", this%fid, .true.)
        ! Leo el numero de instrucciones
        read(this%fid, *) this%num
        ! Adjudico las etiquetas numericas
        allocate( this%ilabels(this%num) )
        ! Leo las etiquetas numericas
        read(this%fid, *) this%ilabels
        ! Adjudico memoria para los parametros de las instrucciones
        allocate( this%lista(this%num) )
        ! ---
        ! -----
    end subroutine read_labels
    ! ----------

    ! ----------
    subroutine read_instr_list(this)
        ! -----
        implicit none
        class(instrucciones), intent(inout) :: this
        integer :: i
        integer :: istat
        ! -----
        ! ---
        do i=1,this%num
            istat = read_from_configfile( this%lista(i), this%fid, get_str_label(this%ilabels(i)) )
        end do
        ! ---
        ! -----
    end subroutine read_instr_list
    ! ----------

    ! ----------
    subroutine read_instrucciones(this)
        ! -----
        implicit none
        class(instrucciones), intent(inout) :: this
        ! -----
        call read_labels(this)
        call read_instr_list(this)
        ! ---
        close(unit=this%fid)
        ! -----
    end subroutine read_instrucciones
    ! ----------

    ! ----------
    subroutine imprimir(this)
        ! -----
        implicit none
        class(instrucciones), intent(in) :: this
        integer :: i
        ! -----
        write(6,*) "---"
        write(6,*) "Instrucciones"
        write(6,*) "n: ", this%num
        write(6,*) "labels: ", this%ilabels
        write(6,*) "---"
        do i=1,this%num
            write(6,*) get_str_label(this%ilabels(i))
            write(6,*) this%lista(i)%tipo
        end do
        ! -----
    end subroutine
    ! ----------

    ! ----------
    function get_str_label(ilabel) result(slabel)
        ! -----
        implicit none
        integer, intent(in) :: ilabel
        character(len=3) :: slabel
        ! -----
        ! ---
        WRITE(slabel,'(A2,I1)') "* ", ilabel
        ! ---
        ! -----
    end function
    ! ----------

    ! ==========
end module class_instrucciones
! ==========
