
program hello
    USE class_malla_completa
    use class_mallita
    USE Aux
    implicit none
    ! ==========================================================================
    ! ==========================================================================
    ! Current working directory
    CHARACTER(LEN=255) :: cwd
    ! ==========================================================================
    TYPE(MallaCom) :: MC, MC2
    TYPE(MallaSim) :: ms, ms2
    CHARACTER(LEN=120) :: configfile
    integer :: fid_cf
    integer :: fid_lista_mallas
    integer :: nmallas
    integer :: num_instruc
    integer, allocatable :: lista_instruc(:)
    integer :: i_etiqueta
    character(len=3) :: str_etiqueta
    character(len=120) :: str_instruccion
    integer :: j_instr
    ! ----------
    ! variables para instruccion "Intersectar"
    integer :: opcion_archivo
    character(len=120) :: nombre_archivo, nombre_malla
    integer :: num_pasadas
    logical :: period_intersec
    integer :: j_malla
    ! ----------
    ! variables para instruccion "Simplificar"
    ! ----------
    ! variables para instruccion "Equilibrar"
    ! ----------
    CHARACTER(LEN=120) :: filename, mallaname
    integer :: i,j,k,n,m
    INTEGER :: iStat, iStat1, iStat2
    real(8), allocatable :: r1(:,:), Ag(:,:), bg(:), dr(:,:)
    character(len=120) :: formato
    real(8), allocatable :: fuerzas_f(:,:), fuerzas_n(:,:)
    real(8) :: Fmacro(2,2)
    integer :: opcion ! 1=intersectar, 2=simplificar, 3=equilibrio
    ! ==========================================================================

    ! --------------------------------------------------------------------------
    ! Imprimo carpeta de trabajo actual
    print *, "Multiscale Nanofibers Mesh RVE 2.0"
    CALL GETCWD(cwd)
    PRINT *, "Current Working Directory:", cwd
    ! --------------------------------------------------------------------------
    ! Leer ConfigurationFile.txt para obtener instrucciones
    configfile = "ConfigurationFile.txt"
    fid_cf = get_file_unit()
    OPEN(unit=fid_cf, file=trim(configfile), status="old")
    iStat = FindStringInFile("* Numero de acciones", fid_cf, .true.)
    read(fid_cf,*) num_instruc
    allocate( lista_instruc(num_instruc) )
    read(fid_cf,*) lista_instruc
    CLOSE(unit=fid_cf)
    ! --------------------------------------------------------------------------
    ! recorro las etiquetas de las instrucciones a realizar
    do j_instr=1,num_instruc
        ! --------------------------------------------------------------------------
        ! busco cada etiqueta en configfile y leo la string de instruccion (es el identificador: Intersectar, Simplificar o Equilibrar)
        i_etiqueta = lista_instruc(j_instr)
        WRITE(str_etiqueta,'(A2,I1)') "* ", i_etiqueta
        OPEN(unit=fid_cf, file=trim(configfile), status="old")
        iStat = FindStringInFile(str_etiqueta, fid_cf, .true.)
        read(fid_cf,*) str_instruccion
        ! --------------------------------------------------------------------------
        ! luego, segun la instruccion me fijo que hay que hacer
        select case (trim(str_instruccion))
        ! --------------------------------------------------------------------------
        ! INTERSECTAR
        case ("Intersectar")
            ! Leo los parametros de configuracion
            read(fid_cf,*) opcion_archivo, nombre_archivo
            read(fid_cf,*) num_pasadas
            read(fid_cf,*) period_intersec
            ! Comienzo a intersectar
            if (opcion_archivo==1) then
                ! Caso de una sola malla a intersectar
                call main_intersectar(nombre_archivo, num_pasadas, period_intersec)
            elseif (opcion_archivo==2) then
                ! Caso de una lista de mallas a intersectar, la lista se da en un archivo
                fid_lista_mallas = get_file_unit()
                open(unit=fid_lista_mallas, file=trim(nombre_archivo), status="old")
                read(fid_lista_mallas,*) nmallas
                ! Recorro las mallas de la lista, calculando las intersecciones para cada una
                do j_malla=1,nmallas
                    read(fid_lista_mallas,*) nombre_malla
                    write(*,*) "Intersectando malla:"
                    write(*,*) nombre_malla
                    call main_intersectar(nombre_malla, num_pasadas, period_intersec)
                end do
                ! Cierro el archivo de la lista de mallas
                close(unit=fid_lista_mallas)
            else
                ! Si no encontre opcion 1 o 2, entonces hay algun error!!!
                write(*,*) "Error, para Intersectar, opcion_archivo debe ser 1 o 2, y es: ", opcion_archivo
                write(*,*) "En etiqueta: ", str_etiqueta
                write(*,*) "En Instruccion: ", str_instruccion
                stop
            end if
        ! FIN INTERSECTAR
        ! --------------------------------------------------------------------------
        ! SIMPLIFICAR
        case ("Simplificar")
            write(*,*) str_instruccion
        ! FIN SIMPLIFICAR
        ! --------------------------------------------------------------------------
        ! EQUILIBRAR
        case ("Equilibrar")
            write(*,*) str_instruccion
        ! FIN EQUILIBRAR
        ! --------------------------------------------------------------------------
        case default
            write(*,*) "Instruccion desconocida: ", str_instruccion
            write(*,*) "En etiqueta: ", str_etiqueta
            stop
        ! --------------------------------------------------------------------------
        end select
        ! --------------------------------------------------------------------------

        close(unit=fid_cf)
    end do
    ! --------------------------------------------------------------------------



    !    mallaname = "default"
    !    opcion = 1
!    ! ==========
!    select case (opcion)
!    ! ==========
!    case (1)
!    ! ==========
!        call main_intersectar(mallaname)
!    ! ==========
!    case (2)
!    ! ==========
!        call main_simplificar(mallaname)
!    ! ==========
!    case (3)
!    ! ==========
!        call main_equilibrio(mallaname)
!    ! ==========
!    end select
!    ! ==========

end program
! ==========================================================================
! ==========================================================================
! ==========================================================================

! ==========================================================================
subroutine main_intersectar(filename_malla_in, npasadas, periodicidad)
    use class_malla_completa
    implicit none
    CHARACTER(LEN=120), intent(in) :: filename_malla_in
    integer, intent(in) :: npasadas
    logical, intent(in) :: periodicidad
    character(len=120) :: filename_malla_in2, filename_malla_out
    TYPE(MallaCom) :: MC, MC2
    integer :: i
    integer :: iStat1, iStat2

    if (trim(filename_malla_in) == "default") then
        filename_malla_in2 = "Malla.txt"
    else
        filename_malla_in2 = filename_malla_in
    end if

    write(*,*) "Leer malla, intersectar fibras y reescribir:"
    CALL leer_malla(MC, filename_malla_in2)
!    if (MC%sidelen > 99.d0) then
!        write(*,*) "malla con problema"
!    end if
    ! Hago la interseccion muchas veces porque cada vez tengo la limitacion de no cortar al mismo segmento dos veces
    i = 0
    write(*,*) "Intersectando fibras"
    iStat1 = 0
    iStat2 = 0
    write(*,*) mc%nsegs
    DO WHILE (.true.)
        i = i+1
        WRITE(*,'(I4)', ADVANCE='no') i
        CALL intersectar_fibras_3(MC, MC2, .FALSE., periodicidad, iStat1) ! dentro de la misma capa
        MC = MC2
        CALL intersectar_fibras_3(MC, MC2, .TRUE., periodicidad, iStat2) ! con capas adyacentes
        MC = MC2
        write(*,*) mc%nsegs
        IF ( (iStat1 == 1).AND.(iStat2 == 1) )  EXIT
        if (i==npasadas) exit
    END DO
    write(*,*)

    write(*,*) "Escribiendo malla intersectada"
    filename_malla_out = "_i"
    call modify_txt_filename(filename_malla_in2, filename_malla_out)
    CALL escribir_malla(mc, filename_malla_out)
    write(*,*) "Malla intersectada OK"

end subroutine main_intersectar
! ==========================================================================

! ==========================================================================
subroutine main_simplificar(filename_malla_in)
    use class_malla_completa
    use class_mallita
    implicit none
    CHARACTER(LEN=120), intent(in) :: filename_malla_in
    character(len=120) :: filename_malla_in2, filename_malla_out
    type(MallaCom) :: mc
    type(MallaSim) :: ms

    if (trim(filename_malla_in) == "default") then
        filename_malla_in2 = "Malla_i.txt"
    else
        filename_malla_in2 = filename_malla_in
    end if

    write(*,*) "Leer malla intersectada y generar malla simplificada:"
    call leer_malla(mc, filename_malla_in2)
    call Desde_MallaCom(mc, ms, 2, [10.d0, .1d0])

    write(*,*) "Escribiendo mallita"
    filename_malla_out = "_s"
    call modify_txt_filename(filename_malla_in2, filename_malla_out)
    CALL escribir_mallita(ms, filename_malla_out)
    write(*,*) "Malla simplificada OK"

end subroutine main_simplificar
! ==========================================================================

! ==========================================================================
subroutine main_equilibrio(filename_malla_in)
    use class_mallita
    implicit none
    CHARACTER(LEN=120), intent(in) :: filename_malla_in
    character(len=120) :: filename_malla_in2, filename_malla_out
    type(MallaSim) :: ms
    real(8), allocatable :: r1(:,:)
    real(8) :: Fmacro(2,2)
    integer :: n, iStat1

    if (trim(filename_malla_in) == "default") then
        filename_malla_in2 = "Malla_i_s.txt"
    else
        filename_malla_in2 = filename_malla_in
    end if

        write(*,*) "Calculando equilibrio"
        call leer_mallita(ms, filename_malla_in2)
        n = ms%nnods
        allocate( r1(2,n) )
        Fmacro(1,:) = [1.2d0, 0.d0]
        Fmacro(2,:) = [0.0d0, 1.d0]
        call deformar_afin(ms, Fmacro)
        call calcular_equilibrio(ms,r1,10000,1.d-1,iStat1)

        ms%rnods = r1
        filename_malla_out = "_e"
        call modify_txt_filename(filename_malla_in2, filename_malla_out)
        call escribir_mallita(ms, filename_malla_out)
        write(*,*) "Equilibrio calculado OK"
end subroutine main_equilibrio
! ==========================================================================

! ==========================================================================
subroutine modify_txt_filename(txt_filename_in, txt_filename_out)
    ! Agrego caracteres a un archivo txt antes de la extension ".txt"7
    ! El nombre original del archivo viene en txt_filename_in, que se mantiene sin modificar
    ! Los caracteres a agregar van en la variable txt_filename_out
    ! que es donde tambien se guarda el nombre modificado del archivo
    implicit none
    character(len=120), intent(in) :: txt_filename_in
    character(len=120), intent(inout) :: txt_filename_out
    integer :: nchars

    nchars = len(trim(txt_filename_in))
    txt_filename_out = trim(txt_filename_in(1:nchars-4)) // trim(txt_filename_out) // ".txt"

end subroutine modify_txt_filename
! ==========================================================================
