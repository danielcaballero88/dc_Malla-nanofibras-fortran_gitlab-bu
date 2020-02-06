
program hello
    USE class_malla_completa
    use class_mallita
    USE Aux
    use programs
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
    integer :: opcion_F
    character(len=120) :: archivo_F
    integer :: fid_F
    real(8) :: Fmacro(2,2)
    integer :: num_pasos_vibracion
    integer, allocatable :: vec_veces(:)
    real(8), allocatable :: vec_drmags(:)
    integer :: num_F
    integer :: j_F
    character(len=8) :: str_j_F
    ! ----------
    CHARACTER(LEN=120) :: filename, mallaname
    integer :: i,j,k,n,m
    INTEGER :: iStat, iStat1, iStat2
    real(8), allocatable :: r1(:,:), Ag(:,:), bg(:), dr(:,:)
    character(len=120) :: formato
    real(8), allocatable :: fuerzas_f(:,:), fuerzas_n(:,:)
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
            ! Leo los parametros de configuracion
            read(fid_cf,*) opcion_archivo, nombre_archivo
            ! Comienzo a simplificar
            if (opcion_archivo==1) then
                ! Caso de una sola malla
                call main_simplificar(nombre_archivo)
            elseif (opcion_archivo==2) then
                ! Caso de una lista de mallas en un archivo
                fid_lista_mallas = get_file_unit()
                open(unit=fid_lista_mallas, file=trim(nombre_archivo), status="old")
                read(fid_lista_mallas,*) nmallas
                ! Recorro las mallas de la lista
                do j_malla=1,nmallas
                    read(fid_lista_mallas,*) nombre_malla
                    write(*,*) "Simplificando malla:"
                    write(*,*) nombre_malla
                    call main_simplificar(nombre_malla)
                end do
                ! Cierro el archivo de la lista de mallas
                close(unit=fid_lista_mallas)
            else
                ! Si no encontre opcion 1 o 2, entonces hay algun error!!!
                write(*,*) "Error, para Simplificar, opcion_archivo debe ser 1 o 2, y es: ", opcion_archivo
                write(*,*) "En etiqueta: ", str_etiqueta
                write(*,*) "En Instruccion: ", str_instruccion
                stop
            end if
        ! FIN SIMPLIFICAR
        ! --------------------------------------------------------------------------
        ! EQUILIBRAR
        case ("Equilibrar")
            read(fid_cf,*) opcion_archivo, nombre_archivo
            read(fid_cf,*) num_pasos_vibracion
            allocate( vec_veces(num_pasos_vibracion) )
            allocate( vec_drmags(num_pasos_vibracion) )
            read(fid_cf,*) vec_veces
            read(fid_cf,*) vec_drmags
            read(fid_cf,*) opcion_F
            if (opcion_F==1) then
                ! caso de un solo F dado en configfile
                read(fid_cf,*) Fmacro
                call main_equilibrar(nombre_archivo, Fmacro, num_pasos_vibracion, vec_veces, vec_drmags)
            elseif (opcion_F==2) then
                ! caso de un archivo dando muchos F
                read(fid_cf,*) archivo_F
                fid_F = get_file_unit()
                open(unit=fid_F, file=trim(archivo_F), status="old")
                read(fid_F,*) num_F
                do j_F=1,num_F
                    read(fid_F,*) Fmacro
                    write(str_j_F, "(A1, I7.7)") "_", j_F ! 4.4 indica que el campo es de 4 y se usan como minimo 4, entonces imprime los ceros
                    call main_equilibrar(nombre_archivo, Fmacro, num_pasos_vibracion, vec_veces, vec_drmags, str_j_F)
                end do
            else
                write(*,*) "Error en opcion_F, deberia ser 1 o 2, y es: ", opcion_F
                write(*,*) "En etiqueta: ", str_etiqueta
                write(*,*) "En instruccion: ", str_instruccion
                stop
            end if
        close(unit=fid_F)
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
