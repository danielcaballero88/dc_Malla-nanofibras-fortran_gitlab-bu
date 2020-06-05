
program Malla_Nanofibras_Fortran
    USE class_malla_completa
    use class_mallita
    USE Aux
    use programs
    use class_instruccion, only : instruccion, read_from_configfile
    implicit none
    ! ==========================================================================
    ! ==========================================================================
    ! Current working directory
    CHARACTER(LEN=255) :: cwd
    ! ---
    ! Instrucciones
    integer :: num_instruc
    integer, allocatable :: lista_instruc(:)
    type(instruccion), allocatable :: lista_instrucciones(:)
    ! ==========================================================================
    TYPE(MallaCom) :: MC, MC2
    TYPE(MallaSim) :: ms, ms2
    CHARACTER(LEN=120) :: configfile
    integer :: fid_cf
    integer :: fid_lista_mallas
    integer :: nmallas
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
    integer :: numparam
    real(8), allocatable :: param(:)
    ! ----------
    ! variables para instruccion "Equilibrar"
    integer :: opcion_Fmacro
    character(len=120) :: archivo_Fmacro
    integer :: fid_Fmacro
    real(8) :: Fmacro(2,2)
    integer :: num_pasos_vibracion
    integer, allocatable :: vec_veces(:)
    real(8), allocatable :: vec_drmags(:)
    real(8) :: fuerza_ref, fuerza_tol
    integer :: num_Fmacro
    real(8), allocatable :: vec_Fmacro(:,:,:)
    integer :: j_F
    character(len=8) :: str_j_F
    ! variables para instruccion "Traccion"
    real(8) :: delta_t, dot_F11, dot_F22, F11_fin
    CHARACTER(LEN=120) :: filename_curvacon
    integer :: opcion_guardar
    real(8) :: dF_guardar
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
    ! -- NEW --
    allocate( lista_instrucciones(num_instruc) )
    do j_instr=1,num_instruc
        i_etiqueta = lista_instruc(j_instr)
        WRITE(str_etiqueta,'(A2,I1)') "* ", i_etiqueta
        iStat = read_from_configfile(lista_instrucciones(j_instr), fid_cf, str_etiqueta)
    end do
    ! -- --- --
    ! Ya que estamos tambien busco y leo los parametros constitutivos
    iStat = FindStringInFile("* Parametros Constitutivos", fid_cf, .false.)
    if (iStat==0) then
        read(fid_cf,*) numparam
        allocate( param(numparam) )
        read(fid_cf,*) param
    end if
    CLOSE(unit=fid_cf)
    ! --------------------------------------------------------------------------
    ! Recorro las etiquetas de las instrucciones a realizar
    do j_instr=1,num_instruc
        ! --------------------------------------------------------------------------
        ! Busco cada etiqueta en configfile y leo la string de instruccion (es el identificador: Intersectar, Simplificar o Equilibrar)
        i_etiqueta = lista_instruc(j_instr)
        WRITE(str_etiqueta,'(A2,I1)') "* ", i_etiqueta
        OPEN(unit=fid_cf, file=trim(configfile), status="old")
        iStat = FindStringInFile(str_etiqueta, fid_cf, .true.)
        read(fid_cf,*) str_instruccion
        ! --------------------------------------------------------------------------
        ! Luego, segun la instruccion me fijo que hay que hacer
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
                call main_simplificar(nombre_archivo, numparam, param)
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
                    call main_simplificar(nombre_malla, numparam, param)
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
            ! Leo parametros
            read(fid_cf,*) opcion_archivo, nombre_archivo
            read(fid_cf,*) num_pasos_vibracion
            if (allocated(vec_veces)) deallocate(vec_veces)
            if (allocated(vec_drmags)) deallocate(vec_drmags)
            if (num_pasos_vibracion>0) then
                allocate( vec_veces(num_pasos_vibracion) )
                allocate( vec_drmags(num_pasos_vibracion) )
            end if
            read(fid_cf,*) vec_veces
            read(fid_cf,*) vec_drmags
            read(fid_cf,*) fuerza_ref, fuerza_tol
            read(fid_cf,*) opcion_Fmacro
            ! Armo array de deformaciones
            if (opcion_Fmacro==1) then
                ! Caso de un solo Fmacro
                read(fid_cf,*) Fmacro
                allocate( vec_Fmacro(2,2,1) )
                vec_Fmacro(:,:,1) = Fmacro
            elseif (opcion_Fmacro==2) then
                ! Caso de un archivo con una lista de deformaciones
                read(fid_cf,*) archivo_Fmacro
                fid_Fmacro = get_file_unit()
                open(unit=fid_Fmacro, file=trim(archivo_Fmacro), status="old")
                read(fid_Fmacro,*) num_Fmacro
                allocate( vec_Fmacro(2,2,num_Fmacro) )
                do j_F=1,num_Fmacro
                    read(fid_Fmacro,*) vec_Fmacro(:,:,j_F)
                end do
            else
                write(*,*) "Error en opcion_Fmacro, deberia ser 1 o 2, y es: ", opcion_Fmacro
                write(*,*) "En etiqueta: ", str_etiqueta
                write(*,*) "En instruccion: ", str_instruccion
                stop
            end if
            ! Ahora resuelvo el equilibrio para las deformaciones que tengo
            if (opcion_archivo==1) then
                ! Caso de una sola malla
                ! Recorro los Fmacro de la lista para hacer todos los equilibrios
                do j_F=1,num_Fmacro
                    Fmacro = vec_Fmacro(:,:,j_F)
                    write(str_j_F, "(A1, I7.7)") "_", j_F ! 7.7 indica que el campo es de 7 y se usan como minimo 7, entonces imprime los ceros
                    call main_equilibrar(nombre_archivo, numparam, param, Fmacro, num_pasos_vibracion, vec_veces, vec_drmags, fuerza_ref, fuerza_tol, str_j_F)
                end do
            elseif (opcion_archivo==2) then
                ! Caso de una lista de mallas en un archivo
                fid_lista_mallas = get_file_unit()
                open(unit=fid_lista_mallas, file=trim(nombre_archivo), status="old")
                read(fid_lista_mallas,*) nmallas
                ! Recorro las mallas de la lista
                do j_malla=1,nmallas
                    read(fid_lista_mallas,*) nombre_malla
                    write(*,*) "Equilibrando malla:"
                    write(*,*) nombre_malla
                    ! Recorro los Fmacro de la lista para hacer todos los equilibrios
                    do j_F=1,num_Fmacro
                        Fmacro = vec_Fmacro(:,:,j_F)
                        write(str_j_F, "(A1, I7.7)") "_", j_F ! 7.7 indica que el campo es de 7 y se usan como minimo 7, entonces imprime los ceros
                        call main_equilibrar(nombre_archivo, numparam, param, Fmacro, num_pasos_vibracion, vec_veces, vec_drmags, fuerza_ref, fuerza_tol, str_j_F)
                    end do
                end do
                ! Cierro el archivo de la lista de mallas
                close(unit=fid_lista_mallas)
            else
                ! Si no encontre opcion 1 o 2, entonces hay algun error!!!
                write(*,*) "Error, opcion_archivo debe ser 1 o 2, y es: ", opcion_archivo
                write(*,*) "En etiqueta: ", str_etiqueta
                write(*,*) "En Instruccion: ", str_instruccion
                stop
            end if
            ! Cierro archivo de deformaciones
            close(unit=fid_Fmacro)
        ! FIN EQUILIBRAR
        ! --------------------------------------------------------------------------
        ! TRACCION
        case ("Traccion")
            ! Leo parametros
            read(fid_cf,*) opcion_archivo, nombre_archivo
            read(fid_cf,*) num_pasos_vibracion
            if (allocated(vec_veces)) deallocate(vec_veces)
            if (allocated(vec_drmags)) deallocate(vec_drmags)
            if (num_pasos_vibracion>0) then
                allocate( vec_veces(num_pasos_vibracion) )
                allocate( vec_drmags(num_pasos_vibracion) )
            end if
            read(fid_cf,*) vec_veces
            read(fid_cf,*) vec_drmags
            read(fid_cf,*) fuerza_ref, fuerza_tol
            read(fid_cf,*) delta_t, dot_F11, dot_F22, F11_fin
            read(fid_cf,*) filename_curvacon
            read(fid_cf,*) opcion_guardar, dF_guardar
            !
            if (opcion_archivo==1) then
                ! Caso de una sola malla
                ! Recorro los Fmacro de la lista para hacer todos los equilibrios
                call main_traccion(nombre_archivo, numparam, param, num_pasos_vibracion, vec_veces, vec_drmags, fuerza_ref, fuerza_tol, delta_t, dot_F11, dot_F22, F11_fin, filename_curvacon, opcion_guardar, dF_guardar)
            elseif (opcion_archivo==2) then
                ! Caso de una lista de mallas en un archivo
                fid_lista_mallas = get_file_unit()
                open(unit=fid_lista_mallas, file=trim(nombre_archivo), status="old")
                read(fid_lista_mallas,*) nmallas
                ! Recorro las mallas de la lista
                do j_malla=1,nmallas
                    read(fid_lista_mallas,*) nombre_malla
                    write(*,*) "Traccion sobre malla:"
                    write(*,*) nombre_malla
                    ! Recorro los Fmacro de la lista para hacer todos los equilibrios
                    call main_traccion(nombre_malla, numparam, param, num_pasos_vibracion, vec_veces, vec_drmags, fuerza_ref, fuerza_tol, delta_t, dot_F11, dot_F22, F11_fin, filename_curvacon, opcion_guardar, dF_guardar)
                end do
                ! Cierro el archivo de la lista de mallas
            else
                ! Si no encontre opcion 1 o 2, entonces hay algun error!!!
                write(*,*) "Error, opcion_archivo debe ser 1, y es: ", opcion_archivo
                write(*,*) "En etiqueta: ", str_etiqueta
                write(*,*) "En Instruccion: ", str_instruccion
                stop
            end if
        ! FIN TRACCION
        ! --------------------------------------------------------------------------
        ! UNIAXIAL
        case ("Uniaxial")
            ! Leo parametros
            read(fid_cf,*) opcion_archivo, nombre_archivo
            read(fid_cf,*) num_pasos_vibracion
            if (allocated(vec_veces)) deallocate(vec_veces)
            if (allocated(vec_drmags)) deallocate(vec_drmags)
            if (num_pasos_vibracion>0) then
                allocate( vec_veces(num_pasos_vibracion) )
                allocate( vec_drmags(num_pasos_vibracion) )
            end if
            read(fid_cf,*) vec_veces
            read(fid_cf,*) vec_drmags
            read(fid_cf,*) fuerza_ref, fuerza_tol
            read(fid_cf,*) delta_t, dot_F11, dot_F22, F11_fin
            read(fid_cf,*) filename_curvacon
            read(fid_cf,*) opcion_guardar, dF_guardar
            !
            if (opcion_archivo==1) then
                ! Caso de una sola malla
                ! Recorro los Fmacro de la lista para hacer todos los equilibrios
                call main_uniaxial(nombre_archivo, numparam, param, num_pasos_vibracion, vec_veces, vec_drmags, fuerza_ref, fuerza_tol, delta_t, dot_F11, dot_F22, F11_fin, filename_curvacon, opcion_guardar, dF_guardar)
!            elseif (opcion_archivo==2) then
!                stop
            else
                ! Si no encontre opcion 1 o 2, entonces hay algun error!!!
                write(*,*) "Error, opcion_archivo debe ser 1, y es: ", opcion_archivo
                write(*,*) "En etiqueta: ", str_etiqueta
                write(*,*) "En Instruccion: ", str_instruccion
                stop
            end if
        ! FIN UNIAXIAL
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


end program Malla_Nanofibras_Fortran
! ==========================================================================
! ==========================================================================
! ==========================================================================
