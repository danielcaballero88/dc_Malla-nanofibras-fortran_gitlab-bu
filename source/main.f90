
program hello
    USE class_malla_completa
    use class_mallita
    USE Aux
    implicit none
    ! ==========================================================================
    CHARACTER(LEN=120) :: inFile
    CHARACTER(LEN=120) :: outFile
    integer :: fid
    ! ==========================================================================
    ! Current working directory
    CHARACTER(LEN=255) :: cwd
    ! ==========================================================================
    TYPE(MallaCom) :: MC, MC2
    TYPE(MallaSim) :: ms, ms2
    integer :: nmallas
    CHARACTER(LEN=120) :: filename, mallaname
    integer :: i,j,k,n,m
    INTEGER :: iStat1, iStat2
    real(8), allocatable :: r1(:,:), Ag(:,:), bg(:), dr(:,:)
    character(len=120) :: formato
    real(8), allocatable :: fuerzas_f(:,:), fuerzas_n(:,:)
    real(8) :: Fmacro(2,2)
    integer :: opcion ! 1=intersectar, 2=simplificar, 3=equilibrio
    ! ==========================================================================

    print *, "Multiscale Nanofibers Mesh RVE 2.0"
    CALL GETCWD(cwd)
    PRINT *, "Current Working Directory:", cwd

!    filename = "000_mallas_filenames.txt"
!    fid = get_file_unit()
!    OPEN(unit=fid, file=trim(filename), status="old")
!    read(fid,*) nmallas
!    CLOSE(unit=fid)

    ! ==========
    opcion = 2
    ! ==========
    select case (opcion)
    ! ==========
    case (1)
    ! ==========
        write(*,*) "Leer malla, intersectar fibras y reescribir:"
        filename = "Malla.txt"
        CALL leer_malla(MC, filename)
        ! Hago la interseccion muchas veces porque cada vez tengo la limitacion de no cortar al mismo segmento dos veces
        i = 0
        write(*,*) "Intersectando fibras"
        DO WHILE (.true.)
            i = i+1
            WRITE(*,'(I4)', ADVANCE='no') i
            CALL intersectar_fibras_3(MC, MC2, .FALSE., iStat1) ! dentro de la misma capa
            MC = MC2
            CALL intersectar_fibras_3(MC, MC2, .TRUE., iStat2) ! con capas adyacentes
            MC = MC2
            IF ( (iStat1 == 1).AND.(iStat2 == 1) )  EXIT
        END DO
        write(*,*)

        write(*,*) "Escribiendo malla intersectada"
        filename = "Malla_i.txt"
        CALL escribir_malla(mc, filename)
        write(*,*) "Malla intersectada OK"
    ! ==========
    case (2)
    ! ==========
        write(*,*) "Leer malla intersectada y generar malla simplificada:"
        filename = "Malla_i.txt"
        call leer_malla(mc, filename)
        call Desde_MallaCom(mc, ms, 2, [10.d0, .1d0])
        write(*,*) "Escribiendo mallita"
        filename = "Malla_s.txt"
        CALL escribir_mallita(ms, filename)
        write(*,*) "Malla simplificada OK"
    ! ==========
    case (3)
    ! ==========
        write(*,*) "Calculando equilibrio"
        filename = "Malla_s.txt"
        call leer_mallita(ms, filename)
        n = ms%nnods
        allocate( r1(2,n) )
        Fmacro(1,:) = [1.2d0, 0.d0]
        Fmacro(2,:) = [0.0d0, 1.d0]
        call deformar_afin(ms, Fmacro)
        call calcular_equilibrio(ms,r1,10000,1.d-1,iStat1)

        ms%rnods = r1

        filename = "Malla_sd.txt"
        call escribir_mallita(ms, filename)
        write(*,*) "Equilibrio calculado OK"
    ! ==========
    end select
    ! ==========

end program


