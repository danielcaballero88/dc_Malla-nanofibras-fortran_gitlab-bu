module programs

use Aux
USE class_malla_completa, ONLY : MallaCom
use class_mallita, only : MallaSim

contains

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
subroutine main_simplificar(filename_malla_in, nparamcon, paramcon)
    use class_malla_completa
    use class_mallita
    implicit none
    CHARACTER(LEN=120), intent(in) :: filename_malla_in
    integer, intent(in) :: nparamcon
    real(8), intent(in) :: paramcon(nparamcon)
    character(len=120) :: filename_malla_in2, filename_malla_out
    type(MallaCom) :: mc
    type(MallaSim) :: ms

    if (trim(filename_malla_in) == "default") then
        filename_malla_in2 = "Malla_i.txt"
    else
        filename_malla_in2 = "_i"
        call modify_txt_filename(filename_malla_in, filename_malla_in2)
    end if


    write(*,*) "Leer malla intersectada y generar malla simplificada:"
    call leer_malla(mc, filename_malla_in2)
    call Desde_MallaCom(mc, ms, nparamcon, paramcon)

    write(*,*) "Escribiendo mallita"
    filename_malla_out = "_s"
    call modify_txt_filename(filename_malla_in2, filename_malla_out)
    CALL escribir_mallita(ms, filename_malla_out)
    write(*,*) "Malla simplificada OK"

end subroutine main_simplificar
! ==========================================================================

! ==========================================================================
subroutine main_equilibrar(filename_malla_in, nparcon, parcon, Fmacro, num_pasos, lista_veces, lista_drmags, str_num_output_opt)
    ! Calcula el equilibrio elastico de una malla dado un tensor F macroscopico (Fmacro)
    ! ----------
    use class_mallita
    implicit none
    ! ----------
    CHARACTER(LEN=120), intent(in) :: filename_malla_in
    integer, intent(in) :: nparcon
    real(8), intent(in) :: parcon(nparcon)
    real(8), intent(in) :: Fmacro(2,2)
    integer, intent(in) :: num_pasos
    integer, intent(in) :: lista_veces(num_pasos)
    real(8), intent(in) :: lista_drmags(num_pasos)
    character(len=8), intent(in), optional :: str_num_output_opt
    ! ----------
    character(len=120) :: filename_malla_in2, filename_malla_out, aux_string
    type(MallaSim) :: ms
!    real(8), allocatable :: r1(:,:)
    integer :: n, iStat1
    ! ----------

    if (trim(filename_malla_in) == "default") then
        filename_malla_in2 = "Malla_i_s.txt"
    else
        filename_malla_in2 = "_i_s"
        call modify_txt_filename(filename_malla_in, filename_malla_in2)
    end if

    write(*,*) "Calculando equilibrio"
    call leer_mallita(ms, filename_malla_in2, nparcon, parcon)
    n = ms%nnods
!    allocate( r1(2,n) )
!    Fmacro(1,:) = [1.2d0, 0.d0]
!    Fmacro(2,:) = [0.0d0, 1.d0]
    call calcular_equilibrio_vibracional(ms, num_pasos, lista_veces, lista_drmags, Fmacro)
!        call calcular_equilibrio(ms,r1,10000,1.d-1,iStat1)

!    ms%rnods = r1
    filename_malla_out = "_e"
    call modify_txt_filename(filename_malla_in2, filename_malla_out)
    if (present(str_num_output_opt)) then
        aux_string = str_num_output_opt
        call modify_txt_filename(filename_malla_out, aux_string)
        filename_malla_out = aux_string
    end if
    call escribir_mallita(ms, filename_malla_out)
    write(*,*) "Equilibrio calculado OK"

    ! ----------
end subroutine main_equilibrar
! ==========================================================================


! ==========================================================================
subroutine main_traccion(filename_malla_in, nparcon, parcon, num_pasos, lista_veces, lista_drmags, dtime, dotF11, dotF22, F11fin, filename_curva, opcion_save, dF_save)
    ! Simula un ensayo de traccion con un esquema explicito
    ! imponiendo tasas de deformacion axial y transversal
    ! ----------
    use class_mallita
    implicit none
    ! ----------
    CHARACTER(LEN=120), intent(in) :: filename_malla_in
    integer, intent(in) :: nparcon
    real(8), intent(in) :: parcon(nparcon)
    integer, intent(in) :: num_pasos
    integer, intent(in) :: lista_veces(num_pasos)
    real(8), intent(in) :: lista_drmags(num_pasos)
    real(8), intent(in) :: dtime
    real(8), intent(in) :: dotF11
    real(8), intent(in) :: dotF22
    real(8), intent(in) :: F11fin
    CHARACTER(LEN=120), intent(in) :: filename_curva
    integer, intent(in) :: opcion_save
    real(8), intent(in) :: dF_save
    ! ----------
    character(len=120) :: filename_malla_in2, filename_malla_out
    real(8) :: F11ini
    integer :: fid_curva
    real(8) :: time
    real(8) :: Fmacro(2,2)
    real(8) :: lamps
    type(MallaSim) :: ms
    integer :: nsaves
    real(8), allocatable :: lista_saves_F(:)
    logical, allocatable :: lista_saves_if(:)
    logical :: listo_saves = .false.
    integer :: isave
    ! ----------
    integer :: f
    ! ----------

    write(*,*) "Empezando Traccion"
    ! Preparo la lista de saves si es que hay
    if (opcion_save==1) then
        nsaves =  int((F11fin-1.0d0) / dF_save) + 1
        allocate( lista_saves_F(nsaves) )
        allocate( lista_saves_if(nsaves) )
        do isave=1,nsaves
            lista_saves_F(isave) = 1. + dfloat(isave - 1)*dF_save
        end do
        lista_saves_if = .false.
    end if
    ! Leo la malla
    write(*,*) "Leyendo Malla:"
    if (trim(filename_malla_in) == "default") then
        filename_malla_in2 = "Malla_i_s.txt"
    else
!        filename_malla_in2 = "_i_s"
!        call modify_txt_filename(filename_malla_in, filename_malla_in2)
        filename_malla_in2 = trim(filename_malla_in)
    end if
    write(*,*) "archivo: ", filename_malla_in2
    call leer_mallita(ms, filename_malla_in2, nparcon, parcon)
    if (ms%status_deformed) then
        ! Si la malla esta previamente deformada, empiezo a trabajar desde alli
        Fmacro = ms%Fmacro
        F11ini = Fmacro(1,1)
        lista_saves_if = (F11ini > lista_saves_F)
        isave = count(lista_saves_if) + 1
        ! Abro un archivo viejo para continuar la curva constitutiva
        fid_curva = get_file_unit()
        open(unit=fid_curva, file=trim(filename_curva), status="old", position="append", action="write")
    else
        ! Si la malla esta virge, empiezo desde deformacion nula
        Fmacro = reshape(source=[1.d0, 0.d0, 0.d0, 1.d0], shape=shape(Fmacro))
        isave = 1
        ! Abro un archivo nuevo para escribir la curva constitutiva
        fid_curva = get_file_unit()
        open(unit=fid_curva, file=trim(filename_curva), status="replace")
    end if

    ! Comienzo esquema temporal
    time = 0.d0
    do while ( Fmacro(1,1) .le. F11fin )
        ! Calculo el equilibrio de la malla para el Fmacro dado en este paso de tiempo
        call calcular_equilibrio_vibracional(ms, num_pasos, lista_veces, lista_drmags, Fmacro)
        call escribir_mallita(ms, "malla_test.txt")
        ! Guardo en archivo la informacion constitutiva para este step
        write(*,"(2E20.8E4)") Fmacro(1,1), ms%Tmacro(1,1)
        write(fid_curva,"(8E20.8E4)") Fmacro, ms%Tmacro
        ! Calculo plasticidad y/o rotura de fibras
        call calcular_plasticidad_rotura(ms, dtime)
        ! Me fijo si guardo la malla o no
        if (.not. listo_saves) then
            if ( dabs( Fmacro(1,1) - lista_saves_F(isave) ) < dotF11*dtime ) then
                write(filename_malla_out,"(A6, I4.4)") "_save_", isave
                call modify_txt_filename(filename_malla_in2, filename_malla_out)
                call escribir_mallita(ms, filename_malla_out)
                isave = isave + 1
                if (isave > nsaves) listo_saves = .true.
            end if
        end if
        ! Incremento tiempo y deformacion para siguiente paso
        Fmacro(1,1) = Fmacro(1,1) + dotF11*dtime
        Fmacro(2,2) = Fmacro(2,2) + dotF22*dtime
        time = time + dtime
    end do
    ! Cierro el archivo donde guarda la curva constitutiva
    close(unit=fid_curva)

    ! ----------
end subroutine main_traccion
! ==========================================================================


end module
