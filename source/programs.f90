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
        filename_malla_in2 = "_i"
        call modify_txt_filename(filename_malla_in, filename_malla_in2)
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
subroutine main_equilibrar(filename_malla_in, Fmacro, num_pasos, lista_veces, lista_drmags, str_num_output_opt)
    ! Calcula el equilibrio elastico de una malla dado un tensor F macroscopico (Fmacro)
    ! ----------
    use class_mallita
    implicit none
    ! ----------
    CHARACTER(LEN=120), intent(in) :: filename_malla_in
    real(8), intent(in) :: Fmacro(2,2)
    integer, intent(in) :: num_pasos
    integer, intent(in) :: lista_veces(num_pasos)
    real(8), intent(in) :: lista_drmags(num_pasos)
    character(len=8), intent(in), optional :: str_num_output_opt
    ! ----------
    character(len=120) :: filename_malla_in2, filename_malla_out, aux_string
    type(MallaSim) :: ms
    real(8), allocatable :: r1(:,:)
    integer :: n, iStat1
    ! ----------

    if (trim(filename_malla_in) == "default") then
        filename_malla_in2 = "Malla_i_s.txt"
    else
        filename_malla_in2 = "_i_s"
        call modify_txt_filename(filename_malla_in, filename_malla_in2)
    end if

    write(*,*) "Calculando equilibrio"
    call leer_mallita(ms, filename_malla_in2)
    n = ms%nnods
    allocate( r1(2,n) )
!    Fmacro(1,:) = [1.2d0, 0.d0]
!    Fmacro(2,:) = [0.0d0, 1.d0]
    call deformar_afin(ms, Fmacro, r1)
    call calcular_equilibrio_vibracional(ms, num_pasos, lista_veces, lista_drmags, r1)
!        call calcular_equilibrio(ms,r1,10000,1.d-1,iStat1)

    ms%rnods = r1
    filename_malla_out = "_e"
    call modify_txt_filename(filename_malla_in2, filename_malla_out)
    if (present(str_num_output_opt)) then
        aux_string = str_num_output_opt
        call modify_txt_filename(filename_malla_out, aux_string)
        filename_malla_out = aux_string
    end if
    call escribir_mallita(ms, filename_malla_out)
    write(*,*) "Equilibrio calculado OK"
end subroutine main_equilibrar
! ==========================================================================



end module
