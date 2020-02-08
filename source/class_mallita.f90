! ==============================================================================
MODULE class_mallita
! ==============================================================================
    USE Aux
    USE class_malla_completa, ONLY : MallaCom
    IMPLICIT NONE

    real(8), parameter :: delta = 1.d-4
    real(8), parameter :: delta21 = 1.d0 / (2.d0 * delta)
    real(8), parameter :: deltax(2) = delta * [1.d0, 0.d0]
    real(8), parameter :: deltay(2) = delta * [0.d0, 1.d0]

    !PRIVATE
    !PUBLIC :: MallaCom
    !PUBLIC :: leer_malla

    TYPE MallaSim
        ! diemensiones generales
        REAL(8) :: sidelen
        real(8) :: diamed
        integer :: ncapas
        ! nodos
        INTEGER :: nnods
        INTEGER, ALLOCATABLE :: tipos(:)
        REAL(8), ALLOCATABLE :: rnods0(:,:), rnods(:,:)
        LOGICAL, ALLOCATABLE :: mf(:), mi(:)
        ! fibras
        INTEGER :: nfibs
        INTEGER, ALLOCATABLE :: fibs(:,:)
        REAL(8), ALLOCATABLE :: letes0(:)
        REAL(8), ALLOCATABLE :: lamsr(:)
        REAL(8), ALLOCATABLE :: diams(:)
        ! parametros constitutivos
        INTEGER :: nparam
        REAL(8), ALLOCATABLE :: param(:) ! por ahora son los mismos para todas las fibras
        ! informacion de deformacion y tension
        logical :: status_deformed = .false.
        real(8) :: Fmacro(2,2)
        real(8), allocatable :: lams(:)
        real(8), allocatable :: tens(:)
        real(8) :: Tmacro(2,2)
    CONTAINS
        !
    END TYPE MallaSim

! ==============================================================================
CONTAINS
! ==============================================================================


! ================================================================================
! ================================================================================
SUBROUTINE Desde_MallaCom(macom, masim, nparcon, parcon)
    ! ----------
    IMPLICIT NONE
    TYPE(MallaCom), INTENT(IN) :: macom
    TYPE(MallaSim), INTENT(OUT) :: masim
    integer, intent(in) :: nparcon
    real(8), intent(in) :: parcon(nparcon)
    ! ----------
    INTEGER :: i,j,k ! contadors
    INTEGER :: f,s,n,n1,n2 ! indices de fibras, segmentos, nodos
    INTEGER :: nins_f ! para contar el numero de intersecciones por fibra
    integer :: nfibs ! numero de fibras en la malla simplificada
    integer :: nnods ! numero de nodos en la malla simplificada
    integer :: ifib, inod
    integer, allocatable :: oldnods(:) ! array para saber como se conectan los nodos de macom y de masim
    real(8) :: r1(2), r2(2), dr(2), loco0_s, loco0, lete0
    logical :: cond
    ! ----------

    ! ----------
    ! parametros
    masim%sidelen = macom%sidelen
    masim%diamed = macom%diamed
    masim%ncapas = macom%ncaps
    ! ----------
    write(*,*) "Contando el numero de nodos y fibras en la malla simplificada"
    ! primero voy a contar el numero de intersecciones en cada fibra de la malla completa
    ! para saber el numero de fibras y de nodos en la malla simplificada
    nfibs = 0
    nnods = 0
    allocate( oldnods(macom%nnods) ) ! aca voy chequeando si a cada nodo ya lo considere o si es una nueva interseccion
    oldnods = 0
    do f=1,macom%nfibs
        nins_f = 0
        ! el primer nodo de cada fibra siempre es un nodo nuevo de masim
        nnods = nnods + 1
        do j=1,macom%fibsne(f)-1 ! el ultimo no lo tengo en cuenta porque el ultimo nodo siempre es frontera, nunca interseccion
            s = macom%fibsje(macom%fibsie(f)-1+j)
            n1 = macom%segs(1,s)
            n2 = macom%segs(2,s)
            if (macom%tipos(n2) == 2) THEN
                nins_f = nins_f + 1
                ! chequeo si ya esta considerado o aun no
                if (oldnods(n2)==0) then
                    ! este nodo no lo sume aun a los nodos de masim
                    oldnods(n2) = 1
                    nnods = nnods + 1
!                else
!                    ! este nodo ya lo tuve en cuenta, no hace falta hacer nada
                end if
            end if
        end do
        ! sumo las fibras de masim que salen de esta fibra de macom segun las intersecciones que tuvo
        nfibs = nfibs + 1 + nins_f
        ! y sumo el ultimo nodo de la fibra que siempre es un nodo nuevo de masim
        nnods = nnods + 1
    end do
    write(*,*) "nnods: ", nnods, " nfibs: ", nfibs
    ! ----------
    ! ahora que se cuantas fibras y nodos voy a tener, adjudico memoria
    masim%nnods = nnods
    masim%nfibs = nfibs
    allocate( masim%tipos(nnods) )
    allocate( masim%rnods0(2,nnods) )
    allocate( masim%mf(nnods) )
    allocate( masim%mi(nnods) )
    allocate( masim%fibs(2,nfibs) )
    allocate( masim%letes0(nfibs) )
    allocate( masim%lamsr(nfibs) )
    allocate( masim%diams(nfibs) )
!    allocate( masim%param(10,nfibs) )
    ! ----------
    write(*,*) "Construyendo malla simplificada"
    ! ahora voy llenando los valores
    ifib = 0
    inod = 0
    oldnods = 0
    do f=1,macom%nfibs
        ! siempre que empieza una fibra de macom, empieza una nueva fibra de masim
        ifib = ifib + 1
        loco0 = 0.d0
        ! identifico el primer nodo de la fibra antes de hacer un loop por los demas segmentos
        s = macom%fibsje(macom%fibsie(f))
        n1 = macom%segs(1,s)
        ! este nodo siempre es un nodo nuevo de masim, entonces
        ! sumo este nodo a los nodos de masim y a la conectividad de la nueva fibra de masim
        inod = inod + 1
        masim%tipos(inod) = 1 ! frontera
        masim%rnods0(:,inod) = macom%rnods(:,n1)
        oldnods(n1) = inod ! guardo el indice que este nodo tendra en masim
        masim%fibs(1,ifib) = inod
        masim%diams(ifib) = macom%fibsds(f)
        ! ahora recorro todos los segmentos para ir chequeando todos los n1, asi comprendo todos los nodos de la fibra
        do j=1,macom%fibsne(f)
            ! identifico el segmento y sus nodos en macom
            s = macom%fibsje(macom%fibsie(f)-1+j)
            n1 = macom%segs(1,s)
            n2 = macom%segs(2,s)
            ! calculo long de seg para sumar a la long de contorno (loco0)
            r1 = macom%rnods(:,n1)
            r2 = macom%rnods(:,n2)
            dr = r2-r1
            loco0_s = dsqrt(sum(dr*dr))
            cond = iguales(loco0_s,0.d0)
            if (cond) then
                write(*,*) "longitud de segmento nula"
                stop
            end if
            loco0 = loco0 + loco0_s
            ! me fijo si llegue al final de la fibra de masim (si encuentro una interseccion o una frontera)
            if ((macom%tipos(n2) == 2).or.(macom%tipos(n2)==1)) THEN
                ! me tengo que fijar si este nodo que encontre ya lo tengo presente en oldnods
                if (oldnods(n2)==0) then
                    ! si no esta presente tengo que agregarlo a los nodos de masim
                    inod = inod + 1 !nuevo nodo de masim
                    masim%tipos(inod) = macom%tipos(n2) ! interseccion
                    masim%rnods0(:,inod) = macom%rnods(:,n2)
                    oldnods(n2) = inod
                end if
                ! estuviera ya presente o no, lo agrego a la conectividad de la fibra
                masim%fibs(2,ifib) =  oldnods(n2)
                ! calculo la longitud extremo-extremo de la fibra sim (lete0)
                r1 = masim%rnods0(:,masim%fibs(1,ifib))
                r2 = masim%rnods0(:,masim%fibs(2,ifib))
                dr = r2-r1
                lete0 = dsqrt(sum(dr*dr))
                masim%letes0(ifib) = lete0
                masim%lamsr(ifib) = loco0 / lete0
                ! si no llegue a un nodo frontera, comienzo una nueva fibra
                if (macom%tipos(n2)==2) then
                    ifib = ifib + 1
                    loco0 = 0.d0
                    masim%fibs(1,ifib) = oldnods(n2)
                    masim%diams(ifib) = macom%fibsds(f)
                end if
            end if
        end do
        nfibs = nfibs + 1 + nins_f
        nnods = nnods + 2 + nins_f
    end do
    ! ----------
    ! ya tengo los nodos (tipos y coordenadas), ahora hago las masks
    masim%mf = masim%tipos == 1
    masim%mi = masim%tipos == 2
    ! y copio rnods0 a rnods
    allocate( masim%rnods(2,nnods) )
    masim%rnods = masim%rnods0
    ! ----------
    ! parametros constitutivos
    allocate( masim%param(nparcon) )
    masim%nparam = nparcon
    masim%param = parcon
    ! ----------
    write(*,*) "Malla simplificada lista"

    ! ----------
END SUBROUTINE Desde_MallaCom
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
SUBROUTINE leer_mallita(masim, nomarch)
    IMPLICIT NONE
    CHARACTER(LEN=120), INTENT(IN) :: nomarch
    TYPE(MallaSim), INTENT(OUT) :: masim
    !
    INTEGER :: fid
    integer :: iStatus
    INTEGER :: i,j,k
    integer :: tipo
    real(8) :: r0(2), r(2)
    real(8) :: diam, lete0, lamr
    integer :: n0, n1
    CHARACTER(LEN=120) :: formato

    ! ----------
    fid = get_file_unit()
    OPEN(UNIT=fid, FILE=TRIM(nomarch), STATUS="OLD")
    ! ----------
    ! Parametros
    iStatus = FindStringInFile("*parametros", fid, .TRUE.)
    READ(fid,*) masim%sidelen
    READ(fid,*) masim%diamed
    read(fid,*) masim%ncapas
    READ(fid,*) masim%nparam
    allocate( masim%param(masim%nparam) )
    read(fid,*) masim%param
    ! ----------
    ! Nodos
    iStatus = FindStringInFile("*coordenadas", fid, .TRUE.)
    READ(fid,*) masim%nnods
    ALLOCATE( masim%rnods0(2,masim%nnods) )
    allocate( masim%rnods(2,masim%nnods) )
    ALLOCATE( masim%tipos(masim%nnods) )
    DO i=1,masim%nnods
        READ(fid,*) j, tipo, r0, r
        masim%rnods0(:,i) = r0
        masim%rnods(:,i) = r
        masim%tipos(i) = tipo
    END DO
    ! ----------
    ! Fibras
    iStatus = FindStringInFile("*fibras", fid, .TRUE.)
    READ(fid,*) masim%nfibs
    allocate( masim%fibs(2,masim%nfibs) )
    allocate( masim%diams(masim%nfibs) )
    allocate( masim%letes0(masim%nfibs) )
    allocate( masim%lamsr(masim%nfibs) )
    DO i=1,masim%nfibs
        READ(fid,*) j, diam, lete0, lamr, n0, n1
        masim%diams(i) = diam
        masim%letes0(i) = lete0
        masim%lamsr(i) = lamr
        masim%fibs(1,i) = n0+1 !+1 porque vengo de python que tiene base 0
        masim%fibs(2,i) = n1+1
    END DO
    ! ----------
    CLOSE(fid)
    ! ----------


END SUBROUTINE leer_mallita
! ================================================================================
! ================================================================================



! ================================================================================
! ================================================================================
SUBROUTINE escribir_mallita(masim, nomarch)
    IMPLICIT NONE
    CHARACTER(LEN=120), INTENT(IN) :: nomarch
    TYPE(MallaSim), INTENT(IN) :: masim
    !
    INTEGER :: fid
    INTEGER :: i
    CHARACTER(LEN=120) :: formato

    ! ----------
    fid = get_file_unit()
    OPEN(UNIT=fid, FILE=TRIM(nomarch), STATUS="REPLACE")
    ! ----------
    ! Parametros
    WRITE(fid,'(A11)') "*Parametros"
    WRITE(fid,'(E20.8E4)') masim%sidelen
    WRITE(fid,'(E20.8E4)') masim%diamed
    WRITE(fid,'(I10)') masim%ncapas
    WRITE(fid,'(I10)') masim%nparam
    WRITE(formato,'(A1,I0,A8)') "(", masim%nparam, "E20.8E4)"
    WRITE(fid,formato) masim%param
    ! ----------
    ! Deformacion
    if (masim%status_deformed) then
        write(fid,'(A12)') "*Deformacion"
        WRITE(fid,'(4E20.8E4)') masim%Fmacro
        WRITE(fid,'(4E20.8E4)') masim%Tmacro
    end if
    ! ----------
    ! Nodos
    WRITE(fid,'(A12)') "*Coordenadas"
    WRITE(fid,'(I12)') masim%nnods
    DO i=1,masim%nnods
        WRITE(fid, '(I12,I2,4E20.8E4)') i-1, masim%tipos(i), masim%rnods0(:,i), masim%rnods(:,i)
    END DO
    ! ----------
    ! Segmentos
    WRITE(fid,'(A7)') "*Fibras"
    WRITE(fid, '(I12)') masim%nfibs
    DO i=1,masim%nfibs
        WRITE(fid,'(I12, 3E20.8E4, 2I12)') i-1, masim%diams(i), masim%letes0(i), masim%lamsr(i), masim%fibs(:,i) - 1
    END DO
    ! ----------
    CLOSE(fid)
    ! ----------


END SUBROUTINE escribir_mallita
! ================================================================================
! ================================================================================



! ================================================================================
! ================================================================================
subroutine deformar_afin(masim, FF, r_afin)
    ! Deforma la mallita de manera Afin: r = F r0
    ! input: FF (tensor gradiente de deformaciones 2x2)
    !
    implicit none
    ! ----------
    type(MallaSim), intent(inout) :: masim
    real(8), intent(in) :: FF(2,2) ! tensor de deformaciones
    real(8), intent(out) :: r_afin(2,masim%nnods)
    ! ----------
    integer :: i,j,k
    ! ----------

    ! ----------
    r_afin = 0.d0
    do k=1,masim%nnods
        do i=1,2
            do j=1,2
                r_afin(i,k) = r_afin(i,k) + FF(i,j)*masim%rnods0(j,k)
            end do
        end do
    end do
    ! ----------
!    masim%rnods = r_afin
    ! ----------

    ! ----------
end subroutine deformar_afin
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
subroutine calcular_tension_fibra(masim, f, dr_f, tension, fuerzav)
    ! Calcula la tension ingenieril de una fibra (identificada mediante numeracion f)
    ! y la fuerza como vector (yendo del nodo inicial de la fibra al nodo final)
    ! dr_f es el vector de la fibra desde r_n0 hasta r_n1 (desde nodo inicial hasta nodo final)
    ! ----------
    implicit none
    type(MallaSim), intent(in) :: masim
    integer, intent(in) :: f
    real(8), intent(in) :: dr_f(2)
    real(8), intent(out) :: tension, fuerzav(2)
    ! ----------
    integer :: selector
    integer :: last
    ! ---
    real(8) :: k1, k2, Et, Eb
    real(8) :: diam, lete0, lamr, area
    real(8) :: lete, lam
    real(8) :: fuerza, fuerzar, tensionr
    ! ----------

    ! ----------
    ! Por ahora tengo un param con 3 valores: selector, k1 o Et, k2 o Eb
    selector = nint( masim%param(1) )
    last = 1
    select case (selector)

    case (1) ! se dan k1 y k2 para calcular la fuerza directamente
        k1 = masim%param(last+1)
        k2 = masim%param(last+2)
        diam = masim%diams(f)
        area = pi * diam * diam / 4.d0
        lete0 = masim%letes0(f)
        lamr = masim%lamsr(f)
        lete = dsqrt(sum(dr_f*dr_f))
        lam = lete / lete0
        if ( lam <= lamr ) then
            fuerza = k2*(lam - 1.0d0)
        else
            fuerzar = k2 * (lamr - 1.0d0)
            fuerza = fuerzar + k1 * (lam/lamr - 1.0d0)
        end if
        tension = fuerza / area
        fuerzav = fuerza * dr_f / lete
    case (2) ! se dan Et y Eb para calcular la tension ingenieril
        Et = masim%param(last+1)
        Eb = masim%param(last+2)
        lamr = masim%lamsr(f)
        diam = masim%diams(f)
        lete0 = masim%letes0(f)
        lete = dsqrt(sum(dr_f*dr_f))
        lam = lete / lete0
        if ( lam <= lamr ) then
            tension = Eb*(lam - 1.0d0) * diam
        else
            tensionr = Eb * (lamr - 1.0d0)
            tension = tensionr + Et * (lam/lamr - 1.0d0)
        end if
        fuerza = tension * area
        fuerzav = fuerza * dr_f / lete
    end select

    ! ----------

    ! ----------
end subroutine calcular_tension_fibra
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
subroutine calcular_tensiones_fibras(masim, r1, lams_fibs, tens_fibs)
    ! calcula las fuerzas de las fibras y las fuerzas nodales (resultantes)
    ! se calculan las fuerzas como vectores, entonces las dimensiones son (2,nfibs) y (2,nnods)
    implicit none
    type(MallaSim), intent(in) :: masim
    real(8), intent(inout) :: r1(2,masim%nnods)
    real(8), intent(out) :: lams_fibs(masim%nfibs)
    real(8), intent(out) :: tens_fibs(masim%nfibs)
    integer :: f
    integer :: n1f, n2f
    real(8) :: rn1f(2), rn2f(2), drf(2)
    real(8) :: lamf, tenf, fzafv(2)

    do f=1,masim%nfibs
        ! nodos de la fibra f
        n1f = masim%fibs(1,f)
        n2f = masim%fibs(2,f)
        ! posiciones de esos nodos
        rn1f = r1(:,n1f)
        rn2f = r1(:,n2f)
        ! vector fibra
        drf = rn2f - rn1f
        ! fuerza
        call calcular_tension_fibra(masim, f, drf, tenf, fzafv)
        tens_fibs(f) = tenf
        ! elongacion
        lamf = dsqrt(sum(drf*drf)) / masim%letes0(f)
        lams_fibs(f) = lamf
    end do

end subroutine calcular_tensiones_fibras
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
subroutine calcular_fuerzas(masim, r1, fzas_fibs, fzas_nods)
    ! calcula las fuerzas de las fibras y las fuerzas nodales (resultantes)
    ! se calculan las fuerzas como vectores, entonces las dimensiones son (2,nfibs) y (2,nnods)
    implicit none
    type(MallaSim), intent(in) :: masim
    real(8), intent(inout) :: r1(2,masim%nnods)
    real(8), intent(out) :: fzas_fibs(2,masim%nfibs)
    real(8), intent(out) :: fzas_nods(2,masim%nnods)
    integer :: f
    integer :: n1f, n2f
    real(8) :: rn1f(2), rn2f(2), drf(2)
    real(8) :: tenf, fzafv(2)

    fzas_nods = 0.d0
    do f=1,masim%nfibs
        ! nodos de la fibra f
        n1f = masim%fibs(1,f)
        n2f = masim%fibs(2,f)
        ! posiciones de esos nodos
        rn1f = r1(:,n1f)
        rn2f = r1(:,n2f)
        ! vector fibra
        drf = rn2f - rn1f
        ! fuerza
        call calcular_tension_fibra(masim, f, drf, tenf, fzafv)
        fzas_fibs(:,f) = fzafv
        fzas_nods(:,n1f) = fzas_nods(:,n1f) + fzafv
        fzas_nods(:,n2f) = fzas_nods(:,n2f) - fzafv
    end do

end subroutine calcular_fuerzas
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
subroutine vibrar_malla(masim, nveces, drmag, r1)
    ! Mueve los nodos nveces segun la direccion de la resultante
    ! Todos los nodos se mueven segun el valor drmag
    ! Cada vez se calculan las fuerzas de las fibras y la resultante
    ! Entonces la malla va "vibrando" hasta su posicion de equilibrio
    ! ----------
    implicit none
    ! ----------
    type(MallaSim), intent(inout) :: masim
    integer, intent(in) :: nveces
    real(8), intent(in) :: drmag
    real(8), intent(inout) :: r1(2,masim%nnods)
    ! ----------
    real(8) :: fzas_fibs(2,masim%nfibs)
    real(8) :: fzas_nods(2,masim%nnods)
    integer :: vez, n
    real(8) :: fza_n(2), fza_n_mag, dr_n(2)
    ! ----------

    do vez=1,nveces
        call calcular_fuerzas(masim, r1, fzas_fibs, fzas_nods)
        do n=1,masim%nnods
            if (masim%tipos(n) == 1) then
                ! nodo de frontera = nodo de Dirichlet
                r1(:,n) = r1(:,n) ! sentencia innecesaria pero por ahora la dejo para que quede claro
            else
                fza_n = fzas_nods(:,n)
                fza_n_mag = dsqrt(sum(fza_n*fza_n))
                if (fza_n_mag < 1.d-10) then
                    dr_n = 0.d0
                else
                    dr_n = drmag * fza_N / fza_n_mag
                end if
                r1(:,n) = r1(:,n) + dr_n(:)
            end if
        end do
    end do

    ! ----------
end subroutine vibrar_malla
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
subroutine calcular_equilibrio_vibracional(masim, npasos, vec_veces, vec_drmags, Fmacro)
    ! Forma de calcular el equilibrio moviendo los nodos de forma iterativa segun
    ! la direccion de la resultante en un valor prescripto de desplazamiento (drmag)
    ! Se hace de forma progresiva empezando con drmag grande y terminando con drmag chico
    ! ----------
    implicit none
    ! ----------
    type(MallaSim), intent(inout) :: masim
    integer, intent(in) :: npasos
    integer, intent(in) :: vec_veces(npasos)
    real(8), intent(in) :: vec_drmags(npasos)
    real(8), intent(in) :: Fmacro(2,2)
    ! ----------
    real(8) :: r1(2,masim%nnods)
    integer :: paso
    integer :: nveces
    real(8) :: drmag
    ! ----------
    real(8) :: lams_fibras(masim%nfibs)
    real(8) :: tens_fibras(masim%nfibs)
    real(8) :: Tmacro(2,2)
    ! ----------

    ! deformo afin
    call deformar_afin(masim, Fmacro, r1)
    ! encuentro perturbaciones respecto de la afin haciendo pasos de vibracion
    do paso=1,npasos
        nveces = vec_veces(paso)
        drmag = vec_drmags(paso)
        call vibrar_malla(masim, nveces, drmag, r1)
    end do
    ! recalculo el equilibrio para tener las tensiones
    call calcular_tensiones_fibras(masim, r1, lams_fibras, tens_fibras)
    call homogeneizacion(masim, r1, Tmacro)
    ! guardo la deformacion en la malla
    masim%rnods = r1
    masim%status_deformed = .true.
    masim%Fmacro = Fmacro
    masim%lams = lams_fibras
    masim%tens = tens_fibras
    masim%Tmacro = Tmacro
    ! y voila

    ! ----------
end subroutine calcular_equilibrio_vibracional
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
subroutine homogeneizacion(masim, r1, Tmacro)
    ! Calcula el tensor de tensiones como producto de la homogeneizacion de las
    ! tensiones de las fibras.
    ! ----------
    implicit none
    ! ----------
    type(MallaSim), intent(in) :: masim
    real(8), intent(in) :: r1(2,masim%nnods)
    real(8), intent(out) :: Tmacro(2,2) ! tension macroscopica de Cauchy
    ! ----------
    real(8) :: fzas_fibs(2,masim%nfibs)
    real(8) :: fzas_nods(2,masim%nnods)
    integer :: f
    real(8) :: d_f
    real(8) :: Area_f
    real(8) :: lete0_f
    real(8) :: Vol_f
    integer :: n1_f, n2_f
    real(8) :: r_n1_f(2), r_n2_f(2)
    real(8) :: dr_f(2)
    real(8) :: lete_f
    real(8) :: lam_f
    real(8) :: a_f(2)
    real(8) :: ten_f, vfza_f(2)
    real(8) :: Vol_RVE
    ! ----------

!    call calcular_fuerzas(masim, r1, fzas_fibs, fzas_nods)
    Tmacro = 0.d0
    do f=1,masim%nfibs
        d_f = masim%diams(f)
        Area_f = pi * d_f * d_f / 4.d0
        lete0_f = masim%letes0(f)
        Vol_f = Area_f * lete0_f
        ! nodos de la fibra f
        n1_f = masim%fibs(1,f)
        n2_f = masim%fibs(2,f)
        ! posiciones de esos nodos
        r_n1_f = r1(:,n1_f)
        r_n2_f = r1(:,n2_f)
        ! vector fibra
        dr_f = r_n2_f - r_n1_f
        ! elongacion fibra (extremo-extremo)
        lete_f = dsqrt(sum(dr_f*dr_f))
        lam_f = lete_f / lete0_f
        ! vector orientacion de la fibra
        a_f = dr_f/lete_f * lam_f
        ! tension fibra
        call calcular_tension_fibra(masim, f, dr_f, ten_f, vfza_f)
        ! homogeneizacion (queda poco eficiente, a mejorar luego de que funcione si hace falta)
        Tmacro(1,1) = Tmacro(1,1) + ten_f/lam_f * a_f(1) * a_f(1) * Vol_f
        Tmacro(1,2) = Tmacro(1,2) + ten_f/lam_f * a_f(1) * a_f(2) * Vol_f
        Tmacro(2,1) = Tmacro(2,1) + ten_f/lam_f * a_f(2) * a_f(1) * Vol_f
        Tmacro(2,2) = Tmacro(2,2) + ten_f/lam_f * a_f(2) * a_f(2) * Vol_f
    end do
    ! termino de acomodar la tension
    Vol_RVE = masim%ncapas * masim%diamed * masim%sidelen * masim%sidelen
    Tmacro = Tmacro / Vol_RVE

    ! ----------
end subroutine homogeneizacion
! ================================================================================
! ================================================================================


!! ================================================================================
!! ================================================================================
!subroutine calcular_equilibrio(masim, r1, maxiter, tol, iStatus)
!    ! subrutina obsoleta para calcular el equilibrio de la malla
!    implicit none
!    ! ----------
!    type(MallaSim), intent(in) :: masim
!    real(8), intent(out) :: r1(2,masim%nnods)
!    integer, intent(in) :: maxiter
!    real(8), intent(in) :: tol
!    integer, intent(out) :: iStatus
!    ! ----------
!    integer :: i,j,k,n,m
!    real(8) :: dr(2,masim%nnods)
!    real(8) :: A(2*masim%nnods,2*masim%nnods), b(2*masim%nnods), x(2*masim%nnods)
!    real(8) :: pivot(2*masim%nnods), residuo(2*masim%nnods), orp
!    integer :: concrit, iters
!    integer :: iStat
!    real(8) :: error_x
!    real(8) :: aux1, aux2, aux3(2*masim%nnods,2*masim%nnods)
!    integer :: iaux1, iaux2, iaux3
!    ! ----------
!
!    ! ----------
!    ! empiezo a resolver con el valor previo de coordenadas de la malla
!    r1 = masim%rnods
!    ! ----------
!    m = 2*masim%nnods
!    write(*,*) "Newton Raphson"
!    aux1 = 1.d-3
!    do i=1,maxiter
!        write(*,*) "iter: ", i
!        write(*,*) "Calculando A b"
!        call calcular_A_b(masim,r1,A,b)
!!        aux3 = 0.d0
!!        do iaux1=1,m
!!            iaux3 = 0
!!            do iaux2=1,m
!!                if (dabs(A(iaux1,iaux2)) > 1.d-6) then
!!                    iaux3 = iaux3 + 1
!!                    aux3(iaux1,iaux3) = A(iaux1,iaux2)
!!                end if
!!            end do
!!        end do
!
!!        do j=1,m
!!            A(j,1:j-1) = A(j,1:j-1)*0.5d0
!!            A(j,j+1:m) = A(j,j+1:m)*0.5d0
!!        end do
!        write(*,*) "Gauss-Seidel"
!!        call seidel(1,m,A,b,1.d0,x,residuo,iters,iStat)
!!        call jacobi(m, A, b, x, iStat)
!!        call gauss_seidel(m, A, b, x, 100, tol*1.d-1, iStat)
!        call diagonal_iterativo(m, A, b, x, iStat)
!!        call directo(m,A,b,x,iStat)
!        error_x = maxval(dabs(x))
!        aux2 = aux1*error_x
!        if (aux2>aux1) aux2=aux1
!        do j=1,m
!            if (dabs(x(j))>aux2) x(j) = dsign(aux2,x(j))
!!            if (dr(1,j)>aux2) dr(1,j) = aux2
!!            if (dr(2,j)>aux2) dr(2,j) = aux2
!        end do
!        dr = reshape(x,shape(dr))
!        r1 = r1 + dr
!        write(*,*) "i:", i, "error_x: ", error_x, "/ tol:", tol
!        if (error_x<tol) then
!            ! convergio
!            write(*,*) "Newton-Raphson convergio"
!            iStatus = 0
!            return
!        end if
!    end do
!    ! ----------
!    ! si llegue hasta aqui llegue a maxiter
!    iStatus = 1
!    write(*,*) "WARNING: maxiter alcanzado en metodo Newton-Raphson"
!    ! ----------
!
!    ! ----------
!end subroutine calcular_equilibrio
!! ================================================================================
!! ================================================================================


!! ================================================================================
!! ================================================================================
!subroutine calcular_A_b(masim,r1, A, b)
!    ! ----------
!    ! Calculo la matriz tangente A y el vector de carga b
!    ! el vector de carga esta dado por las fuerzas que ejercen las fibras en la configuracion dada por r1 (desplazamientos nodales)
!    ! la matriz tangente la calculo numericamente variando las coordenadas de los nodos (para eso hago un loop sobre las fibras)
!    ! ----------
!    implicit none
!    ! ----------
!    type(MallaSim), intent(in) :: masim
!    real(8), intent(in) :: r1(2,masim%nnods) ! coordenadas sobre las cuales voy a calcular la matriz tangente y el vector de cargas
!    real(8), intent(out) :: A(2*masim%nnods,2*masim%nnods)
!    real(8), intent(out) :: b(2*masim%nnods)
!    ! ----------
!    real(8) :: Af(2,2), bf(2) ! matriz y vectores elementales de una fibra que suman a los globales A y b
!    integer :: f, n1f, n2f ! indice de fibra y sus nodos inicial y final
!    real(8) :: rn1f(2), rn2f(2), drf(2) ! posiciones inicial y final de la fibra f, y su vector de fibra
!    real(8) :: fzaf, fzavf(2) ! fuerza escalar y vectorial de la fibra f
!    ! variaciones:
!    real(8) :: rn1f_mx(2), rn1f_px(2), rn1f_my(2), rn1f_py(2)
!    real(8) :: drf_n1mx(2), drf_n1px(2), drf_n1my(2), drf_n1py(2)
!    real(8) :: fzaf_mx, fzavf_mx(2)
!    real(8) :: fzaf_px, fzavf_px(2)
!    real(8) :: fzaf_my, fzavf_my(2)
!    real(8) :: fzaf_py, fzavf_py(2)
!    ! derivada numerica
!    real(8) :: dfza_dx(2), dfza_dy(2)
!    ! fila y columna para ensamblar
!    integer :: row, col
!    !
!    integer :: n ! indice de nodo
!    ! ----------
!
!    ! ----------
!    A = 0.d0
!    b = 0.d0
!    ! ----------
!    do f=1,masim%nfibs
!        ! nodos de la fibra f
!        n1f = masim%fibs(1,f)
!        n2f = masim%fibs(2,f)
!        ! posiciones de esos nodos
!        rn1f = r1(:,n1f)
!        rn2f = r1(:,n2f)
!        ! vector fibra
!        drf = rn2f - rn1f
!        ! y ahora con variaciones
!        rn1f_mx = rn1f - deltax
!        rn1f_px = rn1f + deltax
!        rn1f_my = rn1f - deltay
!        rn1f_py = rn1f + deltay
!        drf_n1mx = rn2f - rn1f_mx
!        drf_n1px = rn2f - rn1f_px
!        drf_n1my = rn2f - rn1f_my
!        drf_n1py = rn2f - rn1f_py
!        ! fuerzas
!        call calcular_fuerza_fibra(masim, f, drf, fzaf, fzavf)
!        call calcular_fuerza_fibra(masim, f, drf_n1mx, fzaf_mx, fzavf_mx)
!        call calcular_fuerza_fibra(masim, f, drf_n1px, fzaf_px, fzavf_px)
!        call calcular_fuerza_fibra(masim, f, drf_n1my, fzaf_my, fzavf_my)
!        call calcular_fuerza_fibra(masim, f, drf_n1py, fzaf_py, fzavf_py)
!        ! ahora la derivada numerica
!        dfza_dx = (fzavf_px - fzavf_mx) * delta21
!        dfza_dy = (fzavf_py - fzavf_my) * delta21
!        ! matriz local y vector local
!        Af(1,1) = dfza_dx(1)
!        Af(2,2) = dfza_dy(2)
!        Af(1,2) = dfza_dx(2)
!        Af(2,1) = Af(1,2)
!        bf(:) = - fzavf
!        ! ----------
!        ! ya tengo la derivada que aporta esta fibra
!        ! ahora ensamblo en matriz y vector globales
!        ! ----------
!        ! primero el vector de cargas
!        row = (n1f-1)*2 + 1 ! apunta a donde se ubica el primer dof del nodo n1f
!        b(row:row+1) = b(row:row+1) + bf
!        row = (n2f-1)*2 + 1 ! apunta a donde se ubica el primer dof del nodo n2f
!        b(row:row+1) = b(row:row+1) - bf
!        ! ahora la matriz  en 4 partes
!        ! n1f n1f
!        row = (n1f-1)*2 + 1
!        col = (n1f-1)*2 + 1
!        A(row:row+1,col:col+1) = A(row:row+1,col:col+1) + Af
!        ! n2f n2f
!        row = (n2f-1)*2 + 1
!        col = (n2f-1)*2 + 1
!        A(row:row+1,col:col+1) = A(row:row+1,col:col+1) + Af
!        ! cruzadas
!        row = (n1f-1)*2 + 1
!        col = (n2f-1)*2 + 1
!        A(row:row+1,col:col+1) = A(row:row+1,col:col+1) - Af
!        !
!        row = (n2f-1)*2 + 1
!        col = (n1f-1)*2 + 1
!        A(row:row+1,col:col+1) = A(row:row+1,col:col+1) - Af
!        ! presto
!    end do
!    ! ----------
!    ! por ultimo las condiciones de dirichlet
!    do n=1,masim%nnods
!        if (masim%tipos(n)==1) then
!            ! nodo frontera es dirichlet
!            ! debo ensamblar la ecuacion: 1.d0 * dr_i = 0.0d0
!            row = (n-1)*2 + 1
!            A(row:row+1,:) = 0.d0 ! hago cero las dos filas del nodo n
!            A(:,row:row+1) = 0.d0 ! y tambien las columnas ya que este dr siempre va a ser cero
!            A(row,row) = 1.d0  ! pongo un uno en la diagonal
!            A(row+1,row+1) = 1.d0  ! idem
!            b(row:row+1) = 0.d0 ! hago cero el vector de carga
!        end if
!    end do
!    ! ----------
!
!    ! ----------
!end subroutine calcular_A_b
!! ================================================================================
!! ================================================================================


! ==============================================================================
END MODULE class_mallita
! ==============================================================================
