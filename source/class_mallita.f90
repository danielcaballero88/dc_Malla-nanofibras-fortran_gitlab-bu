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
        real(8), allocatable :: lamps(:)
        logical, allocatable :: brokens(:)
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
    real(8) :: L_2
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
    allocate( masim%lamps(nfibs) )
    allocate( masim%brokens(nfibs) )
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
    ! tambien agrego informacion sobre deformacion plastica y rotura de fibras
    masim%lamps = 1.d0 ! sin plasticidad
    masim%brokens = .false. ! todas sanas
    ! ----------
    ! ya tengo los nodos (tipos y coordenadas), ahora hago las masks
    masim%mf = masim%tipos == 1
    masim%mi = masim%tipos == 2
    ! Desplazo las coordenadas para dejar el cero en el medio del rve
    L_2 = 0.5d0 * masim%sidelen
    masim%rnods0 = masim%rnods0 - L_2
    ! y copio rnods0 a rnods
    allocate( masim%rnods(2,nnods) )
    masim%rnods = masim%rnods0
    ! ----------
    ! parametros constitutivos
    allocate( masim%param(nparcon) )
    masim%nparam = nparcon
    masim%param = parcon
    ! ----------
    ! ----------
    write(*,*) "Malla simplificada lista"

    ! ----------
END SUBROUTINE Desde_MallaCom
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
SUBROUTINE leer_mallita(masim, nomarch, numparamcon, paramcon)
    IMPLICIT NONE
    CHARACTER(LEN=120), INTENT(IN) :: nomarch
    TYPE(MallaSim), INTENT(OUT) :: masim
    integer, intent(in) :: numparamcon
    real(8), intent(in) :: paramcon(numparamcon)
    !
    INTEGER :: fid
    integer :: iStatus
    INTEGER :: i,j,k
    integer :: tipo
    real(8) :: r0(2), r(2)
    real(8) :: diam, lete0, lamr, lamp
    logical :: broken
    integer :: n0, n1
    real(8) :: L_2
    ! ----------
    fid = get_file_unit()
    OPEN(UNIT=fid, FILE=TRIM(nomarch), STATUS="OLD")
    ! ----------
    ! Parametros
    iStatus = FindStringInFile("*parametros", fid, .TRUE.)
    READ(fid,*) masim%sidelen
    READ(fid,*) masim%diamed
    read(fid,*) masim%ncapas
    if (.false.) then ! Codigo viejo que estoy reemplazando
        READ(fid,*) masim%nparam
        allocate( masim%param(masim%nparam) )
        read(fid,*) masim%param
    end if
    masim%nparam = numparamcon
    allocate( masim%param(masim%nparam) )
    masim%param = paramcon
    ! ----------
    ! Deformacion (puede no estar)
    iStatus = FindStringInFile("*deformacion", fid, .false.)
    if (iStatus == 0) then
        masim%status_deformed = .true.
        read(fid,*) masim%Fmacro
        read(fid,*) masim%Tmacro
    end if
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
    ! Cambio coordenadas para tener el cero en el centro del RVE
    ! OJO por ahora si hay distorsion (F12 o F21) esto no anda tan facil
    L_2 = 0.5d0 * masim%sidelen
    masim%rnods0 = masim%rnods0 - L_2
    if (masim%status_deformed) then
        L_2 = 0.5d0 * masim%sidelen * masim%Fmacro(1,1)
        masim%rnods(1,:) = masim%rnods(1,:) - L_2
        L_2 = 0.5d0 * masim%sidelen * masim%Fmacro(2,2)
        masim%rnods(2,:) = masim%rnods(2,:) - L_2
    else
        masim%rnods = masim%rnods0
    end if
    ! ----------
    ! Fibras
    iStatus = FindStringInFile("*fibras", fid, .TRUE.)
    READ(fid,*) masim%nfibs
    allocate( masim%fibs(2,masim%nfibs) )
    allocate( masim%diams(masim%nfibs) )
    allocate( masim%letes0(masim%nfibs) )
    allocate( masim%lamsr(masim%nfibs) )
    allocate( masim%lamps(masim%nfibs) )
    allocate( masim%brokens(masim%nfibs) )
    DO i=1,masim%nfibs
        READ(fid,*) j, diam, lete0, lamr, lamp, broken, n0, n1
        masim%diams(i) = diam
        masim%letes0(i) = lete0
        masim%lamsr(i) = lamr
        masim%lamps(i) = lamp
        masim%brokens(i) = broken
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
    real(8) :: L_2, r0(2,masim%nnods), r1(2,masim%nnods)
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
        ! Tengo que escribir las coordenadas con el cero en el vertice inferior izquierdo
        L_2 = 0.5d0 * masim%sidelen
        r0 = masim%rnods0 + L_2
        if (masim%status_deformed) then
            L_2 = 0.5d0 * masim%sidelen * masim%Fmacro(1,1)
            r1(1,:) = masim%rnods(1,:) + L_2
            L_2 = 0.5d0 * masim%sidelen * masim%Fmacro(2,2)
            r1(2,:) = masim%rnods(2,:) + L_2
        else
            r1 = r0
        end if
        WRITE(fid, '(I12,I2,4E20.8E4)') i-1, masim%tipos(i), r0(:,i), r1(:,i)
    END DO
    ! ----------
    ! Fibras
    WRITE(fid,'(A7)') "*Fibras"
    WRITE(fid, '(I12)') masim%nfibs
    DO i=1,masim%nfibs
        WRITE(fid,'(I12, 4E20.8E4, L20 ,2I12)') i-1, masim%diams(i), masim%letes0(i), masim%lamsr(i), masim%lamps(i), masim%brokens(i), masim%fibs(:,i) - 1
    END DO
    ! ----------
    CLOSE(fid)
    ! ----------


END SUBROUTINE escribir_mallita
! ================================================================================
! ================================================================================



! ================================================================================
! ================================================================================
subroutine deformar_afin(masim, FF, r1, axis)
    ! Deforma la mallita de manera Afin: r = F r0
    ! input: FF (tensor gradiente de deformaciones 2x2)
    !
    implicit none
    ! ----------
    type(MallaSim), intent(inout) :: masim
    real(8), intent(in) :: FF(2,2) ! tensor de deformaciones
    real(8), intent(inout) :: r1(2,masim%nnods)
    integer, optional, intent(in) :: axis ! eje que se deforma afin: 0= los dos, 1= x, 2=y
    ! ----------
    integer :: i,j,k
    integer :: axisL
    ! ----------

    ! Chequeo si se da eje
    axisL = 0
    if (present(axis)) axisL=axis


    ! ----------
    if (axisL==0) then
        ! Deformo de manera afin ambos ejes
        do k=1,masim%nnods
            r1(:,k) = 0.d0 ! lo hago cero antes de empezar la sumatoria
            do i=1,2
                do j=1,2
                    r1(i,k) = r1(i,k) + FF(i,j)*masim%rnods0(j,k)
                end do
            end do
        end do
    else if (axisL==1) then
        ! Solamente deformo de manera afin en x (y queda igual)
        do k=1,masim%nnods
            i=1
            r1(i,k) = 0.d0
            do j=1,2
                r1(i,k) = r1(i,k) + FF(i,j)*masim%rnods0(j,k)
            end do
        end do
    else if (axisL==2) then
        ! Solamente deformo de manera afin en y (x queda igual)
        do k=1,masim%nnods
            i=2
            r1(i,k) = 0.d0
            do j=1,2
                r1(i,k) = r1(i,k) + FF(i,j)*masim%rnods0(j,k)
            end do
        end do
    else
        print *, "Mal valor de axisL"
        stop
    end if
    ! ----------
!    masim%rnods = r1
    ! ----------

    ! ----------
end subroutine deformar_afin
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
subroutine deformar_afin_frontera(masim, FF, r1)
    ! Toma un array de coordenadas r1
    ! Y a los nodos de la frontera le aplica la deformacion afin
    ! Al resto los deja como estaban
    ! input: FF (tensor gradiente de deformaciones 2x2)
    !
    implicit none
    ! ----------
    type(MallaSim), intent(inout) :: masim
    real(8), intent(in) :: FF(2,2) ! tensor de deformaciones
    real(8), intent(inout) :: r1(2,masim%nnods)
    ! ----------
    integer :: i,j,k
    ! ----------

    ! ----------
    do k=1,masim%nnods
        if (masim%tipos(k)==1) then
            r1(:,k) = 0.d0
            do i=1,2
                do j=1,2
                    r1(i,k) = r1(i,k) + FF(i,j)*masim%rnods0(j,k)
                end do
            end do
        end if
    end do
    ! ----------
!    masim%rnods = r1
    ! ----------

    ! ----------
end subroutine deformar_afin_frontera
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
    real(8), intent(out) :: tension ! tension ingenieril de la fibra (escalar, >0 elongacion, <0 compresion)
    real(8), intent(out) :: fuerzav(2) ! fuerza de la fibra, vector desde nodo inicial hacia nodo final
    ! ----------
    integer :: selector
    integer :: last
    ! ---
    real(8) :: k1 ! Constante elastica (en fuerza) de fibra recta a la traccion
    real(8) :: k2 ! Constante elastica (en fuerza) de fibra enrulada
    real(8) :: Et ! Modulo elastico de fibra recta a la traccion
    real(8) :: Eb_Et ! Eb/Et
    real(8) :: Eb ! Modulo elastico de fibra enrulada
    real(8) :: Ep_Et ! Ep/Et
    real(8) :: Ep ! Modulo elastico para simular plasticidad
    real(8) :: lamp0_lamr ! lamp0/lamr
    real(8) :: lamp0 ! lamp0 = valor de lam al cual empezaria la plasticidad
    real(8) :: diam ! diametro inicial de la fibra
    real(8) :: lete0 ! longitud extremo-extremo inicial de la fibra
    real(8) :: lamr ! valor de reclutamiento de la fibra
    real(8) :: lamp
    real(8) :: lamrp ! valor de reclutamiento aumentado por la plasticidad
    logical :: broken
    real(8) :: area ! area inicial de la fibra (seccion transversal)
    real(8) :: lete ! longitud extremo-extremo actual
    real(8) :: lam ! elongacion extremo-extremo = lete/lete0
    real(8) :: fuerza ! fuerza actual (modulo) de la fibra
    real(8) :: fuerzar ! fuerza de la fibra para la cual se recluta
    real(8) :: tensionr ! tension de la fibra para la cual se recluta
    real(8) :: tensiony ! tension de la fibra para la cual plastifica
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
        Eb_Et = masim%param(last+2)
        Eb = Eb_Et*Et
        lamr = masim%lamsr(f)
        diam = masim%diams(f)
        area = pi * diam * diam / 4.d0
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
    case (3) ! doy Et, Eb y Ep, es una ley trilineal para simular a groso modo una respuesta con deformacion plastica
        Et = masim%param(last+1)
        Eb_Et = masim%param(last+2)
        Eb = Eb_Et*Et
        Ep_Et = masim%param(last+3)
        Ep = Ep_Et*Et
        lamp0_lamr = masim%param(last+4)
        lamr = masim%lamsr(f)
        lamp0 = lamp0_lamr*lamr
        diam = masim%diams(f)
        area = pi * diam * diam / 4.d0
        lete0 = masim%letes0(f)
        lete = dsqrt(sum(dr_f*dr_f))
        lam = lete / lete0
        if ( lam <= lamr ) then
            tension = Eb*(lam - 1.0d0) * diam
        else if (lam <= lamp0 ) then
            tensionr = Eb * (lamr - 1.0d0)
            tension = tensionr + Et * (lam/lamr - 1.0d0)
        else ! simulando plasticidad en la respuesta
            tensionr = Eb * (lamr - 1.0d0)
            tensiony = tensionr + Et * (lamp0/lamr - 1.0d0)
            tension = tensionr + tensiony + Ep*(lam/lamp0 - 1.0d0)
        end if
        fuerza = tension * area
        fuerzav = fuerza * dr_f / lete
    case (4) ! Ecuacion con plasticidad
        Et = masim%param(last+1)
        Eb_Et = masim%param(last+2)
        Eb = Eb_Et*Et
        lamr = masim%lamsr(f)
        lamp = masim%lamps(f)
        lamrp = lamr*lamp
        broken = masim%brokens(f)
        diam = masim%diams(f)
        area = pi * diam * diam / 4.d0
        lete0 = masim%letes0(f)
        lete = dsqrt(sum(dr_f*dr_f))
        lam = lete / lete0
        ! Ahora si calculo la tension
        if (broken) then
            tension = 0.d0
        elseif ( lam <= lamrp ) then
            tension = Eb*(lam/lamp - 1.0d0) * diam
        else
            tensionr = Eb * (lamrp - 1.0d0)
            tension = tensionr + Et * (lam/lamrp - 1.0d0)
        end if
        fuerza = tension * area
        fuerzav = fuerza * dr_f / lete
    case default
        write(*,*) "Ley constitutiva desconocida, selector: ", selector
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
    real(8), intent(in) :: r1(2,masim%nnods)
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
subroutine desplazar_nodos_malla(masim, nveces, drmag, fza_ref, fza_tol, r1, balanced)
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
    real(8), intent(in) :: fza_ref
    real(8), intent(in) :: fza_tol
    real(8), intent(inout) :: r1(2,masim%nnods)
    logical, intent(out) :: balanced(masim%nnods)
    ! ----------
    real(8) :: fzas_fibs(2,masim%nfibs)
    real(8) :: fzas_nods(2,masim%nnods)
    real(8) :: fzas_nods_mags(masim%nnods)
    integer :: vez, n
    real(8) :: fza_n(2), fza_n_mag, dr_n(2)
    ! ----------

    balanced = .false.
    do vez=1,nveces
        call calcular_fuerzas(masim, r1, fzas_fibs, fzas_nods)
        do n=1,masim%nnods
            if (masim%tipos(n) == 1) then
                ! nodo de frontera = nodo de Dirichlet
                r1(:,n) = r1(:,n) ! sentencia innecesaria pero por ahora la dejo para que quede claro
                fzas_nods_mags(n) = 0.d0
                balanced(n) = .true.
            else
                fza_n = fzas_nods(:,n)
                fza_n_mag = dsqrt(sum(fza_n*fza_n))
                fzas_nods_mags(n) = fza_n_mag
                if (fza_n_mag < fza_tol) then
                    dr_n = 0.d0
                    balanced(n) = .true.
                else
                    dr_n = drmag * fza_n / fza_ref
                    balanced(n) = .false.
                end if
                r1(:,n) = r1(:,n) + dr_n(:)
            end if
        end do
        if ( all(balanced) ) then
            write(*,*) "balanced!", vez
            exit
        end if
    end do
    if ( .not. all(balanced) ) then
        write(*,*) "maxres: ", maxval( fzas_nods_mags )
    end if

    ! ----------
end subroutine desplazar_nodos_malla
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
subroutine calcular_equilibrio(masim, npasos, vec_veces, vec_drmags, f_ref, f_tol, Fmacro)
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
    real(8), intent(in) :: f_ref
    real(8), intent(in) :: f_tol
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
    logical :: equilibrio(masim%nnods)
    ! ----------

    ! deformo afin
    if (npasos==0) then
        ! Deformo toda la malla afin antes de comenzar la vibracion
        call deformar_afin(masim, Fmacro, r1)
    else
        ! Deformo solamente la frontera de forma afin
        r1 = masim%rnods
        call deformar_afin_frontera(masim, Fmacro, r1)
!        call deformar_afin(masim, Fmacro, r1, axis=2)
    end if
    ! Comienzo a mover los nodos segun las resultantes nodales
    do paso=1,npasos
        nveces = vec_veces(paso)
        drmag = vec_drmags(paso)
        call desplazar_nodos_malla(masim, nveces, drmag, f_ref, f_tol, r1, equilibrio)
        ! chequeo equilibrio
        if ( all(equilibrio) ) then
            exit
        end if
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
end subroutine calcular_equilibrio
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
subroutine calcular_plasticidad_rotura(masim, dt)
    ! Ante una deformacion dada (masim%rnods) calcula la deformacion plastica
    ! que sufren las fibras. Son dos valores:
    ! dotlamp = tasa de deformacion plastica
    ! broken = indicador de si la fibra se rompe o no
    ! ----------
    implicit none
    ! ----------
    type(MallaSim), intent(inout) :: masim
    real(8), intent(in) :: dt
    ! ----------
    logical :: brokens(masim%nfibs)
    real(8) :: dotlamps(masim%nfibs)
    ! parametros constitutivos de plasticidad y rotura
    real(8) :: doteps0 ! parametro dimensional proporcional de plasticidad
    real(8) :: s0 ! resistencia a la fluencia inicial
    real(8) :: nhard ! parametro de endurecimiento por deformacion plastica
    real(8) :: elonlimit ! limite de rotura por elongacion
    ! resto
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
    real(8) :: lamp_f
    real(8) :: lamr_f
    real(8) :: lamef_f ! es un valor de elongacion efectiva, o elongacion elastica (la que va con la tension) = lam/lamr/lamp
    ! ----------
    real(8) :: s_f
    real(8) :: dotlamp_f
    logical :: broken_f
    ! ----------

    ! Parametros constitutivos importantes para plasticidad y rotura
    doteps0 = masim%param(4)
    s0 = masim%param(5)
    nhard = masim%param(6)
    elonlimit = masim%param(7)
    ! Bucle por las fibras para ir calculando tensiones actuales
    do f=1,masim%nfibs
        ! Guardo informacion en variables mas sencillas
        d_f = masim%diams(f) ! Diametro inicial de la fibra f
        Area_f = pi * d_f * d_f / 4.d0 ! Seccion transversal inicial de la fibra f
        lete0_f = masim%letes0(f) ! Longitud extremo-extremo inicial de la fibra f
        Vol_f = Area_f * lete0_f ! Volumen de la fibra f
        n1_f = masim%fibs(1,f) ! Nodo inicial de la fibra f
        n2_f = masim%fibs(2,f) ! Nodo final de la fibra f
        r_n1_f = masim%rnods(:,n1_f) ! Posicion actual de n1_f
        r_n2_f = masim%rnods(:,n2_f) ! Posicion actual de n2_f
        dr_f = r_n2_f - r_n1_f ! Vector fibra: apunta de nodo inicial a nodo final y su magnitud es la longitud actual
        lete_f = dsqrt(sum(dr_f*dr_f)) ! Longitud extremo-extremo actual
        lam_f = lete_f / lete0_f ! Elongacion extremo-extremo
        lamr_f = masim%lamsr(f)
        a_f = dr_f/lete_f * lam_f ! Vector orientacion de la fibra, su magnitud es la elongacion lam_f
        ! Calculo las tensiones
        call calcular_tension_fibra(masim, f, dr_f, ten_f, vfza_f)
        ! Calculo la plasticidad
        broken_f =  masim%brokens(f)
        lamp_f = masim%lamps(f)
        lamef_f = lam_f / lamp_f / lamr_f
        if (broken_f) then
            ! esta fibra esta rota y no deforma mas, su plasticidad queda constante
            dotlamp_f = 0.d0
        else
            ! esta fibra puede romperse o plastificar
            if (ten_f>elonlimit) then
                ! se rompe
                broken_f = .true.
                dotlamp_f = 0.d0
            else
                ! plastifica poco o mucho
                s_f = s0 * lamp_f**nhard
                dotlamp_f = doteps0 * dsinh(ten_f/s_f)
            end if
        end if
        brokens(f) = broken_f
        dotlamps(f) = dotlamp_f
        lamp_f = lamp_f + dotlamp_f*dt
        if ( (lam_f>1.d0) .and. (lamp_f > lam_f) ) then
            !print *, "lamp_f>lam_f en fibra: ", f, lam_f, lamp_f
            if (lamp_f > elonlimit) then
                print *, "Tac!", lam_f, lamp_f
                !call escribir_mallita(masim, "malla_con_error.txt"//repeat(" ", 120))
                brokens(f) = .true.
                dotlamps(f) = 0.d0 !
                lamp_f = masim%lamps(f) ! no plastifico, porque rompi
            end if
        end if
        ! Actualizo informacion en la propia malla
        masim%brokens(f) = broken_f
        masim%lamps(f) = lamp_f
    end do
    ! Actualizo la propia malla

    ! ----------
end subroutine calcular_plasticidad_rotura
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


! ==============================================================================
END MODULE class_mallita
! ==============================================================================
