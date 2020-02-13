! ==============================================================================
MODULE class_malla_completa
! ==============================================================================
    USE Aux
    IMPLICIT NONE

    !PRIVATE
    !PUBLIC :: MallaCom
    !PUBLIC :: leer_malla

    TYPE MallaCom
        ! diemensiones generales
        REAL(8) :: sidelen
        REAL(8) :: diamed
        real(8) :: volfrac
        real(8) :: lenseg
        real(8) :: themax
        ! nodos
        INTEGER :: nnods
        INTEGER, ALLOCATABLE :: tipos(:)
        REAL(8), ALLOCATABLE :: rnods(:,:)
        ! segmentos
        INTEGER :: nsegs
        INTEGER, ALLOCATABLE :: segs(:,:)
        ! fibras
        INTEGER :: nfibs
        INTEGER :: lenfibsje
        INTEGER, ALLOCATABLE :: fibsne(:), fibsie(:), fibsje(:)
        REAL(8), ALLOCATABLE :: fibsds(:), fibsdls(:), fibsdths(:)
        ! capas
        INTEGER :: ncaps
        INTEGER :: lencapsje
        INTEGER, ALLOCATABLE :: capsne(:), capsie(:), capsje(:)
    CONTAINS
        !
    END TYPE MallaCom

    TYPE Interseccion
        INTEGER :: tipo
        INTEGER :: fibs(2)
        INTEGER :: segs(2)
        INTEGER :: segs_f(2)
        INTEGER :: seg1(2)
        INTEGER :: seg2(2)
        INTEGER :: inewseg1
        INTEGER :: inewseg2
        INTEGER :: inewnod
        REAL(8) :: newnodr(2)
    CONTAINS
        !
    END TYPE Interseccion

! ==============================================================================
CONTAINS
! ==============================================================================


! ================================================================================
! ================================================================================
SUBROUTINE leer_malla(malla, nomarch)
    IMPLICIT NONE
    CHARACTER(LEN=120), INTENT(IN) :: nomarch
    TYPE(MallaCom), INTENT(INOUT) :: malla
    !
    INTEGER :: fid, iStatus
    INTEGER :: i, j, k
    INTEGER :: tipo
    REAL(8) :: r(2)
    INTEGER :: seg(2)
    REAL(8) :: dl, d, dth
    INTEGER :: nfibsegs, lenfibsje
    INTEGER :: ncapfibs, lencapsje
    INTEGER :: ie0, ie1

    ! ----------
    fid = get_file_unit()
    OPEN(UNIT=fid, FILE=TRIM(nomarch), STATUS="OLD")
    ! ----------
    iStatus = FindStringInFile("*parametros", fid, .TRUE.)
    READ(fid,*) malla%sidelen ! viene em micrones
    READ(fid,*) malla%diamed ! viene em micrones
    READ(fid,*) malla%volfrac
    READ(fid,*) malla%lenseg ! viene en micrones
    READ(fid,*) malla%themax ! lo leo en grados
    malla%themax = malla%themax * PI / 180.d0 ! lo paso a radianes
    ! ----------
    iStatus = FindStringInFile("*coordenadas", fid, .TRUE.)
    READ(fid,*) malla%nnods
    ALLOCATE( malla%rnods(2,malla%nnods) )
    ALLOCATE( malla%tipos(malla%nnods) )
    DO i=1,malla%nnods
        READ(fid,*) j, tipo, r ! los leo en micrones
        malla%rnods(:,i) = r
        malla%tipos(i) = tipo
    END DO
    ! ----------
    iStatus = FindStringInFile("*segmentos", fid, .TRUE.)
    READ(fid,*) malla%nsegs
    ALLOCATE( malla%segs(2,malla%nsegs) )
    DO i=1,malla%nsegs
        READ(fid,*) j, seg
        malla%segs(:,i) = seg + 1 ! sumo uno porque vengo de python con base 0
    END DO
    ! ----------
    iStatus = FindStringInFile("*fibras", fid, .TRUE.)
    READ(fid,*) malla%nfibs
    ALLOCATE( malla%fibsdls(malla%nfibs) )
    ALLOCATE( malla%fibsds(malla%nfibs) )
    ALLOCATE( malla%fibsdths(malla%nfibs) )
    ALLOCATE( malla%fibsne(malla%nfibs) )
    ALLOCATE( malla%fibsie(malla%nfibs + 1) )
    ! Recorro dos veces las fibras
    ! la primera guardo el ne y calculo lenfibsje (cant total de segs)
    ! la segunda calculo el ie y guardo el je
    ! notar que el lenje de las fibras deberia ser el numero total de segmentos
    ! porque cada segmento solo pertenece a una fibra, y el total de fibras contienen a todos los segmentos
    lenfibsje = 0
    DO i=1,malla%nfibs
        READ(fid,*) j, dl, d, dth, nfibsegs
        malla%fibsne(i) = nfibsegs
        lenfibsje = lenfibsje + nfibsegs
    END DO
    malla%lenfibsje = lenfibsje
    iStatus = FindStringInFile("*fibras", fid, .TRUE.)
    READ(fid,*) malla%nfibs ! redundante pero necesario leer algo
    ALLOCATE( malla%fibsje(lenfibsje) )
    malla%fibsie(1) = 1
    DO i=1,malla%nfibs
        malla%fibsie(i+1) = malla%fibsie(i) + malla%fibsne(i)
        ie0 = malla%fibsie(i)
        ie1 = malla%fibsie(i+1)
        READ(fid,*) j, dl, d, dth, nfibsegs, malla%fibsje(ie0:ie1-1)
        malla%fibsje(ie0:ie1-1) = malla%fibsje(ie0:ie1-1) + 1 ! necesario porque vengo de python con base 0
        malla%fibsdls(i) = dl
        malla%fibsds(i) = d
        malla%fibsdths(i) = dth
    END DO
    ! ----------
    iStatus = FindStringInFile("*capas", fid, .TRUE.)
    READ(fid,*) malla%ncaps
    ALLOCATE( malla%capsne(malla%ncaps) )
    ALLOCATE( malla%capsie(malla%ncaps + 1) )
    ! Recorro dos veces las capas
    ! la primera guardo el ne y calculo lenje (cant total de fibras)
    ! la segunda calculo el ie y guardo el je
    ! notar que el lenje de las capas deberia ser el numero total de fibras
    ! porque cada fibra solo pertenece a una capa, y el total de capas contienen a todas las fibras
    lencapsje = 0
    DO i=1,malla%ncaps
        READ(fid,*) j, ncapfibs
        malla%capsne(i) = ncapfibs
        lencapsje = lencapsje + ncapfibs
    END DO
    malla%lencapsje = lencapsje
    iStatus = FindStringInFile("*capas", fid, .TRUE.)
    READ(fid,*) malla%ncaps ! redundante pero necesario leer algo
    ALLOCATE( malla%capsje(lencapsje) )
    malla%capsie(1) = 1
    DO i=1,malla%ncaps
        malla%capsie(i+1) = malla%capsie(i) + malla%capsne(i)
        ie0 = malla%capsie(i)
        ie1 = malla%capsie(i+1)
        READ(fid,*) j, ncapfibs, malla%capsje(ie0:ie1-1)
        malla%capsje(ie0:ie1-1) = malla%capsje(ie0:ie1-1) + 1 ! vengo de python con base 0
    END DO
    ! ----------
    CLOSE(fid)
    ! ----------


END SUBROUTINE leer_malla
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
SUBROUTINE escribir_malla(malla, nomarch)
    IMPLICIT NONE
    CHARACTER(LEN=120), INTENT(IN) :: nomarch
    TYPE(MallaCom), INTENT(INOUT) :: malla
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
    WRITE(fid,'(E20.8E4)') malla%sidelen
    WRITE(fid,'(E20.8E4)') malla%diamed
    WRITE(fid,'(E20.8E4)') malla%volfrac
    WRITE(fid,'(E20.8E4)') malla%lenseg
    WRITE(fid,'(E20.8E4)') malla%themax * 180.d0 / PI ! lo guardo en grados
    ! ----------
    ! Nodos
    WRITE(fid,'(A12)') "*Coordenadas"
    WRITE(fid,'(I12)') malla%nnods
    DO i=1,malla%nnods
        WRITE(fid, '(I12,I2,2E20.8E4)') i-1, malla%tipos(i), malla%rnods(:,i)
    END DO
    ! ----------
    ! Segmentos
    WRITE(fid,'(A10)') "*Segmentos"
    WRITE(fid, '(I12)') malla%nsegs
    DO i=1,malla%nsegs
        WRITE(fid,'(I12, 2I12)') i-1, malla%segs(:,i) - 1
    END DO
    ! ----------
    ! Fibras
    WRITE(fid,'(A7)') "*Fibras"
    WRITE(fid,'(I12)') malla%nfibs
    DO i=1,malla%nfibs
        WRITE(formato, '(A18,I0,A4)') "(I12,3E20.8E4,I12,",malla%fibsne(i),"I12)"
        WRITE(fid,formato) i-1, malla%fibsdls(i), malla%fibsds(i), malla%fibsdths(i), malla%fibsne(i), malla%fibsje(malla%fibsie(i):malla%fibsie(i+1)-1) - 1
    END DO
    ! ----------
    ! Capas
    WRITE(fid,'(A6)') "*Capas"
    WRITE(fid,'(I12)') malla%ncaps
    DO i=1,malla%ncaps
        WRITE(formato, '(A6,I0,A4)') "(2I12,",malla%capsne(i),"I12)"
        WRITE(fid,formato) i-1, malla%capsne(i), malla%capsje(malla%capsie(i):malla%capsie(i+1)-1) - 1
    END DO
    ! ----------
    CLOSE(fid)
    ! ----------


END SUBROUTINE escribir_malla
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
SUBROUTINE intersectar_fibras_3(m_in, m_out, label_intercapas, label_periodicidad, iStatus)
    ! voy a recorrer m_in, sus capas y fibras, encontrando intersecciones entre fibras
    ! como resultado la malla de salida es una malla aumentada porque tiene mas nodos y mas segmentos
    ! en primera instancia la voy a sobredimensionar en memoria
    ! luego readjudico la memoria correcta
    IMPLICIT NONE
    TYPE(MallaCom), INTENT(IN) :: m_in
    TYPE(MallaCom), INTENT(OUT) :: m_out
    LOGICAL, INTENT(IN) :: label_intercapas
    logical, intent(in) :: label_periodicidad
    INTEGER, INTENT(OUT) :: iStatus
    !
    INTEGER :: c1, c2 ! numeracion de capas
    INTEGER :: f1_c1, f1 ! numeracion de fibra 1: local dentro de la capa c1 y global
    INTEGER :: s1_f1, s1 ! numeracion de segmento 1: local dentro de la fibra f1 y global
    INTEGER :: f2_c1, f2 ! numeracion de fibra 2: local y global
    INTEGER :: s2_f2, s2 ! numeracion de segmento 2: local y global
    INTEGER :: n11, n12, n21, n22
    REAL(8) :: r11(2), r12(2), r21(2), r22(2)
    INTEGER :: iStat
    REAL(8) :: rin(2)
    !
    INTEGER :: i1fibs1, nfibs1
    INTEGER, ALLOCATABLE :: fibs1(:)
    !
    ! preadjudico las variables con el doble de tamanio para que sobre
    INTEGER :: nins ! numero de intersecciones
    TYPE(Interseccion) :: ins(100000) ! preadjudicado grande para que sobre
    INTEGER :: nfibins(m_in%nfibs) ! aqui guardo la cantidad de intersecciones sobre cada fibra
    INTEGER :: ifibins(10000, m_in%nfibs) ! aqui guardo en que indice local de cada fibra colocar los segmentos nuevos (i_ins, j_fib)
    INTEGER :: jfibins(10000, m_in%nfibs) ! aqui guardo el indice global de cada segmento nuevo de cada fibra (i_ins, j_fib)
    INTEGER :: nsegins(m_in%nsegs) ! cantidad de intersecciones por cada segmento
    INTEGER :: isegins(100, m_in%nsegs) ! nodo viejo para cada interseccion para cada segmento (i_ins, j_seg)
    INTEGER :: jsegins(100, m_in%nsegs) ! nodo nuevo para cada interseccion para cada segmento (i_ins, j_seg)
    INTEGER :: newisegins(m_in%nsegs) ! indice nuevo del ultimo segmento creado debido a una interseccion del segmento i_seg
    !
    INTEGER :: i, j
    INTEGER :: m_in_fie0, m_in_fie1, m_out_fie0, m_out_fie1
    INTEGER :: nins2
    INTEGER :: aux1, aux2, aux3, aux4
    INTEGER :: inews1, inews2, inewn
    !


    ! ----------
    ! inicializo algunos arrays
    nins = 0
    nfibins = 0
    nsegins = 0
    DO i=1,m_in%nsegs
        newisegins(i) = i
    END DO
    ! ----------
    ! recorro las capas e intersecto fibras dentro de cada capa y con capas adyacentes
    DO c1 = 1,m_in%ncaps
        ! recorro las fibras de la capa c1
        DO f1_c1 = 1,m_in%capsne(c1)
            f1 = m_in%capsje( m_in%capsie(c1)-1 + f1_c1 )
            ! recorro los segmentos de la fibra f1
            DO s1_f1 = 1,m_in%fibsne(f1)
                s1 = m_in%fibsje( m_in%fibsie(f1)-1 + s1_f1 )
!                if (s1==1) then
!                    write(*,*) s1
!                end if
                ! si el segmento s1 ya sufrio una interseccion paso al siguiente
                n11 = m_in%segs(1,s1)
                n12 = m_in%segs(2,s1)
                r11 = m_in%rnods(:,n11)
                r12 = m_in%rnods(:,n12)
                ! puedo estar buscando intersecciones dentro de la misma capa
                ! o entre capas
                IF (ALLOCATED(fibs1)) DEALLOCATE( fibs1 )
                IF (label_intercapas) THEN
                    ! recorro las fibras de la capa siguiente
                    ! o sea que si estoy en la ultima capa, no lo hago
                    IF (c1==m_in%ncaps) THEN
                        if (.not. label_periodicidad) CYCLE ! si activo esto la ultima capa no intersecta con la primera
                        i1fibs1 = 1
                        nfibs1 = m_in%capsne(1) ! voy a chequear la ultima capa con la primera
                        allocate( fibs1(nfibs1) )
                        fibs1 = m_in%capsje( m_in%capsie(1) : m_in%capsie(1)-1 + nfibs1 )
                    else
                        i1fibs1 = 1
                        nfibs1 = m_in%capsne(c1+1)
                        ALLOCATE( fibs1(nfibs1) )
                        fibs1 = m_in%capsje( m_in%capsie(c1+1) : m_in%capsie(c1+1)-1 + nfibs1 )
                    END IF
                ELSE
                    ! recorro las fibras restantes de la capa c1 para intersectar
                    i1fibs1 = f1_c1 + 1
                    nfibs1 = m_in%capsne(c1)
                    ALLOCATE( fibs1(nfibs1) )
                    fibs1 = m_in%capsje( m_in%capsie(c1) : m_in%capsie(c1)-1 + nfibs1 )
                END IF
                DO f2_c1 = i1fibs1,nfibs1
                    f2 = fibs1(f2_c1)
                    ! recorro los segmentos de la fibra f2
                    DO s2_f2 = 1,m_in%fibsne(f2)
                        s2 = m_in%fibsje( m_in%fibsie(f2)-1 + s2_f2 )
                        ! si el segmento s2 ya sufrio una interseccion paso al siguiente
                        IF ( (nsegins(s2)>0) .OR. (nsegins(s1)>0) ) THEN
                            CYCLE
                        END IF
                        n21 = m_in%segs(1,s2)
                        n22 = m_in%segs(2,s2)
                        r21 = m_in%rnods(:,n21)
                        r22 = m_in%rnods(:,n22)
!                        WRITE(*,*) s1, s2, n11, n12, n21, n22
!                        IF ( (s1==81).AND.(s2==393) ) THEN
!                            aux1 = 1
!                        END IF
!                        if (s2==472) then
!                            write(*,*) s2
!                        end if
                        CALL intersectar_segmentos(r11, r12, r21, r22, rin, iStat)
                        IF (iStat == 2) THEN ! 2 es interseccion tipo medio-medio
!                            IF ((n11==386).OR.(n21==386)) THEN
!                                iStat = 2 ! solo para debuggear
!                            END IF
                            nins = nins + 1
                            inewn = m_in%nnods + nins
                            inews1 = m_in%nsegs + 2*nins - 1
                            inews2 = inews1 + 1
                            nsegins(s1) = nsegins(s1) + 1
                            nsegins(s2) = nsegins(s2) + 1
                            isegins(nsegins(s1),s1) = inews1
                            jsegins(nsegins(s1),s1) = inewn
                            isegins(nsegins(s2),s2) = inews2
                            jsegins(nsegins(s2),s2) = inewn
                            nfibins(f1) = nfibins(f1) + 1
                            nfibins(f2) = nfibins(f2) + 1
                            ifibins(nfibins(f1),f1) = s1_f1 + 1
                            ifibins(nfibins(f2),f2) = s2_f2 + 1
                            jfibins(nfibins(f1),f1) = inews1
                            jfibins(nfibins(f2),f2) = inews2
                            ins(nins)%tipo = iStat ! ojo que no es el tipo de nodo, sino el tipo de interseccion (aunque coincida que es un 2 de tuje)
                            ins(nins)%fibs = [f1,f2] ! fibras que participan en la interseccion
                            ins(nins)%segs = [s1,s2] ! segmentos que participan en la interseccion
                            ins(nins)%segs_f = [s1_f1, s2_f2] ! indices locales de los segmentos dentro de sus fibras
                            ins(nins)%seg1 = [n11, n12]
                            ins(nins)%seg2 = [n21, n22]
                            ins(nins)%inewseg1 = inews1
                            ins(nins)%inewseg2 = inews2
                            ins(nins)%inewnod = inewn
                            ins(nins)%newnodr = rin
                            newisegins(s1) = inews1 ! en caso de que vuelva a tener alguna interseccion
                            newisegins(s2) = inews2 ! en caso de que vuelva a tener alguna interseccion
                        END IF
                    END DO
                END DO
            END DO
        END DO
    END DO
    ! ----------

    IF (nins==0) THEN
        iStatus = 1
        m_out = m_in
        RETURN
    END IF

    ! ----------
    ! aqui ya tengo todas las intersecciones guardadas, ahora tengo que crear la nueva malla
    ! dimensiones generales
    m_out%sidelen = m_in%sidelen
    m_out%diamed = m_in%diamed
    m_out%volfrac = m_in%volfrac
    m_out%lenseg = m_in%lenseg
    m_out%themax = m_in%themax
    ! ----------
    ! nodos
    m_out%nnods = m_in%nnods + nins ! cada interseccion implica un nuevo nodo
    ALLOCATE( m_out%tipos(m_out%nnods) ) ! adjudico
    ALLOCATE( m_out%rnods(2,m_out%nnods) ) ! adjudico
    m_out%tipos(1:m_in%nnods) = m_in%tipos ! copio los tipos viejos en los tipos nuevos
    m_out%tipos(m_in%nnods+1:m_out%nnods) = 2 ! los tipos nuevos extras son todos tipo interseccion
    m_out%rnods(:,1:m_in%nnods) = m_in%rnods ! copio coordenadas de los nodos viejos en los nodos nuevos
    DO i=1,nins ! copio coordenadas de los nodos interseccion en los nodos nuevos
        m_out%rnods(:,m_in%nnods+i) = ins(i)%newnodr(:)
    END DO
    ! ----------
    ! segmentos
    m_out%nsegs = m_in%nsegs + 2*nins ! cada interseccion implica dos nuevos segmentos
    ALLOCATE( m_out%segs(2,m_out%nsegs) ) ! adjudico
    m_out%segs(:,1:m_in%nsegs) = m_in%segs ! copio los segmentos viejos en los segmentos nuevos
    ! (los segmentos quedan incompletos, hay que modificarlos luego en un bucle sobre las intersecciones)
    DO i=1,nins
        s1 = ins(i)%segs(1)
        s2 = ins(i)%segs(2)
        n11 = m_in%segs(1,s1)
        n12 = m_in%segs(2,s1)
        n21 = m_in%segs(1,s2)
        n22 = m_in%segs(2,s2)
        inews1 = ins(i)%inewseg1
        inews2 = ins(i)%inewseg2
        inewn = ins(i)%inewnod
        ! cambio la conectividad de los dos segmentos viejos
        ! y agrego los dos segmentos nuevos
        m_out%segs(2,s1) = inewn
        m_out%segs(2,s2) = inewn
        m_out%segs(:,inews1) = [inewn, n12]
        m_out%segs(:,inews2) = [inewn, n22]
    END DO
    ! ----------
    ! fibras
    m_out%nfibs = m_in%nfibs
    m_out%lenfibsje = m_in%lenfibsje + 2*nins ! tengo dos nuevos elemento en la conectividad de las fibras por cada interseccion
    ALLOCATE( m_out%fibsds(m_out%nfibs) )
    ALLOCATE( m_out%fibsdls(m_out%nfibs) )
    ALLOCATE( m_out%fibsdths(m_out%nfibs) )
    ALLOCATE( m_out%fibsne(m_out%nfibs) )
    ALLOCATE( m_out%fibsie(m_out%nfibs + 1) )
    ALLOCATE( m_out%fibsje(m_out%lenfibsje) )
    m_out%fibsds = m_in%fibsds
    m_out%fibsdls = m_in%fibsdls
    m_out%fibsdths = m_in%fibsdths
    m_out%fibsie(1) = 1
    ! primero copio la conectividad vieja en la conectividad nueva
    ! teniendo cuidado de donde va cada fibra
    m_out%fibsje = -1 ! inicializo en -1 indicando lugares libres (sin segmento asignado)
    DO i=1,m_in%nfibs
        m_out%fibsne(i) = m_in%fibsne(i) + nfibins(i)
        m_out%fibsie(i+1) = m_out%fibsie(i) + m_out%fibsne(i)
    END DO
    !
    m_out%fibsje(1:m_in%lenfibsje) = m_in%fibsje
    ! ahora divido cada segmento en dos recorriendo el je de las fibras
    nins2 = 0
    DO i=1,m_in%lenfibsje
        s1 = m_in%fibsje(i)
        IF (nsegins(s1)>0) THEN
            ! este segmento fue intersectado, hay que partirlo en dos
            nins2 = nins2 + 1
            inews1 = isegins(1,s1) ! nodo nuevo generado al partir al segmento s1 en dos
            ! le hago lugar al nuevo segmento desplazando la segunda parte de la conectividad hacia adelante
            m_out%fibsje(i+nins2+1:m_in%lenfibsje+nins2) = m_out%fibsje(i+nins2:m_in%lenfibsje+nins2-1)
            ! pongo el segmento nuevo
            m_out%fibsje(i+nins2) = inews1
        END IF
    END DO
    ! ----------
    ! capas
    m_out%ncaps = m_in%ncaps
    m_out%lencapsje = m_in%lencapsje
    ALLOCATE( m_out%capsne(m_out%ncaps) )
    ALLOCATE( m_out%capsie(m_out%ncaps + 1) )
    ALLOCATE( m_out%capsje(m_out%lencapsje) )
    m_out%capsne = m_in%capsne
    m_out%capsie = m_in%capsie
    m_out%capsje = m_in%capsje
    ! ----------

    iStatus = 0


END SUBROUTINE intersectar_fibras_3
! ================================================================================
! ================================================================================



! ================================================================================
! ================================================================================
SUBROUTINE intersectar_fibras_2(m_in, m_out)
    ! voy a recorrer m_in, sus capas y fibras, encontrando intersecciones entre fibras
    ! como resultado la malla de salida es una malla aumentada porque tiene mas nodos y mas segmentos
    ! en primera instancia la voy a sobredimensionar en memoria
    ! luego readjudico la memoria correcta
    IMPLICIT NONE
    TYPE(MallaCom), INTENT(IN) :: m_in
    TYPE(MallaCom), INTENT(OUT) :: m_out
    !
    INTEGER :: c1, c2 ! numeracion de capas
    INTEGER :: f1_c1, f1 ! numeracion de fibra 1: local dentro de la capa c1 y global
    INTEGER :: s1_f1, s1 ! numeracion de segmento 1: local dentro de la fibra f1 y global
    INTEGER :: f2_c1, f2 ! numeracion de fibra 2: local y global
    INTEGER :: s2_f2, s2 ! numeracion de segmento 2: local y global
    INTEGER :: n11, n12, n21, n22
    REAL(8) :: r11(2), r12(2), r21(2), r22(2)
    INTEGER :: iStat
    REAL(8) :: rin(2)
    !
    ! preadjudico las variables con el doble de tamanio para que sobre
    INTEGER :: nins ! numero de intersecciones
    TYPE(Interseccion) :: ins(1000) ! preadjudicado grande para que sobre
    INTEGER :: nfibins(m_in%nfibs) ! aqui guardo la cantidad de intersecciones sobre cada fibra
    INTEGER :: ifibins(100, m_in%nfibs) ! aqui guardo en que indice local de cada fibra colocar los segmentos nuevos (i_ins, j_fib)
    INTEGER :: jfibins(100, m_in%nfibs) ! aqui guardo el indice global de cada segmento nuevo de cada fibra (i_ins, j_fib)
    INTEGER :: nsegins(m_in%nsegs) ! cantidad de intersecciones por cada segmento
    INTEGER :: isegins(10, m_in%nsegs) ! nodo viejo para cada interseccion para cada segmento (i_ins, j_seg)
    INTEGER :: jsegins(10, m_in%nsegs) ! nodo nuevo para cada interseccion para cada segmento (i_ins, j_seg)
    INTEGER :: newisegins(m_in%nsegs) ! indice nuevo del ultimo segmento creado debido a una interseccion del segmento i_seg
    !
    INTEGER :: i, j
    INTEGER :: m_in_fie0, m_in_fie1, m_out_fie0, m_out_fie1
    INTEGER :: nins2
    INTEGER :: aux1, aux2, aux3, aux4
    INTEGER :: inews1, inews2, inewn
    !


    ! ----------
    ! inicializo algunos arrays
    nins = 0
    nfibins = 0
    nsegins = 0
    DO i=1,m_in%nsegs
        newisegins(i) = i
    END DO
    ! ----------
    ! recorro las capas e intersecto fibras dentro de cada capa y con capas adyacentes
    DO c1 = 1,m_in%ncaps
        ! recorro las fibras de la capa c1
        DO f1_c1 = 1,m_in%capsne(c1)
            f1 = m_in%capsje( m_in%capsie(c1)-1 + f1_c1 )
            ! recorro los segmentos de la fibra f1
            DO s1_f1 = 1,m_in%fibsne(f1)
                s1 = m_in%fibsje( m_in%fibsie(f1)-1 + s1_f1 )
                ! si el segmento s1 ya sufrio una interseccion paso al siguiente
                n11 = m_in%segs(1,s1)
                n12 = m_in%segs(2,s1)
                r11 = m_in%rnods(:,n11)
                r12 = m_in%rnods(:,n12)
                ! recorro las fibras restantes de la capa c1 para intersectar
                DO f2_c1 = f1_c1+1,m_in%capsne(c1)
                    f2 = m_in%capsje( m_in%capsie(c1)-1 + f2_c1 )
                    ! recorro los segmentos de la fibra f2
                    DO s2_f2 = 1,m_in%fibsne(f2)
                        s2 = m_in%fibsje( m_in%fibsie(f2)-1 + s2_f2 )
                        ! si el segmento s2 ya sufrio una interseccion paso al siguiente
                        IF ( (nsegins(s2)>0) .OR. (nsegins(s1)>0) ) THEN
                            CYCLE
                        END IF
                        n21 = m_in%segs(1,s2)
                        n22 = m_in%segs(2,s2)
                        r21 = m_in%rnods(:,n21)
                        r22 = m_in%rnods(:,n22)
!                        WRITE(*,*) s1, s2, n11, n12, n21, n22
!                        IF ( (s1==81).AND.(s2==393) ) THEN
!                            aux1 = 1
!                        END IF
                        CALL intersectar_segmentos(r11, r12, r21, r22, rin, iStat)
                        IF (iStat == 2) THEN ! 2 es interseccion tipo medio-medio
!                            IF ((n11==386).OR.(n21==386)) THEN
!                                iStat = 2 ! solo para debuggear
!                            END IF
                            nins = nins + 1
                            inewn = m_in%nnods + nins
                            inews1 = m_in%nsegs + 2*nins - 1
                            inews2 = inews1 + 1
                            nsegins(s1) = nsegins(s1) + 1
                            nsegins(s2) = nsegins(s2) + 1
                            isegins(nsegins(s1),s1) = inews1
                            jsegins(nsegins(s1),s1) = inewn
                            isegins(nsegins(s2),s2) = inews2
                            jsegins(nsegins(s2),s2) = inewn
                            nfibins(f1) = nfibins(f1) + 1
                            nfibins(f2) = nfibins(f2) + 1
                            ifibins(nfibins(f1),f1) = s1_f1 + 1
                            ifibins(nfibins(f2),f2) = s2_f2 + 1
                            jfibins(nfibins(f1),f1) = inews1
                            jfibins(nfibins(f2),f2) = inews2
                            ins(nins)%tipo = iStat ! ojo que no es el tipo de nodo, sino el tipo de interseccion (aunque coincida que es un 2 de tuje)
                            ins(nins)%fibs = [f1,f2] ! fibras que participan en la interseccion
                            ins(nins)%segs = [s1,s2] ! segmentos que participan en la interseccion
                            ins(nins)%segs_f = [s1_f1, s2_f2] ! indices locales de los segmentos dentro de sus fibras
                            ins(nins)%seg1 = [n11, n12]
                            ins(nins)%seg2 = [n21, n22]
                            ins(nins)%inewseg1 = inews1
                            ins(nins)%inewseg2 = inews2
                            ins(nins)%inewnod = inewn
                            ins(nins)%newnodr = rin
                            newisegins(s1) = inews1 ! en caso de que vuelva a tener alguna interseccion
                            newisegins(s2) = inews2 ! en caso de que vuelva a tener alguna interseccion
                        END IF
                    END DO
                END DO
            END DO
        END DO
    END DO
    ! ----------

    ! ----------
    ! aqui ya tengo todas las intersecciones guardadas, ahora tengo que crear la nueva malla
    ! dimensiones generales
    m_out%sidelen = m_in%sidelen
    m_out%diamed = m_in%diamed
    ! ----------
    ! nodos
    m_out%nnods = m_in%nnods + nins ! cada interseccion implica un nuevo nodo
    ALLOCATE( m_out%tipos(m_out%nnods) ) ! adjudico
    ALLOCATE( m_out%rnods(2,m_out%nnods) ) ! adjudico
    m_out%tipos(1:m_in%nnods) = m_in%tipos ! copio los tipos viejos en los tipos nuevos
    m_out%tipos(m_in%nnods+1:m_out%nnods) = 2 ! los tipos nuevos extras son todos tipo interseccion
    m_out%rnods(:,1:m_in%nnods) = m_in%rnods ! copio coordenadas de los nodos viejos en los nodos nuevos
    DO i=1,nins ! copio coordenadas de los nodos interseccion en los nodos nuevos
        m_out%rnods(:,m_in%nnods+i) = ins(i)%newnodr(:)
    END DO
    ! ----------
    ! segmentos
    m_out%nsegs = m_in%nsegs + 2*nins ! cada interseccion implica dos nuevos segmentos
    ALLOCATE( m_out%segs(2,m_out%nsegs) ) ! adjudico
    m_out%segs(:,1:m_in%nsegs) = m_in%segs ! copio los segmentos viejos en los segmentos nuevos
    ! (los segmentos quedan incompletos, hay que modificarlos luego en un bucle sobre las intersecciones)
    DO i=1,nins
        s1 = ins(i)%segs(1)
        s2 = ins(i)%segs(2)
        n11 = m_in%segs(1,s1)
        n12 = m_in%segs(2,s1)
        n21 = m_in%segs(1,s2)
        n22 = m_in%segs(2,s2)
        inews1 = ins(i)%inewseg1
        inews2 = ins(i)%inewseg2
        inewn = ins(i)%inewnod
        ! cambio la conectividad de los dos segmentos viejos
        ! y agrego los dos segmentos nuevos
        m_out%segs(2,s1) = inewn
        m_out%segs(2,s2) = inewn
        m_out%segs(:,inews1) = [inewn, n12]
        m_out%segs(:,inews2) = [inewn, n22]
    END DO
    ! ----------
    ! fibras
    m_out%nfibs = m_in%nfibs
    m_out%lenfibsje = m_in%lenfibsje + 2*nins ! tengo dos nuevos elemento en la conectividad de las fibras por cada interseccion
    ALLOCATE( m_out%fibsds(m_out%nfibs) )
    ALLOCATE( m_out%fibsdls(m_out%nfibs) )
    ALLOCATE( m_out%fibsdths(m_out%nfibs) )
    ALLOCATE( m_out%fibsne(m_out%nfibs) )
    ALLOCATE( m_out%fibsie(m_out%nfibs + 1) )
    ALLOCATE( m_out%fibsje(m_out%lenfibsje) )
    m_out%fibsds = m_in%fibsds
    m_out%fibsdls = m_in%fibsdls
    m_out%fibsdths = m_in%fibsdths
    m_out%fibsie(1) = 1
    ! primero copio la conectividad vieja en la conectividad nueva
    ! teniendo cuidado de donde va cada fibra
    m_out%fibsje = -1 ! inicializo en -1 indicando lugares libres (sin segmento asignado)
    DO i=1,m_in%nfibs
        m_out%fibsne(i) = m_in%fibsne(i) + nfibins(i)
        m_out%fibsie(i+1) = m_out%fibsie(i) + m_out%fibsne(i)
!        ! auxiliares
!        m_in_fie0 = m_in%fibsie(i)
!        m_in_fie1 = m_in%fibsie(i+1)-1
!        m_out_fie0 = m_out%fibsie(i)
!        m_out_fie1 = m_out%fibsie(i+1)-1
!        ! notar que quedan lugares libres (uno por interseccion) en la conectividad nueva para cada fibra que sufrio alguna interseccion:
!        m_out%fibsje(m_out_fie0:m_out_fie0+m_in%fibsne(i)) =  m_in%fibsje(m_in_fie0:m_in_fie1)
    END DO
    !
    m_out%fibsje(1:m_in%lenfibsje) = m_in%fibsje
    ! ahora divido cada segmento en dos recorriendo el je de las fibras
    nins2 = 0
    DO i=1,m_in%lenfibsje
        s1 = m_in%fibsje(i)
        IF (nsegins(s1)>0) THEN
            ! este segmento fue intersectado, hay que partirlo en dos
            nins2 = nins2 + 1
            inews1 = isegins(1,s1) ! nodo nuevo generado al partir al segmento s1 en dos
            ! le hago lugar al nuevo segmento desplazando la segunda parte de la conectividad hacia adelante
            m_out%fibsje(i+nins2+1:m_in%lenfibsje+nins2) = m_out%fibsje(i+nins2:m_in%lenfibsje+nins2-1)
            ! pongo el segmento nuevo
            m_out%fibsje(i+nins2) = inews1
        END IF
    END DO

!    DO i=1,nins
!        f1 = ins(i)%fibs(1)
!        f2 = ins(i)%fibs(2)
!        s1_f1 = ins(i)%segs_f(1)
!        s2_f2 = ins(i)%segs_f(2)
!        inews1 = ins(i)%inewseg1
!        inews2 = ins(i)%inewseg2
!
!    END DO
!    nins2 = 0
!    DO i=1,m_in%nfibs
!        m_out%fibsne(i) = m_in%fibsne(i) + nfibins(i)
!        m_out%fibsie(i+1) = m_out%fibsie(i) + m_out%fibsne(i)
!        ! auxiliares
!        m_in_fie0 = m_in%fibsie(i)
!        m_in_fie1 = m_in%fibsie(i+1)-1
!        m_out_fie0 = m_out%fibsie(i)
!        m_out_fie1 = m_out%fibsie(i+1)-1
!        ! parto de haber copiado la conectividad de la fibra vieja a la nueva (notar que quedan lugares libres al final en la nueva)
!        m_out%fibsje(m_out_fie0:m_out_fie0+m_in%fibsne(i)) =  m_in%fibsje(m_in_fie0:m_in_fie1)
!        ! ahora si tengo segmentos nuevos a mechar, hago un bucle
!        IF (nfibins(i)>0) THEN
!            DO j=1,nfibins(i)
!                nins2 = nins2 + 1
!                ! para cada interseccion, debo:
!                ! 1: hacer lugar para el segmento nuevo moviendo los segmentos subsiguientes hacia adelante
!                ! 2: agregar el segmento nuevo
!                ! si esta bien hecho deberia llegar a completar toda la conectividad nueva
!                ! ifibins(j,i) es la posicion local de la interseccion local j de la fibra global i
!                ! jfibins(j,i) es el indice global del segmento nuevo que va en la posicion local j de la fibra i
!                aux1 = m_out_fie0 - 1 + ifibins(j,i) +j-1 ! posicion local del nuevo nodo
!                aux2 = m_out_fie0 + m_in%fibsne(i)-1 +j-1 ! posicion local del ultimo segmento en conectividad nueva
!                m_out%fibsje(aux1+1:aux2+1) = m_out%fibsje(aux1:aux2)
!                m_out%fibsje(aux1) = jfibins(j,i)
!            END DO
!        END IF
!    END DO
    ! ----------
    ! capas
    m_out%ncaps = m_in%ncaps
    m_out%lencapsje = m_in%lencapsje
    ALLOCATE( m_out%capsne(m_out%ncaps) )
    ALLOCATE( m_out%capsie(m_out%ncaps + 1) )
    ALLOCATE( m_out%capsje(m_out%lencapsje) )
    m_out%capsne = m_in%capsne
    m_out%capsie = m_in%capsie
    m_out%capsje = m_in%capsje
    ! ----------


END SUBROUTINE intersectar_fibras_2
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
SUBROUTINE intersectar_fibras(m_in, m_out)
    ! voy a recorrer m_in, sus capas y fibras, encontrando intersecciones entre fibras
    ! como resultado la malla de salida es una malla aumentada porque tiene mas nodos y mas segmentos
    ! en primera instancia la voy a sobredimensionar en memoria
    ! luego readjudico la memoria correcta
    IMPLICIT NONE
    TYPE(MallaCom), INTENT(IN) :: m_in
    TYPE(MallaCom), INTENT(OUT) :: m_out
    !
    INTEGER :: c1, c2 ! numeracion de capas
    INTEGER :: f1_c1, f1 ! numeracion de fibra 1: local dentro de la capa c1 y global
    INTEGER :: s1_f1, s1 ! numeracion de segmento 1: local dentro de la fibra f1 y global
    INTEGER :: f2_c1, f2 ! numeracion de fibra 2: local y global
    INTEGER :: s2_f2, s2 ! numeracion de segmento 2: local y global
    INTEGER :: n11, n12, n21, n22
    REAL(8) :: r11(2), r12(2), r21(2), r22(2)
    INTEGER :: iStat
    REAL(8) :: rin(2)
    !
    ! preadjudico las variables con el doble de tamanio para que sobre
    INTEGER :: nins ! numero de intersecciones
    TYPE(Interseccion) :: ins(1000) ! preadjudicado grande para que sobre
    INTEGER :: nfibins(m_in%nfibs) ! aqui guardo la cantidad de intersecciones sobre cada fibra
    INTEGER :: ifibins(100, m_in%nfibs) ! aqui guardo en que indice local de cada fibra colocar los segmentos nuevos
    INTEGER :: jfibins(100, m_in%nfibs) ! aqui guardo el indice global de cada segmento nuevo de cada fibra
    !
    INTEGER :: i, j
    INTEGER :: m_in_fie0, m_in_fie1, m_out_fie0, m_out_fie1
    INTEGER :: nins2
    INTEGER :: aux1, aux2, aux3, aux4
    INTEGER :: news1, news2, newn
    !


    ! ----------
    nins = 0
    nfibins = 0
    ifibins = 0
    ! ----------
    ! recorro las capas e intersecto fibras dentro de cada capa y con capas adyacentes
    DO c1 = 1,m_in%ncaps
        ! recorro las fibras de la capa c1
        DO f1_c1 = 1,m_in%capsne(c1)
            f1 = m_in%capsje( m_in%capsie(c1)-1 + f1_c1 )
            ! recorro los segmentos de la fibra f1
            DO s1_f1 = 1,m_in%fibsne(f1)
                s1 = m_in%fibsje( m_in%fibsie(f1)-1 + s1_f1 )
                n11 = m_in%segs(1,s1)
                n12 = m_in%segs(2,s1)
                r11 = m_in%rnods(:,n11)
                r12 = m_in%rnods(:,n12)
                ! recorro las fibras restantes de la capa c1 para intersectar
                DO f2_c1 = f1_c1+1,m_in%capsne(c1)
                    f2 = m_in%capsje( m_in%capsie(c1)-1 + f2_c1 )
                    ! recorro los segmentos de la fibra f2
                    DO s2_f2 = 1,m_in%fibsne(f2)
                        s2 = m_in%fibsje( m_in%fibsie(f2)-1 + s2_f2 )
                        n21 = m_in%segs(1,s2)
                        n22 = m_in%segs(2,s2)
                        r21 = m_in%rnods(:,n21)
                        r22 = m_in%rnods(:,n22)
                        CALL intersectar_segmentos(r11, r12, r21, r22, rin, iStat)
                        IF (iStat == 2) THEN ! 2 es interseccion tipo medio-medio
                            nins = nins + 1
                            ins(nins)%fibs = [f1,f2]
                            ins(nins)%segs = [s1,s2]
                            ins(nins)%inewnod = m_in%nnods + nins
                            ins(nins)%newnodr = rin
                            ins(nins)%inewseg1 = m_in%nsegs + 2*nins - 1
                            ins(nins)%inewseg2 = m_in%nsegs + 2*nins
                            ins(nins)%tipo = iStat ! ojo que no es el tipo de nodo, sino el tipo de interseccion (aunque coincida que es un 2 de tuje)
                            nfibins(f1) = nfibins(f1) + 1
                            nfibins(f2) = nfibins(f2) + 1
                            ifibins(nfibins(f1),f1) = s1_f1 + 1
                            ifibins(nfibins(f2),f2) = s2_f2 + 1
                            jfibins(nfibins(f1),f1) = m_in%nsegs + 2*nins - 1
                            jfibins(nfibins(f2),f2) = m_in%nsegs + 2*nins
                        END IF
                    END DO
                END DO
            END DO
        END DO
    END DO
    ! ----------

    ! ----------
    ! aqui ya tengo todas las intersecciones guardadas, ahora tengo que crear la nueva malla
    ! dimensiones generales
    m_out%sidelen = m_in%sidelen
    m_out%diamed = m_in%diamed
    ! ----------
    ! nodos
    m_out%nnods = m_in%nnods + nins ! cada interseccion implica un nuevo nodo
    ALLOCATE( m_out%tipos(m_out%nnods) ) ! adjudico
    ALLOCATE( m_out%rnods(2,m_out%nnods) ) ! adjudico
    m_out%tipos(1:m_in%nnods) = m_in%tipos ! copio los tipos viejos en los tipos nuevos
    m_out%tipos(m_in%nnods+1:m_out%nnods) = 2 ! los tipos nuevos extras son todos tipo interseccion
    m_out%rnods(:,1:m_in%nnods) = m_in%rnods ! copio coordenadas de los nodos viejos en los nodos nuevos
    DO i=1,nins ! copio coordenadas de los nodos interseccion en los nodos nuevos
        m_out%rnods(:,m_in%nnods+i) = ins(i)%newnodr(:)
    END DO
    ! ----------
    ! segmentos
    m_out%nsegs = m_in%nsegs + 2*nins ! cada interseccion implica dos nuevos segmentos
    ALLOCATE( m_out%segs(2,m_out%nsegs) ) ! adjudico
    m_out%segs(:,1:m_in%nsegs) = m_in%segs ! copio los segmentos viejos en los segmentos nuevos
    ! (los segmentos quedan incompletos, hay que modificarlos luego en un bucle sobre las intersecciones)
    DO i=1,nins
        s1 = ins(i)%segs(1)
        s2 = ins(i)%segs(2)
        n11 = m_in%segs(1,s1)
        n12 = m_in%segs(2,s1)
        n21 = m_in%segs(1,s2)
        n22 = m_in%segs(2,s2)
        news1 = ins(i)%inewseg1
        news2 = ins(i)%inewseg2
        newn = ins(i)%inewnod
        ! cambio la conectividad de los dos segmentos viejos
        ! y agrego los dos segmentos nuevos
        m_out%segs(2,s1) = newn
        m_out%segs(2,s2) = newn
        m_out%segs(:,news1) = [newn, n12]
        m_out%segs(:,news2) = [newn, n22]
    END DO
    ! ----------
    ! fibras
    m_out%nfibs = m_in%nfibs
    m_out%lenfibsje = m_in%lenfibsje + 2*nins ! tengo dos nuevos elemento en la conectividad de las fibras por cada interseccion
    ALLOCATE( m_out%fibsds(m_out%nfibs) )
    ALLOCATE( m_out%fibsdls(m_out%nfibs) )
    ALLOCATE( m_out%fibsdths(m_out%nfibs) )
    ALLOCATE( m_out%fibsne(m_out%nfibs) )
    ALLOCATE( m_out%fibsie(m_out%nfibs + 1) )
    ALLOCATE( m_out%fibsje(m_out%lenfibsje) )
    m_out%fibsds = m_in%fibsds
    m_out%fibsdls = m_in%fibsdls
    m_out%fibsdths = m_in%fibsdths
    m_out%fibsie(1) = 1
    nins2 = 0
    DO i=1,m_in%nfibs
        m_out%fibsne(i) = m_in%fibsne(i) + nfibins(i)
        m_out%fibsie(i+1) = m_out%fibsie(i) + m_out%fibsne(i)
        ! auxiliares
        m_in_fie0 = m_in%fibsie(i)
        m_in_fie1 = m_in%fibsie(i+1)-1
        m_out_fie0 = m_out%fibsie(i)
        m_out_fie1 = m_out%fibsie(i+1)-1
        ! parto de haber copiado la conectividad de la fibra vieja a la nueva (notar que quedan lugares libres al final en la nueva)
        m_out%fibsje(m_out_fie0:m_out_fie0+m_in%fibsne(i)) =  m_in%fibsje(m_in_fie0:m_in_fie1)
        ! ahora si tengo segmentos nuevos a mechar, hago un bucle
        IF (nfibins(i)>0) THEN
            DO j=1,nfibins(i)
                nins2 = nins2 + 1
                ! para cada interseccion, debo:
                ! 1: hacer lugar para el segmento nuevo moviendo los segmentos subsiguientes hacia adelante
                ! 2: agregar el segmento nuevo
                ! si esta bien hecho deberia llegar a completar toda la conectividad nueva
                ! ifibins(j,i) es la posicion local de la interseccion local j de la fibra global i
                ! jfibins(j,i) es el indice global del segmento nuevo que va en la posicion local j de la fibra i
                aux1 = m_out_fie0 - 1 + ifibins(j,i) +j-1 ! posicion local del nuevo nodo
                aux2 = m_out_fie0 + m_in%fibsne(i)-1 +j-1 ! posicion local del ultimo segmento en conectividad nueva
                m_out%fibsje(aux1+1:aux2+1) = m_out%fibsje(aux1:aux2)
                m_out%fibsje(aux1) = jfibins(j,i)
            END DO
        END IF
    END DO
    ! ----------
    ! capas
    m_out%ncaps = m_in%ncaps
    m_out%lencapsje = m_in%lencapsje
    ALLOCATE( m_out%capsne(m_out%ncaps) )
    ALLOCATE( m_out%capsie(m_out%ncaps + 1) )
    ALLOCATE( m_out%capsje(m_out%lencapsje) )
    m_out%capsne = m_in%capsne
    m_out%capsie = m_in%capsie
    m_out%capsje = m_in%capsje
    ! ----------


END SUBROUTINE intersectar_fibras
! ================================================================================
! ================================================================================


! ==============================================================================
END MODULE class_malla_completa
! ==============================================================================
