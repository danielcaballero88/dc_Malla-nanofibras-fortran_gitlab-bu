MODULE Aux


real(8), parameter :: pi = 3.1415926536
REAL(8), PARAMETER :: pi2 = pi*2.0d0, pi15 = pi*1.5d0, pi05 = pi*0.5d0
real(8), parameter :: e = 2.7182818285


CONTAINS

! ================================================================================
! ================================================================================
FUNCTION FindStringInFile(Str, ioUnit, Mandatory) RESULT (iError)
    ! ================================================================================
    ! Busca un String en un archivo (STR), sino lo encuentra rebovina el archivo
    ! y pone iError < 0 como indicador de no hallazgo
    ! Str: String to find, ioUnit: Unit assigned to Input File; iError: Error Status variable
    ! ================================================================================
    IMPLICIT NONE
    ! ================================================================================
    ! Parameters
    LOGICAL,PARAMETER  :: Segui=.True.
    ! ================================================================================
    ! Arguments
    CHARACTER(*), INTENT(IN) :: Str
    INTEGER, INTENT (IN) :: ioUnit
    LOGICAL, OPTIONAL, INTENT(IN) :: Mandatory
    ! ================================================================================
    ! Locals
    LOGICAL :: MandatoryL
    CHARACTER(LEN=120) :: DummyString
    CHARACTER(LEN=120) :: upStr
    INTEGER :: iError
    INTEGER :: Leng
    ! ================================================================================

    ! ================================================================================
    IF ( PRESENT(Mandatory) ) THEN
        MandatoryL = Mandatory
    ELSE
        MandatoryL = .FALSE.
    END IF
    ! ================================================================================
    iError=0
    Leng = LEN_TRIM(Str)
    upStr = Upper_Case(Str)       ! Line added by NB. Converts Str to Upper Case string
    ! ================================================================================
    REWIND(ioUnit)
    Search_Loop: DO WHILE (segui)
        READ(ioUnit,'(1A120)',IOSTAT=iError) DummyString
        DummyString = Upper_Case(DummyString)   ! line added by NB
        !       if (iError==59) backspace(ioUnit)
        IF (iError.lt.0) THEN
            REWIND (ioUnit)
            EXIT Search_Loop
        END IF
        IF ( DummyString(1:1)    /=    '*'      ) CYCLE Search_Loop
        IF ( DummyString(1:Leng) == upStr(1:Leng) ) EXIT Search_Loop
    END DO Search_Loop
    ! ================================================================================
    IF (MandatoryL) THEN
        IF ( .NOT. ( iError == 0 ) ) THEN
            WRITE(*,*) upStr, 'NOT FOUND'
            STOP
        END IF
    END IF
    ! ================================================================================

END FUNCTION FindStringInFile
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
FUNCTION Upper_Case(string)

    !-------------------------------------------------------------------------------------------
    !! The Upper_Case function converts the lower case characters of a string to upper case one.
    !! Use this function in order to achieve case-insensitive: all character variables used
    !! are pre-processed by Uppper_Case function before these variables are used. So the users
    !! can call functions without pey attention of case of the keywords passed to the functions.
    !-------------------------------------------------------------------------------------------

    IMPLICIT NONE

    !-------------------------------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(IN):: string     ! string to be converted
    CHARACTER(LEN=LEN(string))::   Upper_Case ! converted string
    INTEGER::                      n1         ! characters counter
    !-------------------------------------------------------------------------------------------

    Upper_Case = string
    DO n1=1,LEN(string)
        SELECT CASE(ICHAR(string(n1:n1)))
            CASE(97:122)
                Upper_Case(n1:n1)=CHAR(ICHAR(string(n1:n1))-32) ! Upper case conversion
        END SELECT
    ENDDO
    RETURN

    !-----------------------------------------------------------------------------------------------
END FUNCTION Upper_Case
! ================================================================================
! ================================================================================



! ================================================================================
! ================================================================================
FUNCTION get_file_unit() RESULT(lu)
    ! ================================================================================
    ! get_file_unit returns a unit number that is not in use
    ! ================================================================================
    INTEGER :: lu_max,  lu, checkIOSTAT
    LOGICAL   checkOPENED
    INTEGER, PARAMETER :: m = 99
    !
    DO lu = m,1,-1
        INQUIRE (UNIT=lu, OPENED=checkOPENED, IOSTAT=checkIOSTAT)
        IF (checkIOSTAT.ne.0) CYCLE
        IF (.NOT.checkOPENED) EXIT
    END DO
    !

    RETURN
END FUNCTION get_file_unit
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
SUBROUTINE angulo_de_segmento(r0, r1, th, iError)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: r0(2), r1(2)
    REAL(8), INTENT(OUT) :: th ! angulo
    INTEGER, INTENT(OUT) :: iError ! -1: No hay interseccion, 0: hay interseccion en el medio, 1: hay interseccion en un punto extremo de uno de los segmentos
    !
    REAL(8) :: dr(2), dx, dy


    iError = 0
    ! ----------
    dr = r1 - r0
    dx = dr(1)
    dy = dr(2)
    IF (iguales(dx,0.d0)) THEN ! segmento vertical
        IF (iguales(dy,0.d0)) THEN ! segmento nulo
            WRITE(*,*) "Error, segmento nulo!"
            STOP
        ELSE IF (dy>0.d0) THEN ! vertical hacia arriba
            th = pi05
            RETURN
        ELSE ! vertical hacia abajo
            th = pi15
            RETURN
        END IF
    ELSE IF (iguales(dy,0.d0)) THEN ! segmento horizontal
        IF (dx>0.d0) THEN ! segmento hacia derecha
            th = 0.d0
            RETURN
        ELSE ! segmento hacia izquierda
            th = pi
            RETURN
        END IF
    ELSE ! segmento oblicuo
        IF (dx<0) THEN ! segundo o tercer cuadrante
            th = pi + DATAN(dy/dx)
            RETURN
        ELSE IF (dy>0) THEN ! primer cuadrante
            th = DATAN(dy/dx)
            RETURN
        ELSE ! cuarto cuadrante
            th = pi2 + DATAN(dy/dx)
            RETURN
        END IF
    END IF
    ! ----------
    th = 0.d0
    iError = 1 ! si llegamos hasta aca estamos mal
    WRITE(*,*) "Error calculando el angulo de un segmento!"
    STOP


END SUBROUTINE angulo_de_segmento
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
SUBROUTINE intersectar_segmentos(rs1n1, rs1n2, rs2n1, rs2n2, r_in, iStatus)
    ! Busca interseccion entre dos segmentos
    ! devuelve iStatus >= 0 si hay interseccion, -1 si no hay interseccion
    ! r_in son las coordenadas de la interseccion
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: rs1n1(2), rs1n2(2), rs2n1(2), rs2n2(2)
    REAL(8), INTENT(OUT) :: r_in(2)
    INTEGER, INTENT(OUT) :: iStatus ! -1: No hay interseccion, 0: hay interseccion en el medio, 1: hay interseccion en un punto extremo de uno de los segmentos
    !
    INTEGER :: ierr
    REAL(8) :: th1, th2
    LOGICAL :: paralelos
    REAL(8) :: bot1, rig1, top1, lef1, bot2, rig2, top2, lef2
    REAL(8) :: maxbot, minrig, mintop, maxlef
    LOGICAL :: lejos
    REAL(8) :: th2rel
    LOGICAL :: m2relinf
    REAL(8) :: cos1, sin1
    REAL(8) :: A1(2,2) ! matriz de cambio de base con th1
    REAL(8) :: s12(2), s21(2), s22(2), s_in ! posiciones relativas en el sistema intrinseco al segmento 1 (s11 = [0,0])
    REAL(8) :: m2rel
    LOGICAL :: inex11, inex12, inex21, inex22
    INTEGER :: inex1, inex2 ! interseccion en extremo


    ! ----------
    ! chequeo paralelismo (no habria interseccion)
    CALL angulo_de_segmento(rs1n1, rs1n2, th1, ierr)
    CALL angulo_de_segmento(rs2n1, rs2n2, th2, ierr)
    paralelos = iguales(th1, th2, pi*1.d-5)
    IF (paralelos) THEN
        iStatus = -1
        RETURN
    END IF
    ! ----------
    ! chequeo lejania (no habria interseccion)
    bot1 = MIN(rs1n1(2),rs1n2(2))
    rig1 = MAX(rs1n1(1),rs1n2(1))
    top1 = MAX(rs1n1(2),rs1n2(2))
    lef1 = MIN(rs1n1(1),rs1n2(1))
    bot2 = MIN(rs2n1(2),rs2n2(2))
    rig2 = MAX(rs2n1(1),rs2n2(1))
    top2 = MAX(rs2n1(2),rs2n2(2))
    lef2 = MIN(rs2n1(1),rs2n2(1))
    ! ----------
    maxbot = MAX(bot1,bot2)
    minrig = MIN(rig1,rig2)
    mintop = MIN(top1,top2)
    maxlef = MAX(lef1,lef2)
    ! ----------
    lejos = ( ((maxlef-minrig)>1.0d-8) .OR. ((maxbot-mintop)>1.0d-8) )
    IF (lejos) THEN
        iStatus = -1
        RETURN
    END IF
    ! ----------
    ! sistema de referencia relativo, intrinseco al segmento 1
    th2rel = th2 - th1
    DO WHILE ( (th2rel<0.d0) .OR. (th2rel>=pi2) )
        IF (th2rel<0.d0) THEN
            th2rel = th2rel + pi2
        ELSE IF (th2rel>=pi2) THEN
            th2rel = th2rel - pi2
        END IF
    END DO
    ! me fijo que en el sistema intrinsieco al segmento 1, el segmento 2 no sea vertical
    m2relinf = ( (iguales(th2rel,pi05,pi*1.d-5)) .OR. (iguales(th2rel,pi15,pi*1.d-5)) )
    ! coseno y seno de th1
    cos1 = DCOS(th1)
    sin1 = DSIN(th1)
    A1 = RESHAPE( [cos1, -sin1, sin1, cos1], [2,2] )
    s12 = MATMUL(A1,rs1n2-rs1n1)
    s21 = MATMUL(A1,rs2n1-rs1n1)
    s22 = MATMUL(A1,rs2n2-rs1n1)
    ! ----------
    ! chequeo que el segmento 2 corte al eje del segmento 1
    IF ( NINT(SIGN(1.0d0,s21(2))) == NINT(SIGN(1.0d0,s22(2))) ) THEN
        iStatus = -1
        RETURN
    END IF
    ! ----------
    ! calculo a que distancia sobre el eje del segmento 1 se produce la posible interseccion
    IF (m2relinf) THEN
        s_in = s21(1)
    ELSE
        m2rel = DTAN(th2rel)
        s_in = s21(1) - s21(2)/m2rel
    END IF
    ! ----------
    ! ahora chequeo si la interseccion se produjo por fuera del segmento 1
    IF ((s_in<0.d0) .OR. (s_in>s12(1))) THEN
        iStatus = -1
        RETURN
    END IF
    ! ----------
    ! HAY INTERSECCION
    ! me fijo si ocurre en el medio de los segmentos, o si coincide con los nodos ya existentes (extremos de segmentos)
    inex11 = iguales(s_in,0.d0)
    inex12 = iguales(s_in,s12(1))
    inex21 = iguales(s21(2),0.d0)
    inex22 = iguales(s22(2),0.d0)
    ! HAY VARIAS POSIBILIDADES
    ! 4: extremo-extremo
    IF ((inex11 .OR. inex12) .AND. (inex21 .OR. inex22)) THEN
!        WRITE(*,*) " Interseccion tipo extremo-extremo. Por ahora no esta implementada."
        iStatus = -1
        RETURN
    END IF
    ! 3: extremo-medio
    IF ((inex11 .OR. inex12) .OR. (inex21 .OR. inex22)) THEN
!        WRITE(*,*) " Interseccion tipo extremo-medio. Por ahora no esta implementada."
        iStatus = -1
        RETURN
    END IF
    ! 2: medio-medio
    ! chequeo redundante
    IF ( (s_in>0.d0) .AND. (s_in<s12(1)) ) THEN
        iStatus = 2
        r_in = rs1n1 + s_in * [cos1, sin1]
        RETURN
    END IF
    ! ----------
    ! si llegue hasta aca hay algo mal
    WRITE(*,*) "Error en subrutina calcular interseccion entre segmentos!"
    STOP
    ! ----------



END SUBROUTINE intersectar_segmentos
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
FUNCTION iguales(x,y,tol_in) RESULT(igual)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x,y
    REAL(8), OPTIONAL, INTENT(IN) :: tol_in
    LOGICAL :: igual
    REAL(8) :: tol

    IF (PRESENT(tol_in)) THEN
        tol = tol_in
    ELSE
        tol = 1.d-6
    END IF

    igual = DABS(x-y)<tol

END FUNCTION
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
subroutine diagonal_iterativo(n, A, b, x, iStatus)
    ! ----------
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: A(n,n), b(n)
    real(8), intent(out) :: x(n)
    integer, intent(out) :: iStatus
    ! ----------
    integer, parameter :: maxiter = 100
    real(8), parameter :: tol = 1.d-6
    ! ----------
    real(8) :: x1(n)
    real(8) :: tmp
    integer i,j,k
    real(8) :: error_dx
    ! ----------

    ! ----------
    ! valor semilla
    x(:) = 0.d0
    x1(:) = x(:)
    ! ----------
    ! iteraciones
    do k=1,maxiter
        ! hago la sumatoria de A(i,j)*x(j) exceptuando la diagonal
        do i=1,n
            ! calculo cada valor de x
            x(i) = b(i) / A(i,i)
        end do
        ! chequeo error y convergencia
        error_dx = maxval(dabs(x-x1))
        if (error_dx<tol) then
            iStatus = 0
            return
        end if
        ! actualizo el valor de la iteracion anterior
        x1(:) = x(:)
    end do
    ! ----------
    ! si llegue hasta aca es porque llegue a maxiter
    write(*,*) "WARNING: maxiter alcanzado en metodo diagonal_iterativo"
    iStatus = 1
    ! ----------

    ! ----------
end subroutine diagonal_iterativo
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
subroutine jacobi(n, A, b, x, iStatus)
    ! ----------
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: A(n,n), b(n)
    real(8), intent(out) :: x(n)
    integer, intent(out) :: iStatus
    ! ----------
    integer, parameter :: maxiter = 100
    real(8), parameter :: tol = 1.d-6
    ! ----------
    real(8) :: x1(n)
    real(8) :: tmp
    integer i,j,k
    real(8) :: error_dx
    ! ----------

    ! ----------
    ! valor semilla
    x(:) = 0.d0
    x1(:) = x(:)
    ! ----------
    ! iteraciones de gauss-seidel
    do k=1,maxiter
        ! hago la sumatoria de A(i,j)*x(j) exceptuando la diagonal
        do i=1,n
            tmp = 0
            do j=1,n
                if (j<i) then
                    tmp = tmp + A(i,j)*x1(j)
                else if (j>i) then
                    tmp = tmp + A(i,j)*x1(j)
                end if
            end do
            ! calculo cada valor de x
            x(i) = (b(i)-tmp) / A(i,i)
        end do
        ! chequeo error y convergencia
        error_dx = maxval(dabs(x-x1))
        if (error_dx<tol) then
            iStatus = 0
            return
        end if
        ! actualizo el valor de la iteracion anterior
        x1(:) = x(:)
    end do
    ! ----------
    ! si llegue hasta aca es porque llegue a maxiter
    write(*,*) "WARNING: maxiter alcanzado en metodo jacobi"
    iStatus = 1
    ! ----------

    ! ----------
end subroutine jacobi
! ================================================================================
! ================================================================================


! ================================================================================
! ================================================================================
subroutine gauss_seidel(n, A, b, x, maxiter, tol, iStatus)
    ! ----------
    implicit none
    integer, intent(in) :: n
    real(8), intent(inout) :: A(n,n), b(n)
    real(8), intent(out) :: x(n)
    integer, intent(out) :: iStatus
    integer, intent(in) :: maxiter
    real(8), intent(in) :: tol
    ! ----------
    real(8) :: x1(n)
    real(8) :: tmp, mintmp, maxtmp
    integer i,j,k
    real(8) :: error_dx
    ! ----------

    ! ----------
    ! valor semilla
    x(:) = 0.d0
    x1(:) = x(:)
    ! ----------
    ! divido cada fila de matriz y vector por el valor de la diagonal
    ! ojo que no deberia ser un valor muy pequeno porque error
!    mintmp = dabs( A(1,1) )
!    maxtmp = dabs( A(1,1) )
!    do i=1,n
!        tmp = 1.d0 / dabs(A(i,i))
!        A(i,:) = A(i,:) * tmp
!        b(i) = b(i) * tmp
!        if (tmp<mintmp) mintmp = tmp
!        if (tmp>maxtmp) maxtmp = tmp
!    end do
    ! ----------
    ! iteraciones de gauss-seidel
    do k=1,maxiter
        ! hago la sumatoria de A(i,j)*x(j) exceptuando la diagonal
        do i=1,n
            tmp = 0
            do j=1,n
                if (j<i) then
                    tmp = tmp + A(i,j)*x(j)
                else if (j>i) then
                    tmp = tmp + A(i,j)*x1(j)
                end if
            end do
            ! calculo cada valor de x
            x(i) = (b(i)-tmp) / A(i,i)
        end do
        ! chequeo error y convergencia
        error_dx = maxval(dabs(x-x1))
        write(*,*) "k=", k, "error_dx=", error_dx, "/ tol=", tol
        if (error_dx<tol) then
            write(*,*) "Gauss-Seidel converge"
            iStatus = 0
            return
        end if
        ! actualizo el valor de la iteracion anterior
        x1(:) = x(:)
    end do
    ! ----------
    ! si llegue hasta aca es porque llegue a maxiter
    write(*,*) "WARNING: maxiter alcanzado en metodo gauss_seidel: ", maxiter
    iStatus = 1
    ! ----------

    ! ----------
end subroutine gauss_seidel
! ================================================================================
! ================================================================================

! ================================================================================
! ================================================================================
subroutine directo(n, A, b, x, iStatus)
    ! ----------
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: A(n,n), b(n)
    real(8), intent(out) :: x(n)
    integer, intent(out) :: iStatus
    ! ----------
    integer, parameter :: maxiter = 100
    real(8), parameter :: tol = 1.d-6
    ! ----------
    real(8) :: A1(n,n)
    real(8) :: residuo(n)
    real(8) :: tmp
    integer i,j,k
    real(8) :: error_dx
    ! ----------


    ! invierto la matriz A
    A1 = matinv(A)
    ! calculo x haciendo x=A1 b
    do i=1,n
        x(i) = 0.d0
        do j=1,n
            x(i) = x(i) + A1(i,j)*b(j)
        end do
    end do
    iStatus = 0

    ! ----------
end subroutine directo
! ================================================================================
! ================================================================================

function matinv(A) result(Ainv)
  real(8), dimension(:,:), intent(in) :: A
  real(8), dimension(size(A,1),size(A,2)) :: Ainv

  real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function matinv



END MODULE Aux
