program jcb
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: pi = acos(-1.0_dp)               ! Literalmente pi
    real(dp), parameter :: sigma = 5.67e-8_dp               ! Stefan-Boltzmann [W/m²K⁴]
    real(dp), parameter :: AU_km = 1.496e8_dp               ! 1 UA en km
    real(dp), parameter :: Rsun_km = 695700.0_dp            ! 1 Radio solar en km
    real(dp), parameter :: Rsun_to_AU = Rsun_km / AU_km     ! Conversión Rsun a UA
    real(dp), parameter :: L_sun = 3.828e26_dp              ! Luminosidad solar [W]
    real(dp), parameter :: F_sun = 1361.0_dp                ! Irradiancia solar a 1 UA [W/m²]
    real(dp), parameter :: T_sun = 5778.0_dp                ! Temperatura solar [K]

    real(dp), dimension(3) :: E_1, E_2, CM                  ! Vectores de posición [UA]
    real(dp) :: m1, Re1, L_1, Temp1, m2, Re2, L_2, Temp2    ! Masas [Msun], Radios [Rsun], Temperaturas [K], Luminosidades [Lsun]
    real(dp) :: masa_binaria, Dis                           ! Masa total [Msun], Distancia entre estrellas [UA]
    real(dp) :: d, R_out, R_in                              ! Resolución [UA], Radios del disco [UA]
    real(dp) :: i, j, k, h                                  ! Posiciones en la malla [UA], Distancia radial [UA]
    real(dp) :: flujo_en_p, temperatura_en_p
    integer :: zona
    integer, parameter :: archivo = 10

    ! Datos de las estrellas
    E_1 = [0.0_dp, 0.0_dp, 0.0_dp]            ! Posición inicial [UA]
    Re1 = 1.0_dp                              ! Radio en [Rsun]
    Temp1 = 5778.0_dp                         ! Temperatura [K]
    m1 = 1.0_dp                               ! Masa [Msun]

    E_2 = [0.4_dp, 0.0_dp, 0.0_dp]            ! Posición [UA]
    Re2 = 1.0_dp                              ! Radio [Rsun]
    Temp2 = 5778.0_dp                         ! Temperatura [K]
    m2 = 1.0_dp                               ! Masa [Msun]

    ! Calcular centro de masa y distancia entre estrellas [UA]
    Dis = sqrt((E_2(1)-E_1(1))**2 + (E_2(2)-E_1(2))**2 + (E_2(3)-E_1(3))**2)
    write(*,*) "Distancia entre estrellas [UA]:", Dis
    masa_binaria = m1 + m2
    CM = (m1*E_1 + m2*E_2)/masa_binaria

    ! Convertir coordenadas al sistema centrado en el CM [UA]
    E_1 = E_1 - CM
    E_2 = E_2 - CM

    ! Parámetros de la malla [UA]
    d = 5000000.0_dp / AU_km                  ! Resolución [UA] (5000000 km -> UA)
    R_in = Dis * 3.0                         ! Radio interno [UA]
    R_out = 15.0_dp                            ! Radio externo [UA]

    ! Cálculo de luminosidades en [Lsun]
    L_1 = (Re1)**2 * (Temp1/T_sun)**4
    L_2 = (Re2)**2 * (Temp2/T_sun)**4

    open(unit=archivo, file="jcb.dat", status="replace")
    
    ! Escribir datos de estrellas y CM (convertidos a km para consistencia con output anterior)
    write(archivo, '(a,3f25.12)') '#POSSTAR', E_1 * AU_km
    write(archivo, '(a, f12.2, 1X, f12.2)') '#STAR', Re1 * Rsun_km, Temp1
    write(archivo, '(a,3f25.12)') '#POSSTAR', E_2 * AU_km
    write(archivo, '(a, f12.2, 1X, f12.2)') '#STAR', Re2 * Rsun_km, Temp2
    write(archivo, '(a,f25.12)') '#DistanciaEstrellas', Dis * AU_km
    write(archivo, '(a,3f25.12)') '#CM  ', 0.0_dp, 0.0_dp, 0.0_dp
    write(10,'(a, 1X, f16.2, 1X, f16.2, 1X, f12.2)') '#DSIZE', R_in * AU_km, R_out * AU_km, d * AU_km
    
    ! Generar malla centrada en el CM [UA]
    i = -R_out
    k = 0.0_dp
    do while (i <= R_out)
        j = -R_out
        do while (j <= R_out)
            j = j + d
            h = sqrt(i**2 + j**2 + k**2)
            if (h >= R_in .and. h <= R_out) then
                call determinar_zona_cartesiana(i, j, k, E_1, E_2, Re1, Re2, zona)
                flujo_en_p = calcular_flujo(i, j, k, E_1, E_2, zona, L_1, L_2)
                temperatura_en_p = calcular_temperatura(flujo_en_p, sigma)
                write(10,'(F20.6, 1X, F20.6, 1X, I4, 1X, F20.6)') i * AU_km, j * AU_km, zona, temperatura_en_p 
            end if
        end do
        i = i + d
    end do

    close(archivo)
    
contains

    subroutine calcular_distancias(x, y, z, E1, E2, r1, r2)
        real(dp), intent(in) :: x, y, z
        real(dp), dimension(3), intent(in) :: E1, E2
        real(dp), intent(out) :: r1, r2
        r1 = sqrt( (x - E1(1))**2 + (y - E1(2))**2 + (z - E1(3))**2 )
        r2 = sqrt( (x - E2(1))**2 + (y - E2(2))**2 + (z - E2(3))**2 )
    end subroutine calcular_distancias

    subroutine determinar_zona_cartesiana(x, y, z, E1, E2, Re1, Re2, zona)
        real(dp), intent(in) :: x, y, z
        real(dp), dimension(3), intent(in) :: E1, E2
        real(dp), intent(in) :: Re1, Re2  ! En [Rsun]
        integer, intent(out) :: zona
        real(dp) :: r1, r2, theta1, theta2, alpha, arg_alpha, dot_product_val
        real(dp) :: Re1_au, Re2_au

        ! Convertir radios estelares a [UA]
        Re1_au = Re1 * Rsun_to_AU
        Re2_au = Re2 * Rsun_to_AU

        call calcular_distancias(x, y, z, E1, E2, r1, r2)
    
        ! Cálculo de radios angulares
        theta1 = merge(atan(Re1_au/r1), 0.0_dp, r1 > 0.0_dp)
        theta2 = merge(atan(Re2_au/r2), 0.0_dp, r2 > 0.0_dp)
    
        ! Cálculo del ángulo entre vectores
        dot_product_val = (E1(1)-x)*(E2(1)-x) + (E1(2)-y)*(E2(2)-y) + (E1(3)-z)*(E2(3)-z)
        if (r1 > 0.0_dp .and. r2 > 0.0_dp) then
            arg_alpha = dot_product_val / (r1 * r2)
            arg_alpha = max(min(arg_alpha, 1.0_dp), -1.0_dp)
            alpha = acos(arg_alpha)
        else
            alpha = 0.0_dp
        end if
    
        ! Lógica de zonas
        if (alpha <= abs(theta1 - theta2)) then
            zona = 2    ! Zona de sombra total
        else if (alpha < (theta1 + theta2)) then
            zona = 1    ! Zona de penumbra
        else
            zona = 0    ! Zona iluminada
        end if
    end subroutine determinar_zona_cartesiana

    subroutine calcular_alpha(x, y, z, E1, E2, alpha)
        real(dp), intent(in) :: x, y, z
        real(dp), dimension(3), intent(in) :: E1, E2
        real(dp), intent(out) :: alpha
        real(dp) :: dot_product_val, r1, r2
        call calcular_distancias(x, y, z, E1, E2, r1, r2)
        if (r1 > 0.0_dp .and. r2 > 0.0_dp) then
            dot_product_val = (E1(1)-x)*(E2(1)-x) + (E1(2)-y)*(E2(2)-y) + (E1(3)-z)*(E2(3)-z)
            alpha = acos(dot_product_val / (r1 * r2))
        else
            alpha = 0.0_dp
        end if
    end subroutine calcular_alpha

    subroutine calcular_fraccion_bloqueo(theta1, theta2, alpha, f)
        real(dp), intent(in) :: theta1, theta2, alpha
        real(dp), intent(out) :: f
        real(dp) :: d, Ri, Rj, theta_i, theta_j, area_tri, arg1, arg2
    
        Ri = theta1
        Rj = theta2
        d = alpha
    
        if (d >= Ri + Rj .or. d <= abs(Ri - Rj)) then
            f = 0.0_dp
            return
        end if
    
        arg1 = (cos(d) - cos(Ri)*cos(Rj)) / (sin(Ri)*sin(Rj))
        arg1 = max(min(arg1, 1.0_dp), -1.0_dp)
        theta_i = acos(arg1)
    
        arg2 = (cos(Rj) - cos(Ri)*cos(d)) / (sin(Ri)*sin(d))
        arg2 = max(min(arg2, 1.0_dp), -1.0_dp)
        theta_j = acos(arg2)
    
        area_tri = 2.0_dp*(theta_i + theta_j - pi)*Rj**2
        f = area_tri/(2.0_dp*pi*Rj**2*(1.0_dp - cos(Rj)))
        f = max(min(f, 1.0_dp), 0.0_dp)
    end subroutine calcular_fraccion_bloqueo

    function calcular_flujo(x, y, z, E1, E2, zona, L_1, L_2) result(f_t)
        real(dp), intent(in) :: x, y, z, L_1, L_2  ! L_1, L_2 en [Lsun]
        real(dp), dimension(3), intent(in) :: E1, E2
        integer, intent(in) :: zona
        real(dp) :: f_t, r1, r2, alpha, f
        real(dp) :: theta1, theta2
        real(dp) :: gamma1, gamma2
        real(dp) :: Re1_au, Re2_au
        real(dp), parameter :: Alb = 0.2_dp
        real(dp), parameter :: haltura = 0.05e8_dp / AU_km  ! Altura del disco [UA] (5e6 km -> UA)
        
        call calcular_distancias(x, y, z, E1, E2, r1, r2)
        
        ! Factores de corrección geométrica (ángulo de incidencia)
        if (r1 > 0.0_dp) then
            gamma1 = haltura / sqrt(r1**2 + haltura**2)
        else
            gamma1 = 0.0_dp
        end if
        if (r2 > 0.0_dp) then
            gamma2 = haltura / sqrt(r2**2 + haltura**2)
        else
            gamma2 = 0.0_dp
        end if
        
        ! Convertir radios estelares a [UA] para cálculo de radios angulares
        Re1_au = Re1 * Rsun_to_AU
        Re2_au = Re2 * Rsun_to_AU
        theta1 = merge(atan(Re1_au/r1), 0.0_dp, r1 > 0.0_dp)
        theta2 = merge(atan(Re2_au/r2), 0.0_dp, r2 > 0.0_dp)
        
        call calcular_alpha(x, y, z, E1, E2, alpha)
        
        select case(zona)
            case(0)  ! Zona completamente iluminada
                f_t = (L_1 * F_sun * (1-Alb) * gamma1) / (r1**2) + &
                      (L_2 * F_sun * (1-Alb) * gamma2) / (r2**2)
                
            case(1)  ! Zona de penumbra
                call calcular_fraccion_bloqueo(theta1, theta2, alpha, f)
                if (theta1 > theta2) then
                    f_t = (L_1 * F_sun * (1-Alb) * gamma1) / (r1**2) + &
                          (1-f) * (L_2 * F_sun * (1-Alb) * gamma2) / (r2**2)
                else
                    f_t = (L_2 * F_sun * (1-Alb) * gamma2) / (r2**2) + &
                          (1-f) * (L_1 * F_sun * (1-Alb) * gamma1) / (r1**2)
                end if
                
            case(2)  ! Zona de umbra
                if (theta1 > theta2) then
                    f_t = (L_1 * F_sun * (1-Alb) * gamma1) / (r1**2)
                else
                    f_t = (L_2 * F_sun * (1-Alb) * gamma2) / (r2**2)
                end if
                
            case default
                f_t = 0.0_dp
        end select
        
    end function calcular_flujo

    function calcular_temperatura(flujo_en_p, sigma) result(temp)
        real(dp), intent(in) :: flujo_en_p, sigma
        real(dp) :: temp
        temp = (flujo_en_p / sigma)**(0.25_dp)
    end function calcular_temperatura

end program jcb