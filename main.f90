program main
    use fonction
    implicit none

    !déclaration des variables
    
    integer :: unit_num, N, i, j, k, cas
    real(PR) :: t, delta_x, delta_t, x, bord_gauche, bord_droit, t_real
    real(PR) :: h, u, v, a, b, h_moy, a_moy, u_moy, alpha, alpha_max, c, ca, err
    real(PR), allocatable,dimension(:) :: f_ui, f_ui1
    real(PR), allocatable, dimension(:,:) :: UN, UN1, F
    
    
    !instructions
    unit_num = 10

    ! Ouvrir le fichier de paramètres
    OPEN(UNIT=unit_num, FILE="param.txt", STATUS="OLD", ACTION="READ")
    ! Lire les valeurs du fichier
    READ(unit_num, *) N
    READ(unit_num, *) t
    READ(unit_num, *) bord_gauche
    READ(unit_num, *) bord_droit
    READ(unit_num, *) cas
    ! Fermer le fichier
    CLOSE(unit_num)
    alpha_max = 0
    !défition du pas en espace
    delta_x = (bord_droit-bord_gauche)/N
    !intitialisation des variables
    allocate(UN(5,N),UN1(5,N),F(5,N+1))
    !structure du code 
    OPEN(UNIT=unit_num+1, FILE="result.txt", STATUS="OLD", ACTION="WRITE")
    !initialisation de U dans toutes les mailles du domaine
    do i = 1, N
        ! calcul de x au centre d'une maille pour initialiser les vecteurs UN et F(U)
        x = ((2*i-1)*delta_x)/2._PR
        call initialisation(x, h, u, v, a, b, cas)
        UN(1,i) = h
        UN(2,i) = h*u
        UN(3,i) = h*v
        UN(4,i) = h*a
        UN(5,i) = h*b
    end do
    do i = 1, N-1
        ! !h_moy = (UN(1,i) + UN(1,i+1)) * (1._PR/2)
        ! h_moy = max(UN(1,i),UN(1,i+1))
        ! !u_moy = (UN(2,i)/UN(1,i) + UN(2,i+1)/UN(1,i+1)) * (1._PR/2)
        ! u_moy = max(UN(2,i)/UN(1,i),UN(2,i+1)/UN(1,i+1))
        ! !a_moy = (UN(4,i)/UN(1,i) + UN(4,i+1)/UN(1,i+1)) * (1._PR/2)
        ! a_moy = max(UN(4,i)/UN(1,i),UN(4,i+1)/UN(1,i+1))
        ! alpha = max(u_moy, u_moy - abs(a_moy), u_moy + abs(a_moy), &
        ! u_moy + sqrt(a_moy**2 + g*h_moy), &
        ! u_moy - sqrt(a_moy**2 + g*h_moy))
        !print *, "alpha : " , alpha, " arète : ", i


        !calcul des valeurs propres
        h_moy = max(UN(1,i),UN(1,i+1))
        u_moy = max(UN(2,i)/UN(1,i),UN(2,i+1)/UN(1,i+1))
        a_moy = max(UN(4,i)/UN(1,i),UN(4,i+1)/UN(1,i+1))
        c = h_moy*sqrt(a_moy**2+g*h_moy)
        ca = h*abs(a)
        alpha = max(u_moy - c/h_moy, u - ca/h_moy, u_moy, u_moy + ca/h_moy, u + c/h_moy)
        if (alpha >= alpha_max) then 
            alpha_max = alpha
            !print *, "alpha_max : ", alpha_max
        end if
    end do
    !boucle en temps
    t_real = 0.0_PR
    k = 0 
    do while(t_real < t .and. k < 100000000)
        !on impose un flux nul au bord
        F(:,1) = 0.0
        ! F(1,1) = 1
        ! F(2,1) = 0.2
        ! F(3,1) = 0.7
        ! F(4,1) = 0.5
        ! F(5,1) = 0.4
        F(:,N+1) = 0.0
        ! F(1,N+1) = 1
        ! F(2,N+1) = 0.2
        ! F(3,N+1) = 0.7
        ! F(4,N+1) = 0.5
        ! F(5,N+1) = 0.4
        !parcours de toutes les arêtes central
        do i = 1, N-1
            !calcul de flux phyique en Ui et Ui+1
            allocate(f_ui(5), f_ui1(5))
            call flux_physique(f_ui, UN, i)
            call flux_physique(f_ui1, UN, i+1)
            !calcul des valeur moyenne de u, a et h pour majorer la valeur propre
            ! !h_moy = (UN(1,i) + UN(1,i+1)) * (1._PR/2)
            ! h_moy = max(UN(1,i),UN(1,i+1))
            ! !u_moy = (UN(2,i)/UN(1,i) + UN(2,i+1)/UN(1,i+1)) * (1._PR/2)
            ! u_moy = max(UN(2,i)/UN(1,i),UN(2,i+1)/UN(1,i+1))
            ! !a_moy = (UN(4,i)/UN(1,i) + UN(4,i+1)/UN(1,i+1)) * (1._PR/2)
            ! a_moy = max(UN(4,i)/UN(1,i),UN(4,i+1)/UN(1,i+1))
            ! alpha = max(u_moy, u_moy - abs(a_moy), u_moy + abs(a_moy), &
            ! u_moy + sqrt(a_moy**2 + g*h_moy), &
            ! u_moy - sqrt(a_moy**2 + g*h_moy))
            !calcul des valeurs propres
            h_moy = max(UN(1,i),UN(1,i+1))
            u_moy = max(UN(2,i)/UN(1,i),UN(2,i+1)/UN(1,i+1))
            a_moy = max(UN(4,i)/UN(1,i),UN(4,i+1)/UN(1,i+1))
            c = h_moy*sqrt(a_moy**2+g*h_moy)
            ca = h*abs(a)
            alpha = max(u_moy - c/h_moy, u - ca/h_moy, u_moy, u_moy + ca/h_moy, u + c/h_moy)
            
            if (alpha >= alpha_max) then 
                alpha_max = alpha
                !print *, alpha_max
            end if
            !flux de Rusanov : 1/2*(F(UNi)+F(UNi+1)) - b/2*(UNi+1 - UNi)
            do j = 1,5
                F(j,1) = f_ui(j)
                F(j,i+1) = (1._PR/2._PR) * (f_ui(j) + f_ui1(j)) - alpha * (1._PR/2._PR) * (UN(j,i+1) - UN(j,i))
                F(j,N+1) = f_ui1(j)
            end do
            deallocate(f_ui,f_ui1) 
        end do
        !calcul du pas de temps
        delta_t = 0.5_PR*delta_x/(alpha_max)
        !Mise à jour de la contribution du flux dans chaque maille/ mise à jour de la solution
        do i = 2, N-1
            do j = 1, 5
                UN1(j,1) = UN(j,1)
                UN1(j,i) = UN(j,i) - delta_t/delta_x * (F(j,i+1) - F(j,i))
                UN1(j,N) = UN(j,N)
                UN(j,i) = UN1(j,i)
            end do  
        end do

        t_real = t_real +delta_t
        k = k + 1 
        !print *, t_real
    end do

    do i = 1, N
        x = ((2*i-1)*delta_x)/2._PR
        WRITE(unit_num+1, '(F8.3, 1X, F8.3, 1X, F8.3, 1X, F8.3, 1X, F8.3, 1X, F8.3)') &
        x, UN(1,i), UN(2,i)/UN(1,i), &
        UN(3,i)/UN(1,i), UN(5,i)/UN(1,i), &
        UN(4,i)/UN(1,i)
    end do 

    !calcul de l'erreur
    call erreur(UN, "solution_cas_1.txt", err)




    deallocate(UN,UN1)
    deallocate(F)
end program