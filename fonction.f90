module fonction
    implicit none
    integer, parameter :: PR = 8
    real(PR), parameter :: g = 9.81
  contains
    ! Fonction qui calcule la somme de deux nombres
    subroutine initialisation(x, h, u, v, a, b, cas)
      real(PR), intent(inout) :: x, h, u, v, a, b
      integer, intent(in) :: cas
      if ( cas == 0 ) then
        h = 1
        u = 0.2_PR
        v = 0.7_PR
        a = 0.5_PR
        b = 0.4_PR
      end if
      if ( cas == 1 ) then
        if (x <= 0.5) then
            h = 1
            u = 0.2_PR
            v = 0.7_PR
            a = 0.5_PR
            b = 0.4_PR
        else if (0.5< x .and. x <= 1) then
            h = 0.5_PR
            u = -0.1_PR
            v = 0.3_PR
            a = 1._PR
            b = 0.1_PR
        end if 
      else if ( cas == 2) then
        if (x <= 0.5) then
            h = 1.4_PR
            u = 0.2_PR
            v = 0.6_PR
            a = 1.0_PR
            b = 0.4_PR
        else if (0.5< x .and. x <= 1) then
            h = 0.2_PR
            u = -0.1_PR
            v = 0.3_PR
            a = 1.2_PR
            b = 0.1_PR
        end if
      else if ( cas == 3 ) then 
        if (x <= 0.5) then
            h = 2.0_PR
            u = 1.0_PR
            v = 2.5_PR
            a = 0.8_PR
            b = 0.4_PR
        else if (0.5< x .and. x <= 1) then
            h = 0.0_PR
            u = 0.0_PR
            v = 0.0_PR
            a = 0.0_PR
            b = 0.0_PR
        end if
      end if
    end subroutine initialisation

    subroutine majoration_valeur_propre()

    end subroutine majoration_valeur_propre
  
    subroutine flux_physique(fui, UN, i)
        integer, intent(in) :: i
        real(PR), dimension(5), intent(inout) :: fui
        real(PR), dimension(:,:), intent(in) :: UN

        fui(1) = UN(2,i)
        fui(2) = UN(2,i)**2/UN(1,i) + g*UN(1,i)**2*(1._PR/2)-UN(4,i)**2/UN(1,i)
        fui(3) = UN(2,i)*UN(3,i)/UN(1,i) - UN(4,i)*UN(5,i)/UN(1,i)
        fui(4) = 0
        fui(5) = UN(2,i)*UN(5,i)/UN(1,i) - UN(3,i)*UN(4,i)/UN(1,i)
    end subroutine flux_physique

    subroutine erreur(UN, file, err)
      ! Arguments
      real(PR), dimension(:,:), intent(in) :: UN      ! Solution actuelle
      character(len=*), intent(in) :: file             ! Nom du fichier de référence
      real(PR), intent(out) :: err                     ! Erreur globale (norme L2)
  
      ! Variables locales
      integer :: i, j, num_ref_points, N
      real(PR), dimension(:), allocatable :: ref_x
      real(PR), dimension(:,:), allocatable :: ref_U
      real(PR) :: delta_x, x, error_L2, diff1, diff2, diff3, diff4, diff5
  
      ! Taille de la solution actuelle
      N = size(UN, 2)
      delta_x = 1.0_PR / N  ! Espacement des points dans UN
  
      ! Lecture du fichier de référence
      OPEN(UNIT=20, FILE=file, STATUS="OLD", ACTION="READ")
      read(20, *) num_ref_points  ! Nombre de points dans le fichier de référence
      allocate(ref_x(num_ref_points))
      allocate(ref_U(5,num_ref_points))
      do i = 1, num_ref_points
          read(20, *) ref_x(i), ref_U(1,i), ref_U(2,i), ref_U(3,i), ref_U(4,i), ref_U(5,i)
      end do
      close(20)
 
      !Calcul des erreurs L2
      error_L2 = 0.0_PR
      do i=1, N
        do j = 1, num_ref_points-1
          x = ((2*i-1)*delta_x)/2._PR
          if ( x >= ref_x(j) .and. x < ref_x(j+1) ) then
            ! Calcul des différences individuelles
            diff1 = (UN(1, i) - ref_U(1, j))**2
            diff2 = (UN(2, i) / UN(1, i) - ref_U(2, j))**2
            diff3 = (UN(3, i) / UN(1, i) - ref_U(3, j))**2
            diff4 = (UN(5, i) / UN(1, i) - ref_U(4, j))**2
            diff5 = (UN(4, i) / UN(1, i) - ref_U(5, j))**2
            !print *, diff1, diff2, diff3, diff4, diff5
            ! Addition des contributions
            error_L2 = error_L2 + diff1 + diff2 + diff3 + diff4 + diff5

              !print *, error_L2 + (UN(k,i)-ref_U(k,j))**2 
          end if
        end do 
      end do

      error_L2 = 1._PR/N*sqrt(error_L2 / (5.0_PR * N))

      ! Retourner l'erreur comme norme L2
      err = error_L2
  
      ! Optionnel : imprimer l'erreur pour contrôle
      write(*, *) "Erreur L2 :", error_L2
  
      ! Libération de la mémoire
      deallocate(ref_x)
      deallocate(ref_U)
  
  end subroutine erreur
  



end module fonction
  