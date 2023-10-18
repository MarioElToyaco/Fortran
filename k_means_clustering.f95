program k_means

    implicit none
  
    integer :: n_points, n_clusters, i, j, iter
    real :: distance, min_distance
    real, parameter :: tolerance = 1e-5
    real, dimension(:, :), allocatable :: points, centroids, centroids_old 
    integer, dimension(:), allocatable :: cluster_assignment
    integer, dimension(:,:), allocatable :: cluster_count
    character(len=100) :: filename, fileout, answer
    integer, parameter :: max_iter = 1000
    integer :: status

    print *, "Nombre del archivo (con extensión):"
    read *, filename
    n_points=0
    n_clusters=0

    call read_constants(filename, n_points, n_clusters)
    
    allocate(points(n_points,3))

    call read_points(filename, n_points, points)
   
    allocate(cluster_assignment(n_points))
    allocate(cluster_count(n_clusters,1))
    allocate(centroids(n_clusters,3))
    allocate(centroids_old(n_clusters,3))
    
    do i = 1, n_clusters
        centroids(i, :) = 0
    end do

    iter = 0
    do while (iter < max_iter)

        cluster_count(:,:)=0
        iter = iter + 1

        do i = 1, n_points
            min_distance = 1e20
            do j = 1, n_clusters
                distance = euclidean_distance(points(i, :), centroids(j, :))
                if (distance < min_distance) then
                    min_distance = distance
                    cluster_assignment(i) = j
                end if
            end do
            cluster_count(cluster_assignment(i),:)=cluster_count(cluster_assignment(i),:) + 1
        end do

        centroids_old = centroids

        call update_centroids(points, centroids, cluster_assignment, n_clusters, n_points)

        if (maxval(abs(centroids - centroids_old)) < tolerance) exit

    end do
    
    print *, "¿Quieres guardar el resultado en un archivo?: (yes or no)"
    read *, answer
    
    if (answer=="no") then
        call print_results(n_clusters, n_points, points, centroids, cluster_assignment, cluster_count)
    end if
    
    if (answer=="yes") then
        fileout = "results.out"
        call save_file(n_clusters, n_points, points, centroids, cluster_assignment, cluster_count, fileout)
    end if
    
    deallocate(points)
    deallocate(cluster_assignment)
    deallocate(cluster_count)
    deallocate(centroids)
    deallocate(centroids_old)    
    

contains

  subroutine read_constants(filename,n_points,n_clusters)

	  character(len=100), intent(in) :: filename
	  character(len=100) :: line
	  integer, intent(inout) :: n_points, n_clusters
	  integer :: status

        open(10, file=trim(filename), status='old', action='read', iostat=status)
        if (status /= 0) then
            print *, "Error al abrir el archivo"
            stop
        end if
        
        read(10, *, iostat=status) n_clusters

        do 
            read(10, '(A)', iostat=status) line
            if (status /= 0) exit
            n_points = n_points + 1
        end do
    
        close(10)

  end subroutine read_constants


  subroutine read_points(filename,n_points,points)

	  character(len=100), intent(in) :: filename
	  integer, intent(in) :: n_points
	  real, dimension(:, :), intent(inout) :: points
	  integer :: status, dummy, i

        open(10, file=trim(filename), status='old', action='read', iostat=status)
        read(10, *, iostat=status) dummy
        do i=1 , n_points
            read(10, *, iostat=status) points(i,:)
            if (status /= 0) exit
        end do
        close(10)

  end subroutine read_points


  function euclidean_distance(p1, p2) result(distance)

    real, dimension(:), intent(in) :: p1, p2
	  integer :: distance
	    
        distance = sqrt(sum((p1 - p2)**2))

  end function euclidean_distance

  subroutine update_centroids(points, centroids, cluster_assignment, n_clusters, n_points)
        real, dimension(:,:), intent(in) :: points
        real, dimension(:,:), intent(inout) :: centroids
        integer, dimension(:), intent(in) :: cluster_assignment
	      integer, intent(in) :: n_clusters, n_points
        integer :: i, j
        integer, dimension(n_clusters) :: cluster_number

        cluster_number = 0
        centroids = 0.0

        do i = 1, n_points
            j = cluster_assignment(i)
            centroids(j, :) = centroids(j, :) + points(i, :)
            cluster_number(j) = cluster_number(j) + 1
        end do

        do i = 1, n_clusters
            if (cluster_number(i) > 0) then
		           centroids(i, :) = centroids(i, :) / cluster_number(i)
	          end if
        end do

  end subroutine update_centroids

  subroutine print_results(n_clusters, n_points, points, centroids, cluster_assignment, cluster_count)

	  integer, intent(in) :: n_clusters, n_points
	  real, dimension(:, :), intent(in) :: points, centroids
    integer, dimension(:), intent(in) :: cluster_assignment
    integer, dimension(:,:), intent(in) :: cluster_count
	  integer :: i
	    
	  print *, "Clusters:", n_clusters
        do i = 1, n_clusters
            print *, i, cluster_count(i, :), int(centroids(i, :))
        end do
	    
        print *, "Puntos:", n_points
        do i = 1, n_points
            print *, int(points(i, :)), "  ", cluster_assignment(i)
        end do

  end subroutine print_results
    
  subroutine save_file(n_clusters, n_points, points, centroids, cluster_assignment, cluster_count, filename)
	    
	  integer, intent(in) :: n_clusters, n_points
	  real, dimension(:, :), intent(in) :: points, centroids
    integer, dimension(:), intent(in) :: cluster_assignment
    integer, dimension(:,:), intent(in) :: cluster_count
    character(len=100), intent(in) :: filename
	  integer :: i
	    
        open(unit=10, file=filename, status='replace')
    
        write(10,*) "Clusters:", n_clusters
        do i = 1, n_clusters
            write(10,*) i, cluster_count(i, :), int(centroids(i, :))
        end do
    
        write(10,*) "Puntos:", n_points
        do i = 1, n_points
            write(10,*) int(points(i, :)), "  ", cluster_assignment(i)
        end do
    
        close(10)
    
  end subroutine save_file
  
end program k_means
