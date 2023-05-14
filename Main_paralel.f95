program levialdi
    
    use mpi_f08
    implicit none

    integer :: my_rank, n_ranks, rank
    integer :: i, j, n_filas, n_col, pasos, par_fila, par_col, k, f, s, q
    integer :: ib, ie, jb, je, root, fila, col,north_rank,south_rank,west_rank,east_rank
    logical,dimension(:,:),allocatable:: matriz, matriz_nueva, matriz_temp, matriz_completa, matriz_falsos
    integer,dimension(:,:),allocatable:: matriz_res, matriz_completa_res
    logical,dimension(:,:,:),allocatable::map_history
    logical::iguales
    integer, dimension(:,:,:), allocatable :: resultado
    type(MPI_Datatype) :: a_short_row, a_long_row, a_col
    type(MPI_Status) :: status

    call MPI_Init()
    call MPI_Comm_rank( MPI_COMM_WORLD, my_rank )
    call MPI_Comm_size( MPI_COMM_WORLD, n_ranks )
    
    root=0
    if (my_rank == root) then
        read *, n_filas, n_col, par_fila, par_col
    end if

    call MPI_Bcast(n_filas, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(n_col, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(par_fila, 1, MPI_INTEGER, root, MPI_COMM_WORLD)
    call MPI_Bcast(par_col, 1, MPI_INTEGER, root, MPI_COMM_WORLD)

    call get_coords(my_rank, par_fila, par_col, fila, col )
    north_rank = get_rank(fila - 1, col, par_fila, par_col)
    south_rank = get_rank(fila + 1, col, par_fila, par_col)
    west_rank  = get_rank(fila, col - 1, par_fila, par_col)
    east_rank  = get_rank(fila, col + 1, par_fila, par_col)

    call partition( fila, par_fila, n_filas, ib, ie )
    call partition( col, par_col, n_col, jb, je )

    allocate(matriz(ib - 1:ie + 1, jb - 1:je + 1))
    allocate(matriz_nueva(ib - 1:ie + 1, jb - 1:je + 1))
    allocate(matriz_temp(ib - 1:ie + 1, jb - 1:je + 1))
    allocate(matriz_falsos(ib - 1:ie + 1, jb - 1:je + 1))
    allocate(matriz_res(ib - 1:ie + 1, jb - 1:je + 1))
    matriz_falsos=.False.
    matriz_res=0

    allocate(matriz_completa(n_filas+2,n_col+2))
    allocate(matriz_completa_res(n_filas+2,n_col+2))
    allocate(map_history(300, n_filas+2, n_col+2))
    allocate(resultado(300, n_filas+2, n_col+2))

    call MPI_Type_contiguous( ie - ib + 1, MPI_LOGICAL, a_col )
    call MPI_Type_commit(a_col)
    block
        type(MPI_Datatype) :: a_tmp_row, another_tmp_row
        integer(kind=MPI_ADDRESS_KIND) :: lb, real_extent

        call MPI_Type_vector(je - jb + 3, 1, ie - ib + 3, MPI_LOGICAL, a_tmp_row)
        call MPI_Type_get_extent( MPI_LOGICAL, lb, real_extent )
        call MPI_Type_create_resized( a_tmp_row, lb, real_extent, a_long_row )
        call MPI_Type_commit( a_long_row )

        call MPI_Type_vector(je - jb + 1, 1, ie - ib + 3, MPI_LOGICAL, another_tmp_row)
        call MPI_Type_get_extent( MPI_LOGICAL, lb, real_extent )
        call MPI_Type_create_resized( another_tmp_row, lb, real_extent, a_short_row )
        call MPI_Type_commit( a_short_row )
    end block

    call read_map(matriz, n_filas, n_col)
    call add_borders(matriz)
    call update_borders( matriz )

    do k=1,1000
        call save_map( matriz, n_filas, n_col, matriz_completa)
        call save_all_maps(matriz_completa,map_history,k)
        call reduccion(matriz, matriz_nueva, matriz_temp, ie, je, ib, jb)
        call add_borders(matriz_nueva)
        call update_borders(matriz_nueva)
        iguales=world_is_still(matriz_nueva,matriz_falsos)
        call MPI_Bcast(iguales, 1, MPI_LOGICAL, root, MPI_COMM_WORLD)
        if (iguales) exit
        matriz_temp = matriz  
        matriz = matriz_nueva  
        matriz_nueva = matriz_temp
    end do
    write(*,*)"Número de pasos:",k


    f=1
    s=0
    do q = k, 1, -1
        call save_map_res(matriz_res, n_filas, n_col, matriz_completa_res)
        call save_all_maps_res(matriz_completa_res,resultado,q)
        call aumento( map_history, resultado, matriz_res, ie, je, ib, jb, q, f, s)
        call add_borders_res(matriz_res)
        call update_borders_res(matriz_res)
    end do
    write(*,*) "Número de componentes:", f-1

    call print_map(matriz_res, n_filas, n_col)

    deallocate(matriz)
    deallocate(matriz_nueva)
    deallocate(matriz_temp)
    deallocate(matriz_completa)
    deallocate(matriz_falsos)
    deallocate(resultado)
    deallocate(map_history)
    deallocate(matriz_res)
    deallocate(matriz_completa_res)
    
    call MPI_Type_free(a_col)
    call MPI_Type_free(a_short_row)
    call MPI_Type_free(a_long_row)
    call MPI_Finalize()

contains

    subroutine reduccion(matriz, matriz_nueva, matriz_temp, w, h, l, m)
        logical, dimension(:, :), intent(inout) :: matriz_nueva, matriz, matriz_temp
        integer, intent(in) :: w, h, l, m

        do i=l, w
            do j=m, h
                if ((matriz(i-1,j) .eqv. .False.) .and. (matriz(i,j-1) .eqv. .False.) .and. (matriz(i-1,j-1) .eqv. .False.)) then
                    matriz_nueva(i,j)  = .False.
                else
                    matriz_nueva(i,j) = matriz(i,j)
                end if
                if ((matriz(i-1,j) .eqv. .True.) .and. (matriz(i,j-1) .eqv. .True.)) then
                    matriz_nueva(i,j) = .True.
                end if
            end do
        end do
        
    end subroutine reduccion
    
    subroutine aumento( map, resultado, matriz, w, h, l, m, q, f, s)
        logical, dimension(:,:,:), intent(in) :: map
        integer, dimension(:,:,:), intent(in) :: resultado
        integer, dimension(:,:), intent(inout) :: matriz
        integer, intent(in) :: w, h, l, m, q
        integer,intent(inout) :: f, s
        integer :: i, j
            do i=l, w
                do j=m, h
                    if (.not. map(q,i,j) .eqv. map(q+1,i,j)) then 
                        if ( map(q,i,j) .eqv. .True.) then
                            if(map(q,i-1,j) .eqv. .True. .and. resultado(q+1,i-1,j)/=0) then
                                s=resultado(q+1,i-1,j)
                            end if
                            if(map(q,i,j-1) .eqv. .True. .and. resultado(q+1,i,j-1)/=0) then
                                s=resultado(q+1,i,j-1)
                            end if
                            if(map(q,i-1,j-1) .eqv. .True. .and. resultado(q+1,i-1,j-1)/=0) then
                                s=resultado(q+1,i-1,j-1)
                            end if
                            if(map(q,i+1,j) .eqv. .True. .and. resultado(q+1,i+1,j)/=0) then
                                s=resultado(q+1,i+1,j)
                            end if
                            if(map(q,i,j+1) .eqv. .True. .and. resultado(q+1,i,j+1)/=0) then
                                s=resultado(q+1,i,j+1)
                            end if
                            if(map(q,i+1,j+1) .eqv. .True. .and. resultado(q+1,i+1,j+1)/=0) then
                                s=resultado(q+1,i+1,j+1)
                            end if
                            if(map(q,i+1,j-1) .eqv. .True. .and. resultado(q+1,i+1,j-1)/=0) then
                                s=resultado(q+1,i+1,j-1)
                            end if
                            if(map(q,i-1,j+1) .eqv. .True. .and. resultado(q+1,i-1,j+1)/=0) then
                                s=resultado(q+1,i-1,j+1)
                            end if
                            matriz(i,j)=s
                            if( map(q,i-1,j+1) .eqv. .False. .and. map(q,i+1,j-1) .eqv. .False. ) then
                                if(map(q,i-1,j) .eqv. .False. .and. map(q,i,j-1) .eqv. .False. ) then
                                    if( map(q,i,j+1) .eqv. .False. .and. map(q,i+1,j) .eqv. .False.) then
                                        if(map(q,i+1,j+1) .eqv. .False. .and. map(q,i-1,j-1) .eqv. .False.) then
                                            matriz(i,j)=f
                                            f=f+1
                                        end if
                                    end if
                                end if
                            end if
                        end if
                        if ( map(q,i,j) .eqv. .False.) then
                            matriz(i,j)=0
                        end if
                    end if
                    if (map(q,i,j) .eqv. map(q+1,i,j)) then  
                        matriz(i,j)=resultado(q+1,i,j)
                    end if
                end do
            end do
    end subroutine aumento
   

    subroutine read_map( map, h, w )
        logical, dimension(:, :), intent(inout) :: map
        integer, intent(in) :: h, w
        character(len=:), allocatable :: linea
        logical, allocatable :: temp(:)

        integer :: i, j, rb, re, cb, ce
       
        if (my_rank == root) then
            block
                integer :: current_row 
                integer :: current_col
                integer :: dst
                allocate(character(len=w) :: linea)
                do current_row = 0, par_fila - 1
                    call partition(current_row, par_fila, h, rb, re)
                    do i = rb, re
                        read *, linea
                        do current_col = 0, par_col - 1
                            call partition(current_col, par_col, w, cb, ce)
                            dst = get_rank(current_row, current_col, par_fila, par_col)
                            allocate(temp(ce - cb + 1))
                            do j = cb, ce 
                                select case (linea(j:j))
                                case ('X')
                                    temp(j - cb + 1) = .True.
                                case ('.')
                                    temp(j - cb + 1) = .False.
                                case default
                                    stop "read_map: wrong input character `" // linea(j:j) // "`"
                                end select
                            end do
                            if (dst == root) then
                                map(i, cb : ce)  = temp
                            else
                                call MPI_Send( temp, ce - cb + 1, MPI_LOGICAL, dst, 0,  MPI_COMM_WORLD )
                            end if
                            if (allocated( temp )) deallocate(temp)
                        end do
                    end do
                end do
                if (allocated( linea )) deallocate(linea)
            end block
        else
            do i = ib, ie
                call MPI_Recv(map(i,jb), 1, a_short_row, root, 0, MPI_COMM_WORLD, status)
            end do
        end if
    end subroutine read_map

    subroutine save_map( map, h, w, big_map)
        logical, dimension(:, :), intent(in) :: map
        logical, dimension(:, :), intent(inout) :: big_map
        integer, intent(in) :: h, w
        logical, allocatable :: temp(:)

        integer :: i, j, rb, re, cb, ce
       
        if (my_rank == root) then
            block
                integer :: current_row 
                integer :: current_col
                integer :: dst
                logical,dimension(:),allocatable :: linea
                allocate(linea(1:w))
                do current_row = 0, par_fila - 1
                    call partition(current_row, par_fila, h, rb, re)
                    do i = rb, re
                        do current_col = 0, par_col - 1
                            call partition(current_col, par_col, w, cb, ce)
                            dst = get_rank(current_row, current_col, par_fila, par_col)
                            allocate(temp(ce - cb + 1))
                            if (dst == root) then
                                do j = cb, ce
                                    linea(j:j) = map(i,j)
                                end do
                            else
                                call MPI_Recv( temp, ce - cb + 1, MPI_LOGICAL, dst, 0,  MPI_COMM_WORLD, status )
                                do j = cb, ce
                                    linea(j:j) = temp(j-cb+1)
                                end do
                            end if
                            if (allocated( temp )) deallocate(temp)
                        end do
                        big_map(i,:) = linea
                    end do
                end do
                if (allocated( linea )) deallocate(linea)
            end block
        else
            do i = ib, ie
                call MPI_Send(map(i,jb), 1, a_short_row, root, 0, MPI_COMM_WORLD)
            end do
        end if
    end subroutine save_map


    subroutine save_map_res( map, h, w, big_map)
        integer, dimension(:, :), intent(in) :: map
        integer, dimension(:, :), intent(inout) :: big_map
        integer, intent(in) :: h, w
        integer, allocatable :: temp(:)

        integer :: i, j, rb, re, cb, ce
       
        if (my_rank == root) then
            block
                integer :: current_row 
                integer :: current_col
                integer :: dst
                integer,dimension(:),allocatable :: linea
                allocate(linea(1:w))
                do current_row = 0, par_fila - 1
                    call partition(current_row, par_fila, h, rb, re)
                    do i = rb, re
                        do current_col = 0, par_col - 1
                            call partition(current_col, par_col, w, cb, ce)
                            dst = get_rank(current_row, current_col, par_fila, par_col)
                            allocate(temp(ce - cb + 1))
                            if (dst == root) then
                                do j = cb, ce
                                    linea(j:j) = map(i,j)
                                end do
                            else
                                call MPI_Recv( temp, ce - cb + 1, MPI_INTEGER, dst, 0,  MPI_COMM_WORLD, status )
                                do j = cb, ce
                                    linea(j:j) = temp(j-cb+1)
                                end do
                            end if
                            if (allocated( temp )) deallocate(temp)
                        end do
                        big_map(i,:) = linea
                    end do
                end do
                if (allocated( linea )) deallocate(linea)
            end block
        else
            do i = ib, ie
                call MPI_Send(map(i,jb), 1, a_short_row, root, 0, MPI_COMM_WORLD)
            end do
        end if
    end subroutine save_map_res


    subroutine partition( id, n_ids, size, b, e )
        integer, intent(in)    :: id, n_ids, size
        integer, intent(inout) :: b, e
        integer :: remainder, quotient

        remainder = modulo( size, n_ids )
        quotient  = (size - remainder) / n_ids
        b = 1 + quotient * (id    ) + min( remainder, id     )
        e =     quotient * (id + 1) + min( remainder, id + 1 )
    end subroutine partition


    subroutine get_coords( rank, n_rows, n_cols, row, col )
        integer, intent(in)    :: rank, n_rows, n_cols
        integer, intent(inout) :: row, col

        row = modulo(rank, n_rows)
        col = (rank - row) / n_rows
        if (0 <= col .and. col < n_cols) then
            return
        else
            print "(a, 2(i0, a))", "get_coords: rank ", rank, &
                " is outside the column range [0, ", n_cols, ")."
            call MPI_Abort( MPI_COMM_WORLD, MPI_ERR_TOPOLOGY )
        end if
    end subroutine get_coords

    integer function get_rank( row, col, n_rows, n_cols )
        integer, intent(in) :: row, col, n_rows, n_cols

        if (      0 <= col .and. col < n_cols &
            .and. 0 <= row .and. row < n_rows) then
                get_rank = row + col * n_rows
        else
            get_rank = MPI_PROC_NULL
        end if
    end function get_rank

    subroutine add_borders(map)

        logical, dimension(:, :), intent(inout) :: map

        map(ib:ie, jb - 1) = .False.
        map(ib:ie, je + 1) = .False.
        map(ib - 1, jb - 1:je + 1) = .False.
        map(ie + 1, jb - 1:je + 1) = .False.

    end subroutine add_borders

    subroutine add_borders_res(map)

        integer, dimension(:, :), intent(inout) :: map

        map(ib:ie, jb - 1) = 0
        map(ib:ie, je + 1) = 0
        map(ib - 1, jb - 1:je + 1) = 0
        map(ie + 1, jb - 1:je + 1) = 0

    end subroutine add_borders_res

    subroutine update_borders( map )
        logical, dimension(:, :), intent(inout) :: map

        call MPI_Sendrecv( map(ib , jb),     1, a_col, west_rank, 3, &
                           map(ib , je + 1), 1, a_col, east_rank, 3, MPI_COMM_WORLD, status)
        call MPI_Sendrecv( map(ib , je),     1, a_col, east_rank, 4, &
                           map(ib , jb - 1), 1, a_col, west_rank, 4, MPI_COMM_WORLD, status)

        call MPI_Sendrecv( map(ib,jb-1),      1, a_long_row, north_rank, 1, &
                           map(ie + 1,jb-1),  1, a_long_row, south_rank, 1, MPI_COMM_WORLD, status)
        call MPI_Sendrecv( map(ie,jb-1),      1, a_long_row, south_rank, 2, &
                           map(ib-1, jb-1),   1, a_long_row, north_rank, 2, MPI_COMM_WORLD, status)

    end subroutine update_borders

    subroutine update_borders_res( map )
        integer, dimension(:, :), intent(inout) :: map

        call MPI_Sendrecv( map(ib , jb),     1, a_col, west_rank, 3, &
                           map(ib , je + 1), 1, a_col, east_rank, 3, MPI_COMM_WORLD, status)
        call MPI_Sendrecv( map(ib , je),     1, a_col, east_rank, 4, &
                           map(ib , jb - 1), 1, a_col, west_rank, 4, MPI_COMM_WORLD, status)

        call MPI_Sendrecv( map(ib,jb-1),      1, a_long_row, north_rank, 1, &
                           map(ie + 1,jb-1),  1, a_long_row, south_rank, 1, MPI_COMM_WORLD, status)
        call MPI_Sendrecv( map(ie,jb-1),      1, a_long_row, south_rank, 2, &
                           map(ib-1, jb-1),   1, a_long_row, north_rank, 2, MPI_COMM_WORLD, status)

    end subroutine update_borders_res

    subroutine save_all_maps(map,all_maps,k)

        logical, dimension(:, :), intent(inout) :: map
        logical, dimension(:, :, :), intent(inout) :: all_maps
        integer, intent(in) :: k

        map(2:n_filas+1, 2:n_col+1) = map
        map(:,1) = .False.
        map(:,n_col+2) = .False.
        map(1,:) = .False.
        map(n_filas+2,:) = .False.
        all_maps(k, :, :) = map

    end subroutine save_all_maps

    subroutine save_all_maps_res(map,all_maps,k)

        integer, dimension(:, :), intent(inout) :: map
        integer, dimension(:, :, :), intent(inout) :: all_maps
        integer, intent(in) :: k

        map(2:n_filas+1, 2:n_col+1) = map
        map(:,1) = 0
        map(:,n_col+2) = 0
        map(1,:) = 0
        map(n_filas+2,:) = 0
        all_maps(k, :, :) = map

    end subroutine save_all_maps_res

    logical function world_is_still( map, falsos )
        logical, dimension(:, :), intent(in) :: map, falsos
        logical:: comp
        logical, dimension(n_ranks) :: all_is_still
    
        comp = all( map .eqv. falsos )
        call MPI_Gather(comp, 1, MPI_LOGICAL, &
                        all_is_still,   1, MPI_LOGICAL, root, MPI_COMM_WORLD)
        if (my_rank == root) then
            comp = all ( all_is_still .eqv. .True.)
        end if

    end function world_is_still
    
    subroutine print_map( map, h, w)
        integer, dimension(:, :), intent(in) :: map
        integer, intent(in) :: h, w
        integer, allocatable :: temp(:)

        integer :: i, j, rb, re, cb, ce
       
        if (my_rank == root) then
            block
                integer :: current_row 
                integer :: current_col
                integer :: dst
                integer,dimension(:),allocatable :: linea
                allocate(linea(1:w))
                do current_row = 0, par_fila - 1
                    call partition(current_row, par_fila, h, rb, re)
                    do i = rb, re
                        do current_col = 0, par_col - 1
                            call partition(current_col, par_col, w, cb, ce)
                            dst = get_rank(current_row, current_col, par_fila, par_col)
                            allocate(temp(ce - cb + 1))
                            if (dst == root) then
                                do j = cb, ce
                                    linea(j:j) = map(i,j)
                                end do
                            else
                                call MPI_Recv( temp, ce - cb + 1, MPI_INTEGER, dst, 0,  MPI_COMM_WORLD, status )
                                do j = cb, ce
                                    linea(j:j) = temp(j-cb+1)
                                end do
                            end if
                            if (allocated( temp )) deallocate(temp)
                        end do
                        print "(I4)", linea
                    end do
                end do
                if (allocated( linea )) deallocate(linea)
            end block
        else
            do i = ib, ie
                call MPI_Send(map(i,jb), 1, a_short_row, root, 0, MPI_COMM_WORLD)
            end do
        end if
    end subroutine print_map

end program levialdi
