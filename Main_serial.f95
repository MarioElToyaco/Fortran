program levialdi
    
    implicit none
    integer :: i, j, n_filas, n_col, pasos
    character(len=1),dimension(:,:),allocatable::matriz_sin, matriz, matriz_nueva, matriz_temp,matriz_puntos
    character(len=1),dimension(:,:,:),allocatable::map_history
    character(len=50) :: nombre_archivo
    logical::iguales
    integer, dimension(:,:,:), allocatable :: resultado
    
    write(*,*) 'Nombre del archivo (con extensión):'
    read(*,*) nombre_archivo
    
    open(unit=10, file=nombre_archivo, status='old', action='read')
    read(10,*) n_filas, n_col
    close(10)

    call asignar_matriz(matriz_sin,n_filas,n_col)
    call asignar_matriz(matriz,n_filas+2,n_col+2)
    call asignar_matriz(matriz_nueva,n_filas+2,n_col+2)
    call asignar_matriz(matriz_temp,n_filas+2,n_col+2)
    call asignar_matriz(matriz_puntos,n_filas+2,n_col+2)
    call asignar_mapa(map_history,300, n_filas+2, n_col+2)
    call asignar_resultado(resultado,300, n_filas+2, n_col+2)
    call read_map( matriz_sin, n_filas, n_col )

    matriz(2:n_filas+1, 2:n_col+1) = matriz_sin
    matriz(:,1) = '.'
    matriz(:,n_col+2) = '.'
    matriz(1,:) = '.'
    matriz(n_filas+2,:) = '.'
    
    call reduccion( matriz, matriz_nueva,map_history, n_filas, n_col, pasos,matriz_temp)

    call aumento( map_history, resultado, n_filas, n_col, pasos )

    do i = 2, n_filas+1
        do j = 2, n_col+1
            write(*, "(I4)", advance='no') resultado(1,i,j)
        end do
        write(*,*)
    end do
    
    deallocate(matriz)
    deallocate(matriz_sin)
    deallocate(matriz_nueva)
    deallocate(matriz_temp)
    deallocate(matriz_puntos)
    deallocate(resultado)
    deallocate(map_history)
    

contains

    subroutine asignar_matriz(matriz, filas, columnas)

        character(len=1), dimension(:,:), allocatable, intent(out) :: matriz
        integer, intent(in) :: filas, columnas
   
        allocate(matriz(filas, columnas))
        matriz = '.'
   
    end subroutine
    
    subroutine asignar_mapa(matriz,pasos, filas, columnas)

        character(len=1), dimension(:,:,:), allocatable, intent(out) :: matriz
        integer, intent(in) ::pasos, filas, columnas
   
        allocate(matriz(pasos,filas, columnas))
        matriz = '.'
   
    end subroutine
    
    subroutine asignar_resultado(matriz,pasos, filas, columnas)

        integer, dimension(:,:,:), allocatable, intent(out) :: matriz
        integer, intent(in) ::pasos, filas, columnas
   
        allocate(matriz(pasos,filas, columnas))
        matriz = 0
   
    end subroutine

    subroutine reduccion(matriz, matriz_nueva, map, w, h, pasos,matriz_temp)
        character(len=1), dimension(:, :), intent(inout) :: matriz_nueva, matriz,matriz_temp
        character(len=1), dimension(:,:,:), intent(inout) :: map
        integer, intent(in) :: w, h
        integer, intent(inout) :: pasos
        integer :: i, j, k
        matriz_nueva=matriz
        do k=1,10000000
            map(k, :, :) = matriz_nueva
            do i = 1, w+2
                do j = 1, h+2
                    if ( matriz(i-1,j) == '.' .and. matriz(i,j-1) == '.' .and. matriz(i-1,j-1) == '.') then
                        matriz_nueva(i,j)  = '.'
                    else
                        matriz_nueva(i,j) = matriz(i,j)
                    end if
                    if ( matriz(i-1,j) == 'X' .and. matriz(i,j-1) == 'X') then
                        matriz_nueva(i,j) = 'X'
                    end if
                end do
            end do
            matriz_temp = matriz  
            matriz = matriz_nueva  
            matriz_nueva = matriz_temp
            iguales = all(matriz_nueva == matriz_puntos)
            if(iguales) exit
        end do
        pasos=k
        write(*,*)"Número de pasos:",pasos
    end subroutine reduccion
    
    subroutine aumento( map, resultado, w, h, t )
        character(len=1), dimension(:,:,:), intent(in) :: map
        integer, dimension(:,:,:), intent(inout) :: resultado
        integer, intent(in) :: h, w, t
        integer :: i, j, q, l, s
        l=1
        s=0
        
        do q = t, 1, -1
            do i = 1, w+2
                do j = 1, h+2
                    if (.not. map(q,i,j)==map(q+1,i,j)) then 
                        if ( map(q,i,j)=='X') then
                            if(map(q,i-1,j)=='X' .and. resultado(q+1,i-1,j)/=0) then
                                s=resultado(q+1,i-1,j)
                            end if
                            if(map(q,i,j-1)=='X' .and. resultado(q+1,i,j-1)/=0) then
                                s=resultado(q+1,i,j-1)
                            end if
                            if(map(q,i-1,j-1)=='X' .and. resultado(q+1,i-1,j-1)/=0) then
                                s=resultado(q+1,i-1,j-1)
                            end if
                            if(map(q,i+1,j)=='X' .and. resultado(q+1,i+1,j)/=0) then
                                s=resultado(q+1,i+1,j)
                            end if
                            if(map(q,i,j+1)=='X' .and. resultado(q+1,i,j+1)/=0) then
                                s=resultado(q+1,i,j+1)
                            end if
                            if(map(q,i+1,j+1)=='X' .and. resultado(q+1,i+1,j+1)/=0) then
                                s=resultado(q+1,i+1,j+1)
                            end if
                            if(map(q,i+1,j-1)=='X' .and. resultado(q+1,i+1,j-1)/=0) then
                                s=resultado(q+1,i+1,j-1)
                            end if
                            if(map(q,i-1,j+1)=='X' .and. resultado(q+1,i-1,j+1)/=0) then
                                s=resultado(q+1,i-1,j+1)
                            end if
                            resultado(q,i,j)=s
                            if( map(q,i-1,j+1)=='.' .and. map(q,i+1,j-1)=='.' ) then
                                if(map(q,i-1,j)=='.' .and. map(q,i,j-1)=='.' ) then
                                    if( map(q,i,j+1)=='.' .and. map(q,i+1,j)=='.') then
                                        if(map(q,i+1,j+1)=='.' .and. map(q,i-1,j-1)=='.') then
                                            resultado(q,i,j)=l
                                            l=l+1
                                        end if
                                    end if
                                end if
                            end if
                        end if
                        if ( map(q,i,j)=='.') then
                            resultado(q,i,j)=0
                        end if
                    end if
                    if (map(q,i,j)==map(q+1,i,j)) then  
                        resultado(q,i,j)=resultado(q+1,i,j)
                    end if
                end do
            end do
        end do
        write(*,*) "Número de componentes:", l-1
    end subroutine aumento
    
    
    subroutine read_map( map, h, w )
        character(len=1), dimension(:, :), intent(inout) :: map
        integer, intent(in) :: h, w
        character(len=:), allocatable :: linea
        integer :: i, j
        
        allocate(character(len=w) :: linea)
        open(unit=10, file=nombre_archivo, status='old', action='read')
        read(10,*) n_filas, n_col

        do i = 1, n_filas
            read(10, '(A)') linea
            do j = 1, n_col
                matriz_sin(i, j) = linea(j:j)
            end do
        end do
        if (allocated( linea )) deallocate(linea)
        close(10)
    end subroutine read_map
    
end program levialdi
