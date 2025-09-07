program hash3_coarray
  use iso_fortran_env
  implicit none
  character(len=256) :: input_file, arg
  integer :: len, iperf, ios
  integer(8) :: n
  character(:), allocatable :: data_
  integer(1), allocatable :: data(:)
  integer(4), allocatable :: freq(:)[:]
  integer :: i, j, k, m, is_eq, ind, sh, sh1
  integer :: numimgs, thisimg
  integer(8) :: t1, t2, rate
  logical :: perf_collect

  call get_command_argument(1, input_file)
  call get_command_argument(2, arg)
  read(arg, *) len
  call get_command_argument(3, arg)
  read(arg, *) iperf
  perf_collect = (iperf /= 0)

  inquire(file=trim(input_file), size=n)
  if (n == -1) stop 'File size error'

  allocate(character(n) :: data_)
  open(10, file=trim(input_file), access='stream', form='unformatted', status='old', iostat=ios)
  if (ios /= 0) stop 'Open error'
  read(10, iostat=ios) data_
  if (ios /= 0) stop 'Read error'
  close(10)

  n = len_trim(data_)
  allocate(data(n))

  do i = 1, n
    select case(data_(i:i))
    case('A'); data(i) = 0
    case('C'); data(i) = 1
    case('G'); data(i) = 2
    case('T'); data(i) = 3
    case default; data(i) = -1 ! error
    end select
  end do
  deallocate(data_)

  m = n - len + 1
  if (m < 1) stop 'Invalid length'


  allocate(freq(1:m)[*])

  numimgs = num_images()
  thisimg = this_image()


  call system_clock(t1, rate)

  do i = thisimg, m, numimgs
    block
      integer :: shift(0:63)
      shift = len - 2
      do j = 2, len - 1
        ind = int(data(i + j - 2), 4) * 16 + int(data(i + j - 1), 4) * 4 + int(data(i + j), 4)
        shift(ind) = len - 1 - j
      end do
      ind = int(data(i + len - 3), 4) * 16 + int(data(i + len - 2), 4) * 4 + int(data(i + len - 1), 4)
      sh1 = shift(ind)
      shift(ind) = 0
      if (sh1 == 0) sh1 = 1

      freq(i) = 0
      j = len - 1
      do while (j < n)
        sh = 1
        do while (sh /= 0 .and. j < n)
          ind = int(data(j - 2), 4) * 16 + int(data(j - 1), 4) * 4 + int(data(j), 4)
          sh = shift(ind)
          j = j + sh
        end do
        if (j < n) then
          is_eq = 1
          do k = 1, len
            if (data(i + k - 1) /= data(j - len + 1 + k - 1)) then
              is_eq = 0
              exit
            end if
          end do
          freq(i) = freq(i) + is_eq
          j = j + sh1
        else
          exit
        end if
      end do
    end block
  end do

  sync all

  call system_clock(t2)

  if (perf_collect .and. thisimg == 1) then
    write(*, *) n, len, real(t2 - t1, kind=8) / real(rate, kind=8)
  end if

end program hash3_coarray