      program debug_program
      use debugging
      implicit none
      logical :: debug,printer
      call init_debug(debug,printer) !Initialize debug, printer
      open(10,file="fortran/temp/debug.txt",status='unknown')
      write(10,*) debug
      close(10)
      open(10,file="fortran/temp/printer.txt",status='unknown')
      write(10,*) printer
      close(10)
      end
