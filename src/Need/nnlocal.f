      program nnlocal
c--- This is an entry point into nnlocal (usually called by nnlocal program)    
      implicit none
      character*72 inputfile,workdir
      call determinefilenames(inputfile,workdir)
      call nnlocal_main(inputfile,workdir)
      end
