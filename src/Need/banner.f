      block data codeversion_data
      implicit none
      include 'codeversion.f'
      data codeversion/'beta'/
      end

      subroutine banner
      implicit none
      include 'codeversion.f'
      write(6,*) '************************** version  '//codeversion// '**************************'
      write(6,*) '*                                                                  *'
      write(6,*) '* ##    ## ##    ## ##        #######   ######     ###    ##       *'
      write(6,*) '* ###   ## ###   ## ##       ##     ## ##    ##   ## ##   ##       *'
      write(6,*) '* ####  ## ####  ## ##       ##     ## ##        ##   ##  ##       *'
      write(6,*) '* ## ## ## ## ## ## ##       ##     ## ##       ##     ## ##       *'
      write(6,*) '* ##  #### ##  #### ##       ##     ## ##       ######### ##       *'
      write(6,*) '* ##   ### ##   ### ##       ##     ## ##    ## ##     ## ##       *'
      write(6,*) '* ##    ## ##    ## ########  #######   ######  ##     ## ######## *'
      write(6,*) '*                                                                  *'
      write(6,*) '*********************** December 20th, 2024 ************************'
      write(6,*) '*                                                                  *'
      write(6,*) '* Authors:                                                         *'
      write(6,*) '* Vittorio Del Duca <Vittorio.Del.Duca@cern.ch>                    *'
      write(6,*) '* Claude Duhr <cduhr@uni-bonn.de>                                  *'
      write(6,*) '* Levente Fekeshazy <levente.fekeshazy@desy.de>                    *'
      write(6,*) '* Flavio Guadagni <guadagni.flavio@gmail.com>                      *'
      write(6,*) '* Pooja Mukherjee <pooja.mukherjee@desy.de>                        *'
      write(6,*) '* Gabor Somogyi <somogyi.gabor@wigner.hun-ren.hu>                  *'
      write(6,*) '* Sam Van Thurenhout <sam.van.thurenhout@wigner.hun-ren.hu>        *'
      write(6,*) '* Francesco Tramontano <francesco.tramontano@unina.it>             *'
      write(6,*) '*                                                                  *'
      write(6,*) '* https://github.com/nnlocal/nnlocal.git                           *'
      write(6,*) '*                                                                  *'
      write(6,*) '********************************************************************'
      return
      end








