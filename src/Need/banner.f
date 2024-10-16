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
      write(6,*) '*********************** September 10th, 2024 ***********************'
      write(6,*) '*                                                                  *'
      write(6,*) '* Authors:                                                         *'
      write(6,*) '* Vittorio Del Duca <Vittorio.Del.Duca@cern.ch>                    *'
      write(6,*) '* Claude Duhr <cduhr@uni-bonn.de>                                  *'
      write(6,*) '* Flavio Guadagni <guadagni.flavio@gmail.com>                      *'
      write(6,*) '* Pooja Mukherjee <pmukherj@uni-bonn.de>                           *'
      write(6,*) '* Gabor Somogyi <gabor.somogyi@cern.ch>                            *'
      write(6,*) '* Sam Van Thurenhout Sam <sam.van.thurenhout@wigner.hu>            *'
      write(6,*) '* Francesco Tramontano <francesco.tramontano@unina.it>             *'
      write(6,*) '*                                                                  *'
      write(6,*) '* https://github.com/nnlocal/nnlocal.git                           *'
      write(6,*) '*                                                                  *'
      write(6,*) '********************************************************************'
      return
      end








