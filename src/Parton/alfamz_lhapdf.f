      double precision function alphas(q,amz,nloop)
      include 'blha.f'
      double precision alphasmz,alphasPDF
      double precision q, amz
      integer nloop

      if (useblha.ne.0) then
         alphas = alphasPDF(q)
      else
         alphas = alphasmz(q,amz,nloop)
      endif
      return
      end
      
      
