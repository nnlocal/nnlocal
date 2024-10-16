      character*72 function checkpath(inpstr)
c--- take an input string (path to PDF file) and strip the leading
c--- directory if the parameter "vanillafiles" is true
      implicit none
      character inpstr*(*)
      checkpath=inpstr
      return
      end
      
