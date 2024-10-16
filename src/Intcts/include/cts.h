      real(ki) na2srsqq(0:2,0:2,-4:0), na2srsgg(0:2,0:2,-4:0),
     1         na2carssrsqq(0:2,0:2,-4:0), na2carssrsgg(0:2,0:2,-4:0),
     2         na2carbssrsgg(0:2,0:2,-4:0),
     3         na2carsggg(2,2,-4:0), na2carsqqg(2,2,-4:0), na2carsqqq(2,2,-4:0), na2carsggq(2,2,-4:0), 
     4         na2carsqqp(2,2,-4:0), na2carsqgg(2,2,-4:0), na2carsgqq(2,2,-4:0), na2carsqpp(2,2,-4:0),
     5         na2carbsqgqg(2,2,-4:0), na2carbsqggg(2,2,-4:0), na2carbsqggq(2,2,-4:0), na2carbsqgqq(2,2,-4:0),
     7         na2carbsggqg(2,2,-4:0), na2carbsgggg(2,2,-4:0), na2carbsgggq(2,2,-4:0), na2carbsggqq(2,2,-4:0),
     8         na2carbsgqqg(2,2,-4:0), na2carbsgqgg(2,2,-4:0), na2carbsgqgq(2,2,-4:0), na2carbsgqqq(2,2,-4:0),
     9         na2carbsqqqg(2,2,-4:0), na2carbsqqgg(2,2,-4:0), na2carbsqqgq(2,2,-4:0), na2carbsqqqq(2,2,-4:0)
      common/a2cts/na2srsqq, na2srsgg, na2carssrsqq, na2carssrsgg, na2carbssrsgg,
     1  na2carsggg, na2carsqqg, na2carsqqq, na2carsggq, 
     2  na2carsqqp, na2carsqgg, na2carsgqq, na2carsqpp,
     3  na2carbsqgqg, na2carbsqggg, na2carbsqggq, na2carbsqgqq,
     4  na2carbsggqg, na2carbsgggg, na2carbsgggq, na2carbsggqq,
     5  na2carbsgqqg, na2carbsgqgg, na2carbsgqgq, na2carbsgqqq,
     6  na2carbsqqqg, na2carbsqqgg, na2carbsqqgq, na2carbsqqqq

      real(ki) na1carsr0g(2,2,-2:2),
     1         na1ncar0qg(2,2,-2:2), na1ncar0gg(2,2,-2:2), na1ncar0gq(2,2,-2:2), na1ncar0qq(2,2,-2:2),
     2         na1sr0g(2,2,-2:2)
      common/a1cts/na1carsr0g,
     1         na1ncar0qg, na1ncar0gg, na1ncar0gq, na1ncar0qq,
     2         na1sr0g
      real(ki) nrva1carsr1g(0:2,0:2,-4:0),
     1         nrva1ncar1qg(2,2,-4:0), nrva1ncar1gg(2,2,-4:0), nrva1ncar1gq(2,2,-4:0), nrva1ncar1qq(2,2,-4:0),
     2         nrva1sr1g(0:2,0:2,-4:0)
      common/rva1cts/nrva1carsr1g,
     1         nrva1ncar1qg, nrva1ncar1gg, nrva1ncar1gq, nrva1ncar1qq,
     2         nrva1sr1g
      real(ki) na1c1carqgqg(2,2,-2:2), na1c1carqggg(2,2,-2:2), na1c1carqggq(2,2,-2:2), na1c1carqgqq(2,2,-2:2),
     1           na1c1carggqg(2,2,-2:2), na1c1cargggg(2,2,-2:2), na1c1cargggq(2,2,-2:2), na1c1carggqq(2,2,-2:2),
     2           na1c1cargqqg(2,2,-2:2), na1c1cargqgg(2,2,-2:2), na1c1cargqgq(2,2,-2:2), na1c1cargqqq(2,2,-2:2),
     3           na1c1carqqqg(2,2,-2:2), na1c1carqqgg(2,2,-2:2), na1c1carqqgq(2,2,-2:2), na1c1carqqqq(2,2,-2:2),
     4           na1c1cbrqgqg(2,2,-2:2), na1c1cbrqggg(2,2,-2:2), na1c1cbrqggq(2,2,-2:2), na1c1cbrqgqq(2,2,-2:2),
     5           na1c1cbrggqg(2,2,-2:2), na1c1cbrgggg(2,2,-2:2), na1c1cbrgggq(2,2,-2:2), na1c1cbrggqq(2,2,-2:2),
     6           na1c1cbrgqqg(2,2,-2:2), na1c1cbrgqgg(2,2,-2:2), na1c1cbrgqgq(2,2,-2:2), na1c1cbrgqqq(2,2,-2:2),
     7           na1c1cbrqqqg(2,2,-2:2), na1c1cbrqqgg(2,2,-2:2), na1c1cbrqqgq(2,2,-2:2), na1c1cbrqqqq(2,2,-2:2)
      common/a1c1cts/ na1c1carqgqg, na1c1carqggg, na1c1carqggq, na1c1carqgqq,
     1           na1c1carggqg, na1c1cargggg, na1c1cargggq, na1c1carggqq,
     2           na1c1cargqqg, na1c1cargqgg, na1c1cargqgq, na1c1cargqqq,
     3           na1c1carqqqg, na1c1carqqgg, na1c1carqqgq, na1c1carqqqq,
     4           na1c1cbrqgqg, na1c1cbrqggg, na1c1cbrqggq, na1c1cbrqgqq,
     5           na1c1cbrggqg, na1c1cbrgggg, na1c1cbrgggq, na1c1cbrggqq,
     6           na1c1cbrgqqg, na1c1cbrgqgg, na1c1cbrgqgq, na1c1cbrgqqq,
     7           na1c1cbrqqqg, na1c1cbrqqgg, na1c1cbrqqgq, na1c1cbrqqqq

      real(ki)
     1  na1a1carcasggg(2,2,-4:0), na1a1carcasggq(2,2,-4:0), na1a1carcasgqg(2,2,-4:0), na1a1carcasgqq(2,2,-4:0), 
     2  na1a1carcasqgg(2,2,-4:0), na1a1carcasqgq(2,2,-4:0), na1a1carcasqqg(2,2,-4:0), na1a1carcasqqq(2,2,-4:0),
     3  na1a1cbrcasqgqg(2,2,-4:0), na1a1cbrcasqggg(2,2,-4:0), na1a1cbrcasqggq(2,2,-4:0), na1a1cbrcasqgqq(2,2,-4:0),
     4  na1a1cbrcasggqg(2,2,-4:0), na1a1cbrcasgggg(2,2,-4:0), na1a1cbrcasgggq(2,2,-4:0), na1a1cbrcasggqq(2,2,-4:0),
     5  na1a1cbrcasgqqg(2,2,-4:0), na1a1cbrcasgqgg(2,2,-4:0), na1a1cbrcasgqgq(2,2,-4:0), na1a1cbrcasgqqq(2,2,-4:0),
     6  na1a1cbrcasqqqg(2,2,-4:0), na1a1cbrcasqqgg(2,2,-4:0), na1a1cbrcasqqgq(2,2,-4:0), na1a1cbrcasqqqq(2,2,-4:0),
     7  na1a1crscasggg(2,2,-4:0), na1a1crscasgqg(2,2,-4:0), na1a1crscasgqq(2,2,-4:0), 
     8  na1a1crscasqgg(2,2,-4:0), na1a1crscasqqg(2,2,-4:0), na1a1crscasqqq(2,2,-4:0),
     9  na1a1cassscarsragg(0:2,0:2,-4:0),
     1  na1a1cassssragg(0:2,0:2,-4:0),
     2  na1a1sscarsragg(0:2,0:2,-4:0),
     3  na1a1sssrgg(0:2,0:2,-4:0),
     4  na1a1crssscarggg(0:2,0:2,-4:0), na1a1crssscargqg(0:2,0:2,-4:0), na1a1crssscarqgg(0:2,0:2,-4:0), na1a1crssscarqqg(0:2,0:2,-4:0),
     5  na1a1sscarggg(0:2,0:2,-4:0), na1a1sscargqg(0:2,0:2,-4:0), na1a1sscarqgg(0:2,0:2,-4:0), na1a1sscarqqg(0:2,0:2,-4:0)

	
      common/a1a1cts/
     1  na1a1carcasggg, na1a1carcasggq, na1a1carcasgqg, na1a1carcasgqq, 
     2  na1a1carcasqgg, na1a1carcasqgq, na1a1carcasqqg, na1a1carcasqqq,
     3  na1a1cbrcasqgqg, na1a1cbrcasqggg, na1a1cbrcasqggq, na1a1cbrcasqgqq,
     4  na1a1cbrcasggqg, na1a1cbrcasgggg, na1a1cbrcasgggq, na1a1cbrcasggqq,
     5  na1a1cbrcasgqqg, na1a1cbrcasgqgg, na1a1cbrcasgqgq, na1a1cbrcasgqqq,
     6  na1a1cbrcasqqqg, na1a1cbrcasqqgg, na1a1cbrcasqqgq, na1a1cbrcasqqqq,
     7  na1a1crscasggg, na1a1crscasgqg, na1a1crscasgqq, 
     8  na1a1crscasqgg, na1a1crscasqqg, na1a1crscasqqq,
     9  na1a1cassscarsragg,
     1  na1a1cassssragg,
     2  na1a1sscarsragg,
     3  na1a1sssrgg,
     4  na1a1crssscarggg, na1a1crssscargqg, na1a1crssscarqgg, na1a1crssscarqqg,
     5  na1a1sscarggg, na1a1sscargqg, na1a1sscarqgg, na1a1sscarqqg


      real(ki)
     1  na12carssrscrsqq(0:2,0:2,-4:0), na12carssrscrsgg(0:2,0:2,-4:0),
     2  na12carssrsssgg(0:2,0:2,-4:0),
     3  na12srscasssgg(0:2,0:2,-4:0),
     4  na12srsssgg(0:2,0:2,-4:0),
     5  na12carssrscasssgg(0:2,0:2,-4:0),
     6  na12carscasssggg(0:2,0:2,-4:0), na12carscasssgqg(0:2,0:2,-4:0), 
     7  na12carscasssqgg(0:2,0:2,-4:0), na12carscasssqqg(0:2,0:2,-4:0),
     8  na12carscrsggg(2,2,-4:0), na12carscrsgqg(2,2,-4:0), na12carscrsgqq(2,2,-4:0), 
     9  na12carscrsqgg(2,2,-4:0), na12carscrsqqg(2,2,-4:0), na12carscrsqqq(2,2,-4:0),
     1  na12carsssggg(0:2,0:2,-4:0), na12carsssgqg(0:2,0:2,-4:0), na12carsssqgg(0:2,0:2,-4:0),
     2  na12carsssqqg(0:2,0:2,-4:0),
     3  na12srscrsgg(0:2,0:2,-4:0), na12srscrsqq(0:2,0:2,-4:0),
     4  na12carscrsssggg(0:2,0:2,-4:0), na12carscrsssgqg(0:2,0:2,-4:0), na12carscrsssqgg(0:2,0:2,-4:0),
     5  na12carscrsssqqg(0:2,0:2,-4:0),
     6  na12carscasggg(2,2,-4:0), na12carscasggq(2,2,-4:0), na12carscasgqg(2,2,-4:0), na12carscasgqq(2,2,-4:0),
     7  na12carscasqgg(2,2,-4:0), na12carscasqgq(2,2,-4:0), na12carscasqqg(2,2,-4:0), na12carscasqqq(2,2,-4:0),
     8  na12casbrcasgggg(2,2,-4:0), na12casbrcasgqgg(2,2,-4:0), na12casbrcasqggg(2,2,-4:0), na12casbrcasqqgg(2,2,-4:0),
     9  na12casbrcasgggq(2,2,-4:0), na12casbrcasgqgq(2,2,-4:0), na12casbrcasqggq(2,2,-4:0), na12casbrcasqqgq(2,2,-4:0),
     1  na12casbrcasggqg(2,2,-4:0), na12casbrcasgqqg(2,2,-4:0), na12casbrcasqgqg(2,2,-4:0), na12casbrcasqqqg(2,2,-4:0),
     2  na12casbrcasggqq(2,2,-4:0), na12casbrcasgqqq(2,2,-4:0), na12casbrcasqgqq(2,2,-4:0), na12casbrcasqqqq(2,2,-4:0)

      common/a12cts/
     1  na12carssrscrsqq, na12carssrscrsgg,
     2  na12carssrsssgg,
     3  na12srscasssgg,
     4  na12srsssgg,
     5  na12carssrscasssgg,
     6  na12carscasssggg, na12carscasssgqg, 
     7  na12carscasssqgg, na12carscasssqqg,
     8  na12carscrsggg, na12carscrsgqg, na12carscrsgqq, 
     9  na12carscrsqgg, na12carscrsqqg, na12carscrsqqq,
     1  na12carsssggg, na12carsssgqg, na12carsssqgg,
     2  na12carsssqqg,
     3  na12srscrsgg, na12srscrsqq,
     4  na12carscrsssggg, na12carscrsssgqg, na12carscrsssqgg,
     5  na12carscrsssqqg,
     6  na12carscasggg, na12carscasggq, na12carscasgqg, na12carscasgqq,
     7  na12carscasqgg, na12carscasqgq, na12carscasqqg, na12carscasqqq,
     8  na12casbrcasgggg, na12casbrcasgqgg, na12casbrcasqggg, na12casbrcasqqgg,
     9  na12casbrcasgggq, na12casbrcasgqgq, na12casbrcasqggq, na12casbrcasqqgq,
     1  na12casbrcasggqg, na12casbrcasgqqg, na12casbrcasqgqg, na12casbrcasqqqg,
     2  na12casbrcasggqq, na12casbrcasgqqq, na12casbrcasqgqq, na12casbrcasqqqq


	
