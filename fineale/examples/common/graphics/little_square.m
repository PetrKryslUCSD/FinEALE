
function  litt

        [fens,fes]= T3_block(1.0,2.0,4,5,[]);
        
vtk_export_mesh ('theVTKFile',fes.conn,fens.xyz3,5,[]);
if paraview()
    paraview('theVTKFile.vtk')
end
  
end