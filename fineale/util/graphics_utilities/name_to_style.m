function styl= name_to_style(description)
switch description
    
    case  {'ref','reference'};
        styl='k-';
           
    case  'QT10MS';
        styl='rh-';
           
    case  'H8MSGS';
        styl='kx-';
            
    case  'H8MSGSO';
        styl='ko-';
           
    case  'H8MSPES';
        styl='rp--';
    
    case  'H8MSPSS';
        styl='m*-';
                  
    case  '1pt-H8-nodalstab';
        styl='rx-';
        
    case  'THEX';
        styl='r<-';
    case  'THEX-L';
        styl='g<:';
    case  'THEX-Q5';
        styl='b<--';
    case  'THEX-Q10';
        styl='k<-.';
        
    case  'NICE-H8';
        styl='mh-';
        
    case  'H8';
        styl='rh-';
    case  'H8-SRI';
        styl='rs-';
    case  'H8-GSRI';
        styl='ms--';
    case  'H8-Bbar';
        styl='b*-';
    case  'H8-Bbar-ISO';
        styl='bx:';
     case  'H8-BbarX';
        styl='bh-.';
     case  'H8-3FS';
        styl='kh-.';
        
    case  'H20';
        styl='bo-';
    case  'H20R';
        styl='ro--';
    case  'H20-GSRI';
        styl='bo--';
    case  'H20-Bbar';
        styl='co-.';
    case  'H20-BbarX';
        styl='ch-.';
    case  'H20-3FS';
        styl='ch-';
    case  'H20-Bbar-ISO';
        styl='cx:';
       
    case  'H27';
        styl='cd-';
    case  'H27-GSRI';
        styl='cd--';
    case  'H27-Bbar';
        styl='cd-.';
    case  'H27-BbarX';
        styl='ch-.';
    case  'H27-3FS';
        styl='cd-';
        
    case  'H64';
        styl='kp-';
    case  'H64-GSRI';
        styl='kp--';
    case  'H64-Bbar';
        styl='kp-.';
        
    case  'T10';
        styl='k^--';
    case  'T10-GSRI';
        styl='kv-';
    case  'T10-Bbar';
        styl='kv-.';
        
    
    case  'C3D8';
        styl='ks--';
    case  'C3D8H';
        styl='kd-.';
    case  'C3D8I';
        styl='kd-';
    case  'C3D8IH';
        styl='kd:';
    case  'C3D8RH';
        styl='kd--';
            
    case  'C3D20';
        styl='ms--';
    case  'C3D20H';
        styl='md-.';
    case  'C3D20R';
        styl='md-';
    case  'C3D20RH';
        styl='md--';
             
    otherwise
        styl='kx-';
        
end
end
