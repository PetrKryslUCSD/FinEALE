function  styl =el_to_style( description)
switch description
    case 'H64';
        styl='k*--';
    case 'H27';
        styl='ks--';
    case 'T10';% tetrahedron
        styl='k^--';
    case 'T10-SRI';
        styl='kv-';
    case 'H20R';
        styl='ro--';
    case 'H8';
        styl='md--';
    case 'H8-SRI';
        styl='md-';
    otherwise
        styl='k-';
end