def calc_f_stat(RSS_R,RSS_U,n,pu,pr) :
    f_GC = ( (RSS_R-RSS_U)/ (RSS_U)) * ((n-pu-1)/(pu-pr));
    return f_GC    
