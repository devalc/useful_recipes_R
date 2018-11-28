cfs.to.mm = function(cfs, km2){
    # ft3/sec * 1m3/35.314667ft3 * 1/km2 * 86400sec/1day * 1km2/1000000m2 * 1000mm/1m
    mm = cfs * (1/km2) * 2.446576
    return(mm)
}

# Test
# cfs.to.mm(1, 0.00404686)