mm.to.cfs = function(mm, km2){
  # mm/day * km2 * (1m/1000mm) * (1000000m2/1km2) * (1day/86400s) * (35.314667ft3/1m3)
  # Simplified constant: (1/1000) * 1000000 * (1/86400) * 35.314667
  # = 1000 * 35.314667 / 86400 = 35314.667 / 86400 = 0.4087569
  cfs = mm * km2 * 0.4087569
  return(cfs)
}
