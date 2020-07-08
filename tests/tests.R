# 
#
# differents R tests
#
# 

# folder list test 

d = list.dirs("test", recursive = FALSE)
d

expanded = NULL

for(dir in d) {
  dd = list.dirs(dir, recursive = FALSE)
  cat(dd)
  expanded = append(expanded, dd)
}

f = list.files("test", recursive = TRUE, include.dirs = TRUE)
f

