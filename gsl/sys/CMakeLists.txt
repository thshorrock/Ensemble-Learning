
#build the gsl sys

project(gslsys)

set(gslsys_SRC  minmax.c prec.c hypot.c log1p.c expm1.c coerce.c invhyp.c pow_int.c infnan.c fdiv.c fcmp.c ldfrexp.c)

INCLUDE_DIRECTORIES( ${gsl_INCLUDE_DIRS} "../")
ADD_LIBRARY(gslsys ${gslsys_SRC})
install (TARGETS gslsys DESTINATION lib)