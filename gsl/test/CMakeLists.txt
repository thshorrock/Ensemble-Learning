
#build the gsl test

project(gsltest)
INCLUDE_DIRECTORIES( ${gsl_INCLUDE_DIRS} "../")
ADD_LIBRARY(gsltest results.c)
install (TARGETS gsltest DESTINATION lib)
install (FILES gsl_test.h DESTINATION include)