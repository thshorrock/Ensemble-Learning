#Build a subset of gsl

project(gsl)

# configure a header file to pass some of the CMake settings
# to the source code
 message ("config file at ${PROJECT_SOURCE_DIR}/config.h.in")
configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )
INCLUDE (CheckIncludeFiles)
include (CheckFunctionExists)
check_function_exists (acosh HAVE_DECL_ACOSH)
check_function_exists (asinh HAVE_DECL_ASINH)
check_function_exists (atanh HAVE_DECL_ATANH)
check_function_exists (expm1 HAVE_DECL_EXPM1)
check_function_exists (feenableexcept HAVE_DECL_FEENABLEEXCEPT)
check_function_exists (fesettrapenable HAVE_DECL_FESETTRAPENABLE)
check_function_exists (finite HAVE_DECL_FINITE)
check_function_exists (frexp HAVE_DECL_FREXP)
check_function_exists (hypot HAVE_DECL_HYPOT)
check_function_exists (isfinite HAVE_DECL_ISFINITE)
check_function_exists (isinf HAVE_DECL_ISINF)
check_function_exists (isnan HAVE_DECL_ISNAN)
check_function_exists (ldexp HAVE_DECL_LDEXP)
check_function_exists (log1p HAVE_DECL_LOG1P)

CHECK_INCLUDE_FILES(dlfcn.h HAVE_DLFCN_H)
check_function_exists (memcpy HAVE_MEMCPY)
check_function_exists (memmove HAVE_MEMMOVE)
CHECK_INCLUDE_FILES(memory.h HAVE_MEMORY_H)
#CHECK_INCLUDE_FILES( HAVE_PRINTF_LONGDOUBLE)
CHECK_INCLUDE_FILES(stdint.h HAVE_STDINT_H)
CHECK_INCLUDE_FILES(stdlib.h HAVE_STDLIB_H)
check_function_exists (strdup HAVE_STRDUP)
CHECK_INCLUDE_FILES(strings.h HAVE_STRINGS_H)
CHECK_INCLUDE_FILES(string.h HAVE_STRING_H)
check_function_exists (strtol HAVE_STRTOL)
check_function_exists (strtoul HAVE_STRTOUL)
CHECK_INCLUDE_FILES(sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(unistd.h HAVE_UNISTD_H)

add_subdirectory(sys)
add_subdirectory(test)