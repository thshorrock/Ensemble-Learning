
set(MARCH_FLAGS)  
if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX ) 
  #let gcc work it out
  set(MARCH_FLAGS "-march=native")
endif()