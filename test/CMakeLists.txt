add_executable(test
  algorithm_gmres.cpp
  create_mat_r.cpp
  create_mat_r.h
)

target_link_libraries(test
${CMAKE_PROJECT_NAME}
)