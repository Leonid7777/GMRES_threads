target_sources(${CMAKE_PROJECT_NAME}
  PRIVATE
  mat_methods.cpp
  vec_methods.cpp
)

find_package(Threads REQUIRED)
target_link_libraries(${CMAKE_PROJECT_NAME} Threads::Threads)