cmake_minimum_required(VERSION 2.8.3)

project(topp)

set(CMAKE_CXX_FLAGS "-std=c++14")
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -fPIC")

find_package(Eigen3 REQUIRED)
#find_package(ompl REQUIRED)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  geometry_msgs
  sensor_msgs
  visualization_msgs
)

include_directories(
    ${catkin_INCLUDE_DIRS}
    #${OMPL_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    include
)

catkin_package()

add_executable(topp_test src/test.cpp)

target_link_libraries(topp_test
  #${OMPL_LIBRARIES}
  ${catkin_LIBRARIES}
)


add_executable(cubic_test src/cubic_test.cpp)

target_link_libraries(cubic_test
  #${OMPL_LIBRARIES}
  ${catkin_LIBRARIES}
)


