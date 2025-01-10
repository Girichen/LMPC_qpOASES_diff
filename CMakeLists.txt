cmake_minimum_required(VERSION 3.0.2)
project(mpc)

find_package(OsqpEigen REQUIRED)
find_package(Eigen3 REQUIRED)
set(QPOASES_INLUDE_DIRS /usr/local/include/qpOASES)
find_package(Boost REQUIRED COMPONENTS system thread)
find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  geometry_msgs
  tf2
  tf2_ros
  message_generation
)

 add_message_files(
   FILES
   TrajectoryPointMsg.msg
   TrajectoryMsg.msg

 )

 generate_messages(
   DEPENDENCIES
   geometry_msgs
   std_msgs
 )

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES mpc
  CATKIN_DEPENDS roscpp rospy std_msgs geometry_msgs message_runtime
#  DEPENDS system_lib 
)


include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/qpOASES/include
)
include_directories(${OsqpEigen_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})


add_executable(mpc src/mpc_node.cpp)

 add_dependencies(mpc ${PROJECT_NAME}_generate_messages_cpp)


 target_link_libraries(mpc
  ${catkin_LIBRARIES}
  ${OsqpEigen_LIBRARIES}
  Eigen3::Eigen
  Boost::system Boost::thread
  ${CMAKE_CURRENT_SOURCE_DIR}/qpOASES/build/libs/libqpOASES.a
 )
