cmake_minimum_required(VERSION 3.5)

project(ADC2017)

find_package(ATK)
find_package(Boost REQUIRED unit_test_framework system)

add_definitions("-DBOOST_ALL_NO_LIB")

FILE(GLOB SRC *.cpp *.h)

include_directories(${Boost_INCLUDE_DIR})
include_directories(${ATK_INCLUDE_DIRS})

add_executable(ATKADC ${SRC})
target_link_libraries(ATKADC ${ATK_EQ_LIBRARY} ${ATK_MOCK_LIBRARY} ${ATK_TOOLS_LIBRARY} ${ATK_CORE_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
