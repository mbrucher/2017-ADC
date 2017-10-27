cmake_minimum_required(VERSION 3.5)

project(ADC2017)

find_package(ATK)
find_package(Boost REQUIRED unit_test_framework system)

FILE(GLOB SRC *.cpp *.h)

add_executable(ATKADC ${SRC})
link_libraries(ATKADC ATKEQ ATKTools ATKMock ATKCore ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
