cmake_minimum_required(VERSION 3.18)

# Project name
project(pyrbd_plusplus)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add subdirectory for core components
add_subdirectory(pyrbd_plusplus/_core/)