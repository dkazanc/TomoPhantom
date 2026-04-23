cmake_minimum_required(VERSION 3.10)

project(tomophantom)
set(TOMOPHANTOM_VERSION "$ENV{TOMOPHANTOM_VERSION}" CACHE INTERNAL "TomoPhantom Library version" FORCE)
add_subdirectory(Core)
