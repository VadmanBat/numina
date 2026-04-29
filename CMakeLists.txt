cmake_minimum_required(VERSION 3.28)
project(numina VERSION 1.0 LANGUAGES CXX)

# ========== Настройка C++ ==========
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(STATIC_BUILD OFF)

# Подключаем заголовки Eigen
set(EIGEN3_INCLUDE_DIR "C:/programming/libraries/eigen-5.0.0")
include_directories(${EIGEN3_INCLUDE_DIR})

# ========== Источники и исполняемый файл ==========
add_executable(${PROJECT_NAME}
        main.cpp
        #numina/classes/poly-solver/tests/poly-solver-test.cpp

        numina/classes/polynomial/src/core/remove-leading-zeros.cpp
        numina/classes/polynomial/src/operators/analytical.cpp
        numina/classes/polynomial/src/operators/arithmetic.cpp
        numina/classes/polynomial/src/core/assignment.cpp
        numina/classes/polynomial/src/operators/bool-ops.cpp
        numina/classes/polynomial/src/operators/io-ops.cpp
        numina/classes/polynomial/src/operators/scalar-ops.cpp
        numina/classes/polynomial/src/roots/multiplicity.cpp
        numina/classes/polynomial/src/roots/root.cpp
        numina/classes/polynomial/src/compose.cpp
        numina/classes/polynomial/src/divmod.cpp
        numina/classes/polynomial/src/gcd.cpp

        numina/classes/poly-solver/src/core/deflate.cpp
        numina/classes/poly-solver/src/core/prepare.cpp
        numina/classes/poly-solver/src/methods/laguerre.cpp
        numina/classes/poly-solver/src/methods/newton.cpp
        numina/classes/poly-solver/src/solve/poly-solver.cpp
        numina/classes/poly-solver/src/solve/solve-api.cpp

        numina/classes/quadratic-deriv-poly/src/quadratic-deriv-poly.cpp
)

# ========== Красивые include ==========
target_include_directories(${PROJECT_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/numina/classes/laplace-expression/include
        ${CMAKE_CURRENT_SOURCE_DIR}/numina/classes/poly-solver/include
        ${CMAKE_CURRENT_SOURCE_DIR}/numina/classes/polynomial/include
        ${CMAKE_CURRENT_SOURCE_DIR}/numina/classes/quadratic-deriv-poly/include
        ${CMAKE_CURRENT_SOURCE_DIR}/numina/classes/transfer-function/include
)
