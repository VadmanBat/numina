cmake_minimum_required(VERSION 3.28)
project(numina VERSION 1.0 LANGUAGES CXX)

# ========== Настройка C++ ==========
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(STATIC_BUILD OFF)

# ========== Источники и исполняемый файл ==========
add_executable(${PROJECT_NAME}
        main.cpp

        numina/classes/polynomial/src/operators/analytical.cpp
        numina/classes/polynomial/src/operators/assignment.cpp
        numina/classes/polynomial/src/operators/arithmetic.cpp
        numina/classes/polynomial/src/operators/scalar-ops.cpp
        numina/classes/polynomial/src/divmod.cpp

        numina/classes/polynomial/src/internal/remove-leading-zeros.cpp
        numina/classes/polynomial/src/operators/bool-ops.cpp
        numina/classes/polynomial/src/operators/io-ops.cpp
        numina/classes/polynomial/src/gcd.cpp
        numina/classes/polynomial/src/compose.cpp
        numina/classes/polynomial/src/roots/root.cpp
        numina/classes/poly-solver/src/internal/calculate.cpp
        numina/classes/poly-solver/src/solve.cpp
        numina/classes/poly-solver/src/internal/deflate.cpp
        numina/classes/poly-solver/src/internal/state.cpp
        numina/classes/poly-solver/src/internal/solve.cpp
        numina/classes/polynomial/src/roots/multiplicity.cpp
)

# ========== Красивые include ==========
target_include_directories(${PROJECT_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/numina/classes/polynomial/include
        ${CMAKE_CURRENT_SOURCE_DIR}/numina/classes/poly-solver/include
        ${CMAKE_CURRENT_SOURCE_DIR}/numina/classes/transfer-function/include
)