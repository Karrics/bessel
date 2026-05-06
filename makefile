# Компилятор: clang++ для macOS, g++ для Linux
CXX ?= clang++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
LDFLAGS ?=

TARGET = bessel
SRC = bessel.cpp
OBJ = $(SRC:.cpp=.o)

# GSL
GSL_CFLAGS := $(shell pkg-config --cflags gsl 2>/dev/null)
GSL_LIBS   := $(shell pkg-config --libs gsl 2>/dev/null)

# OpenMP (macOS + libomp)
USE_OMP ?= 1
ifeq ($(USE_OMP),1)
  ifeq ($(shell uname -s),Darwin)
    # Apple Silicon
    ifneq ($(wildcard /opt/homebrew/opt/libomp/lib/libomp.dylib),)
      OMP_CFLAGS = -Xclang -fopenmp -I/opt/homebrew/opt/libomp/include
      OMP_LIBS   = -lomp -L/opt/homebrew/opt/libomp/lib
      LDFLAGS   += -Wl,-rpath,/opt/homebrew/opt/libomp/lib
    # Intel Mac
    else ifneq ($(wildcard /usr/local/opt/libomp/lib/libomp.dylib),)
      OMP_CFLAGS = -Xclang -fopenmp -I/usr/local/opt/libomp/include
      OMP_LIBS   = -lomp -L/usr/local/opt/libomp/lib
      LDFLAGS   += -Wl,-rpath,/usr/local/opt/libomp/lib
    endif
  else
    # Linux
    OMP_CFLAGS = -fopenmp
    OMP_LIBS   = -fopenmp
  endif
else
  OMP_CFLAGS =
  OMP_LIBS   =
endif

CXXFLAGS += $(GSL_CFLAGS) $(OMP_CFLAGS)
LDLIBS   += $(GSL_LIBS) $(OMP_LIBS)

.PHONY: all clean run debug help

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET) *.csv

run: $(TARGET)
	./$(TARGET)

debug: CXXFLAGS += -g -O0 -DDEBUG
debug: clean $(TARGET)

help:
	@echo "📦 Цели:"
	@echo "  make            - Сборка с OpenMP (по умолчанию)"
	@echo "  make USE_OMP=0  - Сборка без параллелизма"
	@echo "  make run        - Запуск"
	@echo "  make debug      - Отладочная сборка"
	@echo "  make clean      - Очистка"
	@echo "  make CXX=g++    - Использовать g++ вместо clang++"