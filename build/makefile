

# Compiler
CXX := g++
CXXFLAGS := -Wall -Wextra -std=c++26

# libs
LIBS := 
LDFLAGS += $(LIBS)

# Build mode
MODE ?= normal
RUN := .\

ifeq ($(MODE),debug)
    CXXFLAGS += -ggdb -O0
    RUN := gdb 
else ifeq ($(MODE),normal)
    CXXFLAGS += -O1
    RUN := .\\
else ifeq ($(MODE),release)
    CXXFLAGS += -O3 -march=native -funroll-loops -flto \
                -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIE -DNDEBUG
    LDFLAGS  += -flto
    RUN := .\\
else
    $(error Unknown build mode '$(MODE)'. Valid options: normal, debug, release)
endif

STATIC ?= false
ifeq ($(STATIC),true)
    LDFLAGS += -static -static-libgcc -static-libstdc++
endif

# Directories
SRC_DIR := ../src
BUILD_DIR := .
OBJ_DIR := $(BUILD_DIR)/.obj
TARGET := CUXMP$(MODE)

INCLUDES += -I../src -I../include -IC:/msys64/ucrt64/include -I C:/msys64/ucrt64/include -I C:/msys64/ucrt64/include/c++/12.2.0 -I C:/msys64/ucrt64/include/c++/12.2.0/x86_64-w64-mingw32 -I C:/msys64/ucrt64/x86_64-w64-mingw32/include
CXXFLAGS += $(INCLUDES)

# Find all source files
SRCS := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*/*.cpp) $(wildcard $(SRC_DIR)/*/*/*.cpp)

# Map source files to object files in OBJ_DIR
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

# ANSI color codes
GREEN := \033[32m
RESET := \033[0m

# Default target
all: $(OBJ_DIR) $(TARGET)

# Link executable
$(TARGET): $(OBJS)
	@echo -e "$(GREEN)Linking $(TARGET)...$(RESET)"
	$(CXX) $^ $(LDFLAGS) -o $@

# Compile object files and generate dependencies
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	@echo -e "$(GREEN)Compiling $<...$(RESET)"
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Include dependency files
-include $(DEPS)

# Ensure build directories exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Clean
clean:
	@echo -e "$(GREEN)Cleaning...$(RESET)"
	rm -rf $(OBJ_DIR) $(TARGET)

run: all
	@echo -e "$(GREEN)Running $(TARGET)...$(RESET)"
	$(RUN)$(TARGET)

.PHONY: all clean run
