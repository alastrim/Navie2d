SOURCES += \
        main.cpp \
        misc.cpp \
    printers.cpp \
    grid.cpp \
    matvec.cpp \
    discrete_function.cpp

HEADERS += \
    misc.h \
    printers.h \
    grid.h \
    matvec.h \
    discrete_function.h

QMAKE_CXXFLAGS += -llaspack -lxc -g -W -Wall -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -fstack-protector-all -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wcast-align -Wno-long-long -Wmissing-declarations -O2 -ffast-math -lm -std=c++14
QMAKE_LFLAGS += -llaspack -lxc
LIBS += -llaspack -lxc
