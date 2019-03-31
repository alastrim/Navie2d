SOURCES += \
        main.cpp \
        misc.cpp \
    printers.cpp \
    dicrete_function.cpp \
    grid.cpp \
    matvec.cpp

HEADERS += \
    misc.h \
    printers.h \
    dicrete_function.h \
    grid.h \
    matvec.h

QMAKE_CXXFLAGS += -llaspack -lxc -g -W -Wall -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -fstack-protector-all -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wcast-align -Wno-long-long -Wmissing-declarations -O2 -ffast-math -lm -std=c++14
QMAKE_LFLAGS += -llaspack -lxc
LIBS += -llaspack -lxc
